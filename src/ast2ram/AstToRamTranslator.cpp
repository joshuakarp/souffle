/*
 * Souffle - A Datalog Compiler
 * Copyright (c) 2013, 2015, Oracle and/or its affiliates. All rights reserved
 * Licensed under the Universal Permissive License v 1.0 as shown at:
 * - https://opensource.org/licenses/UPL
 * - <souffle root>/licenses/SOUFFLE-UPL.txt
 */

/************************************************************************
 *
 * @file AstToRamTranslator.cpp
 *
 * Translator from AST to RAM structures.
 *
 ***********************************************************************/

#include "ast2ram/AstToRamTranslator.h"
#include "Global.h"
#include "LogStatement.h"
#include "ast/Aggregator.h"
#include "ast/Argument.h"
#include "ast/Atom.h"
#include "ast/BinaryConstraint.h"
#include "ast/BranchInit.h"
#include "ast/Clause.h"
#include "ast/Constant.h"
#include "ast/Constraint.h"
#include "ast/Directive.h"
#include "ast/Functor.h"
#include "ast/Literal.h"
#include "ast/Negation.h"
#include "ast/NilConstant.h"
#include "ast/Node.h"
#include "ast/NumericConstant.h"
#include "ast/Program.h"
#include "ast/QualifiedName.h"
#include "ast/RecordInit.h"
#include "ast/Relation.h"
#include "ast/StringConstant.h"
#include "ast/SubroutineArgument.h"
#include "ast/TranslationUnit.h"
#include "ast/UnnamedVariable.h"
#include "ast/Variable.h"
#include "ast/analysis/AuxArity.h"
#include "ast/analysis/Functor.h"
#include "ast/analysis/IOType.h"
#include "ast/analysis/PolymorphicObjects.h"
#include "ast/analysis/RecursiveClauses.h"
#include "ast/analysis/RelationDetailCache.h"
#include "ast/analysis/RelationSchedule.h"
#include "ast/analysis/SCCGraph.h"
#include "ast/analysis/SumTypeBranches.h"
#include "ast/analysis/TopologicallySortedSCCGraph.h"
#include "ast/analysis/TypeEnvironment.h"
#include "ast/analysis/TypeSystem.h"
#include "ast/transform/Transformer.h"
#include "ast/utility/NodeMapper.h"
#include "ast/utility/SipsMetric.h"
#include "ast/utility/Utils.h"
#include "ast/utility/Visitor.h"
#include "ast2ram/ClauseTranslator.h"
#include "ast2ram/ConstraintTranslator.h"
#include "ast2ram/Location.h"
#include "ast2ram/ProvenanceClauseTranslator.h"
#include "ast2ram/ValueIndex.h"
#include "ast2ram/ValueTranslator.h"
#include "ram/Call.h"
#include "ram/Clear.h"
#include "ram/Condition.h"
#include "ram/Conjunction.h"
#include "ram/Constraint.h"
#include "ram/DebugInfo.h"
#include "ram/EmptinessCheck.h"
#include "ram/ExistenceCheck.h"
#include "ram/Exit.h"
#include "ram/Expression.h"
#include "ram/Extend.h"
#include "ram/Filter.h"
#include "ram/FloatConstant.h"
#include "ram/IO.h"
#include "ram/LogRelationTimer.h"
#include "ram/LogSize.h"
#include "ram/LogTimer.h"
#include "ram/Loop.h"
#include "ram/Negation.h"
#include "ram/Parallel.h"
#include "ram/Program.h"
#include "ram/Project.h"
#include "ram/Query.h"
#include "ram/Relation.h"
#include "ram/RelationSize.h"
#include "ram/Scan.h"
#include "ram/Sequence.h"
#include "ram/SignedConstant.h"
#include "ram/Statement.h"
#include "ram/SubroutineArgument.h"
#include "ram/SubroutineReturn.h"
#include "ram/Swap.h"
#include "ram/TranslationUnit.h"
#include "ram/TupleElement.h"
#include "ram/UndefValue.h"
#include "ram/UnsignedConstant.h"
#include "ram/utility/Utils.h"
#include "reports/DebugReport.h"
#include "reports/ErrorReport.h"
#include "souffle/BinaryConstraintOps.h"
#include "souffle/SymbolTable.h"
#include "souffle/TypeAttribute.h"
#include "souffle/utility/ContainerUtil.h"
#include "souffle/utility/FunctionalUtil.h"
#include "souffle/utility/MiscUtil.h"
#include "souffle/utility/StringUtil.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace souffle::ast2ram {

AstToRamTranslator::AstToRamTranslator() = default;
AstToRamTranslator::~AstToRamTranslator() = default;

/** append statement to a list of statements */
inline void appendStmt(VecOwn<ram::Statement>& stmtList, Own<ram::Statement> stmt) {
    if (stmt) {
        stmtList.push_back(std::move(stmt));
    }
}

Own<ram::TupleElement> AstToRamTranslator::makeRamTupleElement(const Location& loc) {
    return mk<ram::TupleElement>(loc.identifier, loc.element);
}

size_t AstToRamTranslator::getEvaluationArity(const ast::Atom* atom) const {
    std::string relName = atom->getQualifiedName().toString();
    if (isPrefix("@info_", relName)) return 0;

    // Get the original relation name
    if (isPrefix("@delta_", relName)) {
        relName = stripPrefix("@delta_", relName);
    } else if (isPrefix("@new_", relName)) {
        relName = stripPrefix("@new_", relName);
    }

    const auto* originalRelation = relDetail->getRelation(ast::QualifiedName(relName));
    return auxArityAnalysis->getArity(originalRelation);
}

std::vector<std::map<std::string, std::string>> AstToRamTranslator::getInputDirectives(
        const ast::Relation* rel) {
    std::vector<std::map<std::string, std::string>> inputDirectives;

    for (const auto* load : program->getDirectives()) {
        if (load->getQualifiedName() != rel->getQualifiedName() ||
                load->getType() != ast::DirectiveType::input) {
            continue;
        }

        std::map<std::string, std::string> directives;
        for (const auto& currentPair : load->getParameters()) {
            directives.insert(std::make_pair(currentPair.first, unescape(currentPair.second)));
        }
        inputDirectives.push_back(directives);
    }

    if (inputDirectives.empty()) {
        inputDirectives.emplace_back();
    }

    return inputDirectives;
}

std::vector<std::map<std::string, std::string>> AstToRamTranslator::getOutputDirectives(
        const ast::Relation* rel) {
    std::vector<std::map<std::string, std::string>> outputDirectives;

    for (const auto* store : program->getDirectives()) {
        if (store->getQualifiedName() != rel->getQualifiedName() ||
                (store->getType() != ast::DirectiveType::printsize &&
                        store->getType() != ast::DirectiveType::output)) {
            continue;
        }

        std::map<std::string, std::string> directives;
        for (const auto& currentPair : store->getParameters()) {
            directives.insert(std::make_pair(currentPair.first, unescape(currentPair.second)));
        }
        outputDirectives.push_back(directives);
    }

    if (outputDirectives.empty()) {
        outputDirectives.emplace_back();
    }

    return outputDirectives;
}

std::string AstToRamTranslator::getConcreteRelationName(const ast::Atom* atom) {
    return getRelationName(atom->getQualifiedName());
}

std::string AstToRamTranslator::getConcreteRelationName(
        const ast::Relation* rel, const std::string relationNamePrefix) {
    return relationNamePrefix + getRelationName(rel->getQualifiedName());
}

std::string AstToRamTranslator::getDeltaRelationName(const ast::Relation* rel) {
    return getConcreteRelationName(rel, "@delta_");
}

std::string AstToRamTranslator::getNewRelationName(const ast::Relation* rel) {
    return getConcreteRelationName(rel, "@new_");
}

Own<ram::Expression> AstToRamTranslator::translateValue(const ast::Argument* arg, const ValueIndex& index) {
    if (arg == nullptr) return nullptr;
    return ValueTranslator::translate(*this, index, getSymbolTable(), *arg);
}

SymbolTable& AstToRamTranslator::getSymbolTable() {
    static SymbolTable symbolTable;
    return symbolTable;
}

Own<ram::Condition> AstToRamTranslator::translateConstraint(
        const ast::Literal* lit, const ValueIndex& index) {
    assert(lit != nullptr && "literal should be defined");
    return ConstraintTranslator::translate(*this, index, *lit);
}

RamDomain AstToRamTranslator::getConstantRamRepresentation(const ast::Constant& constant) {
    if (auto strConstant = dynamic_cast<const ast::StringConstant*>(&constant)) {
        return getSymbolTable().lookup(strConstant->getConstant());
    } else if (isA<ast::NilConstant>(&constant)) {
        return 0;
    } else if (auto* numConstant = dynamic_cast<const ast::NumericConstant*>(&constant)) {
        assert(numConstant->getFinalType().has_value() && "constant should have valid type");
        switch (numConstant->getFinalType().value()) {
            case ast::NumericConstant::Type::Int:
                return RamSignedFromString(numConstant->getConstant(), nullptr, 0);
            case ast::NumericConstant::Type::Uint:
                return RamUnsignedFromString(numConstant->getConstant(), nullptr, 0);
            case ast::NumericConstant::Type::Float: return RamFloatFromString(numConstant->getConstant());
        }
    }

    fatal("unaccounted-for constant");
}

Own<ram::Expression> AstToRamTranslator::translateConstant(ast::Constant const& c) {
    auto const rawConstant = getConstantRamRepresentation(c);
    if (auto* const c_num = dynamic_cast<const ast::NumericConstant*>(&c)) {
        switch (c_num->getFinalType().value()) {
            case ast::NumericConstant::Type::Int: return mk<ram::SignedConstant>(rawConstant);
            case ast::NumericConstant::Type::Uint: return mk<ram::UnsignedConstant>(rawConstant);
            case ast::NumericConstant::Type::Float: return mk<ram::FloatConstant>(rawConstant);
        }
        fatal("unaccounted-for constant");
    }
    return mk<ram::SignedConstant>(rawConstant);
}

/** generate RAM code for a non-recursive relation */
Own<ram::Statement> AstToRamTranslator::translateNonRecursiveRelation(const ast::Relation& rel) {
    /* start with an empty sequence */
    VecOwn<ram::Statement> res;

    std::string relName = getConcreteRelationName(&rel);

    /* iterate over all clauses that belong to the relation */
    for (ast::Clause* clause : relDetail->getClauses(rel.getQualifiedName())) {
        // skip recursive rules
        if (recursiveClauses->recursive(clause)) {
            continue;
        }

        // translate clause
        Own<ram::Statement> rule = ClauseTranslator(*this).translateClause(*clause, *clause);

        // add logging
        if (Global::config().has("profile")) {
            const std::string& relationName = toString(rel.getQualifiedName());
            const auto& srcLocation = clause->getSrcLoc();
            const std::string clauseText = stringify(toString(*clause));
            const std::string logTimerStatement =
                    LogStatement::tNonrecursiveRule(relationName, srcLocation, clauseText);
            const std::string logSizeStatement =
                    LogStatement::nNonrecursiveRule(relationName, srcLocation, clauseText);
            rule = mk<ram::LogRelationTimer>(std::move(rule), logTimerStatement, relName);
        }

        // add debug info
        std::ostringstream ds;
        ds << toString(*clause) << "\nin file ";
        ds << clause->getSrcLoc();
        rule = mk<ram::DebugInfo>(std::move(rule), ds.str());

        // add rule to result
        appendStmt(res, std::move(rule));
    }

    // add logging for entire relation
    if (Global::config().has("profile")) {
        const std::string& relationName = toString(rel.getQualifiedName());
        const auto& srcLocation = rel.getSrcLoc();
        const std::string logSizeStatement = LogStatement::nNonrecursiveRelation(relationName, srcLocation);

        // add timer if we did any work
        if (!res.empty()) {
            const std::string logTimerStatement =
                    LogStatement::tNonrecursiveRelation(relationName, srcLocation);
            auto newStmt =
                    mk<ram::LogRelationTimer>(mk<ram::Sequence>(std::move(res)), logTimerStatement, relName);
            res.clear();
            appendStmt(res, std::move(newStmt));
        } else {
            // add table size printer
            appendStmt(res, mk<ram::LogSize>(relName, logSizeStatement));
        }
    }

    // done
    return mk<ram::Sequence>(std::move(res));
}

/**
 * A utility function assigning names to unnamed variables such that enclosing
 * constructs may be cloned without losing the variable-identity.
 */
void AstToRamTranslator::nameUnnamedVariables(ast::Clause* clause) {
    // the node mapper conducting the actual renaming
    struct Instantiator : public ast::NodeMapper {
        mutable int counter = 0;

        Instantiator() = default;

        Own<ast::Node> operator()(Own<ast::Node> node) const override {
            // apply recursive
            node->apply(*this);

            // replace unknown variables
            if (dynamic_cast<ast::UnnamedVariable*>(node.get()) != nullptr) {
                auto name = " _unnamed_var" + toString(++counter);
                return mk<ast::Variable>(name);
            }

            // otherwise nothing
            return node;
        }
    };

    // name all variables in the atoms
    Instantiator init;
    for (auto& atom : ast::getBodyLiterals<ast::Atom>(*clause)) {
        atom->apply(init);
    }
}

/** converts the given relation identifier into a relation name */
std::string AstToRamTranslator::getRelationName(const ast::QualifiedName& id) {
    return toString(join(id.getQualifiers(), "."));
}

Own<ram::Sequence> AstToRamTranslator::translateSCC(size_t scc, size_t idx) {
    // make a new ram statement for the current SCC
    VecOwn<ram::Statement> current;

    // find out if the current SCC is recursive
    const auto& isRecursive = sccGraph->isRecursive(scc);

    // make variables for particular sets of relations contained within the current SCC, and,
    // predecessors and successor SCCs thereof
    const auto& allInterns = sccGraph->getInternalRelations(scc);
    const auto& internIns = sccGraph->getInternalInputRelations(scc);
    const auto& internOuts = sccGraph->getInternalOutputRelations(scc);

    // make a variable for all relations that are expired at the current SCC
    const auto& internExps = relationSchedule->schedule().at(idx).expired();

    // load all internal input relations from the facts dir with a .facts extension
    for (const auto& relation : internIns) {
        makeRamLoad(current, relation);
    }

    // compute the relations themselves
    Own<ram::Statement> bodyStatement =
            (!isRecursive) ? translateNonRecursiveRelation(*((const ast::Relation*)*allInterns.begin()))
                           : translateRecursiveRelation(allInterns);
    appendStmt(current, std::move(bodyStatement));

    // store all internal output relations to the output dir with a .csv extension
    for (const auto& relation : internOuts) {
        makeRamStore(current, relation);
    }

    // if provenance is disabled, drop all relations expired as per the topological order
    if (!Global::config().has("provenance")) {
        for (const auto& relation : internExps) {
            makeRamClear(current, relation);
        }
    }

    return mk<ram::Sequence>(std::move(current));
}

/** generate RAM code for recursive relations in a strongly-connected component */
Own<ram::Statement> AstToRamTranslator::translateRecursiveRelation(
        const std::set<const ast::Relation*>& scc) {
    // initialize sections
    VecOwn<ram::Statement> preamble;
    VecOwn<ram::Statement> updateTable;
    VecOwn<ram::Statement> postamble;

    auto genMerge = [&](const ast::Relation* rel, const std::string& destRel,
                            const std::string& srcRel) -> Own<ram::Statement> {
        VecOwn<ram::Expression> values;
        if (rel->getArity() == 0) {
            return mk<ram::Query>(mk<ram::Filter>(mk<ram::Negation>(mk<ram::EmptinessCheck>(srcRel)),
                    mk<ram::Project>(destRel, std::move(values))));
        }
        for (std::size_t i = 0; i < rel->getArity(); i++) {
            values.push_back(mk<ram::TupleElement>(0, i));
        }
        auto stmt = mk<ram::Query>(mk<ram::Scan>(srcRel, 0, mk<ram::Project>(destRel, std::move(values))));
        if (rel->getRepresentation() == RelationRepresentation::EQREL) {
            return mk<ram::Sequence>(mk<ram::Extend>(destRel, srcRel), std::move(stmt));
        }
        return stmt;
    };

    // --- create preamble ---

    /* Compute non-recursive clauses for relations in scc and push
       the results in their delta tables. */
    for (const ast::Relation* rel : scc) {
        /* create update statements for fixpoint (even iteration) */
        Own<ram::Statement> updateRelTable =
                mk<ram::Sequence>(genMerge(rel, getConcreteRelationName(rel), getNewRelationName(rel)),
                        mk<ram::Swap>(getDeltaRelationName(rel), getNewRelationName(rel)),
                        mk<ram::Clear>(getNewRelationName(rel)));

        /* measure update time for each relation */
        if (Global::config().has("profile")) {
            updateRelTable = mk<ram::LogRelationTimer>(std::move(updateRelTable),
                    LogStatement::cRecursiveRelation(toString(rel->getQualifiedName()), rel->getSrcLoc()),
                    getNewRelationName(rel));
        }

        /* drop temporary tables after recursion */
        appendStmt(postamble, mk<ram::Clear>(getDeltaRelationName(rel)));
        appendStmt(postamble, mk<ram::Clear>(getNewRelationName(rel)));

        /* Generate code for non-recursive part of relation */
        /* Generate merge operation for temp tables */
        appendStmt(preamble, translateNonRecursiveRelation(*rel));
        appendStmt(preamble, genMerge(rel, getDeltaRelationName(rel), getConcreteRelationName(rel)));

        /* Add update operations of relations to parallel statements */
        appendStmt(updateTable, std::move(updateRelTable));
    }

    // --- build main loop ---

    VecOwn<ram::Statement> loopSeq;

    // create a utility to check SCC membership
    auto isInSameSCC = [&](const ast::Relation* rel) {
        return std::find(scc.begin(), scc.end(), rel) != scc.end();
    };

    /* Compute temp for the current tables */
    for (const ast::Relation* rel : scc) {
        VecOwn<ram::Statement> loopRelSeq;

        /* Find clauses for relation rel */
        for (const auto& cl : relDetail->getClauses(rel->getQualifiedName())) {
            // skip non-recursive clauses
            if (!recursiveClauses->recursive(cl)) {
                continue;
            }

            // each recursive rule results in several operations
            int version = 0;
            const auto& atoms = ast::getBodyLiterals<ast::Atom>(*cl);
            for (size_t j = 0; j < atoms.size(); ++j) {
                const ast::Atom* atom = atoms[j];
                const ast::Relation* atomRelation = getAtomRelation(atom, program);

                // only interested in atoms within the same SCC
                if (!isInSameSCC(atomRelation)) {
                    continue;
                }

                // modify the processed rule to use delta relation and write to new relation
                auto r1 = souffle::clone(cl);
                r1->getHead()->setQualifiedName(getNewRelationName(rel));
                ast::getBodyLiterals<ast::Atom>(*r1)[j]->setQualifiedName(getDeltaRelationName(atomRelation));
                if (Global::config().has("provenance")) {
                    r1->addToBody(mk<ast::ProvenanceNegation>(souffle::clone(cl->getHead())));
                } else if (r1->getHead()->getArity() > 0) {
                    r1->addToBody(mk<ast::Negation>(souffle::clone(cl->getHead())));
                }

                // replace wildcards with variables (reduces indices when wildcards are used in recursive
                // atoms)
                nameUnnamedVariables(r1.get());

                // reduce R to P ...
                for (size_t k = j + 1; k < atoms.size(); k++) {
                    if (isInSameSCC(getAtomRelation(atoms[k], program))) {
                        auto cur = souffle::clone(ast::getBodyLiterals<ast::Atom>(*r1)[k]);
                        cur->setQualifiedName(getDeltaRelationName(getAtomRelation(atoms[k], program)));
                        r1->addToBody(mk<ast::Negation>(std::move(cur)));
                    }
                }

                Own<ram::Statement> rule = ClauseTranslator(*this).translateClause(*r1, *cl, version);

                /* add logging */
                if (Global::config().has("profile")) {
                    const std::string& relationName = toString(rel->getQualifiedName());
                    const auto& srcLocation = cl->getSrcLoc();
                    const std::string clauseText = stringify(toString(*cl));
                    const std::string logTimerStatement =
                            LogStatement::tRecursiveRule(relationName, version, srcLocation, clauseText);
                    const std::string logSizeStatement =
                            LogStatement::nRecursiveRule(relationName, version, srcLocation, clauseText);
                    rule = mk<ram::LogRelationTimer>(
                            std::move(rule), logTimerStatement, getNewRelationName(rel));
                }

                // add debug info
                std::ostringstream ds;
                ds << toString(*cl) << "\nin file ";
                ds << cl->getSrcLoc();
                rule = mk<ram::DebugInfo>(std::move(rule), ds.str());

                // add to loop body
                appendStmt(loopRelSeq, std::move(rule));

                // increment version counter
                version++;
            }

            if (cl->getExecutionPlan() != nullptr) {
                // ensure that all required versions have been created, as expected
                int maxVersion = -1;
                for (auto const& cur : cl->getExecutionPlan()->getOrders()) {
                    maxVersion = std::max(cur.first, maxVersion);
                }
                assert(version > maxVersion && "missing clause versions");
            }
        }

        // if there was no rule, continue
        if (loopRelSeq.empty()) {
            continue;
        }

        // label all versions
        if (Global::config().has("profile")) {
            const std::string& relationName = toString(rel->getQualifiedName());
            const auto& srcLocation = rel->getSrcLoc();
            const std::string logTimerStatement = LogStatement::tRecursiveRelation(relationName, srcLocation);
            const std::string logSizeStatement = LogStatement::nRecursiveRelation(relationName, srcLocation);
            auto newStmt = mk<ram::LogRelationTimer>(
                    mk<ram::Sequence>(std::move(loopRelSeq)), logTimerStatement, getNewRelationName(rel));
            loopRelSeq.clear();
            appendStmt(loopRelSeq, std::move(newStmt));
        }

        /* add rule computations of a relation to parallel statement */
        appendStmt(loopSeq, mk<ram::Sequence>(std::move(loopRelSeq)));
    }
    auto loop = mk<ram::Parallel>(std::move(loopSeq));

    /* construct exit conditions for odd and even iteration */
    auto addCondition = [](Own<ram::Condition>& cond, Own<ram::Condition> clause) {
        cond = ((cond) ? mk<ram::Conjunction>(std::move(cond), std::move(clause)) : std::move(clause));
    };

    Own<ram::Condition> exitCond;
    VecOwn<ram::Statement> exitStmts;
    for (const ast::Relation* rel : scc) {
        addCondition(exitCond, mk<ram::EmptinessCheck>(getNewRelationName(rel)));
        if (ioType->isLimitSize(rel)) {
            Own<ram::Condition> limit = mk<ram::Constraint>(BinaryConstraintOp::GE,
                    mk<ram::RelationSize>(getConcreteRelationName(rel)),
                    mk<ram::SignedConstant>(ioType->getLimitSize(rel)));
            appendStmt(exitStmts, mk<ram::Exit>(std::move(limit)));
        }
    }

    /* construct fixpoint loop  */
    VecOwn<ram::Statement> res;
    if (preamble.size() > 0) {
        appendStmt(res, mk<ram::Sequence>(std::move(preamble)));
    }
    if (!loop->getStatements().empty() && exitCond && updateTable.size() > 0) {
        appendStmt(res,
                mk<ram::Loop>(mk<ram::Sequence>(std::move(loop), mk<ram::Exit>(std::move(exitCond)),
                        mk<ram::Sequence>(std::move(exitStmts)), mk<ram::Sequence>(std::move(updateTable)))));
    }
    if (postamble.size() > 0) {
        appendStmt(res, mk<ram::Sequence>(std::move(postamble)));
    }
    if (res.size() > 0) {
        return mk<ram::Sequence>(std::move(res));
    }

    fatal("Not Implemented");
}

/** make a subroutine to search for subproofs */
Own<ram::Statement> AstToRamTranslator::makeSubproofSubroutine(const ast::Clause& clause) {
    auto intermediateClause = mk<ast::Clause>(souffle::clone(clause.getHead()));

    // create a clone where all the constraints are moved to the end
    for (auto bodyLit : clause.getBodyLiterals()) {
        // first add all the things that are not constraints
        if (!isA<ast::Constraint>(bodyLit)) {
            intermediateClause->addToBody(souffle::clone(bodyLit));
        }
    }

    // now add all constraints
    for (auto bodyLit : ast::getBodyLiterals<ast::Constraint>(clause)) {
        intermediateClause->addToBody(souffle::clone(bodyLit));
    }

    // name unnamed variables
    nameUnnamedVariables(intermediateClause.get());

    // add constraint for each argument in head of atom
    ast::Atom* head = intermediateClause->getHead();
    size_t auxiliaryArity = auxArityAnalysis->getArity(head);
    auto args = head->getArguments();
    for (size_t i = 0; i < head->getArity() - auxiliaryArity; i++) {
        auto arg = args[i];

        if (auto var = dynamic_cast<ast::Variable*>(arg)) {
            // FIXME: float equiv (`FEQ`)
            auto constraint = mk<ast::BinaryConstraint>(
                    BinaryConstraintOp::EQ, souffle::clone(var), mk<ast::SubroutineArgument>(i));
            constraint->setFinalType(BinaryConstraintOp::EQ);
            intermediateClause->addToBody(std::move(constraint));
        } else if (auto func = dynamic_cast<ast::Functor*>(arg)) {
            TypeAttribute returnType;
            if (auto* inf = dynamic_cast<ast::IntrinsicFunctor*>(func)) {
                assert(inf->getFinalReturnType().has_value() && "functor has missing return type");
                returnType = inf->getFinalReturnType().value();
            } else if (auto* udf = dynamic_cast<ast::UserDefinedFunctor*>(func)) {
                assert(udf->getFinalReturnType().has_value() && "functor has missing return type");
                returnType = udf->getFinalReturnType().value();
            } else {
                assert(false && "unexpected functor type");
            }
            auto opEq = returnType == TypeAttribute::Float ? BinaryConstraintOp::FEQ : BinaryConstraintOp::EQ;
            auto constraint =
                    mk<ast::BinaryConstraint>(opEq, souffle::clone(func), mk<ast::SubroutineArgument>(i));
            constraint->setFinalType(opEq);
            intermediateClause->addToBody(std::move(constraint));
        } else if (auto rec = dynamic_cast<ast::RecordInit*>(arg)) {
            auto constraint = mk<ast::BinaryConstraint>(
                    BinaryConstraintOp::EQ, souffle::clone(rec), mk<ast::SubroutineArgument>(i));
            constraint->setFinalType(BinaryConstraintOp::EQ);
            intermediateClause->addToBody(std::move(constraint));
        }
    }

    // index of level argument in argument list
    size_t levelIndex = head->getArguments().size() - auxiliaryArity;

    // add level constraints, i.e., that each body literal has height less than that of the head atom
    const auto& bodyLiterals = intermediateClause->getBodyLiterals();
    for (auto lit : bodyLiterals) {
        if (auto atom = dynamic_cast<ast::Atom*>(lit)) {
            auto arity = atom->getArity();
            auto atomArgs = atom->getArguments();
            // arity - 1 is the level number in body atoms
            auto constraint = mk<ast::BinaryConstraint>(BinaryConstraintOp::LT,
                    souffle::clone(atomArgs[arity - 1]), mk<ast::SubroutineArgument>(levelIndex));
            constraint->setFinalType(BinaryConstraintOp::LT);
            intermediateClause->addToBody(std::move(constraint));
        }
    }
    return ProvenanceClauseTranslator(*this).translateClause(*intermediateClause, clause);
}

/** make a subroutine to search for subproofs for the non-existence of a tuple */
Own<ram::Statement> AstToRamTranslator::makeNegationSubproofSubroutine(const ast::Clause& clause) {
    // TODO (taipan-snake): Currently we only deal with atoms (no constraints or negations or aggregates
    // or anything else...)
    //
    // The resulting subroutine looks something like this:
    // IF (arg(0), arg(1), _, _) IN rel_1:
    //   return 1
    // IF (arg(0), arg(1), _ ,_) NOT IN rel_1:
    //   return 0
    // ...

    // clone clause for mutation, rearranging constraints to be at the end
    auto clauseReplacedAggregates = mk<ast::Clause>(souffle::clone(clause.getHead()));

    // create a clone where all the constraints are moved to the end
    for (auto bodyLit : clause.getBodyLiterals()) {
        // first add all the things that are not constraints
        if (!isA<ast::Constraint>(bodyLit)) {
            clauseReplacedAggregates->addToBody(souffle::clone(bodyLit));
        }
    }

    // now add all constraints
    for (auto bodyLit : ast::getBodyLiterals<ast::Constraint>(clause)) {
        clauseReplacedAggregates->addToBody(souffle::clone(bodyLit));
    }

    int aggNumber = 0;
    struct AggregatesToVariables : public ast::NodeMapper {
        int& aggNumber;

        AggregatesToVariables(int& aggNumber) : aggNumber(aggNumber) {}

        Own<ast::Node> operator()(Own<ast::Node> node) const override {
            if (dynamic_cast<ast::Aggregator*>(node.get()) != nullptr) {
                return mk<ast::Variable>("agg_" + std::to_string(aggNumber++));
            }

            node->apply(*this);
            return node;
        }
    };

    AggregatesToVariables aggToVar(aggNumber);
    clauseReplacedAggregates->apply(aggToVar);

    // build a vector of unique variables
    std::vector<const ast::Variable*> uniqueVariables;

    visitDepthFirst(*clauseReplacedAggregates, [&](const ast::Variable& var) {
        if (var.getName().find("@level_num") == std::string::npos) {
            // use find_if since uniqueVariables stores pointers, and we need to dereference the pointer to
            // check equality
            if (std::find_if(uniqueVariables.begin(), uniqueVariables.end(),
                        [&](const ast::Variable* v) { return *v == var; }) == uniqueVariables.end()) {
                uniqueVariables.push_back(&var);
            }
        }
    });

    // a mapper to replace variables with subroutine arguments
    struct VariablesToArguments : public ast::NodeMapper {
        const std::vector<const ast::Variable*>& uniqueVariables;

        VariablesToArguments(const std::vector<const ast::Variable*>& uniqueVariables)
                : uniqueVariables(uniqueVariables) {}

        Own<ast::Node> operator()(Own<ast::Node> node) const override {
            // replace unknown variables
            if (auto varPtr = dynamic_cast<const ast::Variable*>(node.get())) {
                if (varPtr->getName().find("@level_num") == std::string::npos) {
                    size_t argNum = std::find_if(uniqueVariables.begin(), uniqueVariables.end(),
                                            [&](const ast::Variable* v) { return *v == *varPtr; }) -
                                    uniqueVariables.begin();

                    return mk<ast::SubroutineArgument>(argNum);
                } else {
                    return mk<ast::UnnamedVariable>();
                }
            }

            // apply recursive
            node->apply(*this);

            // otherwise nothing
            return node;
        }
    };

    // the structure of this subroutine is a sequence where each nested statement is a search in each
    // relation
    VecOwn<ram::Statement> searchSequence;

    // make a copy so that when we mutate clause, pointers to objects in newClause are not affected
    auto newClause = souffle::clone(clauseReplacedAggregates);

    // go through each body atom and create a return
    size_t litNumber = 0;
    for (const auto& lit : newClause->getBodyLiterals()) {
        if (auto atom = dynamic_cast<ast::Atom*>(lit)) {
            size_t auxiliaryArity = auxArityAnalysis->getArity(atom);

            auto relName = getConcreteRelationName(atom);
            // construct a query
            VecOwn<ram::Expression> query;

            // translate variables to subroutine arguments
            VariablesToArguments varsToArgs(uniqueVariables);
            atom->apply(varsToArgs);

            auto atomArgs = atom->getArguments();
            // add each value (subroutine argument) to the search query
            for (size_t i = 0; i < atom->getArity() - auxiliaryArity; i++) {
                auto arg = atomArgs[i];
                query.push_back(translateValue(arg, ValueIndex()));
            }

            // fill up query with nullptrs for the provenance columns
            for (size_t i = 0; i < auxiliaryArity; i++) {
                query.push_back(mk<ram::UndefValue>());
            }

            // ensure the length of query tuple is correct
            assert(query.size() == atom->getArity() && "wrong query tuple size");

            // create existence checks to check if the tuple exists or not
            auto existenceCheck = mk<ram::ExistenceCheck>(relName, std::move(query));
            auto negativeExistenceCheck = mk<ram::Negation>(souffle::clone(existenceCheck));

            // return true if the tuple exists
            VecOwn<ram::Expression> returnTrue;
            returnTrue.push_back(mk<ram::SignedConstant>(1));

            // return false if the tuple exists
            VecOwn<ram::Expression> returnFalse;
            returnFalse.push_back(mk<ram::SignedConstant>(0));

            // create a ram::Query to return true/false
            appendStmt(searchSequence, mk<ram::Query>(mk<ram::Filter>(std::move(existenceCheck),
                                               mk<ram::SubroutineReturn>(std::move(returnTrue)))));
            appendStmt(searchSequence, mk<ram::Query>(mk<ram::Filter>(std::move(negativeExistenceCheck),
                                               mk<ram::SubroutineReturn>(std::move(returnFalse)))));
        } else if (auto neg = dynamic_cast<ast::Negation*>(lit)) {
            auto atom = neg->getAtom();

            size_t auxiliaryArity = auxArityAnalysis->getArity(atom);
            auto relName = getConcreteRelationName(atom);
            // construct a query
            VecOwn<ram::Expression> query;

            // translate variables to subroutine arguments
            VariablesToArguments varsToArgs(uniqueVariables);
            atom->apply(varsToArgs);

            auto atomArgs = atom->getArguments();
            // add each value (subroutine argument) to the search query
            for (size_t i = 0; i < atom->getArity() - auxiliaryArity; i++) {
                auto arg = atomArgs[i];
                query.push_back(translateValue(arg, ValueIndex()));
            }

            // fill up query with nullptrs for the provenance columns
            for (size_t i = 0; i < auxiliaryArity; i++) {
                query.push_back(mk<ram::UndefValue>());
            }

            // ensure the length of query tuple is correct
            assert(query.size() == atom->getArity() && "wrong query tuple size");

            // create existence checks to check if the tuple exists or not
            auto existenceCheck = mk<ram::ExistenceCheck>(relName, std::move(query));
            auto negativeExistenceCheck = mk<ram::Negation>(souffle::clone(existenceCheck));

            // return true if the tuple exists
            VecOwn<ram::Expression> returnTrue;
            returnTrue.push_back(mk<ram::SignedConstant>(1));

            // return false if the tuple exists
            VecOwn<ram::Expression> returnFalse;
            returnFalse.push_back(mk<ram::SignedConstant>(0));

            // create a ram::Query to return true/false
            appendStmt(searchSequence, mk<ram::Query>(mk<ram::Filter>(std::move(existenceCheck),
                                               mk<ram::SubroutineReturn>(std::move(returnFalse)))));
            appendStmt(searchSequence, mk<ram::Query>(mk<ram::Filter>(std::move(negativeExistenceCheck),
                                               mk<ram::SubroutineReturn>(std::move(returnTrue)))));

        } else if (auto con = dynamic_cast<ast::Constraint*>(lit)) {
            VariablesToArguments varsToArgs(uniqueVariables);
            con->apply(varsToArgs);

            // translate to a ram::Condition
            auto condition = translateConstraint(con, ValueIndex());
            auto negativeCondition = mk<ram::Negation>(souffle::clone(condition));

            // create a return true value
            VecOwn<ram::Expression> returnTrue;
            returnTrue.push_back(mk<ram::SignedConstant>(1));

            // create a return false value
            VecOwn<ram::Expression> returnFalse;
            returnFalse.push_back(mk<ram::SignedConstant>(0));

            appendStmt(searchSequence, mk<ram::Query>(mk<ram::Filter>(std::move(condition),
                                               mk<ram::SubroutineReturn>(std::move(returnTrue)))));
            appendStmt(searchSequence, mk<ram::Query>(mk<ram::Filter>(std::move(negativeCondition),
                                               mk<ram::SubroutineReturn>(std::move(returnFalse)))));
        }

        litNumber++;
    }

    return mk<ram::Sequence>(std::move(searchSequence));
}

bool AstToRamTranslator::removeADTs(const ast::TranslationUnit& translationUnit) {
    struct ADTsFuneral : public ast::NodeMapper {
        mutable bool changed{false};
        const ast::analysis::SumTypeBranchesAnalysis& sumTypesBranches;

        ADTsFuneral(const ast::TranslationUnit& tu)
                : sumTypesBranches(*tu.getAnalysis<ast::analysis::SumTypeBranchesAnalysis>()) {}

        Own<ast::Node> operator()(Own<ast::Node> node) const override {
            // Rewrite sub-expressions first
            node->apply(*this);

            if (!isA<ast::BranchInit>(node)) {
                return node;
            }

            changed = true;
            auto& adt = *as<ast::BranchInit>(node);
            auto& type = sumTypesBranches.unsafeGetType(adt.getConstructor());
            auto& branches = type.getBranches();

            // Find branch ID.
            ast::analysis::AlgebraicDataType::Branch searchDummy{adt.getConstructor(), {}};
            auto iterToBranch = std::lower_bound(branches.begin(), branches.end(), searchDummy,
                    [](const ast::analysis::AlgebraicDataType::Branch& left,
                            const ast::analysis::AlgebraicDataType::Branch& right) {
                        return left.name < right.name;
                    });

            // Branch id corresponds to the position in lexicographical ordering.
            auto branchID = std::distance(std::begin(branches), iterToBranch);

            if (isADTEnum(type)) {
                auto branchTag = mk<ast::NumericConstant>(branchID);
                branchTag->setFinalType(ast::NumericConstant::Type::Int);
                return branchTag;
            } else {
                // Collect branch arguments
                VecOwn<ast::Argument> branchArguments;
                for (auto* arg : adt.getArguments()) {
                    branchArguments.emplace_back(arg->clone());
                }

                // Branch is stored either as [branch_id, [arguments]]
                // or [branch_id, argument] in case of a single argument.
                auto branchArgs = [&]() -> Own<ast::Argument> {
                    if (branchArguments.size() != 1) {
                        return mk<ast::Argument, ast::RecordInit>(std::move(branchArguments));
                    } else {
                        return std::move(branchArguments.at(0));
                    }
                }();

                // Arguments for the resulting record [branch_id, branch_args].
                VecOwn<ast::Argument> finalRecordArgs;

                auto branchTag = mk<ast::NumericConstant>(branchID);
                branchTag->setFinalType(ast::NumericConstant::Type::Int);
                finalRecordArgs.push_back(std::move(branchTag));
                finalRecordArgs.push_back(std::move(branchArgs));

                return mk<ast::RecordInit>(std::move(finalRecordArgs), adt.getSrcLoc());
            }
        }
    };

    ADTsFuneral mapper(translationUnit);
    translationUnit.getProgram().apply(mapper);
    return mapper.changed;
}

void AstToRamTranslator::makeRamLoad(VecOwn<ram::Statement>& curStmts, const ast::Relation* relation) {
    for (auto directives : getInputDirectives(relation)) {
        Own<ram::Statement> statement = mk<ram::IO>(getConcreteRelationName(relation), directives);
        if (Global::config().has("profile")) {
            const std::string logTimerStatement = LogStatement::tRelationLoadTime(
                    toString(relation->getQualifiedName()), relation->getSrcLoc());
            statement = mk<ram::LogRelationTimer>(
                    std::move(statement), logTimerStatement, getConcreteRelationName(relation));
        }
        appendStmt(curStmts, std::move(statement));
    }
}

void AstToRamTranslator::makeRamStore(VecOwn<ram::Statement>& curStmts, const ast::Relation* relation) {
    for (auto directives : getOutputDirectives(relation)) {
        Own<ram::Statement> statement = mk<ram::IO>(getConcreteRelationName(relation), directives);
        if (Global::config().has("profile")) {
            const std::string logTimerStatement = LogStatement::tRelationSaveTime(
                    toString(relation->getQualifiedName()), relation->getSrcLoc());
            statement = mk<ram::LogRelationTimer>(
                    std::move(statement), logTimerStatement, getConcreteRelationName(relation));
        }
        appendStmt(curStmts, std::move(statement));
    }
}

void AstToRamTranslator::makeRamClear(VecOwn<ram::Statement>& curStmts, const ast::Relation* relation) {
    appendStmt(curStmts, mk<ram::Clear>(getConcreteRelationName(relation)));
}

/** translates the given datalog program into an equivalent RAM program  */
void AstToRamTranslator::translateProgram(const ast::TranslationUnit& translationUnit) {
    // keep track of relevant analyses
    ioType = translationUnit.getAnalysis<ast::analysis::IOTypeAnalysis>();
    typeEnv = &translationUnit.getAnalysis<ast::analysis::TypeEnvironmentAnalysis>()->getTypeEnvironment();
    const auto& sccOrder = *translationUnit.getAnalysis<ast::analysis::TopologicallySortedSCCGraphAnalysis>();
    relationSchedule = translationUnit.getAnalysis<ast::analysis::RelationScheduleAnalysis>();
    sccGraph = translationUnit.getAnalysis<ast::analysis::SCCGraphAnalysis>();
    recursiveClauses = translationUnit.getAnalysis<ast::analysis::RecursiveClausesAnalysis>();
    auxArityAnalysis = translationUnit.getAnalysis<ast::analysis::AuxiliaryArityAnalysis>();
    functorAnalysis = translationUnit.getAnalysis<ast::analysis::FunctorAnalysis>();
    relDetail = translationUnit.getAnalysis<ast::analysis::RelationDetailCacheAnalysis>();
    polyAnalysis = translationUnit.getAnalysis<ast::analysis::PolymorphicObjectsAnalysis>();

    // set up the final fixed types
    // TODO (azreika): should be removed once the translator is refactored to avoid cloning
    visitDepthFirst(*program, [&](const ast::NumericConstant& nc) {
        const_cast<ast::NumericConstant&>(nc).setFinalType(polyAnalysis->getInferredType(&nc));
    });
    visitDepthFirst(*program, [&](const ast::Aggregator& aggr) {
        const_cast<ast::Aggregator&>(aggr).setFinalType(polyAnalysis->getOverloadedOperator(&aggr));
    });
    visitDepthFirst(*program, [&](const ast::BinaryConstraint& bc) {
        const_cast<ast::BinaryConstraint&>(bc).setFinalType(polyAnalysis->getOverloadedOperator(&bc));
    });
    visitDepthFirst(*program, [&](const ast::IntrinsicFunctor& inf) {
        const_cast<ast::IntrinsicFunctor&>(inf).setFinalOpType(polyAnalysis->getOverloadedFunctionOp(&inf));
        const_cast<ast::IntrinsicFunctor&>(inf).setFinalReturnType(functorAnalysis->getReturnType(&inf));
    });
    visitDepthFirst(*program, [&](const ast::UserDefinedFunctor& udf) {
        const_cast<ast::UserDefinedFunctor&>(udf).setFinalReturnType(functorAnalysis->getReturnType(&udf));
    });

    // determine the sips to use
    std::string sipsChosen = "all-bound";
    if (Global::config().has("RamSIPS")) {
        sipsChosen = Global::config().get("RamSIPS");
    }
    sips = ast::SipsMetric::create(sipsChosen, translationUnit);

    // replace ADTs with record representatives
    removeADTs(translationUnit);

    // handle the case of an empty SCC graph
    if (sccGraph->getNumberOfSCCs() == 0) return;

    // create all Ram relations in ramRels
    for (const auto& scc : sccOrder.order()) {
        const auto& isRecursive = sccGraph->isRecursive(scc);
        const auto& allInterns = sccGraph->getInternalRelations(scc);
        for (const auto& rel : allInterns) {
            std::string name = rel->getQualifiedName().toString();
            auto arity = rel->getArity();
            auto auxiliaryArity = auxArityAnalysis->getArity(rel);
            auto representation = rel->getRepresentation();
            const auto& attributes = rel->getAttributes();

            std::vector<std::string> attributeNames;
            std::vector<std::string> attributeTypeQualifiers;
            for (size_t i = 0; i < rel->getArity(); ++i) {
                attributeNames.push_back(attributes[i]->getName());
                if (typeEnv != nullptr) {
                    attributeTypeQualifiers.push_back(
                            getTypeQualifier(typeEnv->getType(attributes[i]->getTypeName())));
                }
            }
            ramRels[name] = mk<ram::Relation>(
                    name, arity, auxiliaryArity, attributeNames, attributeTypeQualifiers, representation);

            // recursive relations also require @delta and @new variants, with the same signature
            if (isRecursive) {
                std::string deltaName = "@delta_" + name;
                std::string newName = "@new_" + name;
                ramRels[deltaName] = mk<ram::Relation>(deltaName, arity, auxiliaryArity, attributeNames,
                        attributeTypeQualifiers, representation);
                ramRels[newName] = mk<ram::Relation>(newName, arity, auxiliaryArity, attributeNames,
                        attributeTypeQualifiers, representation);
            }
        }
    }

    // maintain the index of the SCC within the topological order
    size_t indexOfScc = 0;

    // iterate over each SCC according to the topological order
    for (const auto& scc : sccOrder.order()) {
        // create subroutine for this stratum
        ramSubs["stratum_" + std::to_string(indexOfScc)] = translateSCC(scc, indexOfScc);
        indexOfScc++;
    }

    // invoke all strata
    VecOwn<ram::Statement> res;
    for (size_t i = 0; i < indexOfScc; i++) {
        appendStmt(res, mk<ram::Call>("stratum_" + std::to_string(i)));
    }

    // add main timer if profiling
    if (res.size() > 0 && Global::config().has("profile")) {
        auto newStmt = mk<ram::LogTimer>(mk<ram::Sequence>(std::move(res)), LogStatement::runtime());
        res.clear();
        appendStmt(res, std::move(newStmt));
    }

    // done for main prog
    ramMain = mk<ram::Sequence>(std::move(res));

    // add subroutines for each clause
    if (Global::config().has("provenance")) {
        visitDepthFirst(*program, [&](const ast::Clause& clause) {
            std::stringstream relName;
            relName << clause.getHead()->getQualifiedName();

            // do not add subroutines for info relations or facts
            if (relName.str().find("@info") != std::string::npos || clause.getBodyLiterals().empty()) {
                return;
            }

            std::string subroutineLabel =
                    relName.str() + "_" + std::to_string(getClauseNum(program, &clause)) + "_subproof";
            ramSubs[subroutineLabel] = makeSubproofSubroutine(clause);

            std::string negationSubroutineLabel = relName.str() + "_" +
                                                  std::to_string(getClauseNum(program, &clause)) +
                                                  "_negation_subproof";
            ramSubs[negationSubroutineLabel] = makeNegationSubproofSubroutine(clause);
        });
    }
}

Own<ram::TranslationUnit> AstToRamTranslator::translateUnit(ast::TranslationUnit& tu) {
    auto ram_start = std::chrono::high_resolution_clock::now();
    program = &tu.getProgram();

    translateProgram(tu);

    SymbolTable& symTab = getSymbolTable();
    ErrorReport& errReport = tu.getErrorReport();
    DebugReport& debugReport = tu.getDebugReport();
    VecOwn<ram::Relation> rels;
    for (auto& cur : ramRels) {
        rels.push_back(std::move(cur.second));
    }

    if (ramMain == nullptr) {
        ramMain = mk<ram::Sequence>();
    }
    auto ramProg = mk<ram::Program>(std::move(rels), std::move(ramMain), std::move(ramSubs));

    // add the translated program to the debug report
    if (!Global::config().get("debug-report").empty()) {
        auto ram_end = std::chrono::high_resolution_clock::now();
        std::string runtimeStr =
                "(" + std::to_string(std::chrono::duration<double>(ram_end - ram_start).count()) + "s)";
        std::stringstream ramProgStr;
        ramProgStr << *ramProg;
        debugReport.addSection("ram-program", "RAM Program " + runtimeStr, ramProgStr.str());
    }

    return mk<ram::TranslationUnit>(std::move(ramProg), std::move(symTab), errReport, debugReport);
}

}  // namespace souffle::ast2ram
