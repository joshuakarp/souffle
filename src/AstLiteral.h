/*
 * Souffle - A Datalog Compiler
 * Copyright (c) 2013, Oracle and/or its affiliates. All rights reserved
 * Licensed under the Universal Permissive License v 1.0 as shown at:
 * - https://opensource.org/licenses/UPL
 * - <souffle root>/licenses/SOUFFLE-UPL.txt
 */

/************************************************************************
 *
 * @file AstLiteral.h
 *
 * Define classes for Literals and its subclasses atoms, negated atoms,
 * and binary relations.
 *
 ***********************************************************************/

#pragma once

#include "AstAbstract.h"
#include "AstArgument.h"
#include "AstNode.h"
#include "AstQualifiedName.h"
#include "BinaryConstraintOps.h"
#include "Util.h"

#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "AstAbstract.h"
#include "AstNode.h"

namespace souffle {

class AstRelation;
class AstClause;
class AstProgram;
class AstAtom;

/**
 * Subclass of Literal that represents the use of a relation
 * either in the head or in the body of a Clause, e.g., parent(x,y).
 * The arguments of the atom can be variables or constants.
 */
class AstAtom : public AstLiteral {
public:
    AstAtom(AstQualifiedName name = AstQualifiedName()) : name(std::move(name)) {}

    AstAtom(AstQualifiedName name, std::vector<std::unique_ptr<AstArgument>> args, SrcLocation srcLoc)
            : name(std::move(name)), arguments(std::move(args)) {
        setSrcLoc(srcLoc);
    }

    /** get qualified name */
    const AstQualifiedName& getQualifiedName() const {
        return name;
    }

    /** get arity of the atom */
    size_t getArity() const {
        return arguments.size();
    }

    /** set qualified name */
    void setQualifiedName(const AstQualifiedName& n) {
        name = n;
    }

    /** add argument to the atom */
    void addArgument(std::unique_ptr<AstArgument> arg) {
        arguments.push_back(std::move(arg));
    }

    /** get arguments */
    std::vector<AstArgument*> getArguments() const {
        return toPtrVector(arguments);
    }

    AstAtom* clone() const override {
        auto res = new AstAtom(name);
        res->setSrcLoc(getSrcLoc());
        for (const auto& cur : arguments) {
            res->arguments.emplace_back(cur->clone());
        }
        return res;
    }

    void apply(const AstNodeMapper& map) override {
        for (auto& arg : arguments) {
            arg = map(std::move(arg));
        }
    }

    std::vector<const AstNode*> getChildNodes() const override {
        std::vector<const AstNode*> res;
        for (auto& cur : arguments) {
            res.push_back(cur.get());
        }
        return res;
    }

protected:
    void print(std::ostream& os) const override {
        os << getQualifiedName() << "(";
        os << join(arguments, ",", print_deref<std::unique_ptr<AstArgument>>());
        os << ")";
    }

    bool equal(const AstNode& node) const override {
        const auto& other = static_cast<const AstAtom&>(node);
        return name == other.name && equal_targets(arguments, other.arguments);
    }

    /** name */
    AstQualifiedName name;

    /** arguments */
    std::vector<std::unique_ptr<AstArgument>> arguments;
};

/**
 * Subclass of Literal that represents a negated atom, * e.g., !parent(x,y).
 * A Negated atom occurs in a body of clause and cannot occur in a head of a clause.
 */
class AstNegation : public AstLiteral {
public:
    AstNegation(std::unique_ptr<AstAtom> atom) : atom(std::move(atom)) {}

    /** get negated atom */
    AstAtom* getAtom() const {
        return atom.get();
    }

    AstNegation* clone() const override {
        auto* res = new AstNegation(std::unique_ptr<AstAtom>(atom->clone()));
        res->setSrcLoc(getSrcLoc());
        return res;
    }

    void apply(const AstNodeMapper& map) override {
        atom = map(std::move(atom));
    }

    std::vector<const AstNode*> getChildNodes() const override {
        return {atom.get()};
    }

protected:
    void print(std::ostream& os) const override {
        os << "!" << *atom;
    }

    bool equal(const AstNode& node) const override {
        assert(nullptr != dynamic_cast<const AstNegation*>(&node));
        const auto& other = static_cast<const AstNegation&>(node);
        return equal_ptr(atom, other.atom);
    }

    /** negated atom */
    std::unique_ptr<AstAtom> atom;
};

/**
 * Subclass of Literal that represents a negated atom, * e.g., !parent(x,y).
 * A Negated atom occurs in a body of clause and cannot occur in a head of a clause.
 *
 * Specialised for provenance: used for existence check that tuple doesn't already exist
 */
class AstProvenanceNegation : public AstNegation {
public:
    AstProvenanceNegation(std::unique_ptr<AstAtom> atom) : AstNegation(std::move(atom)) {}

    AstProvenanceNegation* clone() const override {
        auto* res = new AstProvenanceNegation(std::unique_ptr<AstAtom>(atom->clone()));
        res->setSrcLoc(getSrcLoc());
        return res;
    }

protected:
    void print(std::ostream& os) const override {
        os << "prov!" << *atom;
    }
};

/**
 * Boolean Constraint
 *
 * Representing either 'true' or 'false' values
 */
class AstBooleanConstraint : public AstConstraint {
public:
    AstBooleanConstraint(bool truthValue) : truthValue(truthValue) {}

    /** check whether constraint holds */
    bool isTrue() const {
        return truthValue;
    }

    /** set truth value */
    void set(bool value) {
        truthValue = value;
    }

    AstBooleanConstraint* clone() const override {
        auto* res = new AstBooleanConstraint(truthValue);
        res->setSrcLoc(getSrcLoc());
        return res;
    }

protected:
    void print(std::ostream& os) const override {
        os << (truthValue ? "true" : "false");
    }

    bool equal(const AstNode& node) const override {
        assert(nullptr != dynamic_cast<const AstBooleanConstraint*>(&node));
        const auto& other = static_cast<const AstBooleanConstraint&>(node);
        return truthValue == other.truthValue;
    }

    /** truth value */
    bool truthValue;
};

/**
 * Functional Constraint
 *
 * Representing a functional dependency (choice construct)
 * eg. x -> y
 * x uniquely identifies some y value
 * LHS of dependency stored as a vector of pairs: 
 *      - first: pointer to LHS of dependency variable
 *      - second: index position of this LHS variable in relation arguments (0-indexed)
 * Stored in a vector to support n-ary dependencies.
 * eg. (x,y) -> z
 */
class AstFunctionalConstraint : public AstConstraint {
public:
    AstFunctionalConstraint(
            std::vector<std::unique_ptr<AstVariable>> ls, std::unique_ptr<AstVariable> rs)
            : lhs(std::move(ls)), positions(lhs.size(), 0), rhs(std::move(rs)) {}

    // Unsure if the "std::make_pair" is allowed in an initialiser list
    AstFunctionalConstraint(
            std::unique_ptr<AstVariable> ls, std::unique_ptr<AstVariable> rs)
            : positions(1, 0), rhs(std::move(rs)) 
            {
                lhs.push_back(std::move(ls));
            }

    /** get left-hand side of functional constraint */ 
    const AstVariable *getLHS(size_t lhsNum) const {
        return lhs.at(lhsNum).get();
    }

    /** get arity of the functional constraint (i.e. number of source nodes: (x,y)->z has an arity of 2) */
    const size_t getArity() const {
        return lhs.size();
    }

    /** get left-hand side of functional constraint */ 
    const AstVariable *getRHS() const {
        return rhs.get(); 
    }

    std::vector<const AstNode*> getChildNodes() const override {
        std::vector<const AstNode*> res;
        for (auto& cur : lhs) {
            res.push_back(cur.get());
        }
        res.push_back(rhs.get());
        return res;
    }

    AstFunctionalConstraint* clone() const override {
        std::vector<std::unique_ptr<AstVariable>> newLHS;
        // TODO: should include positions to deep clone
        for (size_t i = 0; i < lhs.size(); i++) {
            newLHS.push_back(std::unique_ptr<AstVariable>(lhs.at(i)->clone()));
        }
        auto* res = new AstFunctionalConstraint(
            std::move(newLHS), // parsing newLHS on its own giving "use of deleted function" error
            std::unique_ptr<AstVariable>(rhs->clone()));
        res->setSrcLoc(getSrcLoc());
        return res;
    }

    /* get index position of LHS (source) in relation arguments (0-indexed) */ 
    // int getPosition() const {
    //     return position;
    // }
    size_t getPosition(size_t lhsNum) const {
        return positions.at(lhsNum);
    }

    /* set index position of LHS (source) in relation arguments (0-indexed) */
    // void setPosition (int pos) {
    //     position = pos;
    // }
    /** set index position of LHS (source) in relation arguments (0-indexed) 
     * lhsNum = the LHS source currently wanting to be set
     * i.e. (x,y)->z: to set y, lhsNum = 1
     * pos: the attribute's index position in relation (0-indexed)
    */
    void setPosition (size_t lhsNum, size_t pos) {
        std::cout << "lhsNum (j): " << lhsNum << " pos (i): " << pos << "\n";
        positions.at(lhsNum) = pos;
    }

protected:
    void print(std::ostream& os) const override {
        // if (lhs.size() > 1) {
        //     os << "(";
        // }
        // os << join(lhs, ",", print_deref<std::unique_ptr<AstVariable>>());
        // if (lhs.size() > 1) {
        //     os << ")";
        // }
        // os << "->" << *rhs;
    }

    bool equal(const AstNode& node) const override {
        assert(nullptr != dynamic_cast<const AstFunctionalConstraint*>(&node));
        const auto& other = static_cast<const AstFunctionalConstraint&>(node);
        if (lhs.size() != other.lhs.size()) {
            return false;
        }
        bool lhsMatch = true;
        // Check that each pointer and position match in LHS
        for (size_t i = 0; i < lhs.size(); i++) {
            if (!equal_ptr(lhs.at(i), other.lhs.at(i))) {
                lhsMatch = false;
                break;
            }

            if (positions.at(i) != other.positions.at(i)) {
                lhsMatch = false;
                break;
            }
        }
        return lhsMatch && equal_ptr(rhs,other.rhs);
    }

    /* lhs of functional constraint */ 
    std::vector<std::unique_ptr<AstVariable>> lhs; 

    /** index positions of LHS (source) nodes in relation arguments (0-indexed)
     * eg. A(x,y,z) constrains y->z 
     * Dependency y->z has position 1
     **/
    std::vector<size_t> positions;

    /* rhs of functional constraint */ 
    std::unique_ptr<AstVariable> rhs;
};

/**
 * Subclass of Constraint that represents a binary constraint
 * e.g., x = y.
 */
class AstBinaryConstraint : public AstConstraint {
public:
    AstBinaryConstraint(
            BinaryConstraintOp o, std::unique_ptr<AstArgument> ls, std::unique_ptr<AstArgument> rs)
            : operation(o), lhs(std::move(ls)), rhs(std::move(rs)) {}

    /** get LHS argument */
    AstArgument* getLHS() const {
        return lhs.get();
    }

    /** get RHS argument */
    AstArgument* getRHS() const {
        return rhs.get();
    }

    /** get binary operator */
    BinaryConstraintOp getOperator() const {
        return operation;
    }

    /** set binary operator */
    void setOperator(BinaryConstraintOp op) {
        operation = op;
    }

    AstBinaryConstraint* clone() const override {
        auto* res = new AstBinaryConstraint(operation, std::unique_ptr<AstArgument>(lhs->clone()),
                std::unique_ptr<AstArgument>(rhs->clone()));
        res->setSrcLoc(getSrcLoc());
        return res;
    }

    void apply(const AstNodeMapper& map) override {
        lhs = map(std::move(lhs));
        rhs = map(std::move(rhs));
    }

    std::vector<const AstNode*> getChildNodes() const override {
        return {lhs.get(), rhs.get()};
    }

protected:
    void print(std::ostream& os) const override {
        os << *lhs;
        os << " " << toBinaryConstraintSymbol(operation) << " ";
        os << *rhs;
    }

    bool equal(const AstNode& node) const override {
        assert(nullptr != dynamic_cast<const AstBinaryConstraint*>(&node));
        const auto& other = static_cast<const AstBinaryConstraint&>(node);
        return operation == other.operation && equal_ptr(lhs, other.lhs) && equal_ptr(rhs, other.rhs);
    }

    /** constraint operator */
    BinaryConstraintOp operation;

    /** left-hand side of binary constraint */
    std::unique_ptr<AstArgument> lhs;

    /** right-hand side of binary constraint */
    std::unique_ptr<AstArgument> rhs;
};

}  // end of namespace souffle
