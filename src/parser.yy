/*
 * Souffle - A Datalog Compiler
 * Copyright (c) 2013, 2015, Oracle and/or its affiliates. All rights reserved
 * Licensed under the Universal Permissive License v 1.0 as shown at:
 * - https://opensource.org/licenses/UPL
 * - <souffle root>/licenses/SOUFFLE-UPL.txt
 */

/************************************************************************
 *
 * @file parser.yy
 *
 * @brief Parser for Datalog
 *
 ***********************************************************************/
%skeleton "lalr1.cc"
%require "3.0.2"

%defines
%define parser_class_name {parser}
%define api.token.constructor
%define api.value.type variant
%define parse.assert
%define api.location.type {SrcLocation}

%locations

%define parse.trace
%define parse.error verbose

/* -- Dependencies -- */
%code requires {
    #include "AggregateOp.h"
    #include "AstArgument.h"
    #include "AstClause.h"
    #include "AstComponent.h"
    #include "AstFunctorDeclaration.h"
    #include "AstIO.h"
    #include "AstNode.h"
    #include "AstParserUtils.h"
    #include "AstPragma.h"
    #include "AstProgram.h"
    #include "BinaryConstraintOps.h"
    #include "FunctorOps.h"
    #include "RamTypes.h"

    using namespace souffle;

    namespace souffle {
        class ParserDriver;
    }

    using yyscan_t = void*;

    #define YY_NULLPTR nullptr

    /* Macro to update locations as parsing proceeds */
    #define YYLLOC_DEFAULT(Cur, Rhs, N)                         \
    do {                                                        \
        if (N) {                                                \
            (Cur).start         = YYRHSLOC(Rhs, 1).start;       \
            (Cur).end           = YYRHSLOC(Rhs, N).end;         \
            (Cur).filenames     = YYRHSLOC(Rhs, N).filenames;   \
        } else {                                                \
            (Cur).start         = YYRHSLOC(Rhs, 0).end;         \
            (Cur).end           = YYRHSLOC(Rhs, 0).end;         \
            (Cur).filenames     = YYRHSLOC(Rhs, 0).filenames;   \
        }                                                       \
    } while (0)
}

%code {
    #include "ParserDriver.h"
}

%param { ParserDriver &driver }
%param { yyscan_t yyscanner }

/* -- Tokens -- */
%token END 0                     "end of file"
%token <std::string> STRING      "symbol"
%token <std::string> IDENT       "identifier"
%token <std::string> NUMBER      "number"
%token <std::string> FLOAT       "float"
%token PRAGMA                    "pragma directive"
%token OUTPUT_QUALIFIER          "relation qualifier output"
%token INPUT_QUALIFIER           "relation qualifier input"
%token PRINTSIZE_QUALIFIER       "relation qualifier printsize"
%token BRIE_QUALIFIER            "BRIE datastructure qualifier"
%token BTREE_QUALIFIER           "BTREE datastructure qualifier"
%token EQREL_QUALIFIER           "equivalence relation qualifier"
%token OVERRIDABLE_QUALIFIER     "relation qualifier overidable"
%token INLINE_QUALIFIER          "relation qualifier inline"
%token TMATCH                    "match predicate"
%token TCONTAINS                 "checks whether substring is contained in a string"
%token CAT                       "concatenation of strings"
%token ORD                       "ordinal number of a string"
%token STRLEN                    "length of a string"
%token SUBSTR                    "sub-string of a string"
%token MEAN                      "mean aggregator"
%token MIN                       "min aggregator"
%token MAX                       "max aggregator"
%token COUNT                     "count aggregator"
%token SUM                       "sum aggregator"
%token TRUE                      "true literal constraint"
%token FALSE                     "false literal constraint"
%token PLAN                      "plan keyword"
%token IF                        ":-"
%token DECL                      "relation declaration"
%token FUNCTOR                   "functor declaration"
%token INPUT_DECL                "input directives declaration"
%token OUTPUT_DECL               "output directives declaration"
%token PRINTSIZE_DECL            "printsize directives declaration"
%token OVERRIDE                  "override rules of super-component"
%token TYPE                      "type declaration"
%token COMPONENT                 "component declaration"
%token INSTANTIATE               "component instantiation"
%token NUMBER_TYPE               "numeric type declaration"
%token SYMBOL_TYPE               "symbolic type declaration"
%token TONUMBER                  "convert string to (signed) number"
%token TOSTRING                  "convert number to string"
%token ITOU                      "convert int to unsigned"
%token ITOF                      "convert int to float"
%token UTOI                      "convert unsigned to int"
%token UTOF                      "convert unsigned to float"
%token FTOI                      "convert float to int"
%token FTOU                      "convert float to unsigned"
%token AS                        "type cast"
%token AT                        "@"
%token NIL                       "nil reference"
%token PIPE                      "|"
%token LBRACKET                  "["
%token RBRACKET                  "]"
%token UNDERSCORE                "_"
%token DOLLAR                    "$"
%token PLUS                      "+"
%token MINUS                     "-"
%token EXCLAMATION               "!"
%token LPAREN                    "("
%token RPAREN                    ")"
%token COMMA                     ","
%token COLON                     ":"
%token DOUBLECOLON               "::"
%token SEMICOLON                 ";"
%token DOT                       "."
%token EQUALS                    "="
%token STAR                      "*"
%token SLASH                     "/"
%token CARET                     "^"
%token PERCENT                   "%"
%token LBRACE                    "{"
%token RBRACE                    "}"
%token SUBTYPE                   "<:"
%token LT                        "<"
%token GT                        ">"
%token RIGHTARROW                "->"
%token LE                        "<="
%token GE                        ">="
%token NE                        "!="
%token BW_AND                    "band"
%token BW_OR                     "bor"
%token BW_XOR                    "bxor"
%token BW_SHIFT_L                "bshl"
%token BW_SHIFT_R                "bshr"
%token BW_SHIFT_R_UNSIGNED       "bshru"
%token BW_NOT                    "bnot"
%token L_AND                     "land"
%token L_OR                      "lor"
%token L_NOT                     "lnot"

/* -- Non-Terminal Types -- */
%type <AstAtom *>                           atom
%type <AstArgument *>                       arg
%type <RuleBody *>                          body
%type <AstComponentType *>                  comp_type
%type <AstComponentInit *>                  comp_init
%type <AstComponent *>                      component
%type <AstComponent *>                      component_body
%type <AstComponent *>                      component_head
%type <RuleBody *>                          conjunction
%type <AstConstraint *>                     constraint
%type <RuleBody *>                          disjunction
%type <AstExecutionOrder *>                 exec_order_list
%type <AstExecutionPlan *>                  exec_plan
%type <AstExecutionPlan *>                  exec_plan_list
%type <AstClause *>                         fact
%type <AstFunctorDeclaration *>             functor_decl
%type <TypeAttribute>                       predefined_type
%type <std::vector<AstAtom *>>              head
%type <std::vector<std::string>>            identifier
%type <std::vector<AstIO *>>                io_directive_list
%type <std::vector<AstIO *>>                io_relation_list
%type <std::string>                         kvp_value
%type <std::vector<AstIO *>>                io_head
%type <std::vector<AstArgument *>>          non_empty_arg_list
%type <std::vector<AstAttribute *>>         non_empty_attributes
%type <AstExecutionOrder *>                 non_empty_exec_order_list
%type <std::vector<TypeAttribute>>          non_empty_functor_arg_type_list
%type <std::vector<std::pair
            <std::string, std::string>>>    non_empty_key_value_pairs
%type <AstRecordType *>                     non_empty_record_type_list
%type <AstPragma *>                         pragma
%type <std::set<RelationTag>>               relation_tags
%type <std::vector<AstRelation *>>          relation_decl
%type <std::vector<AstRelation *>>          relation_list
%type <std::vector<AstClause *>>            rule
%type <std::vector<AstClause *>>            rule_def
%type <RuleBody *>                          term
%type <AstType *>                           type
%type <std::vector<AstQualifiedName>>      type_params
%type <std::vector<AstQualifiedName>>      type_param_list
%type <AstUnionType *>                      union_type_list

/* -- Destructors -- */
%destructor { delete $$; }                                  atom
%destructor { delete $$; }                                  arg
%destructor { delete $$; }                                  body
%destructor { delete $$; }                                  comp_type
%destructor { delete $$; }                                  comp_init
%destructor { delete $$; }                                  component_body
%destructor { delete $$; }                                  component_head
%destructor { delete $$; }                                  conjunction
%destructor { delete $$; }                                  constraint
%destructor { delete $$; }                                  disjunction
%destructor { delete $$; }                                  exec_order_list
%destructor { delete $$; }                                  exec_plan
%destructor { delete $$; }                                  exec_plan_list
%destructor { delete $$; }                                  fact
%destructor { delete $$; }                                  functor_decl
%destructor { }                                             predefined_type
%destructor { for (auto* cur : $$) { delete cur; } }        head
%destructor { for (auto* cur : $$) { delete cur; } }        io_directive_list
%destructor { for (auto* cur : $$) { delete cur; } }        io_relation_list
%destructor { for (auto* cur : $$) { delete cur; } }        non_empty_arg_list
%destructor { for (auto* cur : $$) { delete cur; } }        non_empty_attributes
%destructor { delete $$; }                                  non_empty_exec_order_list
%destructor { }                                             non_empty_functor_arg_type_list
%destructor { }                                             non_empty_key_value_pairs
%destructor { delete $$; }                                  non_empty_record_type_list
%destructor { delete $$; }                                  pragma
%destructor { }                                             relation_tags
%destructor { for (auto* cur : $$) { delete cur; } }        relation_decl
%destructor { for (auto* cur : $$) { delete cur; } }        relation_list
%destructor { for (auto* cur : $$) { delete cur; } }        rule
%destructor { for (auto* cur : $$) { delete cur; } }        rule_def
%destructor { for (auto* cur : $$) { delete cur; } }        io_head
%destructor { delete $$; }                                  term
%destructor { delete $$; }                                  type
%destructor { }                                             type_params
%destructor { }                                             type_param_list
%destructor { delete $$; }                                  union_type_list

/* -- Operator precedence -- */
%left L_OR
%left L_AND
%left BW_OR
%left BW_XOR
%left BW_AND
%left BW_SHIFT_L BW_SHIFT_R BW_SHIFT_R_UNSIGNED
%left PLUS MINUS
%left STAR SLASH PERCENT
%precedence NEG BW_NOT L_NOT
%right CARET

/* -- Grammar -- */
%%

%start program;

/* Program */
program
  : unit
  ;

/* Top-level statement */
unit
  : unit type {
        driver.addType(std::unique_ptr<AstType>($type));

        $type = nullptr;
    }
  | unit functor_decl {
        driver.addFunctorDeclaration(std::unique_ptr<AstFunctorDeclaration>($functor_decl));

        $functor_decl = nullptr;
    }
  | unit relation_decl {
        for (auto* cur : $relation_decl) {
            if (cur->hasQualifier(RelationQualifier::INPUT)) {
                auto load = std::make_unique<AstIO>();
                load->setSrcLoc(cur->getSrcLoc());
                load->setType(AstIO::InputIO);
                load->setQualifiedName(cur->getQualifiedName());
                driver.addIO(std::move(load));
            }
            if (cur->hasQualifier(RelationQualifier::OUTPUT)) {
                auto store = std::make_unique<AstIO>();
                store->setSrcLoc(cur->getSrcLoc());
                store->setType(AstIO::OutputIO);
                store->setQualifiedName(cur->getQualifiedName());
                driver.addIO(std::move(store));
            }
            if (cur->hasQualifier(RelationQualifier::PRINTSIZE)) {
                auto printSize = std::make_unique<AstIO>();
                printSize->setSrcLoc(cur->getSrcLoc());
                printSize->setType(AstIO::PrintsizeIO);
                printSize->setQualifiedName(cur->getQualifiedName());
                driver.addIO(std::move(printSize));
            }
            driver.addRelation(std::unique_ptr<AstRelation>(cur));
        }

        $relation_decl.clear();
    }
  | unit io_head {
        for (auto* cur : $io_head) {
            driver.addIO(std::unique_ptr<AstIO>(cur));
        }

        $io_head.clear();
    }
  | unit fact {
        driver.addClause(std::unique_ptr<AstClause>($fact));

        $fact = nullptr;
    }
  | unit rule {
        for (auto* cur : $rule) {
            driver.addClause(std::unique_ptr<AstClause>(cur));
        }

        $rule.clear();
    }
  | unit component {
        driver.addComponent(std::unique_ptr<AstComponent>($component));

        $component = nullptr;
    }
  | unit comp_init {
        driver.addInstantiation(std::unique_ptr<AstComponentInit>($comp_init));

        $comp_init = nullptr;
    }
  | unit pragma {
        driver.addPragma(std::unique_ptr<AstPragma>($pragma));

        $pragma = nullptr;
    }
  | %empty {
    }
  ;

/**
 * Identifiers
 */

identifier
  : IDENT {
        $$.push_back($IDENT);
    }
    /* TODO (azreika): in next version: DOT -> DOUBLECOLON */
  | identifier[curr_identifier] DOT IDENT {
        $$ = $curr_identifier;
        $$.push_back($IDENT);

        $curr_identifier.clear();
    }
  ;

/**
 * Types
 */

/* Type declarations */
type
  : TYPE IDENT SUBTYPE predefined_type {
        $$ = new AstSubsetType($IDENT, $predefined_type);
        $$->setSrcLoc(@$);
    }
  | NUMBER_TYPE IDENT {
        $$ = new AstSubsetType($IDENT, TypeAttribute::Signed);
        $$->setSrcLoc(@$);

        driver.warning(@1, "Deprecated type declaration used");
    }
  | SYMBOL_TYPE IDENT {
        $$ = new AstSubsetType($IDENT, TypeAttribute::Symbol);
        $$->setSrcLoc(@$);

        driver.warning(@1, "Deprecated type declaration used");
    }
  | TYPE IDENT {
        $$ = new AstSubsetType($IDENT, TypeAttribute::Symbol);
        $$->setSrcLoc(@$);

        driver.warning(@1, "Deprecated type declaration used");
    }
  | TYPE IDENT EQUALS union_type_list {
        $$ = $union_type_list;
        $$->setQualifiedName($IDENT);
        $$->setSrcLoc(@$);

        $union_type_list = nullptr;
    }
  | TYPE IDENT EQUALS LBRACKET RBRACKET {
        $$ = new AstRecordType();
        $$->setQualifiedName($IDENT);
        $$->setSrcLoc(@$);
    }
  | TYPE IDENT EQUALS LBRACKET non_empty_record_type_list RBRACKET {
        $$ = $non_empty_record_type_list;
        $$->setQualifiedName($IDENT);
        $$->setSrcLoc(@$);

        $non_empty_record_type_list = nullptr;
    }
  ;

/* Record type argument declarations */
non_empty_record_type_list
  : IDENT COLON identifier {
        $$ = new AstRecordType();
        $$->add($IDENT, $identifier);
        $$->setSrcLoc(@$);

        $identifier.clear();
    }
  | non_empty_record_type_list[curr_record] COMMA IDENT COLON identifier {
        $$ = $curr_record;
        $$->add($IDENT, $identifier);

        $curr_record = nullptr;
        $identifier.clear();
    }
  ;

/* Union type argument declarations */
union_type_list
  : identifier {
        $$ = new AstUnionType();
        $$->add($identifier);

        $identifier.clear();
    }
  | union_type_list[curr_union] PIPE identifier {
        $$ = $curr_union;
        $$->add($identifier);

        $identifier.clear();
        $curr_union = nullptr;
    }
  ;

/**
 * Relations
 */

/* Relation declaration */
relation_decl
  : DECL relation_list LPAREN RPAREN relation_tags {
        for (auto* rel : $relation_list) {
            for (auto tag : $relation_tags) {
                if (isRelationQualifierTag(tag)) {
                    rel->addQualifier(getRelationQualifierFromTag(tag));
                } else if (isRelationRepresentationTag(tag)) {
                    rel->setRepresentation(getRelationRepresentationFromTag(tag));
                } else {
                    assert(false && "unhandled tag");
                }
            }
        }
        $$ = $relation_list;

        $relation_list.clear();
    }
  | DECL relation_list LPAREN non_empty_attributes RPAREN relation_tags {
        for (auto* rel : $relation_list) {
            for (auto tag : $relation_tags) {
                if (isRelationQualifierTag(tag)) {
                    rel->addQualifier(getRelationQualifierFromTag(tag));
                } else if (isRelationRepresentationTag(tag)) {
                    rel->setRepresentation(getRelationRepresentationFromTag(tag));
                } else {
                    assert(false && "unhandled tag");
                }
            }
            for (auto* attr : $non_empty_attributes) {
                rel->addAttribute(std::unique_ptr<AstAttribute>(attr->clone()));
            }
        }
        $$ = $relation_list;

        $relation_list.clear();
    }
  ;

/* List of relation names to declare */
relation_list
  : IDENT {
        auto* rel = new AstRelation();
        rel->setQualifiedName($IDENT);
        rel->setSrcLoc(@$);

        $$.push_back(rel);
    }
  | relation_list[curr_list] COMMA IDENT {
        auto* rel = new AstRelation();
        rel->setQualifiedName($IDENT);
        rel->setSrcLoc(@IDENT);

        $$ = $curr_list;
        $$.push_back(rel);

        $curr_list.clear();
    }
  ;

/* Attribute definition of a relation */
non_empty_attributes
  : IDENT COLON identifier {
        auto attr = new AstAttribute($IDENT, $identifier);
        attr->setSrcLoc(@identifier);

        $$.push_back(attr);

        $identifier.clear();
    }
  | non_empty_attributes[curr_list] COMMA IDENT COLON identifier {
        auto attr = new AstAttribute($IDENT, $identifier);
        attr->setSrcLoc(@identifier);

        $$ = $curr_list;
        $$.push_back(attr);

        $curr_list.clear();
        $identifier.clear();
    }
  ;

/* Relation tags */
relation_tags
  : relation_tags OUTPUT_QUALIFIER {
        driver.warning(@2, "Deprecated output qualifier used");
        if ($1.find(RelationTag::OUTPUT) != $1.end())
            driver.error(@2, "output qualifier already set");
        $1.insert(RelationTag::OUTPUT);
        $$ = $1;
    }
  | relation_tags INPUT_QUALIFIER {
        driver.warning(@2, "Deprecated input qualifier was used");
        if ($1.find(RelationTag::INPUT) != $1.end())
            driver.error(@2, "input qualifier already set");
        $1.insert(RelationTag::INPUT);
        $$ = $1;
    }
  | relation_tags PRINTSIZE_QUALIFIER {
        driver.warning(@2, "Deprecated printsize qualifier was used");
        if ($1.find(RelationTag::PRINTSIZE) != $1.end())
            driver.error(@2, "printsize qualifier already set");
        $1.insert(RelationTag::PRINTSIZE);
        $$ = $1;
    }
  | relation_tags OVERRIDABLE_QUALIFIER {
        if ($1.find(RelationTag::OVERRIDABLE) != $1.end())
            driver.error(@2, "overridable qualifier already set");
        $1.insert(RelationTag::OVERRIDABLE);
        $$ = $1;
    }
  | relation_tags INLINE_QUALIFIER {
        if ($1.find(RelationTag::INLINE) != $1.end())
            driver.error(@2, "inline qualifier already set");
        $1.insert(RelationTag::INLINE);
        $$ = $1;
    }
  | relation_tags BRIE_QUALIFIER {
        if ($1.find(RelationTag::BRIE) != $1.end() ||
            $1.find(RelationTag::BTREE) != $1.end() ||
            $1.find(RelationTag::EQREL) != $1.end())
                driver.error(@2, "btree/brie/eqrel qualifier already set");
        $1.insert(RelationTag::BRIE);
        $$ = $1;
    }
  | relation_tags BTREE_QUALIFIER {
        if ($1.find(RelationTag::BRIE) != $1.end() ||
            $1.find(RelationTag::BTREE) != $1.end() ||
            $1.find(RelationTag::EQREL) != $1.end())
                driver.error(@2, "btree/brie/eqrel qualifier already set");
        $1.insert(RelationTag::BTREE);
        $$ = $1;
    }
  | relation_tags EQREL_QUALIFIER {
        if ($1.find(RelationTag::BRIE) != $1.end() ||
            $1.find(RelationTag::BTREE) != $1.end() ||
            $1.find(RelationTag::EQREL) != $1.end())
                driver.error(@2, "btree/brie/eqrel qualifier already set");
        $1.insert(RelationTag::EQREL);
        $$ = $1;
    }
  | %empty {
        $$ = std::set<RelationTag>();
    }
  ;

/**
 * Datalog Rule Structure
 */

/* Fact */
fact
  : atom DOT {
        $$ = new AstClause(std::unique_ptr<AstAtom>($atom), {}, {});
        $$->setSrcLoc(@$);

        $atom = nullptr;
    }
  ;

/* Rule */
rule
  : rule_def {
        $$ = $rule_def;

        $rule_def.clear();
    }
  | rule[nested_rule] exec_plan {
        $$ = $nested_rule;
        for (auto* rule : $$) {
            rule->setExecutionPlan(std::unique_ptr<AstExecutionPlan>($exec_plan->clone()));
        }

        $nested_rule.clear();
    }
  ;

/* Rule definition */
rule_def
  : head IF body DOT {
        auto heads = $head;
        auto bodies = $body->toClauseBodies();

        for (const auto* head : heads) {
            for (const auto* body : bodies) {
                AstClause* cur = body->clone();
                cur->setHead(std::unique_ptr<AstAtom>(head->clone()));
                cur->setSrcLoc(@$);
                $$.push_back(cur);
            }
        }

        for (auto* body : bodies) {
            delete body;
        }
    }
  ;

/* Rule head */
head
  : atom {
        $$.push_back($atom);

        $atom = nullptr;
    }
  | head[curr_head] COMMA atom {
        $$ = $curr_head;
        $$.push_back($atom);

        $curr_head.clear();
        $atom = nullptr;
    }
  ;

/* Rule body */
body
  : disjunction {
        $$ = $disjunction;

        $disjunction = nullptr;
    }
  ;

/* Rule body disjunction */
disjunction
  : conjunction {
        $$ = $conjunction;

        $conjunction = nullptr;
    }
  | disjunction[curr_disjunction] SEMICOLON conjunction {
        $$ = $curr_disjunction;
        $$->disjunct(std::move(*$conjunction));

        $curr_disjunction = nullptr;
    }
  ;

/* Rule body conjunction */
conjunction
  : term {
        $$ = $term;

        $term = nullptr;
    }
  | conjunction[curr_conjunction] COMMA term {
        $$ = $curr_conjunction;
        $$->conjunct(std::move(*$term));

        $curr_conjunction = nullptr;
    }
  ;

/* Rule execution plan */
exec_plan
  : PLAN exec_plan_list {
        $$ = $exec_plan_list;

        $exec_plan_list = nullptr;
    }
  ;

/* Rule execution plan list */
exec_plan_list
  : NUMBER COLON LPAREN exec_order_list RPAREN {
        $exec_order_list->setSrcLoc(@LPAREN);
        $$ = new AstExecutionPlan();
        $$->setOrderFor(RamSignedFromString($NUMBER), std::unique_ptr<AstExecutionOrder>($exec_order_list));

        $exec_order_list = nullptr;
    }
  | exec_plan_list[curr_list] COMMA NUMBER COLON LPAREN exec_order_list RPAREN {
        $exec_order_list->setSrcLoc(@LPAREN);
        $$ = $curr_list;
        $$->setOrderFor(RamSignedFromString($NUMBER), std::unique_ptr<AstExecutionOrder>($exec_order_list));

        $curr_list = nullptr;
        $exec_order_list = nullptr;
    }
  ;

/* Rule execution order */
exec_order_list
  : non_empty_exec_order_list {
        $$ = $non_empty_exec_order_list;

        $non_empty_exec_order_list = nullptr;
    }
  | %empty {
        $$ = new AstExecutionOrder();
    }
  ;
non_empty_exec_order_list
  : NUMBER {
        $$ = new AstExecutionOrder();
        $$->appendAtomIndex(RamSignedFromString($NUMBER));
    }
  | non_empty_exec_order_list[curr_list] COMMA NUMBER {
        $$ = $curr_list;
        $$->appendAtomIndex(RamSignedFromString($NUMBER));

        $curr_list = nullptr;
    }
  ;

/**
 * Terms in Rule Bodies
 */

/* Rule body term */
term
  : atom {
        $$ = new RuleBody(RuleBody::atom($atom));

        $atom = nullptr;
    }
  | constraint {
        $$ = new RuleBody(RuleBody::constraint($constraint));

        $constraint = nullptr;
    }
  | EXCLAMATION term[nested_term] {
        $$ = $nested_term;
        $$->negate();

        $nested_term = nullptr;
    }
  | LPAREN disjunction RPAREN {
        $$ = $disjunction;

        $disjunction = nullptr;
    }
  ;

/* Rule body atom */
atom
  : identifier LPAREN non_empty_arg_list RPAREN {
        $$ = new AstAtom();

        for (auto* arg : $non_empty_arg_list) {
            $$->addArgument(std::unique_ptr<AstArgument>(arg));
        }

        $$->setQualifiedName($identifier);
        $$->setSrcLoc(@$);

        $identifier.clear();
        $non_empty_arg_list.clear();
    }
  | identifier LPAREN RPAREN {
        $$ = new AstAtom();
        $$->setQualifiedName($identifier);
        $$->setSrcLoc(@$);

        $identifier.clear();
    }
  ;

/* Rule literal constraints */
constraint
    /* binary infix constraints */
  : arg[left] LT arg[right] {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::LT,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] GT arg[right] {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::GT,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] LE arg[right] {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::LE,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] GE arg[right] {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::GE,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] EQUALS arg[right] {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::EQ,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] NE arg[right] {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::NE,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }

    /* binary prefix constraints */
  | TMATCH LPAREN arg[left] COMMA arg[right] RPAREN {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::MATCH,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | TCONTAINS LPAREN arg[left] COMMA arg[right] RPAREN {
        $$ = new AstBinaryConstraint(BinaryConstraintOp::CONTAINS,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }

    /* choice constraint */ 
  | IDENT RIGHTARROW IDENT  {
        $$ = new AstBooleanConstraint(true);
        $$->setSrcLoc(@$);
  }

    /* zero-arity constraints */
  | TRUE {
        $$ = new AstBooleanConstraint(true);
        $$->setSrcLoc(@$);
    }
  | FALSE {
        $$ = new AstBooleanConstraint(false);
        $$->setSrcLoc(@$);
    }
  ;

/* Argument list */
non_empty_arg_list
  : arg {
        $$.push_back($arg);

        $arg = nullptr;
    }
  | non_empty_arg_list[curr_arg_list] COMMA arg {
        $$ = $curr_arg_list;
        $$.push_back($arg);

        $curr_arg_list.clear();
        $arg = nullptr;
    }
  ;

/* Atom argument */
arg
  : STRING {
        $$ = new AstStringConstant($STRING);
        $$->setSrcLoc(@$);
    }
  | FLOAT {
        $$ = new AstNumericConstant($FLOAT, AstNumericConstant::Type::Float);
        $$->setSrcLoc(@$);
    }
  | NUMBER {
        $$ = new AstNumericConstant($NUMBER);
        $$->setSrcLoc(@$);
    }
  | UNDERSCORE {
        $$ = new AstUnnamedVariable();
        $$->setSrcLoc(@$);
    }
  | DOLLAR {
        $$ = new AstCounter();
        $$->setSrcLoc(@$);
    }
  | IDENT {
        $$ = new AstVariable($IDENT);
        $$->setSrcLoc(@$);
    }
  | LPAREN arg[nested_arg] RPAREN {
        $$ = $nested_arg;

        $nested_arg = nullptr;
    }

    /* type-cast */
  | AS LPAREN arg[nested_arg] COMMA identifier RPAREN {
        $$ = new AstTypeCast(std::unique_ptr<AstArgument>($nested_arg), $identifier);
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
        $identifier.clear();
    }

    /* record constructor */
  | NIL {
        $$ = new AstNilConstant();
        $$->setSrcLoc(@$);
    }
  /* TODO (azreika): in next version: prepend records with identifiers */
  | LBRACKET RBRACKET {
        $$ = new AstRecordInit();
        $$->setSrcLoc(@$);
    }
  | LBRACKET non_empty_arg_list RBRACKET {
        auto record = new AstRecordInit();

        for (auto* arg : $non_empty_arg_list) {
            record->addArgument(std::unique_ptr<AstArgument>(arg));
        }

        $$ = record;
        $$->setSrcLoc(@$);

        $non_empty_arg_list.clear();
    }

    /* user-defined functor */
  | AT IDENT LPAREN RPAREN {
        auto functor = new AstUserDefinedFunctor($IDENT);
        $$ = functor;
        $$->setSrcLoc(@$);
    }
  | AT IDENT LPAREN non_empty_arg_list RPAREN {
        auto functor = new AstUserDefinedFunctor($IDENT);

        for (auto* arg : $non_empty_arg_list) {
            functor->addArgument(std::unique_ptr<AstArgument>(arg));
        }

        $$ = functor;
        $$->setSrcLoc(@$);

        $non_empty_arg_list.clear();
    }

    /* -- intrinsic functor -- */
    /* unary functors */
  | MINUS arg[nested_arg] %prec NEG {

        // If we have a constant, that is not already negated we create a new constant.
        const auto* asNumeric = dynamic_cast<const AstNumericConstant*>($nested_arg);
        if (asNumeric && !isPrefix("-", asNumeric->getConstant())) {
            $$ = new AstNumericConstant("-" + asNumeric->getConstant(), asNumeric->getType());
            $$->setSrcLoc(@nested_arg);

        // Otherwise, create a functor.
        } else {
            $$ = new AstIntrinsicFunctor(FunctorOp::NEG,
                std::unique_ptr<AstArgument>($nested_arg));
            $nested_arg = nullptr;
            $$->setSrcLoc(@$);
        }
    }
  | BW_NOT arg[nested_arg] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BNOT,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | L_NOT arg [nested_arg] {
        $$ = new AstIntrinsicFunctor(FunctorOp::LNOT,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | ORD LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::ORD,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | STRLEN LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::STRLEN,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | TONUMBER LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::TONUMBER,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | TOSTRING LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::TOSTRING,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | ITOU LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::ITOU,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | ITOF LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::ITOF,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | UTOI LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::UTOI,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | UTOF LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::UTOF,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | FTOI LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::FTOI,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }
  | FTOU LPAREN arg[nested_arg] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::FTOU,
                std::unique_ptr<AstArgument>($nested_arg));
        $$->setSrcLoc(@$);

        $nested_arg = nullptr;
    }

    /* binary infix functors */
  | arg[left] PLUS arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::ADD,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] MINUS arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::SUB,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] STAR arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::MUL,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] SLASH arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::DIV,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] PERCENT arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::MOD,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] CARET arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::EXP,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] BW_SHIFT_L arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BSHIFT_L,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] BW_SHIFT_R arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BSHIFT_R,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] BW_SHIFT_R_UNSIGNED arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BSHIFT_R_UNSIGNED,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] BW_OR arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BOR,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] BW_XOR arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BXOR,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] BW_AND arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::BAND,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] L_OR arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::LOR,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }
  | arg[left] L_AND arg[right] {
        $$ = new AstIntrinsicFunctor(FunctorOp::LAND,
                std::unique_ptr<AstArgument>($left),
                std::unique_ptr<AstArgument>($right));
        $$->setSrcLoc(@$);

        $left = nullptr;
        $right = nullptr;
    }

    /* binary (or more) prefix functors */
  | MAX LPAREN arg[first] COMMA non_empty_arg_list[rest] RPAREN {
        std::vector<std::unique_ptr<AstArgument>> args;
        args.emplace_back($first);

        for (auto* arg : $rest) {
            args.emplace_back(arg);
        }

        $$ = new AstIntrinsicFunctor(FunctorOp::MAX, std::move(args));
        $$->setSrcLoc(@$);

        $first = nullptr;
        $rest.clear();
    }
  | MIN LPAREN arg[first] COMMA non_empty_arg_list[rest] RPAREN {
        std::vector<std::unique_ptr<AstArgument>> args;
        args.emplace_back($first);

        for (auto* arg : $rest) {
            args.emplace_back(arg);
        }

        $$ = new AstIntrinsicFunctor(FunctorOp::MIN, std::move(args));
        $$->setSrcLoc(@$);

        $first = nullptr;
        $rest.clear();
    }
  | CAT LPAREN arg[first] COMMA non_empty_arg_list[rest] RPAREN {
        std::vector<std::unique_ptr<AstArgument>> args;
        args.emplace_back($first);

        for (auto* arg : $rest) {
            args.emplace_back(arg);
        }

        $$ = new AstIntrinsicFunctor(FunctorOp::CAT, std::move(args));
        $$->setSrcLoc(@$);

        $first = nullptr;
        $rest.clear();
    }

    /* ternary functors */
  | SUBSTR LPAREN arg[first] COMMA arg[second] COMMA arg[third] RPAREN {
        $$ = new AstIntrinsicFunctor(FunctorOp::SUBSTR,
                std::unique_ptr<AstArgument>($first),
                std::unique_ptr<AstArgument>($second),
                std::unique_ptr<AstArgument>($third));
        $$->setSrcLoc(@$);

        $first = nullptr;
        $second = nullptr;
        $third = nullptr;
    }

    /* -- aggregators -- */
  | MEAN arg[target_expr] COLON atom {
        auto aggr = new AstAggregator(AggregateOp::MEAN, std::unique_ptr<AstArgument>($target_expr));

        std::vector<std::unique_ptr<AstLiteral>> body;
        body.push_back(std::unique_ptr<AstLiteral>($atom));

        aggr->setBody(std::move(body));

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
        $atom = nullptr;
    }
  | MEAN arg[target_expr] COLON LBRACE body RBRACE {
        auto aggr = new AstAggregator(AggregateOp::MEAN, std::unique_ptr<AstArgument>($target_expr));

        auto bodies = $body->toClauseBodies();

        if (bodies.size() != 1) {
            std::cerr << "ERROR: currently not supporting non-conjunctive aggregation clauses!";
            exit(1);
        }

        std::vector<std::unique_ptr<AstLiteral>> body;
        for (auto& cur : bodies[0]->getBodyLiterals()) {
            body.push_back(std::unique_ptr<AstLiteral>(cur->clone()));
        }
        aggr->setBody(std::move(body));
        delete bodies[0];

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
    }

  | COUNT COLON atom {
        auto aggr = new AstAggregator(AggregateOp::COUNT);

        std::vector<std::unique_ptr<AstLiteral>> body;
        body.push_back(std::unique_ptr<AstLiteral>($atom));

        aggr->setBody(std::move(body));

        $$ = aggr;
        $$->setSrcLoc(@$);

        $atom = nullptr;
    }
  | COUNT COLON LBRACE body RBRACE {
        auto aggr = new AstAggregator(AggregateOp::COUNT);

        auto bodies = $body->toClauseBodies();

        if (bodies.size() != 1) {
            std::cerr << "ERROR: currently not supporting non-conjunctive aggregation clauses!";
            exit(1);
        }

        std::vector<std::unique_ptr<AstLiteral>> body;
        for (auto& cur : bodies[0]->getBodyLiterals()) {
            body.push_back(std::unique_ptr<AstLiteral>(cur->clone()));
        }
        aggr->setBody(std::move(body));
        delete bodies[0];

        $$ = aggr;
        $$->setSrcLoc(@$);
    }

  | SUM arg[target_expr] COLON atom {
        auto aggr = new AstAggregator(AggregateOp::SUM, std::unique_ptr<AstArgument>($target_expr));

        std::vector<std::unique_ptr<AstLiteral>> body;
        body.push_back(std::unique_ptr<AstLiteral>($atom));

        aggr->setBody(std::move(body));

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
        $atom = nullptr;
    }
  | SUM arg[target_expr] COLON LBRACE body RBRACE {
        auto aggr = new AstAggregator(AggregateOp::SUM, std::unique_ptr<AstArgument>($target_expr));

        auto bodies = $body->toClauseBodies();

        if (bodies.size() != 1) {
            std::cerr << "ERROR: currently not supporting non-conjunctive aggregation clauses!";
            exit(1);
        }

        std::vector<std::unique_ptr<AstLiteral>> body;
        for (auto& cur : bodies[0]->getBodyLiterals()) {
            body.push_back(std::unique_ptr<AstLiteral>(cur->clone()));
        }
        aggr->setBody(std::move(body));
        delete bodies[0];

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
    }

  | MIN arg[target_expr] COLON atom {
        auto aggr = new AstAggregator(AggregateOp::MIN, std::unique_ptr<AstArgument>($target_expr));

        std::vector<std::unique_ptr<AstLiteral>> body;
        body.push_back(std::unique_ptr<AstLiteral>($atom));

        aggr->setBody(std::move(body));
        $atom = nullptr;

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
        $atom = nullptr;
    }
  | MIN arg[target_expr] COLON LBRACE body RBRACE {
        auto aggr = new AstAggregator(AggregateOp::MIN, std::unique_ptr<AstArgument>($target_expr));

        auto bodies = $body->toClauseBodies();

        if (bodies.size() != 1) {
            std::cerr << "ERROR: currently not supporting non-conjunctive aggregation clauses!";
            exit(1);
        }

        std::vector<std::unique_ptr<AstLiteral>> body;
        for (auto& cur : bodies[0]->getBodyLiterals()) {
            body.push_back(std::unique_ptr<AstLiteral>(cur->clone()));
        }
        aggr->setBody(std::move(body));
        delete bodies[0];

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
    }

  | MAX arg[target_expr] COLON atom {
        auto aggr = new AstAggregator(AggregateOp::MAX, std::unique_ptr<AstArgument>($target_expr));

        std::vector<std::unique_ptr<AstLiteral>> body;
        body.push_back(std::unique_ptr<AstLiteral>($atom));

        aggr->setBody(std::move(body));

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
        $atom = nullptr;
    }
  | MAX arg[target_expr] COLON LBRACE body RBRACE {
        auto aggr = new AstAggregator(AggregateOp::MAX, std::unique_ptr<AstArgument>($target_expr));

        auto bodies = $body->toClauseBodies();

        if (bodies.size() != 1) {
            std::cerr << "ERROR: currently not supporting non-conjunctive aggregation clauses!";
            exit(1);
        }

        std::vector<std::unique_ptr<AstLiteral>> body;
        for (auto& cur : bodies[0]->getBodyLiterals()) {
            body.push_back(std::unique_ptr<AstLiteral>(cur->clone()));
        }
        aggr->setBody(std::move(body));
        delete bodies[0];

        $$ = aggr;
        $$->setSrcLoc(@$);

        $target_expr = nullptr;
    }
  ;

/**
 * Components
 */

/* Component */
component
  : component_head LBRACE component_body RBRACE {
        $$ = $component_body;
        auto* type = $component_head->getComponentType()->clone();
        $$->setComponentType(std::unique_ptr<AstComponentType>(type));
        $$->copyBaseComponents($component_head);
        $$->setSrcLoc(@$);

        $component_body = nullptr;
    }
  ;

/* Component head */
component_head
  : COMPONENT comp_type {
        $$ = new AstComponent();
        $$->setComponentType(std::unique_ptr<AstComponentType>($comp_type));

        $comp_type = nullptr;
    }
  | component_head[comp] COLON comp_type {
        $$ = $comp;
        $$->addBaseComponent(std::unique_ptr<AstComponentType>($comp_type));

        $comp = nullptr;
        $comp_type = nullptr;
    }
  | component_head[comp] COMMA comp_type {
        $$ = $comp;
        $$->addBaseComponent(std::unique_ptr<AstComponentType>($comp_type));

        $comp = nullptr;
        $comp_type = nullptr;
    }
  ;

/* Component type */
comp_type
  : IDENT type_params {
        $$ = new AstComponentType($IDENT, $type_params);

        $type_params.clear();
    }
  ;

/* Component type parameters */
type_params
  : LT type_param_list GT {
        $$ = $type_param_list;

        $type_param_list.clear();
    }
  | %empty {
        $$ = std::vector<AstQualifiedName>();
    }
  ;

/* Component type parameter list */
type_param_list
  : IDENT {
        $$.push_back($IDENT);
    }
  | type_param_list[curr_list] COMMA IDENT {
        $$ = $curr_list;
        $$.push_back($IDENT);

        $curr_list.clear();
    }
  ;

/* Component body */
component_body
  : component_body[comp] type {
        $$ = $comp;
        $$->addType(std::unique_ptr<AstType>($type));

        $comp = nullptr;
        $type = nullptr;
    }
  | component_body[comp] relation_decl {
        $$ = $comp;
        for (auto* rel : $relation_decl) {
            if (rel->hasQualifier(RelationQualifier::INPUT)) {
                auto load = std::make_unique<AstIO>();
                load->setSrcLoc(rel->getSrcLoc());
                load->setType(AstIO::InputIO);
                load->setQualifiedName(rel->getQualifiedName());
                driver.addIO(std::move(load));
            }
            if (rel->hasQualifier(RelationQualifier::OUTPUT)) {
                auto store = std::make_unique<AstIO>();
                store->setSrcLoc(rel->getSrcLoc());
                store->setType(AstIO::OutputIO);
                store->setQualifiedName(rel->getQualifiedName());
                driver.addIO(std::move(store));
            }
            if (rel->hasQualifier(RelationQualifier::PRINTSIZE)) {
                auto printSize = std::make_unique<AstIO>();
                printSize->setSrcLoc(rel->getSrcLoc());
                printSize->setQualifiedName(rel->getQualifiedName());
                printSize->setType(AstIO::PrintsizeIO);
                driver.addIO(std::move(printSize));
            }
            $$->addRelation(std::unique_ptr<AstRelation>(rel));
        }

        $comp = nullptr;
        $relation_decl.clear();
    }
  | component_body[comp] io_head {
        $$ = $comp;
        for (auto* io : $io_head) {
            $$->addIO(std::unique_ptr<AstIO>(io));
        }

        $comp = nullptr;
        $io_head.clear();
    }
  | component_body[comp] fact {
        $$ = $comp;
        $$->addClause(std::unique_ptr<AstClause>($fact));

        $comp = nullptr;
        $fact = nullptr;
    }
  | component_body[comp] rule {
        $$ = $comp;
        for (auto* rule : $rule) {
            $$->addClause(std::unique_ptr<AstClause>(rule));
        }

        $comp = nullptr;
        $rule.clear();
    }
  | component_body[comp] OVERRIDE IDENT {
        $$ = $comp;
        $$->addOverride($IDENT);

        $comp = nullptr;
    }
  | component_body[comp] comp_init {
        $$ = $comp;
        $$->addInstantiation(std::unique_ptr<AstComponentInit>($comp_init));

        $comp = nullptr;
        $comp_init = nullptr;
    }
  | component_body[comp] component {
        $$ = $comp;
        $$->addComponent(std::unique_ptr<AstComponent>($component));

        $comp = nullptr;
        $component = nullptr;
    }
  | %empty {
        $$ = new AstComponent();
    }
  ;

/* Component initialisation */
comp_init
  : INSTANTIATE IDENT EQUALS comp_type {
        $$ = new AstComponentInit();
        $$->setInstanceName($IDENT);
        $$->setComponentType(std::unique_ptr<AstComponentType>($comp_type));
        $$->setSrcLoc(@$);

        $comp_type = nullptr;
    }
  ;

/**
 * User-Defined Functors
 */

/* Functor declaration */
functor_decl
  : FUNCTOR IDENT LPAREN RPAREN COLON predefined_type {
        $$ = new AstFunctorDeclaration($IDENT, {}, $predefined_type);
        $$->setSrcLoc(@$);
    }
  | FUNCTOR IDENT LPAREN non_empty_functor_arg_type_list RPAREN COLON predefined_type {
        auto typesig = $non_empty_functor_arg_type_list;
        $$ = new AstFunctorDeclaration($IDENT, typesig, $predefined_type);
        $$->setSrcLoc(@$);
    }
  ;

/* Functor argument list type */
non_empty_functor_arg_type_list
  : predefined_type {
        $$.push_back($predefined_type);
    }
  | non_empty_functor_arg_type_list[curr_list] COMMA predefined_type {
        $$ = $curr_list;
        $$.push_back($predefined_type);
        $curr_list.clear();
    }
  ;

/* Predefined type */
predefined_type
  : IDENT {
        if ($IDENT == "number") {
            $$ = TypeAttribute::Signed;
        } else if ($IDENT == "symbol") {
            $$ = TypeAttribute::Symbol;
        } else if ($IDENT == "float") {
            $$ = TypeAttribute::Float;
        } else if ($IDENT == "unsigned") {
            $$ = TypeAttribute::Unsigned;
        } else {
            driver.error(@IDENT, "[number | symbol | float | unsigned] identifier expected");
        }
    }
  ;

/**
 * Other Directives
 */

/* Pragma directives */
pragma
  : PRAGMA STRING[key] STRING[value] {
        $$ = new AstPragma($key, $value);
        $$->setSrcLoc(@$);
    }
  | PRAGMA STRING[option] {
        $$ = new AstPragma($option, "");
        $$->setSrcLoc(@$);
    }
  ;

/* io directives */
io_head
  : INPUT_DECL io_directive_list {
        for (const auto* io : $io_directive_list) {
            auto load = new AstIO(*io);
            load->setType(AstIO::InputIO);
            $$.push_back(load);
        }
    }
  | OUTPUT_DECL io_directive_list {
        for (const auto* io : $io_directive_list) {
            auto store = new AstIO(*io);
            store->setType(AstIO::OutputIO);
            $$.push_back(store);
        }
    }
  | PRINTSIZE_DECL io_directive_list {
        for (const auto* io : $io_directive_list) {
            auto printsize = new AstIO(*io);
            printsize->setType(AstIO::PrintsizeIO);
            $$.push_back(printsize);
        }
    }
  ;

/* IO directive list */
io_directive_list
  : io_relation_list {
        $$ = $io_relation_list;

        $io_relation_list.clear();
    }
  | io_relation_list LPAREN RPAREN {
        $$ = $io_relation_list;

        $io_relation_list.clear();
    }
  | io_relation_list LPAREN non_empty_key_value_pairs RPAREN {
        for (auto* io : $io_relation_list) {
            for (const auto& kvp : $non_empty_key_value_pairs) {
                io->addDirective(kvp.first, kvp.second);
            }
        }
        $$ = $io_relation_list;

        $io_relation_list.clear();
        $non_empty_key_value_pairs.clear();
    }
  ;

/* IO relation list */
io_relation_list
  : identifier {
        auto* io = new AstIO();
        io->setQualifiedName($identifier);
        io->setSrcLoc(@identifier);

        $$.push_back(io);

        $identifier.clear();
    }
  | io_relation_list[curr_list] COMMA identifier {
        auto* io = new AstIO();
        io->setQualifiedName($identifier);
        io->setSrcLoc(@identifier);

        $$ = $curr_list;
        $$.push_back(io);

        $curr_list.clear();
        $identifier.clear();
    }
  ;

/* Key-value pairs */
non_empty_key_value_pairs
  : IDENT EQUALS kvp_value {
        $$.push_back(std::make_pair($IDENT, $kvp_value));
    }
  | non_empty_key_value_pairs[curr_io] COMMA IDENT EQUALS kvp_value {
        $$ = $curr_io;
        $$.push_back(std::make_pair($IDENT, $kvp_value));

        $curr_io.clear();
    }
  ;
kvp_value
  : STRING {
        $$ = $STRING;
    }
  | IDENT {
        $$ = $IDENT;
    }
  | TRUE {
        $$ = "true";
    }
  | FALSE {
        $$ = "false";
    }
  ;

%%

void yy::parser::error(const location_type &l, const std::string &m) {
    driver.error(l, m);
}
