#include "exactsolver.hpp"
// #include "scip/scip.h"
// #include "scip/scipdefplugins.h"
#include "gurobi_c++.h"
#include "EvalMaxSAT.h"

void solveGurobiExactIlp(Instance* instance, VertexList& dominatingSet) {
    GRBEnv env;
    GRBModel model(env);

    model.set(GRB_IntParam_LogToConsole, 1);   // Ensure logging is enabled
    model.set(GRB_IntParam_DisplayInterval, 1); // Log progress frequently
    model.set(GRB_DoubleParam_Heuristics, 0.0);  // Disable heuristics entirely
    model.set(GRB_DoubleParam_NoRelHeurTime, 0); // Disable NoRel heuristic

    // model.set(GRB_IntParam_PoolSearchMode, 2); // Store multiple solutions
    // model.set(GRB_DoubleParam_Heuristics, 0.5); // Increase heuristic effort (optional)
    // model.set(GRB_DoubleParam_NoRelHeurTime, 10); // Allow NoRel heuristic extra time




    debug("Solving ILP witn number of nodes", instance->n);

    static vector<GRBVar> varmap(instance->props->n+1);
    vector<GRBVar> SCIPVars;

    for (VD v: MyVertices(*instance->G)) {
        if ((*instance->G)[v].can_be_dominating_set == false) {
            continue;
        }
        GRBVar var;
        var = model.addVar(0, 1, 1, GRB_BINARY);
        varmap[(*instance->G)[v].id] = var;
        SCIPVars.push_back(var);
    }

    for (VD v: MyVertices(*instance->G)) {
        if ((*instance->G)[v].is_dominated) {
            continue;
        }
        GRBLinExpr expr;
        if ((*instance->G)[v].can_be_dominating_set) {
            expr += varmap[(*instance->G)[v].id];
        }
        for (VD w: Neighbors(*instance->G, v)) {
            if ((*instance->G)[w].can_be_dominating_set) {
                expr += varmap[(*instance->G)[w].id];
            }
        }
        model.addConstr(expr >= 1);
    }

    model.optimize();
}

void solveEvalMaxSat(Instance* instance, VertexList& dominatingSet){

    EvalMaxSAT solver;

    debug("Solving ILP witn number of nodes", instance->n);

    static vector<int> varmap(instance->props->n+1);

    for (VD v: MyVertices(*instance->G)) {
        if ((*instance->G)[v].can_be_dominating_set == false) {
            continue;
        }
        auto var = solver.newVar();
        varmap[(*instance->G)[v].id] = var;
        solver.addClause({-var}, 1); // soft clause
    }

    for (VD v: MyVertices(*instance->G)) {
        if ((*instance->G)[v].is_dominated) {
            continue;
        }
        vector<int> clause;
        if ((*instance->G)[v].can_be_dominating_set) {
            clause.push_back(varmap[(*instance->G)[v].id]);
        }
        for (VD w: Neighbors(*instance->G, v)) {
            if ((*instance->G)[w].can_be_dominating_set) {
                clause.push_back(varmap[(*instance->G)[w].id]);
            }
        }
        solver.addClause(clause); //hard clause
    }

    bool solved = solver.solve();
    debug(solved);
    for (VD v: MyVertices(*instance->G)) {
        if ((*instance->G)[v].can_be_dominating_set) {       
            if (solver.getValue(varmap[(*instance->G)[v].id])) {                
                dominatingSet.push_back((*instance->G)[v].id);
            }
        }
    }
}

// SCIP_RETCODE solveScip(Instance* instance, VertexList& dominatingSet) {
//     SCIP *scip = nullptr;
//     SCIP_CALL(SCIPcreate(&scip)); // Creating the SCIP environment

//     /* include default plugins */
//     SCIP_CALL(SCIPincludeDefaultPlugins(scip));
// #ifndef MYLOCAL
//     SCIPsetMessagehdlrQuiet(scip, 1);
// #endif
//     SCIP_CALL(SCIPcreateProbBasic(scip, "dominating set"));
//     SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
//     SCIP_CALL(SCIPsetIntParam(scip, "lp/threads", 1));
//     SCIP_CALL(SCIPsetIntParam(scip, "separating/maxstallroundsroot", -1));
//     SCIP_CALL(SCIPsetIntParam(scip, "separating/zerohalf/freq", 1));
//     SCIP_CALL(SCIPsetIntParam(scip, "separating/gomory/freq", 1));

//     debug("Solving ILP witn number of nodes", instance->n);

//     static vector<SCIP_VAR*> varmap(instance->props->n+1);
//     vector<SCIP_VAR*> SCIPVars;
//     vector<SCIP_CONS*> SCIPCons;
//     #ifdef MYLOCAL
//     #else
//         // TODO No logging
//         string collector;
//         BufferStdout swapstdout(collector);
//     #endif

//     for (VD v: Vertices(*instance->G)) {
//         if ((*instance->G)[v].can_be_dominating_set == false) {
//             continue;
//         }
//         SCIP_VAR* var = nullptr;
//         SCIP_CALL(SCIPcreateVarBasic(scip, &var, "", 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY));
//         SCIP_CALL(SCIPaddVar(scip, var));
//         varmap[(*instance->G)[v].id] = var;
//         SCIPVars.push_back(var);
//     }

//     for (VD v: Vertices(*instance->G)) {
//         if ((*instance->G)[v].is_dominated) {
//             continue;
//         }
//         SCIP_CONS* cons = nullptr;
//         vector<SCIP_VAR*> vars;
//         vector<double> vals;
//         if ((*instance->G)[v].can_be_dominating_set) {
//             vars.push_back(varmap[(*instance->G)[v].id]);
//             vals.push_back(1.0);
//         }
//         for (VD w: Neighbors(*instance->G, v)) {
//             if ((*instance->G)[w].can_be_dominating_set) {
//                 vars.push_back(varmap[(*instance->G)[w].id]);
//                 vals.push_back(1.0);
//             }
//         }
//         SCIP_CALL(SCIPcreateConsBasicSetcover(scip, &cons, "", vars.size(), vars.data()));
//         SCIPCons.push_back(cons);
//         SCIPaddCons(scip, cons);
//         // SCIPreleaseCons(scip, &cons);
//     }

//     #ifdef MYLOCAL
// #else
//     // TODO No logging
//     string collector;
//     BufferStdout swapstdout(collector);
// #endif
//     SCIP_CALL(SCIPsolve(scip));
//     debug("solved");

//     SCIP_SOL *sol = nullptr;
//     sol = SCIPgetBestSol(scip);
//     debug(SCIPgetSolOrigObj(scip, sol));

//     for (auto var: SCIPVars)
//     {
//         SCIP_CALL(SCIPreleaseVar(scip, &var));
//     }
//     for (auto cons: SCIPCons)
//     {
//         SCIP_CALL(SCIPreleaseCons(scip, &cons));
//     }
//     SCIP_CALL(SCIPfree(&scip));
//     return SCIP_OKAY;
// }

// void solveScipExactILP(Instance* instance, VertexList& dominatingSet) {
//     solveScip(instance, dominatingSet);
// }

