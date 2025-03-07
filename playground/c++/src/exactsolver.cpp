#include "exactsolver.hpp"
// #include "scip/scip.h"
// #include "scip/scipdefplugins.h"
#ifdef USE_GUROBI
#include "gurobi_c++.h"
#endif
#include "EvalMaxSAT.h"
#ifdef USE_ORTOOLS
#include "ortools/linear_solver/linear_expr.h"
#include "ortools/linear_solver/linear_solver.h"
#endif

#ifdef USE_GUROBI
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

    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].can_be_dominating_set == false) {
            continue;
        }
        GRBVar var;
        var = model.addVar(0, 1, 1, GRB_BINARY);
        varmap[(*instance->G)[v].id] = var;
        SCIPVars.push_back(var);
    }
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].is_dominated) {
            continue;
        }
        GRBLinExpr expr;
        if ((*instance->G)[v].can_be_dominating_set) {
            expr += varmap[(*instance->G)[v].id];
        }
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            if ((*instance->G)[w].can_be_dominating_set) {
                expr += varmap[(*instance->G)[w].id];
            }
        }
        model.addConstr(expr >= 1);
    }

    model.optimize();
}
#endif

void solveEvalMaxSat(Instance* instance, VertexList& dominatingSet){

    EvalMaxSAT solver;

    debug("Solving ILP witn number of nodes", instance->n);

    static vector<int> varmap(instance->props->n+1);

    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].can_be_dominating_set == false) {
            continue;
        }
        auto var = solver.newVar();
        varmap[(*instance->G)[v].id] = var;
        solver.addClause({-var}, 1); // soft clause
    }

    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].is_dominated) {
            continue;
        }
        vector<int> clause;
        if ((*instance->G)[v].can_be_dominating_set) {
            clause.push_back(varmap[(*instance->G)[v].id]);
        }
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            if ((*instance->G)[w].can_be_dominating_set) {
                clause.push_back(varmap[(*instance->G)[w].id]);
            }
        }
        solver.addClause(clause); //hard clause
    }

    bool solved = solver.solve();
    debug(solved);
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].can_be_dominating_set) {       
            if (solver.getValue(varmap[(*instance->G)[v].id])) {                
                dominatingSet.push_back((*instance->G)[v].id);
            }
        }
    }
}

#ifdef USE_ORTOOLS
void solveCPSat(Instance* instance, VertexList& dominatingSet) {
    using namespace operations_research;
    debug("Solving ILP witn number of nodes", instance->n);

    std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("CP-SAT"));
    if (!solver) {
        debug("solver not available");
        return;
    }



    static vector<MPVariable*> varmap(instance->props->n+1);
    LinearExpr obj;

    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].can_be_dominating_set == false) {
            continue;
        }
        auto var = solver->MakeBoolVar("");
        varmap[(*instance->G)[v].id] = var;
        obj += var;
    }
    MPObjective* const objective = solver->MutableObjective();
    objective->MinimizeLinearExpr(obj);

    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].is_dominated) {
            continue;
        }
        LinearExpr expr;
        if ((*instance->G)[v].can_be_dominating_set) {
            expr += varmap[(*instance->G)[v].id];
        }
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            if ((*instance->G)[w].can_be_dominating_set) {
                expr += varmap[(*instance->G)[w].id];
            }
        }
        solver->MakeRowConstraint(expr >= 1);
    }
    // solver->EnableOutput();
    solver->SetNumThreads(1);
    const MPSolver::ResultStatus result_status = solver->Solve();
    // bool solved = solver.solve();
    // debug(solved);
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if ((*instance->G)[v].can_be_dominating_set) {       
            if (varmap[(*instance->G)[v].id]->solution_value() > 0.5) {                
                dominatingSet.push_back((*instance->G)[v].id);
            }
        }
    }
}
#endif

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

