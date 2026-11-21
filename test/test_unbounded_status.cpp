// test_unbounded_status.cpp - Regression tests for false "optimal" status
// on unbounded/infeasible LPs (GitHub issues #346, #347, #348).
//
// In a values pass (ClpSimplex::primal(1), the default for a cold start),
// an incoming variable with no blocking row was always rejected (flagged)
// without checking whether it could improve the objective without limit.
// Flagging dropped the only dual infeasibility, so the solver went on to
// report "optimal" for unbounded LPs (and even for an infeasible LP whose
// entering column was flagged the same way).
//
// Build:
//   g++ -O2 -std=c++17 -o test_unbounded_status test_unbounded_status.cpp \
//     $(pkg-config --cflags --libs clp)
//
// Run:
//   ./test_unbounded_status

#include "ClpSimplex.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include <cmath>
#include <cstdio>

static int nPass = 0, nFail = 0;

#define CHECK(cond, msg)                    \
  do {                                      \
    if (!(cond)) {                          \
      fprintf(stderr, "  FAIL: %s\n", msg); \
      nFail++;                              \
    } else {                                \
      nPass++;                              \
    }                                       \
  } while (0)

// Run primal both as a values pass (cold-start default) and as a normal
// pass, and check the reported status: 0 optimal, 1 infeasible, 2 unbounded.
static void checkStatus(ClpSimplex &model, int expected, const char *label)
{
  char msg[160];
  ClpSimplex valuesPass(model);
  valuesPass.primal(1);
  snprintf(msg, sizeof(msg), "%s: values-pass primal status %d, expected %d",
    label, valuesPass.status(), expected);
  CHECK(valuesPass.status() == expected, msg);

  ClpSimplex normalPass(model);
  normalPass.primal(0);
  snprintf(msg, sizeof(msg), "%s: normal primal status %d, expected %d",
    label, normalPass.status(), expected);
  CHECK(normalPass.status() == expected, msg);
}

int main()
{
  const double inf = COIN_DBL_MAX;

  {
    // Issue #348: c0 improves the objective (cost -4), has a lower bound
    // (-4), no upper bound and no matrix entries, so the LP is unbounded.
    // -primalS used to report "Optimal" after 0 iterations.
    CoinPackedMatrix matrix(true, 0, 0);
    matrix.setDimensions(0, 2);
    CoinPackedVector row;
    row.insert(1, -3.0); // r0: -3 c1 <= 4
    matrix.appendRow(row);
    double colLower[] = { -4.0, 0.0 };
    double colUpper[] = { inf, 0.0 };
    double objective[] = { -4.0, 0.0 };
    double rowLower[] = { -inf };
    double rowUpper[] = { 4.0 };
    ClpSimplex model;
    model.loadProblem(matrix, colLower, colUpper, objective, rowLower, rowUpper);
    checkStatus(model, 2, "issue #348 (unbounded, empty column)");
  }

  {
    // Issue #347: unbounded 3-variable LP; the default solve reported
    // "Optimal" with objective 0.
    //   min  -2 x + 5 z
    //   s.t. -x + y - 0.5 z >= 0
    //        0 <= -y + z <= 2
    //        2 y - 3 z >= 0
    //        x, y free; z <= 5, no lower bound
    CoinPackedMatrix matrix(true, 0, 0);
    matrix.setDimensions(0, 3);
    CoinPackedVector r0, r1, r2;
    r0.insert(0, -1.0);
    r0.insert(1, 1.0);
    r0.insert(2, -0.5);
    matrix.appendRow(r0);
    r1.insert(1, -1.0);
    r1.insert(2, 1.0);
    matrix.appendRow(r1);
    r2.insert(1, 2.0);
    r2.insert(2, -3.0);
    matrix.appendRow(r2);
    double colLower[] = { -inf, -inf, -inf };
    double colUpper[] = { inf, inf, 5.0 };
    double objective[] = { -2.0, 0.0, 5.0 };
    double rowLower[] = { 0.0, 0.0, 0.0 };
    double rowUpper[] = { inf, 2.0, inf };
    ClpSimplex model;
    model.loadProblem(matrix, colLower, colUpper, objective, rowLower, rowUpper);
    // Direct primal on the raw model is unbounded with or without the fix.
    checkStatus(model, 2, "issue #347 (unbounded)");
    // The reported failure is in the default solve, which presolves first:
    // the presolved model's values-pass primal hit the same bug as #348.
    OsiClpSolverInterface presolved;
    presolved.loadProblem(matrix, colLower, colUpper, objective, rowLower, rowUpper);
    presolved.initialSolve();
    CHECK(presolved.isProvenDualInfeasible() && !presolved.isProvenOptimal(),
      "issue #347: default solve with presolve reports unbounded, not optimal");
  }

  {
    // Issue #346: row r0 has no entries and must equal -5, so the LP is
    // infeasible. The reporter's "Optimal" misdiagnosis was observed on
    // MinGW; on Linux/master this already reports infeasible without the
    // fix, so these checks are a guard against future regressions.
    CoinPackedMatrix matrix(true, 0, 0);
    matrix.setDimensions(0, 2);
    CoinPackedVector r0; // empty equality row
    matrix.appendRow(r0);
    CoinPackedVector r1;
    r1.insert(0, 3.0);
    r1.insert(1, -0.5);
    matrix.appendRow(r1);
    double colLower[] = { -inf, -inf };
    double colUpper[] = { inf, inf };
    double objective[] = { 0.0, -7.0 };
    double rowLower[] = { -5.0, -inf };
    double rowUpper[] = { -5.0, inf };
    ClpSimplex model;
    model.loadProblem(matrix, colLower, colUpper, objective, rowLower, rowUpper);
    checkStatus(model, 1, "issue #346 (infeasible)");
  }

  {
    // Bounded control: same shape as issue #348 but c0 has an upper bound
    // of 10. Optimum is c0 = 10, objective -40; must stay "optimal".
    CoinPackedMatrix matrix(true, 0, 0);
    matrix.setDimensions(0, 2);
    CoinPackedVector row;
    row.insert(1, -3.0);
    matrix.appendRow(row);
    double colLower[] = { -4.0, 0.0 };
    double colUpper[] = { 10.0, 0.0 };
    double objective[] = { -4.0, 0.0 };
    double rowLower[] = { -inf };
    double rowUpper[] = { 4.0 };
    ClpSimplex model;
    model.loadProblem(matrix, colLower, colUpper, objective, rowLower, rowUpper);
    checkStatus(model, 0, "bounded control (optimal)");
    ClpSimplex valuesPass(model);
    valuesPass.primal(1);
    CHECK(std::fabs(valuesPass.objectiveValue() - (-40.0)) < 1e-6,
      "bounded control: objective value -40");
  }

  printf("=== Results: %d passed, %d failed ===\n", nPass, nFail);
  return nFail == 0 ? 0 : 1;
}
