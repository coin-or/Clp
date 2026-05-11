// test_racing_lp.cpp — Tests for OsiClpSolverInterface LP solve correctness
// Tests both normal initialSolve and racing mode on programmatically-generated
// problems: N-Queens LP relaxation, TSP-MTZ LP relaxation, infeasible LP.
//
// Build:
//   g++ -O2 -std=c++17 -o test_racing_lp test_racing_lp.cpp \
//     $(pkg-config --cflags --libs osi-clp) -lpthread
//
// Run:
//   ./test_racing_lp

#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "ClpSimplex.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

static const double kTol = 1e-6;
static int nPass = 0, nFail = 0;

#define CHECK(cond, msg) do { \
  if (!(cond)) { \
    fprintf(stderr, "  FAIL: %s\n", msg); \
    nFail++; \
  } else { \
    nPass++; \
  } \
} while(0)

// Verify LP optimality conditions on a solved OsiSolverInterface.
// Returns number of violations found.
static int checkOptimalityConditions(OsiSolverInterface &si, const char *label)
{
  if (!si.isProvenOptimal())
    return 0; // only check optimal solutions

  const double fTol = 1e-5;
  int violations = 0;
  int nCols = si.getNumCols();
  int nRows = si.getNumRows();
  const double *sol = si.getColSolution();
  const double *rowAct = si.getRowActivity();
  const double *rc = si.getReducedCost();
  const double *dual = si.getRowPrice();
  const double *obj = si.getObjCoefficients();
  const double *colLb = si.getColLower();
  const double *colUb = si.getColUpper();
  const double *rowLb = si.getRowLower();
  const double *rowUb = si.getRowUpper();

  // 1. Primal feasibility: column bounds
  for (int j = 0; j < nCols; j++) {
    if (sol[j] < colLb[j] - fTol || sol[j] > colUb[j] + fTol) {
      if (violations < 3)
        fprintf(stderr, "  %s: col %d sol=%.8g not in [%.8g, %.8g]\n",
          label, j, sol[j], colLb[j], colUb[j]);
      violations++;
    }
  }

  // 2. Primal feasibility: row bounds (Ax within row bounds)
  for (int i = 0; i < nRows; i++) {
    if (rowAct[i] < rowLb[i] - fTol || rowAct[i] > rowUb[i] + fTol) {
      if (violations < 3)
        fprintf(stderr, "  %s: row %d act=%.8g not in [%.8g, %.8g]\n",
          label, i, rowAct[i], rowLb[i], rowUb[i]);
      violations++;
    }
  }

  // 3. Reduced cost sign (dual feasibility for minimization):
  //    If x[j] is at lower bound → rc[j] >= -fTol
  //    If x[j] is at upper bound → rc[j] <= fTol
  //    If x[j] is strictly between bounds → rc[j] ≈ 0
  for (int j = 0; j < nCols; j++) {
    bool atLb = (sol[j] - colLb[j] < fTol);
    bool atUb = (colUb[j] - sol[j] < fTol);
    if (atLb && !atUb && rc[j] < -fTol) {
      if (violations < 3)
        fprintf(stderr, "  %s: col %d at lb, rc=%.8g < 0\n", label, j, rc[j]);
      violations++;
    }
    if (atUb && !atLb && rc[j] > fTol) {
      if (violations < 3)
        fprintf(stderr, "  %s: col %d at ub, rc=%.8g > 0\n", label, j, rc[j]);
      violations++;
    }
    if (!atLb && !atUb && fabs(rc[j]) > fTol) {
      if (violations < 3)
        fprintf(stderr, "  %s: col %d interior (%.8g), rc=%.8g != 0\n",
          label, j, sol[j], rc[j]);
      violations++;
    }
  }

  // 4. Objective consistency: c'x should match getObjValue()
  double computedObj = 0.0;
  for (int j = 0; j < nCols; j++)
    computedObj += obj[j] * sol[j];
  double reportedObj = si.getObjValue();
  if (fabs(computedObj - reportedObj) > fTol * (1.0 + fabs(reportedObj))) {
    fprintf(stderr, "  %s: obj mismatch computed=%.8g reported=%.8g\n",
      label, computedObj, reportedObj);
    violations++;
  }

  // 5. Dual values exist and row slacks are consistent
  (void)dual; // dual values accessed — just verify non-null
  if (!dual) violations++;

  return violations;
}

// ─── N-Queens LP relaxation ──────────────────────────────────────────────────
// Place n queens on an n×n board. Variables x[i][j] ∈ [0,1].
// Constraints: exactly 1 per row, at most 1 per column, at most 1 per diagonal.
// Maximize sum(x[i][j]) — LP relaxation value = n (trivially: 1 per row).
// This creates a non-trivial LP with n² variables and O(n) constraints.
static void buildNQueens(OsiSolverInterface &si, int n)
{
  int nCols = n * n;
  // Objective: maximize sum x[i][j] → minimize -sum x[i][j]
  std::vector<double> obj(nCols, -1.0);
  std::vector<double> colLb(nCols, 0.0);
  std::vector<double> colUb(nCols, 1.0);

  si.loadProblem(CoinPackedMatrix(), colLb.data(), colUb.data(),
    obj.data(), nullptr, nullptr);
  // Re-add columns properly
  si.reset();

  CoinPackedMatrix matrix(false, 0, 0); // row-major
  matrix.setDimensions(0, nCols);
  std::vector<double> rowLb, rowUb;

  // Row constraints: sum_j x[i][j] = 1 for each row i
  for (int i = 0; i < n; i++) {
    CoinPackedVector row;
    for (int j = 0; j < n; j++)
      row.insert(i * n + j, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);
    rowUb.push_back(1.0);
  }

  // Column constraints: sum_i x[i][j] <= 1 for each column j
  for (int j = 0; j < n; j++) {
    CoinPackedVector row;
    for (int i = 0; i < n; i++)
      row.insert(i * n + j, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(-COIN_DBL_MAX);
    rowUb.push_back(1.0);
  }

  // Diagonal constraints (/ direction): sum x[i][j] where i-j=k, <= 1
  for (int k = -(n - 2); k <= n - 2; k++) {
    CoinPackedVector row;
    for (int i = 0; i < n; i++) {
      int j = i - k;
      if (j >= 0 && j < n)
        row.insert(i * n + j, 1.0);
    }
    if (row.getNumElements() > 1) {
      matrix.appendRow(row);
      rowLb.push_back(-COIN_DBL_MAX);
      rowUb.push_back(1.0);
    }
  }

  // Anti-diagonal constraints (\ direction): sum x[i][j] where i+j=k, <= 1
  for (int k = 1; k <= 2 * n - 3; k++) {
    CoinPackedVector row;
    for (int i = 0; i < n; i++) {
      int j = k - i;
      if (j >= 0 && j < n)
        row.insert(i * n + j, 1.0);
    }
    if (row.getNumElements() > 1) {
      matrix.appendRow(row);
      rowLb.push_back(-COIN_DBL_MAX);
      rowUb.push_back(1.0);
    }
  }

  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── TSP-MTZ LP relaxation ───────────────────────────────────────────────────
// Asymmetric TSP with Miller-Tucker-Zemlin subtour elimination.
// n cities, n²-n binary x[i][j] variables + n-1 continuous u[i] variables.
// LP relaxation bound is typically weak but the LP is non-trivial.
// Random distances in [1, 100].
static void buildTSP_MTZ(OsiSolverInterface &si, int n, unsigned seed)
{
  srand(seed);
  int nX = n * (n - 1); // x variables (exclude self-loops)
  int nU = n - 1;        // u variables (u[1]..u[n-1])
  int nCols = nX + nU;

  // Map (i,j) with i!=j to variable index
  auto xIdx = [&](int i, int j) -> int {
    assert(i != j);
    int idx = i * (n - 1) + (j > i ? j - 1 : j);
    return idx;
  };

  // Objective: minimize sum c[i][j] * x[i][j]
  std::vector<double> obj(nCols, 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (i != j)
        obj[xIdx(i, j)] = 1.0 + (rand() % 100);

  std::vector<double> colLb(nCols, 0.0);
  std::vector<double> colUb(nCols, 1.0);
  // u variables: u[i] ∈ [1, n-1]
  for (int i = 0; i < nU; i++) {
    colLb[nX + i] = 1.0;
    colUb[nX + i] = n - 1;
  }

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, nCols);
  std::vector<double> rowLb, rowUb;

  // Assignment constraints: sum_j x[i][j] = 1 for each i (out-degree)
  for (int i = 0; i < n; i++) {
    CoinPackedVector row;
    for (int j = 0; j < n; j++)
      if (i != j)
        row.insert(xIdx(i, j), 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);
    rowUb.push_back(1.0);
  }

  // Assignment constraints: sum_i x[i][j] = 1 for each j (in-degree)
  for (int j = 0; j < n; j++) {
    CoinPackedVector row;
    for (int i = 0; i < n; i++)
      if (i != j)
        row.insert(xIdx(i, j), 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);
    rowUb.push_back(1.0);
  }

  // MTZ subtour elimination: u[i] - u[j] + n*x[i][j] <= n-1
  // for all i,j ∈ {1,...,n-1}, i != j
  for (int i = 1; i < n; i++) {
    for (int j = 1; j < n; j++) {
      if (i != j) {
        CoinPackedVector row;
        row.insert(nX + (i - 1), 1.0);   // u[i]
        row.insert(nX + (j - 1), -1.0);  // -u[j]
        row.insert(xIdx(i, j), (double)n); // n*x[i][j]
        matrix.appendRow(row);
        rowLb.push_back(-COIN_DBL_MAX);
        rowUb.push_back(n - 1);
      }
    }
  }

  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── Infeasible LP ───────────────────────────────────────────────────────────
// x1 + x2 = 1, x1 + x2 = 2 — clearly infeasible.
// Extended to n variables to make it non-trivial for presolve.
static void buildInfeasible(OsiSolverInterface &si, int n)
{
  std::vector<double> obj(n, 1.0);
  std::vector<double> colLb(n, 0.0);
  std::vector<double> colUb(n, 10.0);

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, n);
  std::vector<double> rowLb, rowUb;

  // sum x[i] = 1
  {
    CoinPackedVector row;
    for (int i = 0; i < n; i++)
      row.insert(i, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);
    rowUb.push_back(1.0);
  }
  // sum x[i] = 2 (contradicts above)
  {
    CoinPackedVector row;
    for (int i = 0; i < n; i++)
      row.insert(i, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(2.0);
    rowUb.push_back(2.0);
  }
  // Add some filler constraints to make presolve work harder
  for (int i = 0; i < n - 1; i++) {
    CoinPackedVector row;
    row.insert(i, 1.0);
    row.insert(i + 1, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(-COIN_DBL_MAX);
    rowUb.push_back(5.0);
  }

  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── UFL (Uncapacitated Facility Location) LP relaxation ─────────────────────
// n facilities, m clients. Binary y[i] = facility i open, x[i][j] = client j
// served by facility i. Fixed cost f[i], service cost c[i][j].
// min sum_i f[i]*y[i] + sum_ij c[i][j]*x[i][j]
// s.t. sum_i x[i][j] = 1 for all j (each client served)
//      x[i][j] <= y[i] for all i,j (can only serve from open facility)
//      x[i][j] >= 0, y[i] in [0,1]
static void buildUFL(OsiSolverInterface &si, int n, int m, unsigned seed)
{
  srand(seed);
  int nCols = n + n * m; // y[0..n-1], x[0..n*m-1]

  std::vector<double> obj(nCols, 0.0);
  std::vector<double> colLb(nCols, 0.0);
  std::vector<double> colUb(nCols, 1.0);

  // Fixed costs f[i] in [50, 150]
  for (int i = 0; i < n; i++)
    obj[i] = 50.0 + (rand() % 101);

  // Service costs c[i][j] in [1, 50]
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      obj[n + i * m + j] = 1.0 + (rand() % 50);

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, nCols);
  std::vector<double> rowLb, rowUb;

  // Demand constraints: sum_i x[i][j] = 1 for each client j
  for (int j = 0; j < m; j++) {
    CoinPackedVector row;
    for (int i = 0; i < n; i++)
      row.insert(n + i * m + j, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);
    rowUb.push_back(1.0);
  }

  // Linking constraints: x[i][j] <= y[i]
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      CoinPackedVector row;
      row.insert(n + i * m + j, 1.0);
      row.insert(i, -1.0);
      matrix.appendRow(row);
      rowLb.push_back(-COIN_DBL_MAX);
      rowUb.push_back(0.0);
    }
  }

  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── Set Covering LP relaxation ──────────────────────────────────────────────
// min sum c[j]*x[j], s.t. Ax >= 1 (each row covered), x[j] in [0,1].
// Random 0-1 matrix with ~density nonzeros per row.
static void buildSetCovering(OsiSolverInterface &si, int nRows, int nCols,
  double density, unsigned seed)
{
  srand(seed);
  std::vector<double> obj(nCols);
  std::vector<double> colLb(nCols, 0.0);
  std::vector<double> colUb(nCols, 1.0);
  for (int j = 0; j < nCols; j++)
    obj[j] = 1.0 + (rand() % 20);

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, nCols);
  std::vector<double> rowLb, rowUb;

  for (int i = 0; i < nRows; i++) {
    CoinPackedVector row;
    for (int j = 0; j < nCols; j++)
      if ((rand() % 1000) < (int)(density * 1000))
        row.insert(j, 1.0);
    if (row.getNumElements() == 0)
      row.insert(rand() % nCols, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);        // >= 1 (covering)
    rowUb.push_back(COIN_DBL_MAX);
  }
  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── Set Packing LP relaxation ───────────────────────────────────────────────
// max sum c[j]*x[j], s.t. Ax <= 1 (each row at most one), x[j] in [0,1].
static void buildSetPacking(OsiSolverInterface &si, int nRows, int nCols,
  double density, unsigned seed)
{
  srand(seed);
  std::vector<double> obj(nCols);
  std::vector<double> colLb(nCols, 0.0);
  std::vector<double> colUb(nCols, 1.0);
  for (int j = 0; j < nCols; j++)
    obj[j] = -(1.0 + (rand() % 20)); // negate for maximization via min

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, nCols);
  std::vector<double> rowLb, rowUb;

  for (int i = 0; i < nRows; i++) {
    CoinPackedVector row;
    for (int j = 0; j < nCols; j++)
      if ((rand() % 1000) < (int)(density * 1000))
        row.insert(j, 1.0);
    if (row.getNumElements() == 0)
      row.insert(rand() % nCols, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(-COIN_DBL_MAX);
    rowUb.push_back(1.0);        // <= 1 (packing)
  }
  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── Set Partitioning LP relaxation ──────────────────────────────────────────
// min sum c[j]*x[j], s.t. Ax = 1 (each row exactly one), x[j] in [0,1].
static void buildSetPartitioning(OsiSolverInterface &si, int nRows, int nCols,
  double density, unsigned seed)
{
  srand(seed);
  std::vector<double> obj(nCols);
  std::vector<double> colLb(nCols, 0.0);
  std::vector<double> colUb(nCols, 1.0);
  for (int j = 0; j < nCols; j++)
    obj[j] = 1.0 + (rand() % 20);

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, nCols);
  std::vector<double> rowLb, rowUb;

  // Ensure feasibility: each row must have at least one column covering it
  for (int i = 0; i < nRows; i++) {
    CoinPackedVector row;
    for (int j = 0; j < nCols; j++)
      if ((rand() % 1000) < (int)(density * 1000))
        row.insert(j, 1.0);
    if (row.getNumElements() == 0)
      row.insert(rand() % nCols, 1.0);
    matrix.appendRow(row);
    rowLb.push_back(1.0);        // = 1 (partitioning)
    rowUb.push_back(1.0);
  }
  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── Unbounded LP ────────────────────────────────────────────────────────────
// min -x1 - x2 - ... - xn, x[i] >= 0, no upper bounds, one trivial constraint.
static void buildUnbounded(OsiSolverInterface &si, int n)
{
  std::vector<double> obj(n, -1.0);
  std::vector<double> colLb(n, 0.0);
  std::vector<double> colUb(n, COIN_DBL_MAX);

  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, n);
  std::vector<double> rowLb, rowUb;

  // x[0] - x[1] <= 10 (doesn't bound the objective)
  CoinPackedVector row;
  row.insert(0, 1.0);
  row.insert(1, -1.0);
  matrix.appendRow(row);
  rowLb.push_back(-COIN_DBL_MAX);
  rowUb.push_back(10.0);

  si.loadProblem(matrix, colLb.data(), colUb.data(), obj.data(),
    rowLb.data(), rowUb.data());
}

// ─── Test runner ─────────────────────────────────────────────────────────────

struct TestResult {
  bool optimal;
  bool infeasible;
  bool unbounded;
  double objValue;
  int iterations;
};

static TestResult solveLP(OsiSolverInterface &si, const char *label = "")
{
  si.initialSolve();
  TestResult r;
  r.optimal = si.isProvenOptimal();
  r.infeasible = si.isProvenPrimalInfeasible();
  r.unbounded = si.isProvenDualInfeasible();
  r.objValue = si.getObjValue();
  r.iterations = si.getIterationCount();
  // Verify optimality conditions when optimal
  if (r.optimal) {
    int v = checkOptimalityConditions(si, label);
    if (v > 0)
      fprintf(stderr, "  %s: %d optimality violations!\n", label, v);
    r.optimal = (v == 0); // treat violations as non-optimal
  }
  return r;
}

static void runTest(const char *name, int n,
  void (*builder)(OsiSolverInterface &, int),
  bool expectOptimal, bool expectInfeasible, bool expectUnbounded,
  double expectedObj = 0.0, double objTol = 1e-4)
{
  printf("\n--- %s (n=%d) ---\n", name, n);

  // Normal solve
  OsiClpSolverInterface si1;
  si1.messageHandler()->setLogLevel(0);
  builder(si1, n);
  TestResult r1 = solveLP(si1, "Normal");

  printf("  Normal:  opt=%d inf=%d unb=%d obj=%.6f iters=%d\n",
    r1.optimal, r1.infeasible, r1.unbounded, r1.objValue, r1.iterations);

  CHECK(r1.optimal == expectOptimal, "Normal: optimal status mismatch");
  CHECK(r1.infeasible == expectInfeasible, "Normal: infeasible status mismatch");
  CHECK(r1.unbounded == expectUnbounded, "Normal: unbounded status mismatch");
  if (expectOptimal)
    CHECK(fabs(r1.objValue - expectedObj) < objTol,
      "Normal: objective value mismatch");

  // Racing solve (2 threads)
  OsiClpSolverInterface si2;
  si2.messageHandler()->setLogLevel(0);
  builder(si2, n);
  si2.setRacingLP(2);
  TestResult r2 = solveLP(si2, "Racing2");

  printf("  Racing2: opt=%d inf=%d unb=%d obj=%.6f iters=%d\n",
    r2.optimal, r2.infeasible, r2.unbounded, r2.objValue, r2.iterations);

  CHECK(r2.optimal == expectOptimal, "Racing: optimal status mismatch");
  CHECK(r2.infeasible == expectInfeasible, "Racing: infeasible status mismatch");
  CHECK(r2.unbounded == expectUnbounded, "Racing: unbounded status mismatch");
  if (expectOptimal)
    CHECK(fabs(r2.objValue - expectedObj) < objTol,
      "Racing: objective value mismatch");

  // Cross-check: if both optimal, objectives must match
  if (r1.optimal && r2.optimal)
    CHECK(fabs(r1.objValue - r2.objValue) < objTol,
      "Normal vs Racing: objective mismatch");
}

// Wrapper for TSP-MTZ (needs seed)
static void buildTSP20(OsiSolverInterface &si, int) { buildTSP_MTZ(si, 20, 42); }
static void buildTSP40(OsiSolverInterface &si, int) { buildTSP_MTZ(si, 40, 123); }
static void buildTSP60(OsiSolverInterface &si, int) { buildTSP_MTZ(si, 60, 7); }

int main()
{
  printf("=== OsiClpSolverInterface LP Racing Tests ===\n");

  // N-Queens: LP relaxation value = -n (minimizing -sum)
  runTest("NQueens-8", 8, buildNQueens, true, false, false, -8.0);
  runTest("NQueens-20", 20, buildNQueens, true, false, false, -20.0);
  runTest("NQueens-50", 50, buildNQueens, true, false, false, -50.0);
  runTest("NQueens-100", 100, buildNQueens, true, false, false, -100.0);

  // TSP-MTZ: LP relaxation is feasible, obj > 0
  // We don't know exact obj, just check feasibility and consistency
  {
    printf("\n--- TSP-MTZ-20 ---\n");
    OsiClpSolverInterface si1, si2;
    si1.messageHandler()->setLogLevel(0);
    si2.messageHandler()->setLogLevel(0);
    buildTSP20(si1, 0);
    buildTSP20(si2, 0);
    si2.setRacingLP(2);
    TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
    printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
    printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
    CHECK(r1.optimal, "TSP20 Normal: should be optimal");
    CHECK(r2.optimal, "TSP20 Racing: should be optimal");
    CHECK(fabs(r1.objValue - r2.objValue) < kTol, "TSP20: obj mismatch");
  }
  {
    printf("\n--- TSP-MTZ-40 ---\n");
    OsiClpSolverInterface si1, si2;
    si1.messageHandler()->setLogLevel(0);
    si2.messageHandler()->setLogLevel(0);
    buildTSP40(si1, 0);
    buildTSP40(si2, 0);
    si2.setRacingLP(2);
    TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
    printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
    printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
    CHECK(r1.optimal, "TSP40 Normal: should be optimal");
    CHECK(r2.optimal, "TSP40 Racing: should be optimal");
    CHECK(fabs(r1.objValue - r2.objValue) < kTol, "TSP40: obj mismatch");
  }
  {
    printf("\n--- TSP-MTZ-60 ---\n");
    OsiClpSolverInterface si1, si2;
    si1.messageHandler()->setLogLevel(0);
    si2.messageHandler()->setLogLevel(0);
    buildTSP60(si1, 0);
    buildTSP60(si2, 0);
    si2.setRacingLP(2);
    TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
    printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
    printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
    CHECK(r1.optimal, "TSP60 Normal: should be optimal");
    CHECK(r2.optimal, "TSP60 Racing: should be optimal");
    CHECK(fabs(r1.objValue - r2.objValue) < kTol, "TSP60: obj mismatch");
  }

  // UFL (Uncapacitated Facility Location)
  {
    struct { int n; int m; unsigned seed; const char *name; } uflCases[] = {
      {10, 30, 99, "UFL-10x30"},
      {20, 60, 77, "UFL-20x60"},
      {30, 100, 55, "UFL-30x100"},
      {50, 200, 33, "UFL-50x200"},
    };
    for (auto &c : uflCases) {
      printf("\n--- %s ---\n", c.name);
      OsiClpSolverInterface si1, si2;
      si1.messageHandler()->setLogLevel(0);
      si2.messageHandler()->setLogLevel(0);
      buildUFL(si1, c.n, c.m, c.seed);
      buildUFL(si2, c.n, c.m, c.seed);
      si2.setRacingLP(2);
      TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
      printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
      printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
      CHECK(r1.optimal, (std::string(c.name) + " Normal: should be optimal").c_str());
      CHECK(r2.optimal, (std::string(c.name) + " Racing: should be optimal").c_str());
      CHECK(fabs(r1.objValue - r2.objValue) < kTol,
        (std::string(c.name) + ": obj mismatch").c_str());
    }
  }

  // Infeasible problems
  runTest("Infeasible-10", 10, buildInfeasible, false, true, false);
  runTest("Infeasible-50", 50, buildInfeasible, false, true, false);
  runTest("Infeasible-200", 200, buildInfeasible, false, true, false);

  // Unbounded problems
  runTest("Unbounded-10", 10, buildUnbounded, false, false, true);
  runTest("Unbounded-50", 50, buildUnbounded, false, false, true);

  // Set Covering (min cost to cover all rows)
  {
    struct { int m; int n; double d; unsigned seed; const char *name; } cases[] = {
      {30, 100, 0.15, 11, "SetCover-30x100"},
      {50, 200, 0.10, 22, "SetCover-50x200"},
      {100, 500, 0.08, 33, "SetCover-100x500"},
    };
    for (auto &c : cases) {
      printf("\n--- %s ---\n", c.name);
      OsiClpSolverInterface si1, si2;
      si1.messageHandler()->setLogLevel(0);
      si2.messageHandler()->setLogLevel(0);
      buildSetCovering(si1, c.m, c.n, c.d, c.seed);
      buildSetCovering(si2, c.m, c.n, c.d, c.seed);
      si2.setRacingLP(2);
      TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
      printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
      printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
      CHECK(r1.optimal, (std::string(c.name) + " Normal: optimal").c_str());
      CHECK(r2.optimal, (std::string(c.name) + " Racing: optimal").c_str());
      CHECK(fabs(r1.objValue - r2.objValue) < kTol,
        (std::string(c.name) + ": obj mismatch").c_str());
    }
  }

  // Set Packing (max weight independent set relaxation)
  {
    struct { int m; int n; double d; unsigned seed; const char *name; } cases[] = {
      {40, 120, 0.12, 44, "SetPack-40x120"},
      {80, 300, 0.08, 55, "SetPack-80x300"},
      {150, 600, 0.06, 66, "SetPack-150x600"},
    };
    for (auto &c : cases) {
      printf("\n--- %s ---\n", c.name);
      OsiClpSolverInterface si1, si2;
      si1.messageHandler()->setLogLevel(0);
      si2.messageHandler()->setLogLevel(0);
      buildSetPacking(si1, c.m, c.n, c.d, c.seed);
      buildSetPacking(si2, c.m, c.n, c.d, c.seed);
      si2.setRacingLP(2);
      TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
      printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
      printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
      CHECK(r1.optimal, (std::string(c.name) + " Normal: optimal").c_str());
      CHECK(r2.optimal, (std::string(c.name) + " Racing: optimal").c_str());
      CHECK(fabs(r1.objValue - r2.objValue) < kTol,
        (std::string(c.name) + ": obj mismatch").c_str());
    }
  }

  // Set Partitioning (exact cover relaxation)
  {
    struct { int m; int n; double d; unsigned seed; const char *name; } cases[] = {
      {20, 80, 0.20, 77, "SetPart-20x80"},
      {40, 200, 0.12, 88, "SetPart-40x200"},
      {60, 400, 0.08, 99, "SetPart-60x400"},
    };
    for (auto &c : cases) {
      printf("\n--- %s ---\n", c.name);
      OsiClpSolverInterface si1, si2;
      si1.messageHandler()->setLogLevel(0);
      si2.messageHandler()->setLogLevel(0);
      buildSetPartitioning(si1, c.m, c.n, c.d, c.seed);
      buildSetPartitioning(si2, c.m, c.n, c.d, c.seed);
      si2.setRacingLP(2);
      TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
      printf("  Normal:  opt=%d obj=%.4f iters=%d\n", r1.optimal, r1.objValue, r1.iterations);
      printf("  Racing2: opt=%d obj=%.4f iters=%d\n", r2.optimal, r2.objValue, r2.iterations);
      CHECK(r1.optimal, (std::string(c.name) + " Normal: optimal").c_str());
      CHECK(r2.optimal, (std::string(c.name) + " Racing: optimal").c_str());
      CHECK(fabs(r1.objValue - r2.objValue) < kTol,
        (std::string(c.name) + ": obj mismatch").c_str());
    }
  }

  // Racing=1 (single-thread fallback — no threading overhead)
  {
    printf("\n--- NQueens-50 Racing=1 (single-thread) ---\n");
    OsiClpSolverInterface si1, si2;
    si1.messageHandler()->setLogLevel(0);
    si2.messageHandler()->setLogLevel(0);
    buildNQueens(si1, 50);
    buildNQueens(si2, 50);
    si2.setRacingLP(1);
    TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
    printf("  Normal:  opt=%d obj=%.4f\n", r1.optimal, r1.objValue);
    printf("  Racing1: opt=%d obj=%.4f\n", r2.optimal, r2.objValue);
    CHECK(r1.optimal && r2.optimal, "NQueens50 Racing1: both optimal");
    CHECK(fabs(r1.objValue - r2.objValue) < kTol, "NQueens50 Racing1: obj match");
  }

  // Racing=3 (full portfolio)
  {
    printf("\n--- TSP-MTZ-40 Racing=3 ---\n");
    OsiClpSolverInterface si1, si2;
    si1.messageHandler()->setLogLevel(0);
    si2.messageHandler()->setLogLevel(0);
    buildTSP40(si1, 0);
    buildTSP40(si2, 0);
    si2.setRacingLP(3);
    TestResult r1 = solveLP(si1, "Normal"), r2 = solveLP(si2, "Racing");
    printf("  Normal:  opt=%d obj=%.4f\n", r1.optimal, r1.objValue);
    printf("  Racing3: opt=%d obj=%.4f\n", r2.optimal, r2.objValue);
    CHECK(r1.optimal && r2.optimal, "TSP40 Racing3: both optimal");
    CHECK(fabs(r1.objValue - r2.objValue) < kTol, "TSP40 Racing3: obj match");
  }

  // Infeasible with 3 threads
  {
    printf("\n--- Infeasible-50 Racing=3 ---\n");
    OsiClpSolverInterface si;
    si.messageHandler()->setLogLevel(0);
    buildInfeasible(si, 50);
    si.setRacingLP(3);
    TestResult r = solveLP(si, "Racing3");
    printf("  Racing3: opt=%d inf=%d\n", r.optimal, r.infeasible);
    CHECK(r.infeasible, "Infeasible50 Racing3: detected");
    CHECK(!r.optimal, "Infeasible50 Racing3: not optimal");
  }

  // Unbounded with 3 threads
  {
    printf("\n--- Unbounded-50 Racing=3 ---\n");
    OsiClpSolverInterface si;
    si.messageHandler()->setLogLevel(0);
    buildUnbounded(si, 50);
    si.setRacingLP(3);
    TestResult r = solveLP(si, "Racing3");
    printf("  Racing3: opt=%d unb=%d\n", r.optimal, r.unbounded);
    CHECK(r.unbounded, "Unbounded50 Racing3: detected");
    CHECK(!r.optimal, "Unbounded50 Racing3: not optimal");
  }

  printf("\n=== Results: %d passed, %d failed ===\n", nPass, nFail);
  return nFail > 0 ? 1 : 0;
}
