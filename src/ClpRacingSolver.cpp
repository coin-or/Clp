// Copyright (C) 2026, Haroldo Gambini Santos and others.
// All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "ClpRacingSolver.hpp"
#include "ClpEventHandler.hpp"
#include "ClpOutput.hpp"
#include "CoinTime.hpp"
#include <atomic>
#include <memory>
#include <vector>

// On Windows (MSVC) we use std::thread/std::mutex, which is safe there.
// On Linux/Unix we use POSIX pthreads directly to avoid embedding GCC's
// std::thread weak symbols (_M_thread_deps_never_run, _State_impl vtable)
// into static binaries, which can cause at-exit crashes on certain
// glibc/libstdc++ versions. (std::mutex itself -- used by ClpLpPhaseState,
// see ClpOutput.hpp -- does not carry the same risk; only std::thread does.)
#ifdef _WIN32
#include <mutex>
#include <thread>
#define CLP_RACING_USE_STD_THREAD 1
#else
#include <dlfcn.h>
#include <pthread.h>
#define CLP_RACING_USE_STD_THREAD 0
#endif

// Restrict OpenBLAS to 1 thread per racing thread to avoid thread explosion.
#if !defined(_WIN32)
namespace {
inline void racing_set_openblas_threads(int n)
{
  typedef void (*fn_t)(int);
  static fn_t fn = reinterpret_cast<fn_t>(dlsym(RTLD_DEFAULT, "openblas_set_num_threads"));
  if (fn)
    fn(n);
}
} // namespace
#elif defined(CLP_USE_OPENBLAS)
extern "C" void openblas_set_num_threads(int num_threads);
namespace {
inline void racing_set_openblas_threads(int n) { openblas_set_num_threads(n); }
} // namespace
#else
namespace {
inline void racing_set_openblas_threads(int) {}
} // namespace
#endif

// ─── Racing progress handler ─────────────────────────────────────────────────
// Feeds each racing config's progress into the *same* unified LP-progress
// table (ClpLpTable, see ClpOutput.hpp) that a normal, non-racing LP solve
// uses -- instead of a separate ad hoc table. ClpLpTable::printRow() is
// thread-safe and rate-limits printing *across all racing configs combined*
// via the shared ClpLpPhaseState's timeFreq/lastPrintTime, so only one row
// (from whichever config happens to report once the interval has elapsed)
// is printed per interval -- keeping racing's output just as compact as the
// sequential table, per-config chatter is not multiplied by nThreads.
namespace {

class RacingEventHandler : public ClpEventHandler {
public:
  RacingEventHandler(std::shared_ptr<ClpLpPhaseState> state,
    const std::atomic<bool> *abortFlag, const char *label)
    : ClpEventHandler()
    , state_(state)
    , abortFlag_(abortFlag)
    , label_(label) {}

  RacingEventHandler(const RacingEventHandler &rhs)
    : ClpEventHandler(rhs)
    , state_(rhs.state_)
    , abortFlag_(rhs.abortFlag_)
    , label_(rhs.label_) {}

  ClpEventHandler *clone() const override
  {
    return new RacingEventHandler(*this);
  }

  int event(Event whichEvent) override
  {
    // Check abort
    if (abortFlag_ && abortFlag_->load(std::memory_order_relaxed))
      return 0; // stop

    if (whichEvent != endOfIteration || !model_ || !state_)
      return -1;

    ClpLpTable::printRow(*state_, label_, model_->numberIterations(),
      model_->objectiveValue(), model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities());
    return -1;
  }

private:
  std::shared_ptr<ClpLpPhaseState> state_;
  const std::atomic<bool> *abortFlag_;
  const char *label_;
};

} // namespace

ClpRacingSolver::ClpRacingSolver(ClpSimplex *model, int numThreads)
  : model_(model)
  , numThreads_(numThreads)
{
}

void ClpRacingSolver::addConfig(const ClpSolve &config, ConfigSetupFn setupFn)
{
  configs_.push_back(config);
  setupFns_.push_back(std::move(setupFn));
}

void ClpRacingSolver::addDefaultConfigs(int portfolioSize)
{
  configs_.clear();
  setupFns_.clear();

  int k = portfolioSize > 0 ? portfolioSize : numThreads_;

  // ── Shared config builder helpers ─────────────────────────────────────────

  // dual_pertv72: plain dual simplex + perturbation value 72.
  // Mirrors "-pertValue 72 -lpMethod dual" (no PE steepest pivot).
  // Selected as optimal dual component via exhaustive k-fold search over 12
  // cluster reps (C(12,2)=66 pairs, C(12,3)=220 triples, 10-fold CV).
  auto makeDualPertv72 = []() -> std::pair<ClpSolve, ConfigSetupFn> {
    ClpSolve opts;
    opts.setSolveType(ClpSolve::useDual);
    opts.setPresolveType(ClpSolve::presolveOn);
    opts.setSpecialOption(2, 1);
    ConfigSetupFn fn = [](ClpSimplex *m) {
      m->setPerturbation(72);
    };
    return {opts, fn};
  };

  // primal_idiot50: primal simplex with 50 idiot-crash passes
  auto makePrimalIdiot50 = []() -> std::pair<ClpSolve, ConfigSetupFn> {
    ClpSolve opts;
    opts.setSolveType(ClpSolve::usePrimal);
    opts.setPresolveType(ClpSolve::presolveOn);
    opts.setSpecialOption(1, 2, 50); // idiot, 50 passes
    opts.setSpecialOption(2, 1);
    return {opts, nullptr};
  };

  // primal_sprint: primal simplex with sprint
  auto makePrimalSprint = []() -> std::pair<ClpSolve, ConfigSetupFn> {
    ClpSolve opts;
    opts.setSolveType(ClpSolve::usePrimalorSprint);
    opts.setPresolveType(ClpSolve::presolveOn);
    opts.setSpecialOption(1, 3);
    opts.setSpecialOption(2, 1);
    return {opts, nullptr};
  };

  if (k == 2) {
    // K=2 optimal portfolio (exhaustive k-fold search, 1.47x speedup vs baseline):
    //   dual_pertv72 + primal_idiot50
    auto [o0, f0] = makeDualPertv72();
    auto [o1, f1] = makePrimalIdiot50();
    addConfig(o0, f0);
    addConfig(o1, f1);
  } else {
    // K=3 optimal portfolio (exhaustive k-fold search, 1.59x speedup vs baseline):
    //   dual_pertv72 + primal_idiot50 + primal_sprint
    auto [o0, f0] = makeDualPertv72();
    auto [o1, f1] = makePrimalIdiot50();
    auto [o2, f2] = makePrimalSprint();
    addConfig(o0, f0);
    addConfig(o1, f1);
    addConfig(o2, f2);
  }
}

const char *ClpRacingSolver::winnerName() const
{
  static const char *names[] = { "dual", "primal+idiot", "primal+sprint" };
  if (winnerIndex_ < 0)
    return "";
  return (winnerIndex_ < 3) ? names[winnerIndex_] : "unknown";
}

int ClpRacingSolver::solve()
{
  if (configs_.empty())
    addDefaultConfigs(numThreads_);

  int nConfigs = static_cast<int>(configs_.size());
  int nThreads = numThreads_ > 0 ? numThreads_ : nConfigs;
  if (nThreads > nConfigs)
    nThreads = nConfigs;

  // If only one config, just solve directly (no threading overhead)
  if (nThreads <= 1) {
    model_->initialSolve(configs_[0]);
    if (model_->status() == 0) {
      winnerIndex_ = 0;
      winnerIterations_ = model_->numberIterations();
    }
    return winnerIndex_;
  }

  std::atomic<bool> abortFlag{false};
  std::atomic<int> winner{-1};

  // Clone models for each racing thread; apply per-config setup functions.
  std::vector<ClpSimplex *> clones(nThreads, nullptr);
  for (int i = 0; i < nThreads; i++) {
    clones[i] = new ClpSimplex(*model_);
    if (i < static_cast<int>(setupFns_.size()) && setupFns_[i])
      setupFns_[i](clones[i]);
  }

  // Set up shared progress reporting. Prefer reusing an already-installed
  // ClpLpEventHandler's ClpLpPhaseState (installed by the caller -- e.g.
  // CbcSolver::solveInitialLp() -- for the unified LP progress table) so
  // racing rows land in the exact same Phase/Iter/Objective/Primal inf/
  // Dual inf/Time table as a normal LP solve, and so the caller's own
  // printFinalStatus() call (made after solve() returns) closes the table
  // and prints the final summary line as usual. Falls back to a private
  // table state (same format) when no handler is installed.
  static const char *configLabels[] = {"Dual", "P+Idiot", "Sprint"};
  ClpLpEventHandler *existingHandler
    = dynamic_cast<ClpLpEventHandler *>(model_->eventHandler());
  std::shared_ptr<ClpLpPhaseState> tableState;
  bool ownTableState = false;
  if (existingHandler) {
    tableState = existingHandler->sharedState();
  } else if (model_->logLevel() > 0) {
    tableState = std::make_shared<ClpLpPhaseState>();
    tableState->fp = model_->messageHandler()
      ? model_->messageHandler()->filePointer() : stdout;
    tableState->utf8 = ClpOutput::useUtf8();
    tableState->compact = ClpOutput::useCompact();
    tableState->logLevel = model_->logLevel();
    tableState->timeFreq = 2.0;
    // Shift the local wall-clock reference back by however much overall
    // search time had already elapsed before racing began, so every
    // "now - startTime" computation yields time elapsed since the
    // *overall search* began, not just since this LP race started
    // (matches the same fix applied to the sequential root LP relaxation
    // table in CbcSolver::solveInitialLp()).
    tableState->startTime = CoinWallclockTime() - searchElapsedAtStart_;
    tableState->lastPrintTime = tableState->startTime;
    tableState->title = "LP solve";
    ownTableState = true;
  }

  double startTime = CoinGetTimeOfDay();

#if CLP_RACING_USE_STD_THREAD
  std::vector<std::thread> threads;
  threads.reserve(nThreads);
  for (int i = 0; i < nThreads; i++) {
    threads.emplace_back([i, &clones, &abortFlag, &winner, &tableState, this]() {
      ClpSimplex *clone = clones[i];
      clone->setLogLevel(0);
      racing_set_openblas_threads(1);
      const char *label = (i < 3) ? configLabels[i] : "Config";
      RacingEventHandler handler(tableState, &abortFlag, label);
      clone->passInEventHandler(&handler);
      clone->initialSolve(configs_[i]);
      int st = clone->status();
      if (st == 0 || st == 1 || st == 2) {
        int expected = -1;
        if (winner.compare_exchange_strong(expected, i))
          abortFlag.store(true, std::memory_order_relaxed);
      }
    });
  }
  for (auto &t : threads)
    t.join();
#else
  // Per-thread argument struct (passed to the pthread callback)
  struct ThreadArg {
    ClpSimplex *clone;
    ClpSolve *config;
    std::atomic<bool> *abortFlag;
    std::atomic<int> *winner;
    std::shared_ptr<ClpLpPhaseState> *tableState;
    const char *label;
    int index;
  };

  std::vector<ThreadArg> args(nThreads);
  std::vector<pthread_t> threads(nThreads);

  for (int i = 0; i < nThreads; i++) {
    args[i].clone = clones[i];
    args[i].config = &configs_[i];
    args[i].abortFlag = &abortFlag;
    args[i].winner = &winner;
    args[i].tableState = &tableState;
    args[i].label = (i < 3) ? configLabels[i] : "Config";
    args[i].index = i;

    pthread_create(&threads[i], nullptr, [](void *arg) -> void * {
      ThreadArg *a = static_cast<ThreadArg *>(arg);
      ClpSimplex *clone = a->clone;
      clone->setLogLevel(0);
      racing_set_openblas_threads(1);
      RacingEventHandler handler(*a->tableState, a->abortFlag, a->label);
      clone->passInEventHandler(&handler);
      clone->initialSolve(*a->config);
      int st = clone->status();
      if (st == 0 || st == 1 || st == 2) {
        int expected = -1;
        if (a->winner->compare_exchange_strong(expected, a->index))
          a->abortFlag->store(true, std::memory_order_relaxed);
      }
      return nullptr;
    }, &args[i]);
  }

  for (int i = 0; i < nThreads; i++)
    pthread_join(threads[i], nullptr);
#endif

  double endTime = CoinGetTimeOfDay();
  winnerIndex_ = winner.load();

  if (winnerIndex_ >= 0) {
    ClpSimplex *w = clones[winnerIndex_];
    winnerTime_ = endTime - startTime;
    winnerIterations_ = w->numberIterations();

    // Copy solution back to original model
    int nCols = model_->numberColumns();
    int nRows = model_->numberRows();
    CoinMemcpyN(w->primalColumnSolution(), nCols,
      model_->primalColumnSolution());
    CoinMemcpyN(w->dualColumnSolution(), nCols,
      model_->dualColumnSolution());
    CoinMemcpyN(w->primalRowSolution(), nRows,
      model_->primalRowSolution());
    CoinMemcpyN(w->dualRowSolution(), nRows,
      model_->dualRowSolution());
    CoinMemcpyN(w->statusArray(), nCols + nRows,
      model_->statusArray());
    model_->setObjectiveValue(w->objectiveValue());
    model_->setProblemStatus(w->status());
    model_->setNumberIterations(winnerIterations_);
    model_->setSecondaryStatus(w->secondaryStatus());

    if (tableState) {
      tableState->racingWinner = winnerName();
      if (existingHandler) {
        // Point the caller's ClpLpEventHandler at the original model (now
        // holding the winner's solution/status/iterations) so its later
        // printFinalStatus() call -- made by the caller after solve()
        // returns -- closes the table and prints the final summary line
        // exactly as it would for a non-racing solve.
        existingHandler->setSimplex(model_);
      } else if (ownTableState) {
        // No external handler was installed (e.g. a standalone racing
        // caller with no unified progress table already set up) -- close
        // out the table ourselves using the same shared formatting code.
        ClpLpTable::printFinalStatus(*tableState, model_);
      }
    }
  }

  // Clean up clones
  for (int i = 0; i < nThreads; i++)
    delete clones[i];

  return winnerIndex_;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
