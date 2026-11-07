// Copyright (C) 2024, COIN-OR Foundation
// All Rights Reserved. This code is published under the Eclipse Public License.

#ifndef ClpOutput_HPP
#define ClpOutput_HPP

#include "ClpConfig.h"
#include "ClpSimplex.hpp"
#include "ClpEventHandler.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinTime.hpp"

#include <memory>
#include <mutex>
#include <string>

// ─── ClpProgressEventHandler ──────────────────────────────────────────────────
// Event handler for LP progress output in the Clp standalone solver.
// Installed on the ClpSimplex model before solve and removed afterwards.
// Intercepts endOfIteration to print tabular progress, replacing Clp's built-in
// CLP_SIMPLEX_STATUS chatter.
//
// Presolve support: pass origRows/origCols/origNZ from the original model.
// When the first iteration fires in a smaller (presolved) model, the header
// shows the presolve reduction.
// ──────────────────────────────────────────────────────────────────────────────
class CLPLIB_EXPORT ClpProgressEventHandler : public ClpEventHandler {
public:
  ClpProgressEventHandler(CoinMessageHandler *handler, int logLevel,
    int iterFreq, double timeFreq,
    int origRows, int origCols, int origNZ);

  virtual int event(Event whichEvent) override;
  virtual ClpEventHandler *clone() const override;

  void printFinalStatus() const;

private:
  struct SharedState {
    bool headerPrinted = false;
    double startTime = 0.0;
    double lastPrintTime = 0.0;
    int lastPrintIter = 0;
    int maxIterSeen = 0;
  };
  std::shared_ptr<SharedState> shared_;

  CoinMessageHandler *handler_;
  int logLevel_;
  int iterFreq_;
  double timeFreq_;
  int origRows_;
  int origCols_;
  int origNZ_;

  void printHeader() const;
  void printRow(int iter, double obj, double pInf, double dInf, double elapsed, int algo) const;
};

// ─── ClpLpPhaseState ──────────────────────────────────────────────────────────
// Shared mutable state for ClpLpEventHandler and ClpLpMsgHandler, and for
// ClpLpTable::printRow() callers. All of them hold a std::shared_ptr to the
// same instance.
//
// The table is printed as messages (CLP_TABLE_*), so it shows only when the
// messages are in the en_tabular language (ClpModel::newLanguage()); in us_en
// they print nothing.
// ──────────────────────────────────────────────────────────────────────────────
struct ClpLpPhaseState {
  // ── Table bookkeeping ─────────────────────────────────────────────────────
  bool headerPrinted = false; // table header has been printed
  bool lpStarted = false;     // first LP endOfIteration has fired
  double startTime = 0.0;
  double lastPrintTime = 0.0;
  int lastPrintIter = 0;
  int maxIterSeen = 0;

  bool idiotSeen = false;  // first Idiot row printed
  bool sprintSeen = false; // first Sprint row printed

  // ── Configuration (set at construction, read-only afterwards) ─────────────
  /** Where the table's messages go. If null, ClpLpEventHandler uses its
   *  model's handler and messages; ClpLpTable::printRow() callers without a
   *  model must set them. */
  CoinMessageHandler *handler = nullptr;
  CoinMessages *messages = nullptr;
  FILE *fp = nullptr;
  bool utf8 = false;
  bool compact = true;
  int iterFreq = 0;
  double timeFreq = 5.0;
  int logLevel = 1;
  /*
    Modify this as needed
    at present says if primal or dual
   */
  int modifyMsg = 0; 
  int origRows = 0, origCols = 0;
  std::string title = "LP solve";
  std::string racingWinner; // e.g. "primal+idiot", set by ClpRacingSolver on
                             // success so the closing summary line can report
                             // which config won (empty for a non-racing solve).

  // Serializes access when the state is fed from multiple threads (e.g. LP
  // racing configs reporting progress concurrently via ClpLpTable::printRow()
  // below). Single-threaded callers (the normal ClpLpEventHandler/
  // ClpLpMsgHandler event/print() methods) pay only the (uncontended) lock
  // cost, which is negligible next to an LP iteration.
  mutable std::mutex mu;
};

// ─── ClpLpTable ────────────────────────────────────────────────────────────────
// Free functions that let an external, possibly-multi-threaded progress
// source (currently: ClpRacingSolver) feed rows into the exact same unified
// Phase/Iter/Objective/Primal inf/Dual inf/Time table used by a normal
// (non-racing) LP solve, instead of maintaining a separate ad hoc table.
// ──────────────────────────────────────────────────────────────────────────────
namespace ClpLpTable {

/** Thread-safe: print (or, per the state's timeFreq, rate-limit) a progress
 *  row labelled `phaseLabel` (e.g. a racing config's name such as "Dual" or
 *  "P+Idiot"). Opens the table on the first call from any thread. Only one
 *  row is printed per timeFreq interval *across all callers combined*,
 *  which naturally keeps racing's output as compact as the sequential
 *  table (only the "current leader" reports, one row at a time). */
CLPLIB_EXPORT void printRow(ClpLpPhaseState &state, const char *phaseLabel,
  int iter, double obj, double pInf, double dInf,
  CoinMessageHandler *handler = nullptr, CoinMessages *messages = nullptr);

/** Print the final status line (and close the table if one was opened).
 *  Shared implementation used by both ClpLpEventHandler::printFinalStatus()
 *  and callers (e.g. ClpRacingSolver) that drive their own ClpLpPhaseState
 *  without an installed ClpLpEventHandler. */
CLPLIB_EXPORT void printFinalStatus(ClpLpPhaseState &state, ClpSimplex *model,
  int numInts = 0, int numFrac = 0);

} // namespace ClpLpTable

// ─── ClpLpEventHandler ────────────────────────────────────────────────────────
// ClpEventHandler that drives a unified LP+Idiot+Sprint progress table: it
// decides when to print a row and sends it as a message. Rows for Idiot and
// Sprint (endOfIdiotIteration, endOfSprintPass) show from log level 2.
//
// Install on a ClpSimplex model via passInEventHandler() (Clp clones it),
// with a ClpLpMsgHandler as the model's message handler, and switch the
// model to the tabular style with newLanguage(CoinMessages::en_tabular).
//
// Usage:
//   auto state = std::make_shared<ClpLpPhaseState>();
//   // fill state->fp, utf8, logLevel, iterFreq, timeFreq, title ...
//   ClpLpEventHandler evtH(state);
//   ClpLpMsgHandler   msgH(state);
//   bool msgOldDefault;
//   CoinMessageHandler *savedMsg = model->pushMessageHandler(&msgH, msgOldDefault);
//   model->passInEventHandler(&evtH);
//   model->newLanguage(CoinMessages::en_tabular);
//   auto *evtClone = dynamic_cast<ClpLpEventHandler *>(model->eventHandler());
//   ... model->initialSolve() ...
//   if (evtClone) evtClone->printFinalStatus(numInts, numFrac);
//   model->popMessageHandler(savedMsg, msgOldDefault);
//   ClpEventHandler defEvt; model->passInEventHandler(&defEvt);
//   model->newLanguage(CoinMessages::us_en);
// ──────────────────────────────────────────────────────────────────────────────
class CLPLIB_EXPORT ClpLpEventHandler : public ClpEventHandler {
public:
  explicit ClpLpEventHandler(std::shared_ptr<ClpLpPhaseState> state);

  virtual int event(Event whichEvent) override;
  int eventWithInfo(Event whichEvent, void *info) override;
  virtual ClpEventHandler *clone() const override;

  /// Print presolve stats
  void printPresolveStats();
  /** Print final status line + close the table.
   *  Call once after initialSolve() on the clone returned by eventHandler().
   *  numInts: total MIP integer variables (0 for pure LP).
   *  numFrac: fractional integer variables in LP optimal (0 for pure LP). */
  void printFinalStatus(int numInts = 0, int numFrac = 0);

  /// Whether the progress table has already been opened (Idiot/Sprint started)
  bool tableStarted() const { return s_ && (s_->idiotSeen || s_->sprintSeen || s_->lpStarted); }
  /// File pointer for output
  FILE *fp() const { return s_ ? s_->fp : nullptr; }
  /** Shared table state, so an external progress source (ClpRacingSolver)
   *  can feed rows into the exact same table via ClpLpTable::printRow()
   *  instead of maintaining a separate one. */
  std::shared_ptr<ClpLpPhaseState> sharedState() const { return s_; }
  inline void setModifyMsg(int value)
  {
    s_->modifyMsg = value;
  }
  inline void setTimeFreq(double value)
  {
    s_->timeFreq = value;
  }
private:
  std::shared_ptr<ClpLpPhaseState> s_;

  void openTable();
  void printLpRow(int iter, double obj, double pInf, double dInf, double elapsed);
  void printPassRow(const char *phase, bool &seen, const ClpPassProgress &progress);
};

// ─── ClpLpMsgHandler ──────────────────────────────────────────────────────────
// The message handler for the tabular LP output. It prints every message to
// the state's file at the state's log level, without the Clpnnnn prefix; the
// language of the messages decides what each one looks like, or whether it
// prints at all (en_tabular drops the messages that the table replaces).
// ──────────────────────────────────────────────────────────────────────────────
class CLPLIB_EXPORT ClpLpMsgHandler : public CoinMessageHandler {
public:
  explicit ClpLpMsgHandler(std::shared_ptr<ClpLpPhaseState> state);

  virtual CoinMessageHandler *clone() const override;

private:
  std::shared_ptr<ClpLpPhaseState> s_;
};

// ─── ClpOutput ────────────────────────────────────────────────────────────────
// Centralised clean output for the Clp standalone LP solver.
//
//  • printSolverHeader   — one-line "CLP x.y.z — COIN-OR LP Solver" + args
//  • printProblemSummary — compact problem description + coefficient-range table
//
// All output is written via the provided CoinMessageHandler's filePointer().
// ──────────────────────────────────────────────────────────────────────────────

class CLPLIB_EXPORT ClpOutput {
public:
  /// Print the one-line solver header and (if argc>1) an args line.
  static void printSolverHeader(CoinMessageHandler *handler, int logLevel,
    int argc = 0, const char *const *argv = nullptr);

  /// Print compact problem summary (dimensions + coefficient ranges).
  static void printProblemSummary(CoinMessageHandler *handler,
    const ClpSimplex &model, int logLevel);

  // ── UTF-8 toggle ──────────────────────────────────────────────────────────
  static void setUtf8(bool on) { useUtf8_ = on; }
  static bool useUtf8() { return useUtf8_; }

  // ── Compact table style toggle ────────────────────────────────────────────
  static void setCompact(bool on) { useCompact_ = on; }
  static bool useCompact() { return useCompact_; }

private:
  static bool useUtf8_;
  static bool useCompact_;

  // Helpers
  static void sendLine(CoinMessageHandler *handler, const std::string &line);
  static void sendBlank(CoinMessageHandler *handler);
  static std::string fmtNum(double v);

  struct RangeStat {
    double lo = 1e300, hi = -1e300;
    bool valid = false;
    void update(double v);
    std::string ratio() const;
  };
};

#endif // ClpOutput_HPP
