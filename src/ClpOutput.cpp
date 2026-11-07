// Copyright (C) 2024, COIN-OR Foundation
// All Rights Reserved. This code is published under the Eclipse Public License.

#include "ClpOutput.hpp"
#include "ClpInterior.hpp"
#include "ClpMessage.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <string>

#include "ClpConfig.h"
#include "CoinPackedMatrix.hpp"
#include "CoinTable.hpp"
#include "CoinTime.hpp"

// ─── static members ───────────────────────────────────────────────────────────
bool ClpOutput::useUtf8_ = true;
bool ClpOutput::useCompact_ = true;

// ─── LP progress table layout ─────────────────────────────────────────────────
//
//   Phase   |     Iter |       Objective |   Primal inf |    Dual inf |    Time
//   ─────────┬──────────┬─────────────────┬──────────────┬──────────────┬─────────
//   Dual    |     3123  |     3.709408e+01 |  1.7376e+02  |  2.8598e+14 |   5.0s
//
// Column widths:
static const int LP_W_PHASE = 8;
static const int LP_W_ITER  = 8;
static const int LP_W_OBJ   = 15;
static const int LP_W_PINF  = 12;
static const int LP_W_DINF  = 12;
static const int LP_W_TIME  = 8;

static CoinTable makeLpTable(bool utf8, bool compact)
{
  return CoinTable({
    { "Phase",      LP_W_PHASE, /*leftAlign=*/true },
    { "Iter",       LP_W_ITER  },
    { "Objective",  LP_W_OBJ   },
    { "Primal inf", LP_W_PINF  },
    { "Dual inf",   LP_W_DINF  },
    { "Time",       LP_W_TIME  },
  }, utf8, /*indent=*/2, compact);
}

static const char *clpAlgoShortName(int algo)
{
  switch (algo) {
  case -1: return "Dual";
  case 1: return "Primal";
  case 2: return "Barrier";
  default: return "LP";
  }
}

// ─── ClpProgressEventHandler ──────────────────────────────────────────────────

ClpProgressEventHandler::ClpProgressEventHandler(CoinMessageHandler *handler,
  int logLevel, int iterFreq, double timeFreq,
  int origRows, int origCols, int origNZ)
  : ClpEventHandler()
  , shared_(std::make_shared<SharedState>())
  , handler_(handler)
  , logLevel_(logLevel)
  , iterFreq_(iterFreq)
  , timeFreq_(timeFreq)
  , origRows_(origRows)
  , origCols_(origCols)
  , origNZ_(origNZ)
{
  shared_->startTime = CoinWallclockTime();
}

ClpEventHandler *ClpProgressEventHandler::clone() const
{
  // Copies the shared_ptr — all clones share the same SharedState.
  return new ClpProgressEventHandler(*this);
}

void ClpProgressEventHandler::printHeader() const
{
  if (logLevel_ <= 0 || !handler_ || !model_)
    return;
  FILE *fp = handler_->filePointer();
  if (!fp)
    return;

  const bool u8 = ClpOutput::useUtf8();
  const CoinTable tbl = makeLpTable(u8, ClpOutput::useCompact());
  int presRows = model_->numberRows();
  int presCols = model_->numberColumns();
  if (presRows < origRows_ || presCols < origCols_) {
    fprintf(fp, "  Presolve: %d rows (%+d),  %d cols (%+d)\n",
      presRows, presRows - origRows_,
      presCols, presCols - origCols_);
  }
  fprintf(fp, "\n");

  const char *alg = clpAlgoShortName(model_->algorithm());
  std::string title = std::string("LP solve (") + alg + " simplex)";
  fprintf(fp, "\n%s\n", CoinTable::phaseStart(title, u8).c_str()); // phase start
  fprintf(fp, "%s\n", tbl.sepLine(CoinTable::Top).c_str());       // ─────┬─────
  fprintf(fp, "%s\n", tbl.headerLine().c_str());                  // Phase│Iter│...
  fprintf(fp, "%s\n", tbl.sepLine(CoinTable::Middle).c_str());    // ─────┼─────
  fflush(fp);
}

void ClpProgressEventHandler::printRow(int iter, double obj, double pInf,
  double dInf, double elapsed, int algo) const
{
  if (!handler_)
    return;
  FILE *fp = handler_->filePointer();
  if (!fp)
    return;
  const bool u8 = ClpOutput::useUtf8();
  const bool compact = ClpOutput::useCompact();
  const char *bar = compact ? " " : (u8 ? " \xe2\x94\x82 " : " | ");
  char tbuf[16];
  std::snprintf(tbuf, sizeof(tbuf), "%.1fs", elapsed);
  fprintf(fp, "  %-*s%s%*d%s%*.6e%s%*.4e%s%*.4e%s%*s\n",
    LP_W_PHASE, clpAlgoShortName(algo), bar,
    LP_W_ITER, iter, bar,
    LP_W_OBJ, obj, bar,
    LP_W_PINF, (pInf > 0.0 ? pInf : 0.0), bar,
    LP_W_DINF, (dInf > 0.0 ? dInf : 0.0), bar,
    LP_W_TIME, tbuf);
  fflush(fp);
}

void ClpProgressEventHandler::printFinalStatus() const
{
  if (logLevel_ <= 0 || !handler_ || !model_)
    return;
  if (!shared_->headerPrinted)
    return;
  FILE *fp = handler_->filePointer();
  if (!fp)
    return;

  const bool u8 = ClpOutput::useUtf8();
  const CoinTable tbl = makeLpTable(u8, ClpOutput::useCompact());
  const double elapsed = CoinWallclockTime() - shared_->startTime;
  const int iters = (model_->numberIterations() > shared_->maxIterSeen)
    ? model_->numberIterations()
    : shared_->maxIterSeen;

  fprintf(fp, "%s\n", tbl.sepLine(CoinTable::Bottom).c_str());

  const char *status = "Unknown";
  const int s = model_->status();
  const int ss = model_->secondaryStatus();
  if (s == 0)
    status = "Optimal";
  else if (s == 1) {
    if (ss == 4)
      status = "Iteration limit";
    else if (ss == 9)
      status = "Time limit";
    else
      status = "Infeasible";
  } else if (s == 2)
    status = "Unbounded";
  else if (s == 3)
    status = "Stopped";
  else if (s == 4)
    status = "Numerical error";
  else if (s == 5)
    status = "Event";

  // Build phase-end summary line
  char summary[256];
  std::snprintf(summary, sizeof(summary),
    "%s%sObj: %g   Iters: %d   Time: %.2fs",
    status, CoinTable::dashSep(u8), model_->objectiveValue(), iters, elapsed);
  fprintf(fp, "%s\n", CoinTable::phaseEnd(summary, u8).c_str());
  fflush(fp);
}

int ClpProgressEventHandler::event(Event whichEvent)
{
  if (whichEvent != endOfIteration || !model_)
    return -1;
  if (logLevel_ <= 0 || !handler_)
    return -1;

  if (!shared_->headerPrinted) {
    printHeader();
    shared_->headerPrinted = true;
    const int iter = model_->numberIterations();
    const double elapsed = CoinWallclockTime() - shared_->startTime;
    printRow(iter, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), elapsed, model_->algorithm());
    shared_->lastPrintTime = elapsed;
    shared_->lastPrintIter = iter;
    return -1;
  }

  const int iter = model_->numberIterations();
  const double elapsed = CoinWallclockTime() - shared_->startTime;

  if (iter > shared_->maxIterSeen)
    shared_->maxIterSeen = iter;

  // Detect restart: iteration counter dropped significantly
  const bool isRestart = (shared_->lastPrintIter > 0 && iter < shared_->lastPrintIter - 50);
  if (isRestart) {
    FILE *fp = handler_->filePointer();
    if (fp) {
      const bool u8 = ClpOutput::useUtf8();
      const CoinTable tbl = makeLpTable(u8, ClpOutput::useCompact());
      const int w = tbl.totalWidth();
      const std::string rule = CoinTable::sectionRule("restart", u8, w, w / 3);
      fprintf(fp, "  %s\n", rule.c_str());
      fflush(fp);
    }
    shared_->lastPrintIter = -1; // force immediate row print
    shared_->lastPrintTime = elapsed;
  }

  const bool doIter = (iterFreq_ > 0 && iter - shared_->lastPrintIter >= iterFreq_);
  const bool doTime = (timeFreq_ > 0.0 && elapsed - shared_->lastPrintTime >= timeFreq_);
  const bool doForce = (shared_->lastPrintIter < 0);

  if (doIter || doTime || doForce) {
    printRow(iter, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), elapsed, model_->algorithm());
    shared_->lastPrintTime = elapsed;
    shared_->lastPrintIter = iter;
  }
  return -1;
}

// ─── ClpLpPhaseState helpers ──────────────────────────────────────────────────
//
// The table is printed as messages (CLP_TABLE_*), so its look comes from the
// language of the messages: en_tabular gives the table, and in us_en these
// messages are null and print nothing. The layout that depends on the
// settings (UTF-8, compact) is passed in as values: the rules, the column
// headings and the column separators.

static std::string fmtTime(double t)
{
  char buf[32];
  if (t < 1.0)
    std::snprintf(buf, sizeof(buf), "%.3f", t);
  else if (t < 10.0)
    std::snprintf(buf, sizeof(buf), "%.2f", t);
  else if (t < 100.0)
    std::snprintf(buf, sizeof(buf), "%.1f", t);
  else
    std::snprintf(buf, sizeof(buf), "%.0f", t);
  return buf;
}

// Where the table's messages go: the state's handler and messages if set,
// otherwise the model's.
static CoinMessageHandler *tableHandler(const ClpLpPhaseState &s,
  ClpSimplex *model)
{
  if (s.handler)
    return s.handler;
  return model ? model->messageHandler() : nullptr;
}

static CoinMessages *tableMessages(const ClpLpPhaseState &s, ClpSimplex *model)
{
  if (s.messages)
    return s.messages;
  return model ? model->messagesPointer() : nullptr;
}

static void lpPhaseOpenTable(ClpLpPhaseState &s, CoinMessageHandler *handler,
  CoinMessages *messages)
{
  if (s.headerPrinted || !handler || !messages)
    return;
  const bool u8 = s.utf8;
  std::string lead = "\n";
  if (!s.title.empty())
    lead += CoinTable::phaseStart(s.title, u8) + "\n\n";
  const CoinTable tbl = makeLpTable(s.utf8, s.compact);
  handler->message(CLP_TABLE_HEADER, *messages)
    << lead << tbl.sepLine(CoinTable::Top) << tbl.headerLine()
    << tbl.sepLine(CoinTable::Middle) << CoinMessageEol;
  s.headerPrinted = true;
}

static void lpPhaseRow(ClpLpPhaseState &s, CoinMessageHandler *handler,
  CoinMessages *messages, const char *phase, int iter, double obj,
  double pInf, double dInf, double elapsed)
{
  if (!handler || !messages)
    return;
  const char *bar = s.compact ? " " : (s.utf8 ? " \xe2\x94\x82 " : " | ");
  handler->message(CLP_TABLE_ROW, *messages)
    << phase << bar << iter << bar << obj << bar
    << (pInf > 0.0 ? pInf : 0.0) << bar
    << (dInf > 0.0 ? dInf : 0.0) << bar
    << fmtTime(elapsed) << CoinMessageEol;
}

static const char *lpPhaseName(int algo)
{
  switch (algo) {
  case -1: return "Dual";
  case  1: return "Primal";
  case  2: return "Barrier";
  default: return "LP";
  }
}

// ─── ClpLpEventHandler ────────────────────────────────────────────────────────

ClpLpEventHandler::ClpLpEventHandler(std::shared_ptr<ClpLpPhaseState> state)
  : ClpEventHandler()
  , s_(std::move(state))
{
}

ClpEventHandler *ClpLpEventHandler::clone() const
{
  return new ClpLpEventHandler(*this);
}

void ClpLpEventHandler::openTable()
{
  lpPhaseOpenTable(*s_, tableHandler(*s_, model_), tableMessages(*s_, model_));
}

void ClpLpEventHandler::printLpRow(int iter, double obj, double pInf,
  double dInf, double elapsed)
{
  const int algo = model_ ? model_->algorithm() : 0;
  lpPhaseRow(*s_, tableHandler(*s_, model_), tableMessages(*s_, model_),
    lpPhaseName(algo), iter, obj, pInf, dInf, elapsed);
}

/*
  A row for an Idiot iteration or a Sprint pass. These show from log level 2
  up, one row for each (there are at most a few dozen of them, unlike
  simplex iterations).
*/
void ClpLpEventHandler::printPassRow(const char *phase, bool &seen,
  const ClpPassProgress &progress)
{
  if (s_->logLevel < 2 || !s_->fp)
    return;
  const double now = CoinWallclockTime();
  seen = true;
  openTable();
  lpPhaseRow(*s_, tableHandler(*s_, model_), tableMessages(*s_, model_),
    phase, progress.iteration, progress.objective,
    progress.primalInfeasibility, progress.dualInfeasibility,
    now - s_->startTime);
  s_->lastPrintTime = now;
}

int ClpLpEventHandler::eventWithInfo(Event whichEvent, void *info)
{
  if (!info)
    return -1;
  if (whichEvent == endOfBarrierIteration) {
    // The barrier model has no ClpSimplex to take the handler and messages
    // from, so use its own (it's a copy of the LP model, same language).
    ClpInterior *barrier = static_cast<ClpInterior *>(info);
    ClpLpTable::printRow(*s_, "Barrier", barrier->numberIterations(),
      barrier->objectiveValue(), barrier->sumPrimalInfeasibilities(),
      barrier->sumDualInfeasibilities(),
      s_->handler ? s_->handler : barrier->messageHandler(),
      s_->messages ? s_->messages : barrier->messagesPointer());
  } else if (whichEvent == endOfIdiotIteration) {
    printPassRow("Idiot", s_->idiotSeen,
      *static_cast<const ClpPassProgress *>(info));
  } else if (whichEvent == endOfSprintPass) {
    printPassRow("Sprint", s_->sprintSeen,
      *static_cast<const ClpPassProgress *>(info));
  }
  return -1;
}

int ClpLpEventHandler::event(Event whichEvent)
{
  if (whichEvent != endOfIteration || !model_)
    return -1;
  if (s_->logLevel <= 0 || !s_->fp)
    return -1;

  int iter = model_->numberIterations();
  const double now = CoinWallclockTime();
  const double elapsed = now - s_->startTime;

  if (!s_->lpStarted) {
    // First LP endOfIteration event: open the table if Idiot or Sprint
    // haven't already.
    s_->lpStarted = true;
    openTable();
    // Print the first row with actual iteration count (not fake 0,
    // since this event fires after the first iteration has completed).
    printLpRow(iter, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), elapsed);
    s_->lastPrintTime = now;
    s_->lastPrintIter = iter;
    return -1;
  }

  if (iter > s_->maxIterSeen)
    s_->maxIterSeen = iter;

  // Detect restart (Clp resets iter counter after perturbation / algorithm switch)
  //bool isRestart = (s_->lastPrintIter > 0 && iter < s_->lastPrintIter - 50);
  // Restart if algorithm changed
  bool isRestart = iter&&s_->modifyMsg*model_->algorithm()<0;
  if (isRestart) {
    const CoinTable tbl = makeLpTable(s_->utf8, s_->compact);
    const int w = tbl.totalWidth();
    const std::string rule = CoinTable::sectionRule("restart", s_->utf8, w, w / 3);
    CoinMessageHandler *handler = tableHandler(*s_, model_);
    CoinMessages *messages = tableMessages(*s_, model_);
    if (handler && messages)
      handler->message(CLP_TABLE_RULE, *messages) << rule << CoinMessageEol;
    //s_->lastPrintIter = -1;
    //s_->lastPrintTime = now;
  }
  static int lastIterZ=-1;
  const bool doIter = (s_->iterFreq > 0 && iter - s_->lastPrintIter >= s_->iterFreq);
  const bool doTime = (s_->timeFreq > 0.0 && now - s_->lastPrintTime >= s_->timeFreq);
  const bool doForce = (s_->lastPrintIter < 0);
  // set algorithm
  s_->modifyMsg = model_->algorithm();
  // Output after postsolve gives wrong impression
  // Should say postsolve done
  // This fix is not correct but what can I do
  // it does not matter too much
  iter = std::max(model_->numberIterations(), s_->maxIterSeen);
  // Skip if old style printing
  if ((doIter || doTime || doForce) &&
      (model_->messageHandler()->prefix()&512)==0)  {
    printLpRow(iter, model_->objectiveValue(),
      model_->sumPrimalInfeasibilities(),
      model_->sumDualInfeasibilities(), elapsed);
    s_->lastPrintTime = now;
    s_->lastPrintIter = iter;
    //if (iter<lastIterZ)
    //abort();
    lastIterZ=iter;
  }

  return -1;
}

// Appends "<dashSep>winner: <name>" to a summary buffer when a racing
// winner is recorded, so the closing status line reports which racing
// config won -- instead of relying on a separate, unformatted "LP racing:
// winner=..." print that would otherwise interleave with (and visually
// break) the progress table.
static void appendRacingWinner(char *summary, size_t bufSize,
  const ClpLpPhaseState &s, bool utf8)
{
  if (s.racingWinner.empty())
    return;
  const size_t len = std::strlen(summary);
  if (len >= bufSize)
    return;
  std::snprintf(summary + len, bufSize - len,
    "%swinner: %s", CoinTable::dashSep(utf8), s.racingWinner.c_str());
}

void ClpLpEventHandler::printFinalStatus(int numInts, int numFrac)
{
  ClpLpTable::printFinalStatus(*s_, model_, numInts, numFrac);
}

namespace ClpLpTable {

void printFinalStatus(ClpLpPhaseState &s, ClpSimplex *model,
  int numInts, int numFrac)
{
  if (s.logLevel <= 0 || !s.fp)
    return;
  CoinMessageHandler *handler = tableHandler(s, model);
  CoinMessages *messages = tableMessages(s, model);
  if (!handler || !messages)
    return;

  const bool u8 = s.utf8;
  const double elapsed = CoinWallclockTime() - s.startTime;
  const std::string tStr = fmtTime(elapsed);

  // If no intermediate output was ever produced (e.g., LP solved entirely by
  // presolve with 0 simplex iterations), print a compact single-line summary
  // without opening the full iteration table.
  if (!s.headerPrinted) {
    if (!model)
      return;
    const char *baseStatus = "Unknown";
    const int st = model->status();
    const int ss = model->secondaryStatus();
    if (st == 0)
      baseStatus = "Optimal";
    else if (st == 1) {
      if (ss == 4)
        baseStatus = "Iteration limit";
      else if (ss == 9)
        baseStatus = "Time limit";
      else
        baseStatus = "Infeasible";
    } else if (st == 2)
      baseStatus = "Unbounded";
    else if (st == 3)
      baseStatus = "Stopped";
    else if (st == 4)
      baseStatus = "Numerical difficulties";
    else if (st == 5)
      baseStatus = "Stopped by event";
    const std::string statusStr = (numInts > 0 ? "LP " : "") + std::string(baseStatus);
    char summary[512];
    const int iters = model->numberIterations();
    if (numInts > 0 && numFrac >= 0) {
      const double pct = (numInts > 0) ? 100.0 * numFrac / numInts : 0.0;
      if (iters > 0)
        std::snprintf(summary, sizeof(summary),
          "%s%sFrac: %d/%d (%.1f%%)   Obj: %g   Iters: %d   Time: %ss",
          statusStr.c_str(), CoinTable::dashSep(u8),
          numFrac, numInts, pct, model->objectiveValue(), iters, tStr.c_str());
      else
        std::snprintf(summary, sizeof(summary),
          "%s%sFrac: %d/%d (%.1f%%)   Obj: %g   Time: %ss",
          statusStr.c_str(), CoinTable::dashSep(u8),
          numFrac, numInts, pct, model->objectiveValue(), tStr.c_str());
    } else {
      if (iters > 0)
        std::snprintf(summary, sizeof(summary),
          "%s%sObj: %g   Iters: %d   Time: %ss",
          statusStr.c_str(), CoinTable::dashSep(u8),
          model->objectiveValue(), iters, tStr.c_str());
      else
        std::snprintf(summary, sizeof(summary),
          "%s%sObj: %g   Time: %ss",
          statusStr.c_str(), CoinTable::dashSep(u8),
          model->objectiveValue(), tStr.c_str());
    }
    std::string lead;
    if (!s.title.empty())
      lead = "\n" + CoinTable::phaseStart(s.title, u8) + "\n";
    appendRacingWinner(summary, sizeof(summary), s, u8);
    handler->message(CLP_TABLE_SUMMARY, *messages)
      << lead << CoinTable::phaseEnd(summary, u8) << CoinMessageEol;
    return;
  }

  // Intermediate output was printed: print the final LP row if it was
  // skipped, then close the table.
  const CoinTable tbl = makeLpTable(u8, s.compact);

  if (s.lpStarted && model) {
    const int iters = std::max(model->numberIterations(), s.maxIterSeen);
    if (s.lastPrintIter < iters) {
      const double now = CoinWallclockTime();
      lpPhaseRow(s, handler, messages, lpPhaseName(model->algorithm()),
        iters, model->objectiveValue(),
        model->sumPrimalInfeasibilities(),
        model->sumDualInfeasibilities(), now - s.startTime);
    }
  }

  if (!model) {
    handler->message(CLP_TABLE_END, *messages)
      << tbl.sepLine(CoinTable::Bottom)
      << CoinTable::phaseEnd("LP complete", u8) << CoinMessageEol;
    return;
  }

  // Build status string
  const char *baseStatus = "Unknown";
  const int st = model->status();
  const int ss = model->secondaryStatus();
  if (st == 0)
    baseStatus = "Optimal";
  else if (st == 1) {
    if (ss == 4)
      baseStatus = "Iteration limit";
    else if (ss == 9)
      baseStatus = "Time limit";
    else
      baseStatus = "Infeasible";
  } else if (st == 2)
    baseStatus = "Unbounded";
  else if (st == 3)
    baseStatus = "Stopped";
  else if (st == 4)
    baseStatus = "Numerical difficulties";
  else if (st == 5)
    baseStatus = "Stopped by event";

  const std::string statusStr = (numInts > 0 ? "LP " : "") + std::string(baseStatus);
  const int iters = s.lpStarted
    ? std::max(model->numberIterations(), s.maxIterSeen)
    : 0;

  char summary[512];
  if (numInts > 0 && numFrac >= 0) {
    const double pct = (numInts > 0) ? 100.0 * numFrac / numInts : 0.0;
    std::snprintf(summary, sizeof(summary),
      "%s%sFrac: %d/%d (%.1f%%)   Obj: %g   Iters: %d   Time: %ss",
      statusStr.c_str(), CoinTable::dashSep(u8),
      numFrac, numInts, pct, model->objectiveValue(), iters, tStr.c_str());
  } else {
    std::snprintf(summary, sizeof(summary),
      "%s%sObj: %g   Iters: %d   Time: %ss",
      statusStr.c_str(), CoinTable::dashSep(u8),
      model->objectiveValue(), iters, tStr.c_str());
  }
  appendRacingWinner(summary, sizeof(summary), s, u8);
  handler->message(CLP_TABLE_END, *messages)
    << tbl.sepLine(CoinTable::Bottom) << CoinTable::phaseEnd(summary, u8)
    << CoinMessageEol;
}

void printRow(ClpLpPhaseState &state, const char *phaseLabel, int iter,
  double obj, double pInf, double dInf, CoinMessageHandler *handler,
  CoinMessages *messages)
{
  std::lock_guard<std::mutex> lock(state.mu);
  if (state.logLevel <= 0 || !state.fp)
    return;
  if (!handler)
    handler = state.handler;
  if (!messages)
    messages = state.messages;

  const double now = CoinWallclockTime();
  const double elapsed = now - state.startTime;

  if (!state.lpStarted) {
    // First row from any racing config: open the table (mirrors the
    // sequential path always showing the first LP row) and print it
    // unconditionally, regardless of timeFreq.
    state.lpStarted = true;
    lpPhaseOpenTable(state, handler, messages);
    lpPhaseRow(state, handler, messages, phaseLabel, iter, obj, pInf, dInf,
      elapsed);
    state.lastPrintTime = now;
    state.lastPrintIter = iter;
    return;
  }

  // Rate-limited by timeFreq *across all racing configs combined* -- only
  // whichever config happens to report once the interval has elapsed gets
  // to print, so the table stays as compact as the sequential one instead
  // of every thread printing independently.
  const bool doTime = (state.timeFreq > 0.0 && now - state.lastPrintTime >= state.timeFreq);
  if (doTime) {
    lpPhaseRow(state, handler, messages, phaseLabel, iter, obj, pInf, dInf,
      elapsed);
    state.lastPrintTime = now;
    state.lastPrintIter = iter;
  }
}

} // namespace ClpLpTable

// ─── ClpLpMsgHandler ──────────────────────────────────────────────────────────

ClpLpMsgHandler::ClpLpMsgHandler(std::shared_ptr<ClpLpPhaseState> state)
  : CoinMessageHandler()
  , s_(std::move(state))
{
  // The model's log level: the language decides how each message looks,
  // and its detail level whether it shows.
  setLogLevel(s_->logLevel);
  // Tables read better without the Clpnnnn prefix; set it to show it.
  setPrefix(false);
  if (s_->fp)
    setFilePointer(s_->fp);
}

CoinMessageHandler *ClpLpMsgHandler::clone() const
{
  return new ClpLpMsgHandler(*this);
}

// ─── helpers ──────────────────────────────────────────────────────────────────

void ClpOutput::sendLine(CoinMessageHandler *handler, const std::string &line)
{
  FILE *fp = handler->filePointer();
  if (fp) {
    fprintf(fp, "%s\n", line.c_str());
    fflush(fp);
  }
}

void ClpOutput::sendBlank(CoinMessageHandler *handler)
{
  // CoinMessageHandler suppresses empty strings; write directly.
  FILE *fp = handler->filePointer();
  if (fp)
    fprintf(fp, "\n");
}

// ─── RangeStat ────────────────────────────────────────────────────────────────

void ClpOutput::RangeStat::update(double v)
{
  double a = std::fabs(v);
  if (a == 0.0)
    return;
  valid = true;
  if (a < lo)
    lo = a;
  if (a > hi)
    hi = a;
}

std::string ClpOutput::RangeStat::ratio() const
{
  if (!valid || lo == 0.0)
    return "-";
  double r = hi / lo;
  std::ostringstream os;
  // Use integer-style when the ratio is exact and small
  if (r == std::floor(r) && r < 1e6)
    os << (long long)r;
  else
    os << std::setprecision(4) << r;
  return os.str();
}

// ─── fmtNum ───────────────────────────────────────────────────────────────────

std::string ClpOutput::fmtNum(double v)
{
  if (v == 0.0)
    return "0";
  char buf[64];
  // Use %.5g: compact, switches to scientific when needed.
  std::snprintf(buf, sizeof(buf), "%.5g", v);
  return buf;
}

// ─── printSolverHeader ────────────────────────────────────────────────────────

void ClpOutput::printSolverHeader(CoinMessageHandler *handler, int logLevel,
  int argc, const char *const *argv)
{
  if (logLevel < 1)
    return;

  handler->setPrefix(false);

  std::ostringstream hdr;
  hdr << "CLP " << CLP_VERSION;

#ifdef CLP_GIT_HASH
  hdr << " (git:" << CLP_GIT_HASH << ")";
#endif
  hdr << " \xe2\x80\x94 COIN-OR LP Solver";
  sendLine(handler, hdr.str());

  if (argc > 1 && argv) {
    std::ostringstream args;
    args << "  args:";
    for (int i = 1; i < argc; ++i)
      args << " " << argv[i];
    sendLine(handler, args.str());
  }
}

// ─── printProblemSummary ──────────────────────────────────────────────────────

void ClpOutput::printProblemSummary(CoinMessageHandler *handler,
  const ClpSimplex &model, int logLevel)
{
  if (logLevel < 1)
    return;

  handler->setPrefix(false);

  int nRows = model.numberRows();
  int nCols = model.numberColumns();

  // ── problem name ──────────────────────────────────────────────────────────
  const std::string &name = model.problemName();
  std::string displayName = name.empty() ? "(unnamed)" : name;

  std::ostringstream line;
  line << "Problem: " << displayName;
  sendBlank(handler);
  sendLine(handler, line.str());

  // ── dimensions ────────────────────────────────────────────────────────────
  const CoinPackedMatrix *mat = model.matrix();
  int nNZ = mat ? mat->getNumElements() : 0;

  // Widths for alignment
  int wRows = std::to_string(nRows).size();
  int wCols = std::to_string(nCols).size();
  int wNZ   = std::to_string(nNZ).size();
  int w = std::max({ wRows, wCols, wNZ });

  std::ostringstream dim;
  dim << "  |Rows| = " << std::setw(w) << nRows
      << "   |Cols| = " << std::setw(w) << nCols
      << "   |NZ| = "  << std::setw(w) << nNZ;
  sendLine(handler, dim.str());

  // ── coefficient ranges ────────────────────────────────────────────────────
  RangeStat rMat, rObj, rBnd, rRhs;

  // Matrix elements
  if (mat) {
    const double *elems = mat->getElements();
    int ne = mat->getNumElements();
    for (int i = 0; i < ne; ++i)
      rMat.update(elems[i]);
  }

  // Objective coefficients
  const double *objCoef = model.getObjCoefficients();
  for (int j = 0; j < nCols; ++j)
    rObj.update(objCoef[j]);

  // Column bounds
  const double *lb = model.columnLower();
  const double *ub = model.columnUpper();
  for (int j = 0; j < nCols; ++j) {
    if (lb[j] > -1e20)
      rBnd.update(lb[j]);
    if (ub[j] < 1e20)
      rBnd.update(ub[j]);
  }

  // RHS (row bounds)
  const double *rl = model.rowLower();
  const double *ru = model.rowUpper();
  for (int i = 0; i < nRows; ++i) {
    if (rl[i] > -1e20)
      rRhs.update(rl[i]);
    if (ru[i] < 1e20)
      rRhs.update(ru[i]);
  }

  // ── format ranges table ───────────────────────────────────────────────────
  const char *inSym  = useUtf8_ ? "\xe2\x88\x88" : "in";    // ∈ or "in"
  const char *kapSym = useUtf8_ ? "\xce\xba \xe2\x89\x88" : "ratio ="; // κ ≈ or "ratio ="

  struct Row {
    const char *label;
    const RangeStat &rs;
  };
  Row rows[] = {
    { "Matrix", rMat },
    { "Cost",   rObj },
    { "Bounds", rBnd },
    { "RHS",    rRhs },
  };

  // Compute column widths for alignment
  int wLo = 1, wHi = 1;
  for (auto &r : rows) {
    if (!r.rs.valid)
      continue;
    wLo = std::max(wLo, (int)fmtNum(r.rs.lo).size());
    wHi = std::max(wHi, (int)fmtNum(r.rs.hi).size());
  }

  sendBlank(handler);
  sendLine(handler, "  Coefficient ranges:");

  for (auto &r : rows) {
    std::ostringstream row;
    row << "    " << std::left << std::setw(8) << r.label << std::right;
    if (!r.rs.valid) {
      row << "  " << inSym << " [ -, -]   " << kapSym << " -";
    } else {
      std::string lo = fmtNum(r.rs.lo), hi = fmtNum(r.rs.hi);
      row << "  " << inSym << " [" << std::setw(wLo) << lo
          << ", " << std::setw(wHi) << hi << "]"
          << "   " << kapSym << " " << r.rs.ratio();
    }
    sendLine(handler, row.str());
  }

  sendBlank(handler);
}
// Print presolve stats
void ClpLpEventHandler::printPresolveStats()
{
  // Print presolve stats
  if (model_ && model_->presolveRows() >= 0 && s_->origRows > 0) {
    CoinMessageHandler *handler = tableHandler(*s_, model_);
    CoinMessages *messages = tableMessages(*s_, model_);
    if (handler && messages)
      handler->message(CLP_PRESOLVE_SUMMARY, *messages)
        << s_->origRows << s_->origCols
        << model_->presolveRows() << model_->presolveCols()
        << model_->presolveTime() << CoinMessageEol;
  }
}
