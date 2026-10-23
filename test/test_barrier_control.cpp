// Standalone: c++ -O2 -std=c++11 test_barrier_control.cpp \
//   $(pkg-config --cflags --libs --static clp) -lpthread -ldl -o test_barrier_control
#include "ClpSimplex.hpp"
#include "ClpInterior.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpSolve.hpp"
#include "ClpRacingSolver.hpp"
#include "ClpEventHandler.hpp"
#include "ClpOutput.hpp"
#include "Idiot.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
struct Events {
  std::atomic<bool> stop{false};
  int progress = 0;
};

class Handler : public ClpEventHandler {
public:
  explicit Handler(std::shared_ptr<Events> events) : events_(events) {}
  ClpEventHandler *clone() const override { return new Handler(*this); }
  int event(Event) override { return events_->stop.load() ? 0 : -1; }
  int eventWithInfo(Event event, void *info) override {
    if (event == endOfBarrierIteration && info) {
      const ClpInterior *barrier = static_cast<const ClpInterior *>(info);
      if (!std::isfinite(barrier->objectiveValue()))
        throw std::runtime_error("Non-finite barrier progress");
      ++events_->progress;
    }
    return events_->stop.load() ? 0 : -1;
  }
private:
  std::shared_ptr<Events> events_;
};

struct IdiotEvents {
  int checkpoints = 0;
  int simplexEventsAfterStop = 0;
  bool stopped = false;
};

class SymbolicHandler : public ClpEventHandler {
public:
  SymbolicHandler(std::shared_ptr<int> checkpoints, int stopAt)
    : checkpoints_(checkpoints), stopAt_(stopAt) {}
  ClpEventHandler *clone() const override { return new SymbolicHandler(*this); }
  int event(Event event) override
  {
    if (event == barrierCheckpoint && ++*checkpoints_ == stopAt_)
      return 0;
    return -1;
  }
private:
  std::shared_ptr<int> checkpoints_;
  int stopAt_;
};

class CountedCholesky : public ClpCholeskyBase {
public:
  void seedCount(int count) { sizeFactor_ = count; }
};

class CountHandler : public ClpEventHandler {
public:
  CountHandler(CountedCholesky *cholesky, int count)
    : cholesky_(cholesky), count_(count), seeded_(false) {}
  ClpEventHandler *clone() const override { return new CountHandler(*this); }
  int event(Event event) override
  {
    if (event == barrierCheckpoint && !seeded_) {
      cholesky_->seedCount(count_);
      seeded_ = true;
    }
    return -1;
  }
private:
  CountedCholesky *cholesky_;
  int count_;
  bool seeded_;
};

class IdiotHandler : public ClpEventHandler {
public:
  IdiotHandler(std::shared_ptr<IdiotEvents> events, int mode, int stopAt)
    : events_(events), mode_(mode), stopAt_(stopAt) {}
  ClpEventHandler *clone() const override { return new IdiotHandler(*this); }
  int event(Event event) override
  {
    if (events_->stopped && (event == endOfIteration || event == startOfCrossover))
      ++events_->simplexEventsAfterStop;
    if (events_->stopped && event == presolveAfterFirstSolve && mode_ == 3)
      return 1;
    if (event != idiotCheckpoint)
      return -1;
    if (++events_->checkpoints < stopAt_)
      return -1;
    events_->stopped = true;
    if (mode_ == 0)
      simplex()->setMaximumWallSeconds(0);
    else if (mode_ == 1)
      simplex()->setMaximumSeconds(0);
    else if (mode_ == 2)
      return 0;
    else
      simplex()->setMaximumIterations(0);
    return -1;
  }
private:
  std::shared_ptr<IdiotEvents> events_;
  int mode_;
  int stopAt_;
};

void require(bool condition, const char *message)
{
  if (!condition)
    throw std::runtime_error(message);
}

ClpSimplex example()
{
  const int columns = 192;
  CoinPackedMatrix matrix(false, 0, 0);
  matrix.setDimensions(0, columns);
  std::vector<double> lower(columns, 0), upper(columns, 1), objective(columns);
  std::vector<double> rowLower(128, -COIN_DBL_MAX), rowUpper(128);
  for (int j = 0; j < columns; ++j)
    objective[j] = -1.0 - (j % 13) / 13.0;
  for (int i = 0; i < 128; ++i) {
    CoinPackedVector row;
    double sum = 0;
    for (int j = 0; j < columns; ++j) {
      const double value = 1 + ((i + 7) * (j + 11) % 101);
      row.insert(j, value);
      sum += value;
    }
    matrix.appendRow(row);
    rowUpper[i] = 0.4 * sum;
  }
  ClpSimplex m;
  m.setLogLevel(0);
  m.loadProblem(matrix, lower.data(), upper.data(), objective.data(),
    rowLower.data(), rowUpper.data());
  return m;
}

ClpSolve barrierOptions()
{
  ClpSolve opts;
  opts.setSolveType(ClpSolve::useBarrier);
  opts.setPresolveType(ClpSolve::presolveOff);
  opts.setSpecialOption(4, 0);
  return opts;
}

ClpSolve idiotOptions(bool presolve)
{
  ClpSolve opts;
  opts.setSolveType(ClpSolve::usePrimal);
  opts.setPresolveType(presolve ? ClpSolve::presolveOn : ClpSolve::presolveOff);
  opts.setSpecialOption(1, 2, 50);
  opts.setSpecialOption(2, 1);
  return opts;
}

void testSymbolicCancellation(const ClpSimplex &source)
{
  int totalCheckpoints = 0;
  for (int stopAt = 0; stopAt <= totalCheckpoints; ++stopAt) {
    ClpInterior model(source);
    model.setLogLevel(0);
    Handler orderingHandler(std::make_shared<Events>());
    model.passInEventHandler(&orderingHandler);
    ClpCholeskyBase cholesky;
    cholesky.setGoDense(1.1);
    require(cholesky.order(&model) == 0, "Symbolic test ordering failed");
    auto checkpoints = std::make_shared<int>(0);
    SymbolicHandler handler(checkpoints, stopAt);
    model.passInEventHandler(&handler);
    const int result = cholesky.symbolic();
    if (!stopAt) {
      require(result == 0, "Uninterrupted symbolic factorization failed");
      totalCheckpoints = *checkpoints;
      require(totalCheckpoints >= 3, "Too few symbolic cancellation checkpoints");
    } else {
      require(*checkpoints == stopAt && result != 0 && model.status() == 5,
        "Symbolic cancellation was not preserved");
    }
  }
}

void testOrderingSizeLimit(const ClpSimplex &source)
{
  const int memoryLimit = static_cast<int>(
    std::min(static_cast<double>(COIN_INT_MAX),
      static_cast<double>(4LL * 1024 * 1024 * 1024) / sizeof(CoinBigIndex)));
  // Seed the running count to exercise large-graph limits without allocating one.
  for (int count : {COIN_INT_MAX - 1, memoryLimit - 1}) {
    ClpInterior model(source);
    model.setLogLevel(0);
    CountedCholesky cholesky;
    CountHandler handler(&cholesky, count);
    model.passInEventHandler(&handler);
    require(cholesky.order(&model) != 0, "Oversized ordering graph was accepted");
    require(model.status() != 3 && model.status() != 5,
      "Ordering size limit was misreported as cancellation");
  }
}

void testBarrierFailureCleanup()
{
  const int size = 40000;
  std::vector<CoinBigIndex> starts(size + 1);
  std::vector<int> indices(2 * size);
  std::vector<double> elements(2 * size, 1.0);
  for (int i = 0; i < size; ++i) {
    starts[i] = 2 * i;
    indices[2 * i] = i;
    indices[2 * i + 1] = (i + 1) % size;
  }
  starts[size] = 2 * size;
  std::vector<double> lower(size, 0.0), upper(size, 1.0), objective(size, -1.0);
  for (int presolve = 0; presolve < 2; ++presolve) {
    ClpSimplex model;
    model.setLogLevel(0);
    model.loadProblem(size, size, starts.data(), indices.data(), elements.data(),
      lower.data(), upper.data(), objective.data(),
      lower.data(), upper.data());
    const int scaling = model.scalingFlag();
    const int perturbation = model.perturbation();
    ClpSolve options = barrierOptions();
    options.setPresolveType(presolve ? ClpSolve::presolveOn : ClpSolve::presolveOff);
    options.setSpecialOption(2, 1);
    options.setSpecialOption(4, 1 | 8192); // Dense setup exceeds the existing memory cap.
    bool failed = false;
    try {
      model.initialSolve(options);
    } catch (const CoinError &error) {
      failed = error.message().find("Cholesky setup/factorization failed") != std::string::npos;
      if (!failed)
        throw;
    }
    require(failed && model.status() == 4, "Strict barrier failure was not preserved");
    require(model.scalingFlag() == scaling && model.perturbation() == perturbation,
      "Strict barrier failure did not restore solver settings");
    model.setMaximumIterations(1);
    model.dual();
    require(model.status() == 0 || model.status() == 3,
      "Model could not be reused after strict barrier failure");
  }
}

void testIdiotControl(const ClpSimplex &source, double reference)
{
  for (int presolve = 0; presolve < 2; ++presolve) {
    for (int mode = 0; mode < 4; ++mode) {
      for (int stopAt : {1, 4}) {
        ClpSimplex m(source);
        const int scaling = m.scalingFlag();
        auto events = std::make_shared<IdiotEvents>();
        IdiotHandler handler(events, mode, stopAt);
        m.passInEventHandler(&handler);
        ClpSolve opts = idiotOptions(presolve != 0);
        m.initialSolve(opts);
        require(events->checkpoints >= stopAt, "Idiot checkpoint not reached");
        require(m.status() == (mode == 2 ? 5 : 3), "Incorrect idiot stopped status");
        require(m.secondaryStatus() == (mode == 2 ? ClpEventHandler::idiotCheckpoint
          : (mode == 3 ? 0 : 9)), "Incorrect idiot secondary status");
        require(!events->simplexEventsAfterStop, "Simplex/crossover restarted after idiot stop");
        require(m.scalingFlag() == scaling, "Idiot stop did not restore scaling");
        for (int j = 0; j < m.numberColumns(); ++j)
          require(m.objective()[j] == source.objective()[j]
            && m.columnLower()[j] == source.columnLower()[j]
            && m.columnUpper()[j] == source.columnUpper()[j], "Idiot stop changed the model");
        for (int i = 0; i < m.numberRows(); ++i)
          require(m.rowLower()[i] == source.rowLower()[i]
            && m.rowUpper()[i] == source.rowUpper()[i], "Idiot stop changed row bounds");
        ClpEventHandler clear;
        m.passInEventHandler(&clear);
        m.setMaximumWallSeconds(-1);
        m.setMaximumSeconds(-1);
        m.setMaximumIterations(100000);
        m.dual();
        require(m.status() == 0 && std::fabs(m.objectiveValue() - reference) < 1e-5,
          "Interrupted idiot model could not be reused");
      }
    }
  }
  {
    ClpSimplex m(source);
    auto events = std::make_shared<IdiotEvents>();
    IdiotHandler handler(events, 2, 4);
    m.passInEventHandler(&handler);
    Idiot idiot(m);
    idiot.crash(50, m.messageHandler(), m.messagesPointer());
    require(m.status() == 5 && !events->simplexEventsAfterStop,
      "Direct idiot crash did not skip crossover");
  }
  {
    ClpSimplex m(source);
    ClpRacingSolver racer(&m, 2);
    racer.addConfig(idiotOptions(true));
    ClpSolve dual;
    dual.setSolveType(ClpSolve::useDual);
    dual.setPresolveType(ClpSolve::presolveOff);
    racer.addConfig(dual, [](ClpSimplex *clone) { clone->dual(); });
    require(racer.solve() >= 0 && m.status() == 0
      && std::fabs(m.objectiveValue() - reference) < 1e-5,
      "Idiot-losing race corrupted the winner");
  }
}
}

int main(int argc, char **argv)
{
  try {
    ClpSimplex source = example();
    ClpSolve options = barrierOptions();
    ClpSimplex reference(source);
    reference.dual();
    require(reference.status() == 0, "Reference dual failed");
    testSymbolicCancellation(source);
    testOrderingSizeLimit(source);
    testBarrierFailureCleanup();
    testIdiotControl(source, reference.objectiveValue());
    {
      ClpSimplex m(source);
      auto events = std::make_shared<Events>();
      Handler handler(events);
      m.passInEventHandler(&handler);
      m.initialSolve(options);
      require(m.status() == 0, "Barrier failed");
      require(events->progress > 0, "No typed barrier progress events");
      require(std::fabs(m.objectiveValue() - reference.objectiveValue()) < 1e-5,
        "Barrier objective mismatch");
    }
    for (int mode = 0; mode < 3; ++mode) {
      ClpSimplex m(source);
      auto events = std::make_shared<Events>();
      Handler handler(events);
      m.passInEventHandler(&handler);
      if (mode == 1)
        events->stop = true;
      else if (mode == 0)
        m.setMaximumWallSeconds(0);
      else
        m.setMaximumSeconds(0);
      m.initialSolve(options);
      require(m.status() == (mode == 1 ? 5 : 3), "Incorrect interrupted status");
      m.setMaximumWallSeconds(-1);
      m.setMaximumSeconds(-1);
      events->stop = false;
      m.dual();
      require(m.status() == 0, "Interrupted model could not be reused");
    }
    {
      ClpSimplex m(source);
      FILE *log = std::tmpfile();
      require(log != NULL, "Cannot open progress log");
      auto state = std::make_shared<ClpLpPhaseState>();
      state->fp = log;
      state->logLevel = 1;
      state->timeFreq = 0;
      state->startTime = CoinWallclockTime();
      ClpLpEventHandler handler(state);
      m.passInEventHandler(&handler);
      ClpRacingSolver racer(&m, 2);
      racer.addConfig(barrierOptions());
      ClpSolve stopped;
      stopped.setSolveType(ClpSolve::useDual);
      stopped.setPresolveType(ClpSolve::presolveOff);
      racer.addConfig(stopped, [](ClpSimplex *clone) { clone->setMaximumIterations(0); });
      require(racer.solve() == 0, "Barrier did not win against stopped competitor");
      require(std::string(racer.winnerName()) == "barrier", "Incorrect winner label");
      std::rewind(log);
      char text[1024];
      bool found = false;
      while (std::fgets(text, sizeof(text), log))
        found = found || std::string(text).find("Barrier") != std::string::npos;
      std::fclose(log);
      require(found, "Barrier missing from shared racing log");
      require(m.status() == 0 && m.statusArray(), "Winner solution/basis missing");
    }
    {
      ClpSimplex m(reference);
      ClpRacingSolver racer(&m, 2);
      racer.addConfig(barrierOptions());
      ClpSolve dual;
      dual.setSolveType(ClpSolve::useDual);
      dual.setPresolveType(ClpSolve::presolveOff);
      racer.addConfig(dual);
      require(racer.solve() == 1, "Warm-started dual did not beat barrier");
      require(m.status() == 0 && std::fabs(m.objectiveValue()
        - reference.objectiveValue()) < 1e-5, "Barrier-losing race corrupted winner");
    }
    {
      ClpSimplex m(source);
      ClpRacingSolver racer(&m, 2);
      racer.addConfig(barrierOptions());
      racer.addConfig(barrierOptions());
      require(racer.solve() >= 0 && m.status() == 0, "Two-barrier race failed");
      require(std::fabs(m.objectiveValue() - reference.objectiveValue()) < 1e-5,
        "Racing winner objective mismatch");
    }
    if (argc >= 2) {
      if (argc == 3) {
        require(std::strcmp(argv[2], "--idiot") == 0, "Unknown stress-test option");
        options = idiotOptions(true);
      }
      for (int mode = 0; mode < 3; ++mode) {
        ClpSimplex m;
        m.setLogLevel(0);
        require(m.readMps(argv[1]) == 0, "Cannot read stress fixture");
        auto events = std::make_shared<Events>();
        Handler handler(events);
        m.passInEventHandler(&handler);
        std::thread timer;
        if (mode == 1)
          timer = std::thread([events]() {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            events->stop = true;
          });
        else if (mode == 0)
          m.setMaximumWallSeconds(1);
        else
          m.setMaximumSeconds(1);
        const double start = CoinWallclockTime();
        m.initialSolve(options);
        const double elapsed = CoinWallclockTime() - start;
        if (timer.joinable())
          timer.join();
        std::printf("Stress %s elapsed %.3fs status %d\n",
          mode == 1 ? "cancellation" : (mode == 0 ? "deadline" : "cpu"),
          elapsed, m.status());
        require(m.status() == 0 || m.status() == (mode == 1 ? 5 : 3),
          "Stress solve neither finished nor stopped correctly");
        require(elapsed < 10.0, "Large stress stopping delay (over 10 seconds)");
      }
    }
    std::puts("Barrier and idiot control tests passed");
    return 0;
  } catch (const CoinError &error) {
    std::fprintf(stderr, "FAIL: %s\n", error.message().c_str());
  } catch (const std::exception &error) {
    std::fprintf(stderr, "FAIL: %s\n", error.what());
  }
  return 1;
}
