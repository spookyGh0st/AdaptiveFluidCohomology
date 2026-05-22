#pragma once

#include <chrono>

class Stopwatch {
  public:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;

    // Constructor: takes a reference to where the duration will be stored (in milliseconds)
    explicit Stopwatch(double &outputDuration);

    // Destructor: automatically calculates and stores the elapsed time
    ~Stopwatch();

  private:
    double &output;
    TimePoint startTime;
};
