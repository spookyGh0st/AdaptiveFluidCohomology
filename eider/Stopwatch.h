#pragma once

#include <chrono>

class Stopwatch {
  public:
    Stopwatch();

    // Stop the stopwatch and return elapsed time in seconds
    std::chrono::duration<double> stop();

    // Optional: get elapsed time without stopping
    std::chrono::duration<double> elapsed() const;

  private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::duration<double> elapsed_ms{0.0};
    bool running{false};
};
