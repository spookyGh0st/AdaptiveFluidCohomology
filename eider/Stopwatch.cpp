#include "Stopwatch.h"
Stopwatch::Stopwatch() : start_time(std::chrono::high_resolution_clock::now()), running(true) {}
std::chrono::duration<double> Stopwatch::stop() {
    if (!running) return elapsed_ms;
    auto end_time = std::chrono::high_resolution_clock::now();
    elapsed_ms = std::chrono::duration<double>(end_time - start_time);
    running = false;
    return elapsed_ms;
}
std::chrono::duration<double> Stopwatch::elapsed() const {
    auto now = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(now - start_time);
}
