#include "Stopwatch.h"
Stopwatch::Stopwatch(double &outputDuration) : output(outputDuration), startTime(Clock::now()) {}

Stopwatch::~Stopwatch() {
    auto endTime = Clock::now();
    std::chrono::duration<double, std::milli> duration = endTime - startTime;
    output = duration.count(); // store duration in milliseconds
}
