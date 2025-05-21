#include <iostream>
#include <windows.h>

static constexpr int iterations = 1'000'000;

void init() {
    // TODO
}

void benchmark_with_eigen() {
    for (int i = 0; i < iterations; ++i) {
        // TODO
    }
}

void benchmark_without_eigen() {
    for (int i = 0; i < iterations; ++i) {
        // TODO
    }
}

int main()
{
    init();
    LARGE_INTEGER tick_start, tick_end;

    QueryPerformanceCounter(&tick_start);
    benchmark_with_eigen();
    QueryPerformanceCounter(&tick_end);
    int64_t ticks_with_eigen = tick_end.QuadPart - tick_start.QuadPart;

    QueryPerformanceCounter(&tick_start);
    benchmark_without_eigen();
    QueryPerformanceCounter(&tick_end);
    int64_t ticks_without_eigen = tick_end.QuadPart - tick_start.QuadPart;

    std::cout << "   With Eigen : " << ticks_with_eigen << " ticks\n";
    std::cout << "Without Eigen : " << ticks_without_eigen << " ticks\n";
}
