#include <iostream>
#include <windows.h>
#include <memory>
#include "code.h"

static constexpr int iterations = 1'000'000;

void init(std::unique_ptr<BenchmarkTEST>& b) {
	b->setup();
}

void benchmark_with_eigen(std::unique_ptr<BenchmarkTEST>& b) {
	for (int i = 0; i < iterations; ++i) {
		b->body_with_eigen_operation();
	}
}

void benchmark_without_eigen(std::unique_ptr<BenchmarkTEST>& b) {
	for (int i = 0; i < iterations; ++i) {
		b->body_no_eigen_operation();
	}
}

int main()
{
	std::unique_ptr<BenchmarkTEST> b = std::make_unique<BenchmarkTEST>();

	init(b);
	LARGE_INTEGER tick_start, tick_end;

	QueryPerformanceCounter(&tick_start);
	benchmark_with_eigen(b);
	QueryPerformanceCounter(&tick_end);
	int64_t ticks_with_eigen = tick_end.QuadPart - tick_start.QuadPart;

	QueryPerformanceCounter(&tick_start);
	benchmark_without_eigen(b);
	QueryPerformanceCounter(&tick_end);
	int64_t ticks_without_eigen = tick_end.QuadPart - tick_start.QuadPart;

	std::cout << "   With Eigen : " << ticks_with_eigen << " ticks\n";
	std::cout << "Without Eigen : " << ticks_without_eigen << " ticks\n";
}
