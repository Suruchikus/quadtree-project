#pragma once
#include <ctime>   // clock(), clock_t, CLOCKS_PER_SEC

struct Timer {
    clock_t t0;

    Timer() : t0(clock()) {}

    double ms() const {
        clock_t t1 = clock();
        return 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
    }
};

