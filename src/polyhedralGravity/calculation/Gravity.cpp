#include "Gravity.h"

void Gravity::calculate() {
    SPDLOG_INFO("Calculate...");
}

double Gravity::calculateLNpq(size_t p, size_t q) {
    auto endpoints = _polyhedron.getPolyhedralSegment(p, q);
    //TODO
    return 0.0;
}
