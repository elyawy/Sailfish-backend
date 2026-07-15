#include <iostream>
#include <vector>
#include <cmath>

#include "../../../src/nonRevAccelerator.h"

// Uniform off-diagonal Q for 61 states:
// diagonal = -60, off-diagonal = 1
// Simple, valid, non-reversible (asymmetric by construction once we break one entry)
static std::vector<double> makeQ61() {
    const int N = 61;
    std::vector<double> Q(N * N, 1.0);
    for (int i = 0; i < N; ++i)
        Q[i * N + i] = -60.0;
    // Break one symmetry to make it genuinely non-reversible
    Q[0 * N + 1] = 2.0;
    Q[0 * N + 0] = -61.0;
    return Q;
}

static std::vector<double> makeFreq61() {
    const int N = 61;
    std::vector<double> f(N, 1.0 / N);
    return f;
}

static bool approxEq(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) < tol;
}

static bool testIdentityAtZero61() {
    const int N = 61;
    nonRevAccelerator<N> acc(makeQ61(), makeFreq61());
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (!approxEq(acc.Pij_t(i, j, 0.0), expected, 1e-9)) {
                std::cerr << "FAIL testIdentityAtZero61: P(" << i << "," << j
                          << ",0) = " << acc.Pij_t(i, j, 0.0) << "\n";
                return false;
            }
        }
    std::cout << "PASS testIdentityAtZero61\n";
    return true;
}

static bool testRowSums61() {
    const int N = 61;
    nonRevAccelerator<N> acc(makeQ61(), makeFreq61());
    for (double t : { 0.01, 0.1, 1.0 }) {
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j)
                sum += acc.Pij_t(i, j, t);
            if (!approxEq(sum, 1.0, 1e-9)) {
                std::cerr << "FAIL testRowSums61: row " << i << " at t=" << t
                          << " sums to " << sum << "\n";
                return false;
            }
        }
    }
    std::cout << "PASS testRowSums61\n";
    return true;
}

static bool testNonNegative61() {
    const int N = 61;
    nonRevAccelerator<N> acc(makeQ61(), makeFreq61());
    for (double t : { 0.01, 0.5, 2.0 }) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                double p = acc.Pij_t(i, j, t);
                if (p < -1e-10) {
                    std::cerr << "FAIL testNonNegative61: P(" << i << "," << j
                              << ",t=" << t << ") = " << p << "\n";
                    return false;
                }
            }
    }
    std::cout << "PASS testNonNegative61\n";
    return true;
}

static bool testChapmanKolmogorov61() {
    const int N = 61;
    nonRevAccelerator<N> acc(makeQ61(), makeFreq61());
    const double t1 = 0.3, t2 = 0.7;

    // Only check a sample of (i,j) pairs — full N*N*N is slow
    for (int i : { 0, 10, 30, 60 }) {
        for (int j : { 0, 10, 30, 60 }) {
            double product = 0.0;
            for (int k = 0; k < N; ++k)
                product += acc.Pij_t(i, k, t1) * acc.Pij_t(k, j, t2);
            double direct = acc.Pij_t(i, j, t1 + t2);
            if (!approxEq(product, direct, 1e-7)) {
                std::cerr << "FAIL testChapmanKolmogorov61: P(" << i << "," << j
                          << ") product=" << product << " direct=" << direct << "\n";
                return false;
            }
        }
    }
    std::cout << "PASS testChapmanKolmogorov61\n";
    return true;
}

int main() {
    std::cout << "Testing nonRevAccelerator<61>...\n";
    bool ok = true;
    ok &= testIdentityAtZero61();
    ok &= testRowSums61();
    ok &= testNonNegative61();
    ok &= testChapmanKolmogorov61();
    return ok ? 0 : 1;
}
