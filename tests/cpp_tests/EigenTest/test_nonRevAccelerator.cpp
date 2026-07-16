#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include "../../../src/nonRevAccelerator.h"
#include "../../../libs/Eigen/Eigen/Eigenvalues"

// 4-state cycle Q matrix:
//   A->C, C->G, G->T, T->A, all rate 1
// rows sum to 0, clearly non-reversible (no detailed balance)
static const std::vector<double> Q_CYCLE = {
    -1.0,  1.0,  0.0,  0.0,
     0.0, -1.0,  1.0,  0.0,
     0.0,  0.0, -1.0,  1.0,
     1.0,  0.0,  0.0, -1.0,
};

// Uniform root frequencies (arbitrary for these tests)
static const std::vector<double> FREQ = { 0.25, 0.25, 0.25, 0.25 };

static bool approxEq(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) < tol;
}

// Test 1: P(0) == identity
static bool testIdentityAtZero() {
    nonRevAccelerator<4> acc(Q_CYCLE, FREQ);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (!approxEq(acc.Pij_t(i, j, 0.0), expected, 1e-10)) {
                std::cerr << "FAIL testIdentityAtZero: P(" << i << "," << j << ",0) = "
                          << acc.Pij_t(i, j, 0.0) << "\n";
                return false;
            }
        }
    std::cout << "PASS testIdentityAtZero\n";
    return true;
}

// Test 2: rows of P(t) sum to 1 for several t values
static bool testRowSums() {
    nonRevAccelerator<4> acc(Q_CYCLE, FREQ);
    for (double t : { 0.1, 0.5, 1.0, 5.0, 10.0 }) {
        for (int i = 0; i < 4; ++i) {
            double sum = 0.0;
            for (int j = 0; j < 4; ++j)
                sum += acc.Pij_t(i, j, t);
            if (!approxEq(sum, 1.0, 1e-9)) {
                std::cerr << "FAIL testRowSums: row " << i << " at t=" << t
                          << " sums to " << sum << "\n";
                return false;
            }
        }
    }
    std::cout << "PASS testRowSums\n";
    return true;
}

// Test 3: all entries of P(t) are non-negative
static bool testNonNegative() {
    nonRevAccelerator<4> acc(Q_CYCLE, FREQ);
    for (double t : { 0.01, 0.5, 2.0, 10.0 }) {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) {
                double p = acc.Pij_t(i, j, t);
                if (p < -1e-10) {
                    std::cerr << "FAIL testNonNegative: P(" << i << "," << j
                              << ",t=" << t << ") = " << p << "\n";
                    return false;
                }
            }
    }
    std::cout << "PASS testNonNegative\n";
    return true;
}

// Test 4: asymmetry — P(i,j,t) != P(j,i,t) for off-diagonal, confirming
//         non-reversibility is preserved (not silently symmetrised)
static bool testAsymmetry() {
    nonRevAccelerator<4> acc(Q_CYCLE, FREQ);
    const double t = 0.5;
    // A->C should differ from C->A for this cycle at t=0.5
    double pAC = acc.Pij_t(0, 1, t);
    double pCA = acc.Pij_t(1, 0, t);
    if (approxEq(pAC, pCA, 1e-6)) {
        std::cerr << "FAIL testAsymmetry: P(A,C) == P(C,A) = " << pAC << "\n";
        return false;
    }
    std::cout << "PASS testAsymmetry  [P(A->C)=" << pAC << ", P(C->A)=" << pCA << "]\n";
    return true;
}

// Test 5: known reference value via Chapman-Kolmogorov: P(t1+t2) == P(t1)*P(t2)
static bool testChapmanKolmogorov() {
    nonRevAccelerator<4> acc(Q_CYCLE, FREQ);
    const double t1 = 0.3, t2 = 0.7;

    double P1[4][4], P2[4][4], Psum[4][4];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            P1[i][j]   = acc.Pij_t(i, j, t1);
            P2[i][j]   = acc.Pij_t(i, j, t2);
            Psum[i][j] = acc.Pij_t(i, j, t1 + t2);
        }

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            double product = 0.0;
            for (int k = 0; k < 4; ++k)
                product += P1[i][k] * P2[k][j];
            if (!approxEq(product, Psum[i][j], 1e-8)) {
                std::cerr << "FAIL testChapmanKolmogorov: P(" << i << "," << j
                          << ") product=" << product << " direct=" << Psum[i][j] << "\n";
                return false;
            }
        }
    std::cout << "PASS testChapmanKolmogorov\n";
    return true;
}

// Test 6: stationary distribution recovered from left eigenvector of Q
//   (eigenvalue closest to 0), normalised to sum to 1.
//   For the symmetric cycle Q the stationary is [0.25, 0.25, 0.25, 0.25].
static bool testStationaryFromEigen() {
    using MatNNd = Eigen::Matrix<double, 4, 4>;

    MatNNd Q;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            Q(i, j) = Q_CYCLE[i * 4 + j];

    // Left eigenvectors are right eigenvectors of Q^T
    Eigen::EigenSolver<MatNNd> solver(Q.transpose(), /*computeEigenvectors=*/true);

    // Find eigenvector with eigenvalue closest to 0
    auto eigenvalues = solver.eigenvalues();
    int    zeroIdx   = 0;
    double minAbs    = std::abs(eigenvalues[0]);
    for (int k = 1; k < 4; ++k) {
        double a = std::abs(eigenvalues[k]);
        if (a < minAbs) { minAbs = a; zeroIdx = k; }
    }

    // Extract, take real part (imaginary ~0 for valid Q), normalise
    auto   col = solver.eigenvectors().col(zeroIdx);
    double stationary[4], sum = 0.0;
    for (int i = 0; i < 4; ++i) { stationary[i] = col[i].real(); sum += stationary[i]; }
    for (int i = 0; i < 4; ++i)   stationary[i] /= sum;

    const double expected = 0.25;
    for (int i = 0; i < 4; ++i) {
        if (!approxEq(stationary[i], expected, 1e-9)) {
            std::cerr << "FAIL testStationaryFromEigen: stationary[" << i
                      << "] = " << stationary[i] << ", expected " << expected << "\n";
            return false;
        }
    }
    std::cout << "PASS testStationaryFromEigen\n";
    return true;
}

int main() {
    bool ok = true;
    ok &= testIdentityAtZero();
    ok &= testRowSums();
    ok &= testNonNegative();
    ok &= testAsymmetry();
    ok &= testChapmanKolmogorov();
    ok &= testStationaryFromEigen();
    return ok ? 0 : 1;
}