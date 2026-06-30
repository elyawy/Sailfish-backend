#ifndef ___NONREV_ACCELERATOR
#define ___NONREV_ACCELERATOR

#include <complex>
#include <cmath>
#include <limits>
#include <vector>
#include <cassert>

#include "../libs/Phylolib/includes/pijAccelerator.h"
#include "../libs/Phylolib/includes/definitions.h"

#include "../libs/Eigen/Eigenvalues"

// Non-reversible P(t) accelerator using Eigen decomposition.
//
// Accepts a general (non-reversible) Q matrix and a root frequency vector
// directly — no replacementModel required.
//
// At construction: decomposes Q into V, λ, V⁻¹ such that
//   P(t) = V · diag(exp(λ·t)) · V⁻¹
// where V, λ are complex (general Q has no guarantee of real eigenvalues).
//
// Per distinct t:
//   1. Scale columns of V by exp(λk · t)   (N complex exp calls)
//   2. P = V_scaled · V⁻¹                  (N×N complex matrix multiply)
//   3. Take real part — imaginary parts cancel for valid Q matrices
//
// Eigenvectors/values stored as flat complex arrays (stack-friendly for small N).
// Full P(t) matrix cached; Pij_t(i,j,t) refills only when t changes.
// N is a compile-time constant, enabling loop unrolling and auto-vectorisation.

template<int N>
class nonRevAccelerator : public pijAccelerator {
    using Cplx   = std::complex<double>;
    using MatNNc = Eigen::Matrix<Cplx, N, N>;
    using MatNNd = Eigen::Matrix<double, N, N>;

public:
    // qMatrix:     flat row-major N×N rate matrix (rows sum to 0)
    // frequencies: root frequencies, size N (need not be stationary)
    nonRevAccelerator(const std::vector<double>& qMatrix,
                      const std::vector<double>& frequencies)
        : _cachedT(std::numeric_limits<MDOUBLE>::quiet_NaN())
    {
        assert((int)qMatrix.size()    == N * N);
        assert((int)frequencies.size() == N);

        for (int i = 0; i < N; ++i) {
            _freq[i] = frequencies[i];
            for (int j = 0; j < N; ++j)
                _Q[i][j] = qMatrix[i * N + j];
        }

        buildDecomp(qMatrix);
    }

    nonRevAccelerator(const nonRevAccelerator& other)
        : _cachedT(other._cachedT)
    {
        for (int i = 0; i < N; ++i) {
            _freq[i]    = other._freq[i];
            _lambda[i]  = other._lambda[i];
            for (int j = 0; j < N; ++j) {
                _Q[i][j]    = other._Q[i][j];
                _V[i][j]    = other._V[i][j];
                _Vinv[i][j] = other._Vinv[i][j];
                _mat[i][j]  = other._mat[i][j];
            }
        }
    }

    // --- pijAccelerator interface ---

    inline const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const override {
        if (t != _cachedT) refill(t);
        return _mat[i][j];
    }

    const MDOUBLE Qij(const int i, const int j) const override {
        return _Q[i][j];
    }

    const MDOUBLE freq(const size_t i) const override {
        return _freq[i];
    }

    const size_t alphabetSize() const override { return N; }

    pijAccelerator* clone() const override { return new nonRevAccelerator(*this); }

    // Not used in simulation path — stubs satisfy pure-virtual contract
    const MDOUBLE dPij_dt  (const int, const int, const MDOUBLE) const override { return 0.0; }
    const MDOUBLE d2Pij_dt2(const int, const int, const MDOUBLE) const override { return 0.0; }

    // No underlying replacementModel
    replacementModel* getReplacementModel() const override { return nullptr; }

    virtual ~nonRevAccelerator() = default;

private:
    void buildDecomp(const std::vector<double>& qMatrix)
    {
        MatNNd Q;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                Q(i, j) = qMatrix[i * N + j];

        Eigen::EigenSolver<MatNNd> solver(Q, /*computeEigenvectors=*/true);

        // Eigenvalues and right eigenvector matrix
        auto eigenvalues  = solver.eigenvalues();   // N complex values
        MatNNc V          = solver.eigenvectors();  // columns are right eigenvectors
        MatNNc Vinv       = V.inverse();

        for (int k = 0; k < N; ++k) {
            _lambda[k] = eigenvalues[k];
            for (int i = 0; i < N; ++i) {
                _V[i][k]    = V(i, k);
                _Vinv[k][i] = Vinv(k, i);
            }
        }
    }

    void refill(const MDOUBLE t) const
    {
        // Step 1: compute exp(λk · t) for each eigenvalue
        Cplx expLambda[N];
        for (int k = 0; k < N; ++k)
            expLambda[k] = std::exp(_lambda[k] * t);

        // Step 2: P = V · diag(exp(λ·t)) · V⁻¹, take real part
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                Cplx sum = 0.0;
                for (int k = 0; k < N; ++k)
                    sum += _V[i][k] * expLambda[k] * _Vinv[k][j];
                _mat[i][j] = sum.real();
            }
        }

        _cachedT = t;
    }

    // Root frequencies
    double _freq[N];

    // Raw Q matrix (for Qij queries)
    double _Q[N][N];

    // Eigen decomposition — flat complex arrays
    Cplx _lambda[N];     // eigenvalues
    Cplx _V[N][N];       // right eigenvectors: V[i][k]
    Cplx _Vinv[N][N];    // inverse:            Vinv[k][j]

    // Cached real P(t) and the t it was computed for
    mutable double  _mat[N][N];
    mutable MDOUBLE _cachedT;
};

#endif // ___NONREV_ACCELERATOR
