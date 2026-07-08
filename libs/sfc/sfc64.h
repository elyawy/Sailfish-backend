#pragma once

// SFC64 — header-only C++ wrapper
//
// Core algorithm adapted from the NumPy project:
//   https://github.com/numpy/numpy/blob/main/numpy/random/src/sfc64/sfc64.h
// Original algorithm by Chris Doty-Humphrey.
// NumPy is BSD-licensed; see https://github.com/numpy/numpy/blob/main/LICENSE.txt
//
// This wrapper adds:
//   - Self-contained seeding (no external .c file needed)
//   - Full UniformRandomBitGenerator + standard RNG interface compliance
//     (drop-in for std::mt19937_64, pcg64_fast, etc.)
//   - seed_seq constructor and seed() method for template compatibility

#include <cstdint>
#include <random>   // for std::seed_seq
#include <array>

// ---------------------------------------------------------------------------
// Internal state and core step function (verbatim from NumPy)
// ---------------------------------------------------------------------------

struct sfc64_state {
    uint64_t s[4];
    int      has_uint32;
    uint32_t uinteger;
};

namespace sfc64_detail {

inline uint64_t rotl(uint64_t value, unsigned rot) {
#ifdef _WIN32
    return _rotl64(value, rot);
#else
    return (value << rot) | (value >> ((-rot) & 63));
#endif
}

inline uint64_t next(uint64_t *s) {
    const uint64_t tmp = s[0] + s[1] + s[3]++;
    s[0] = s[1] ^ (s[1] >> 11);
    s[1] = s[2] + (s[2] << 3);
    s[2] = rotl(s[2], 24) + tmp;
    return tmp;
}

} // namespace sfc64_detail

// ---------------------------------------------------------------------------
// C++ wrapper — satisfies UniformRandomBitGenerator and the broader
// named requirement RandomNumberEngine (constructible/seedable the same
// ways as std::mt19937_64 and pcg64_fast)
// ---------------------------------------------------------------------------

class SFC64 {
public:
    using result_type = uint64_t;

    // --- Constructors -------------------------------------------------------

    // Default-constructed with a fixed seed (matches std:: RNG behaviour)
    SFC64() { init(0, 0, 0); }

    // Single integer seed — the form most code uses
    explicit SFC64(uint64_t seed) { init(seed, seed, seed); }

    // seed_seq constructor — lets SimulationContext<SFC64> do things like
    //   std::seed_seq seq{...};  RNG rng(seq);
    // which pcg64_fast and mt19937_64 both support.
    explicit SFC64(std::seed_seq &seq) { seed(seq); }

    // Three-lane constructor — best entropy, avoids the correlated early
    // outputs that result from starting all lanes with the same value.
    SFC64(uint64_t s0, uint64_t s1, uint64_t s2) { init(s0, s1, s2); }

    // --- Seeding ------------------------------------------------------------

    void seed(uint64_t s) { init(s, s, s); }

    void seed(std::seed_seq &seq) {
        // Generate 6 × 32-bit values → 3 × 64-bit lane seeds
        std::array<uint32_t, 6> buf;
        seq.generate(buf.begin(), buf.end());
        uint64_t s0 = (static_cast<uint64_t>(buf[0]) << 32) | buf[1];
        uint64_t s1 = (static_cast<uint64_t>(buf[2]) << 32) | buf[3];
        uint64_t s2 = (static_cast<uint64_t>(buf[4]) << 32) | buf[5];
        init(s0, s1, s2);
    }

    // --- UniformRandomBitGenerator interface --------------------------------

    result_type operator()() { return sfc64_detail::next(state_.s); }

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return UINT64_MAX; }

    // --- Misc ---------------------------------------------------------------

    // Discard n outputs (e.g. for parallel streams)
    void discard(unsigned long long n) {
        for (unsigned long long i = 0; i < n; ++i)
            sfc64_detail::next(state_.s);
    }

    // Direct state access for serialisation / inspection
    const sfc64_state &state() const { return state_; }

    bool operator==(const SFC64 &o) const {
        return state_.s[0] == o.state_.s[0] &&
               state_.s[1] == o.state_.s[1] &&
               state_.s[2] == o.state_.s[2] &&
               state_.s[3] == o.state_.s[3];
    }
    bool operator!=(const SFC64 &o) const { return !(*this == o); }

private:
    sfc64_state state_{};

    void init(uint64_t s0, uint64_t s1, uint64_t s2) {
        state_.s[0]     = s0;
        state_.s[1]     = s1;
        state_.s[2]     = s2;
        state_.s[3]     = 1;   // counter lane — must start non-zero
        state_.has_uint32 = 0;
        state_.uinteger   = 0;
        // 18 warm-up rounds flush the weak early output from low-entropy seeds.
        // This matches NumPy's own initialisation behaviour.
        for (int i = 0; i < 18; ++i)
            sfc64_detail::next(state_.s);
    }
};