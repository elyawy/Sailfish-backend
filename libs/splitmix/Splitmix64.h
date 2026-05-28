#include <cstdint>
#include <cmath>



class Splitmix64 {
public:

    using result_type = uint64_t;
    static constexpr uint64_t min() { return 0; }
    static constexpr uint64_t max() { return UINT64_MAX; }

	Splitmix64() { state = 0; }
	Splitmix64(const uint64_t seed) : state(seed) {}

	void seed(const uint64_t seed) {
		state = seed;
	}

	uint64_t next_int() {
		uint64_t z = ( state += 0x9e3779b97f4a7c15 );
		z = ( z ^ ( z >> 30 ) ) * 0xbf58476d1ce4e5b9;
		z = ( z ^ ( z >> 27 ) ) * 0x94d049bb133111eb;
		return z ^ ( z >> 31 );
	}

	double next_float() {
	    return next_int() / twoPower64;
	}

    uint64_t operator()() {
        return next_int();
    }
    
private:
	uint64_t state;
	const double twoPower64 = pow(2.0, 64);
};