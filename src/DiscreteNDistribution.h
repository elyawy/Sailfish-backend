#pragma once

#include <random>
#include <utility>
#include <vector>
#include <stack>
#include <iostream>
#include <array>
#include <cassert>


template<uint8_t K>
struct SplitRNG {
    // Dynamically asserts any power of two (K=4, 16, 32, 64, 128, 256.)
    static_assert((K > 0) && ((K & (K - 1)) == 0), "K must be a power of two");
    
    static constexpr uint8_t ctz(uint64_t v) {
        uint8_t count = 0;
        while ((v & 1) == 0) { v >>= 1; ++count; }
        return count;
    }

    static constexpr uint64_t int_mask = K - 1;  
    static constexpr uint8_t bits_for_int = ctz(K);  
    static constexpr double divisor = 1ULL << (64 - bits_for_int);  
    
    uint8_t integer; // upgraded for N <= 256
    double uniform;
    
    explicit SplitRNG(uint64_t rng_val) {
        integer = static_cast<uint8_t>(rng_val & int_mask); 
        uniform = (rng_val >> bits_for_int) / divisor;  
    }
};

template<uint8_t N>
class DiscreteNDistribution
{
private:
    std::array<double, N> probabilities_;
    std::array<uint8_t, N> alias_; 

public:
    DiscreteNDistribution<N>(const std::array<double, N> &probabilities, double normalizingFactor=1.0) {
        assert(probabilities.size() == N);

        std::pair<int, double> small_buf[N];
        std::pair<int, double> large_buf[N];
        int small_top = 0, large_top = 0;

        for(int i = 0; i < N; i++) {
            double scaled_prob = N * probabilities[i] * normalizingFactor;
            if (scaled_prob < 1.0) {
                small_buf[small_top++] = {i, scaled_prob};
            } else {
                large_buf[large_top++] = {i, scaled_prob};
            }
        }

        while (small_top > 0 && large_top > 0) {
            auto s = small_buf[--small_top];
            auto l = large_buf[--large_top];

            probabilities_[s.first] = s.second;
            alias_[s.first] = static_cast<uint8_t>(l.first);

            l.second = (l.second + s.second) - 1.0;
            if (l.second < 1.0) {
                small_buf[small_top++] = l;
            } else {
                large_buf[large_top++] = l;
            }
        }

        while (large_top > 0) {
            auto l = large_buf[--large_top];
            probabilities_[l.first] = 1.0;
        }
        while (small_top > 0) {
            auto s = small_buf[--small_top];
            probabilities_[s.first] = 1.0;
        }
    }

    template<typename RngType = std::mt19937_64>
    int drawSample(RngType &rng) const {
        uint64_t full64Bits = rng();
        SplitRNG<N> split(full64Bits);

        uint8_t die_roll = split.integer;  
        double coin_flip = split.uniform;  
        
        int result = alias_[die_roll];
        int keep_original = (coin_flip < probabilities_[die_roll]); 

        int mask = -keep_original; 
        result = (result & ~mask) | (die_roll & mask); 

        return result + 1;
    }

    void printTable() {
        for(auto &i: probabilities_){
            std::cout << i << " ";
        }
        std::cout << "\n";
        for(auto &i: alias_){
            std::cout << i << " ";
        }
        std::cout << "\n";
    }

    std::vector<std::pair<double, int>> getTable() {
        std::vector<std::pair<double, int>> prob_alias_table;
        for(int i=0; i < probabilities_.size(); i++){
            std::pair<double, int> current_item (probabilities_[i], alias_[i]);
            prob_alias_table.push_back(current_item);
        }
        return prob_alias_table;
    }

    double getProb(int i) {
        return probabilities_[i];
    }
};