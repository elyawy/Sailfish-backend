// #include "../../src/encodedSequence.h"
#include <iostream>
#include <climits>

#include <bitset>

int main() {

    unsigned char message[8] = {0,2,3,2,3,3,0,1};

    uint64_t bitRep = 0;

    for (size_t i = 1; i <= 8; i++)
    {
        uint64_t bits(message[i-1]);
        // std::cout << (8*8 - i*5) << "\n";
        // std::cout << (bits << (8*8 - i*5)) << "\n";

        bitRep |= (bits << (8*8 - i*2));
        std::cout << std::bitset<64>(bitRep) << "\n";
    }
    
    bitRep >>= 48;
    bitRep >>= 8;

    uint8_t i = bitRep; 

    std::cout << std::bitset<64>(bitRep) << "\n";
    std::cout << i << "\n";


    return 0;



}