#include "../libs/Phylolib/includes/sequence.h"

#include <climits>

class encodedSequence: public sequence
{
public:

	void push_back(ALPHACHAR p) {
        _remainder[_remainderIndex] = p;
        ++_remainderIndex;

        uint64_t bitSequence = 0;

        if (_remainderIndex == CHAR_BIT) {
            for (size_t i = 0; i < CHAR_BIT; i++) {
                bitSequence |= ( << (_remainderIndex - i - 1))
                
            }
            _remainder.clear()
        }
        _vec.push_back(p);
    }


    ALPHACHAR& operator[](const int i) {
        return _vec[i];
    }
	const ALPHACHAR& operator[](const int pos) const {
        return _vec[pos];
    }

    string toString() const{
        string tmp;
        for (size_t k=0; k < _vec.size() ; ++k ){
            tmp+= _alphabet->fromInt(_vec[k]);
        }
        return tmp;
    }

    string toString(const int pos) const{
        return _alphabet->fromInt(_vec[pos]);
    }
private:
    ALPHACHAR _remainder[CHAR_BIT] = {};
    size_t _remainderIndex = 0;

};

