#include <cstdio>
#include <list>
#include <iterator> 
#include <cstddef> 

class SuperSequence {
public:
    struct columnContainer {
        const size_t position;
        size_t absolutePosition;
        bool isColumn;
    };
    
    using SequenceType = std::list<columnContainer>;

private:
    SequenceType _sequence;
    size_t _randomSequenceCounter;
    size_t _leafNum;
    size_t _numSequences;
    size_t _leafSeqLength;
public:
    SuperSequence(size_t sequenceSize, size_t numSequences) {
        _leafSeqLength = 0;
        _leafNum = 0;
        _numSequences = numSequences;
        for (size_t i = 1; i <= sequenceSize; ++i) {
            columnContainer column = {i, false};
            _sequence.push_back(column);
        }
        _randomSequenceCounter = sequenceSize + 1;
    }

    void referencePosition(SequenceType::iterator position) {
        if (!(*position).isColumn) {
            (*position).isColumn = true;
            ++_leafSeqLength;
        }
    }

    void setAbsolutePositions() {
        size_t i = 0;
        for (auto &column: _sequence) {
            if (!column.isColumn) continue;
            column.absolutePosition = i;
            ++i;
        }
    }

    SequenceType::iterator insertItemAtPosition(SequenceType::iterator position, size_t item, bool isLeaf) {
        // std::cout << "INSERT POS: " << *position << " " << item << "\n";
        // printSequence();
        columnContainer newColumn = {item, false};

        if (isLeaf) {
            newColumn.isColumn = true;
            ++_leafSeqLength;
        }

        return _sequence.insert(position ,newColumn);
    }


    SequenceType::iterator begin() {
        return _sequence.begin();
    }

    SequenceType::iterator end() {
        return _sequence.end();
    }

    size_t size() {
        return _sequence.size();
    }


    size_t getRandomSequencePosition() {
        // std::cout << "get random pos: " << _randomSequenceCounter << "\n";
        return _randomSequenceCounter;
    }


    size_t incrementRandomSequencePosition() {
        size_t positionToReturn = ++_randomSequenceCounter;
        // std::cout << "inceremted random pos: " << positionToReturn << "\n";
        return positionToReturn;
    }

    size_t incrementLeafNum() {
        size_t leafToReturn = ++_leafNum;
        // std::cout << "inceremted random pos: " << positionToReturn << "\n";
        return leafToReturn;
    }

    size_t getNumSequences() {
        return _numSequences;
    }

    size_t getLeafSequenceLength() {
        return _leafSeqLength;
    }





    void printSequence() {
        for (auto &item: _sequence) {
            std::cout << item.position  << " ";
        }
        std::cout << "\n";
    }

    bool checkSequenceValidity() {
        
        for (size_t i = 1; i < _randomSequenceCounter; i++)
        {
            size_t numberOfAppearances = 0;
            for (auto j: _sequence) {
                if (i==j.position) numberOfAppearances++;

            }
            if (numberOfAppearances!=1) {
                std::cout << "position " << i << " appears " << numberOfAppearances << " times\n";
                return false;
            }
        }
        return true;
    }


    ~SuperSequence() {
        _sequence.clear();
    };
};