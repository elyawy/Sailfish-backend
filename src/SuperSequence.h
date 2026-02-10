#include <cstdio>
#include <list>
#include <iterator> 
#include <cstddef> 
#include <limits>

#include "BlockTreeWithRates.h"

template<typename RngType = std::mt19937_64>
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
    std::vector<SequenceType::iterator> _positionToIterator;
    size_t _nextSiteCounter;
    size_t _leafNum;
    size_t _numSequences;
    size_t _msaSeqLength;
    BlockTreeWithRates _blocks;
    RngType & _rng;
    CategorySampler _rateCategorySampler;

public:
    SuperSequence(size_t sequenceSize, size_t numSequences, RngType& rng):
         _rng(rng), _rateCategorySampler({{{1.0}}}, {1.0})  // Single category, no heterogeneity
         {
        _msaSeqLength = 0;
        _leafNum = 0;
        _numSequences = numSequences;
        _positionToIterator.resize(sequenceSize + 1);

        for (size_t i = 1; i <= sequenceSize; ++i) {
            columnContainer column = {i, std::numeric_limits<size_t>::max(), false};
            _sequence.push_back(column);
            _positionToIterator[i] = std::prev(_sequence.end());
        }
        _nextSiteCounter = sequenceSize + 1;
    }

    void referencePosition(SequenceType::iterator position) {
        // if (position->position == 0) return;
        if (!(*position).isColumn) {
            (*position).isColumn = true;
            ++_msaSeqLength;
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

    SequenceType::iterator insertItemAtPosition(SequenceType::iterator position, size_t item, bool isToSave) {
        // std::cout << "INSERT POS: " << *position << " " << item << "\n";
        // printSequence();
        columnContainer newColumn = {item, std::numeric_limits<size_t>::max(), false};

        if (isToSave) {
            newColumn.isColumn = true;
            ++_msaSeqLength;
        }
        auto inserted_iterator = _sequence.insert(position ,newColumn);
        _positionToIterator.push_back(inserted_iterator);

        return inserted_iterator;
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
        return _nextSiteCounter;
    }


    size_t incrementRandomSequencePosition() {
        size_t positionToReturn = ++_nextSiteCounter;
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

    size_t getMsaSequenceLength() {
        return _msaSeqLength;
    }


    SequenceType::iterator getIteratorByPosition(size_t position) {
        return _positionToIterator[position];
    }


    void printSequence() {
        for (auto &item: _sequence) {
            std::cout << item.position  << " ";
        }
        std::cout << "\n";
    }

    bool checkSequenceValidity() {
        
        for (size_t i = 1; i < _nextSiteCounter; i++)
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


    void initBlockTree(size_t seqLength){ _blocks.initTree(seqLength);}
    void logEventInBlockTree(Event ev) { _blocks.handleEvent(ev, _rateCategorySampler, _rng)}

    BlockTree& getBlockTree(){ return _blocks;}

    void initRateSampler(const std::vector<std::vector<MDOUBLE>>& transitionMatrix,
                     const std::vector<MDOUBLE>& stationaryProbs) {
        _rateCategorySampler = CategorySampler(transitionMatrix, stationaryProbs);
    }



    ~SuperSequence() {
        _sequence.clear();
    };
};