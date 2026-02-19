#ifndef _SEQUENCE
#define _SEQUENCE


#include <stddef.h>
#include <iostream>
#include <vector>
#include <list>
#include <numeric>
#include <iterator>

#include "SuperSequence.h"
#include "Event.h"

// TODO: retrieve the rateCategories from the Blocks following the application of the events.
// assign the rateCategories to the column within the SuperSequence.


struct CompressedSequence {
    std::vector<std::pair<size_t, size_t>> runs; // (start_position, length)
    size_t uncompressedSize;
};

template<typename RngType = std::mt19937_64, typename BlockTreeType = BlockTree>
class Sequence
{
    using SuperSeqType = SuperSequence<RngType, BlockTreeType>;
    using iteratorType = typename std::list<typename SuperSeqType::columnContainer>::iterator;
    using SequenceType = std::vector<iteratorType>;

private:
    SuperSeqType* _superSequence;
    bool _isSaveSequence;

    SequenceType _sequence;
    const Sequence* _parent;
    std::vector<size_t> _rateCategories;

public:

    Sequence(SuperSeqType& superSeq, bool isSaveSeq) : 
        _superSequence(&superSeq), _isSaveSequence(isSaveSeq) {}


    Sequence(const CompressedSequence& compressed, SuperSeqType& superSeq) 
        : _superSequence(&superSeq), _isSaveSequence(true) {
        
        _sequence.reserve(compressed.uncompressedSize);
        
        for (const auto& [start, length] : compressed.runs) {
            for (size_t i = 0; i < length; ++i) {
                size_t position = start + i;
                auto it = _superSequence->getIteratorByPosition(position);
                _sequence.push_back(it);
            }
        }
    }

    void initSequence() {
        auto superSeqIterator = _superSequence->begin();

        if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
            _rateCategories.reserve(_superSequence->size());
        }

        while (superSeqIterator != _superSequence->end()) {
            if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
                // Sample and assign rate category to the column
                size_t category = _superSequence->sampleRootCategory();
                (*superSeqIterator).rateCategory = category;
                _rateCategories.push_back(category);
            }

            if (_isSaveSequence) _superSequence->referencePosition(superSeqIterator);
            _sequence.push_back(superSeqIterator);
            superSeqIterator++;
        }
    }


    void generateSequence (const EventSequence &eventlist,const Sequence *parentSeq) {
        _sequence.reserve(parentSeq->_sequence.size());

        if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
            _rateCategories.reserve(parentSeq->_sequence.size());
        }

        size_t position;
        size_t length;
        size_t insertion;
        size_t randomPos = _superSequence->getRandomSequencePosition();
        _parent = (parentSeq);



        // apply events on BlockTree
        // Initialize BlockTree with parent's rate categories
        if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
            const std::vector<size_t>& parentRates = parentSeq->_rateCategories;
            _superSequence->initBlockTree(parentSeq->_sequence.size(), parentRates);
        } else {
            _superSequence->initBlockTree(parentSeq->_sequence.size());
        }


        for (const auto& ev: eventlist) {
            _superSequence->logEventInBlockTree(ev);
        }

        auto& blocks  = _superSequence->getBlockTree();

        for (auto it = blocks.begin(); it != blocks.end(); ++it) {
            position = it.key();
            length = (*it).length;
            insertion = (*it).insertion;


            if (position==0 && length==1 && insertion==0) {
                continue;
            }


            if (position!=0) {
                position--;
            } else {
                length--;
            }

            for (size_t i = 0; i < length; i++) {
                if (_isSaveSequence) {
                    _superSequence->referencePosition(_parent->_sequence[position+i]);
                } 
                _sequence.push_back(_parent->_sequence[position+i]);
                if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
                    _rateCategories.push_back(_parent->_rateCategories[position+i]);
                }
            }
            while (_parent->_sequence.size() == 0) _parent = _parent->_parent;

            auto superSeqIterator = _parent->_sequence[position];
            if (!_sequence.empty()) {
                superSeqIterator = _parent->_sequence[position+length-1];
                superSeqIterator++;
            }
            

            for (size_t i = 0; i < insertion; i++) {
                superSeqIterator = _superSequence->insertItemAtPosition(superSeqIterator, randomPos, _isSaveSequence);
                _sequence.push_back(superSeqIterator);

                if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
                    size_t category = (*it).rateCategories[i];
                    (*superSeqIterator).rateCategory = category;
                    _rateCategories.push_back(category);
                }

                superSeqIterator++;
                randomPos = _superSequence->incrementRandomSequencePosition();
            }
        }

        if (_isSaveSequence) _superSequence->incrementLeafNum();
    }

    SuperSeqType* getSuperSequence() {
        return _superSequence;
    }


    typename SequenceType::iterator begin() {
        return _sequence.begin();
    }

    typename SequenceType::iterator end() {
        return _sequence.end();
    }

    size_t size() {
        return _sequence.size();
    }

    iteratorType getPos(size_t pos) {
        return _sequence[pos];
    }


    void printSequence() {
        for(auto &item: _sequence) {
            std::cout << (*item).position << " ";
        }
        std::cout << "\n";
    }

    bool checkSequenceValidity() {
        size_t maxSeqSize = _superSequence->getRandomSequencePosition();
        for (size_t i = 1; i < maxSeqSize; i++)
        {
            size_t numberOfAppearances = 0;
            for (auto j: _sequence) {
                if (i==(*j).position) numberOfAppearances++;

            }
            if (numberOfAppearances > 1) {
                std::cout << "position " << i << " appears " << numberOfAppearances << " times\n";
                return false;
            }
        }
        return true;
    }



    CompressedSequence compress() const {
        CompressedSequence result;
        result.uncompressedSize = _sequence.size();
        result.runs.reserve(_sequence.size() / 10); // Reserve space assuming average run length of 10
        if (_sequence.empty()) return result;
        
        size_t start = (_sequence[0])->position;
        size_t count = 1;
        
        for (size_t i = 1; i < _sequence.size(); ++i) {
            size_t currentPos = (_sequence[i])->position;
            size_t prevPos = (_sequence[i-1])->position;
            
            if (currentPos == prevPos + 1) {
                // Consecutive, extend current run
                count++;
            } else {
                // Non-consecutive, save current run and start new one
                result.runs.push_back({start, count});
                start = currentPos;
                count = 1;
            }
        }
        
        // Don't forget the last run
        result.runs.push_back({start, count});
        
        return result;
    }


    void clear() {
        _sequence.clear();
    }



    ~Sequence() {
    }
};


#endif 
