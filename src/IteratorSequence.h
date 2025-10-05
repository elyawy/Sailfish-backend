#ifndef _ITERATOR_SEQUENCE
#define _ITERATOR_SEQUENCE

#include <stddef.h>
#include <iostream>
#include <vector>
#include <list>
#include <numeric>
#include <iterator>

#include "FixedList.h"
#include "BlockTree.h"

class IteratorSequence
{
    using iteratorType = FixedList::iterator;
    using SequenceType = std::vector<iteratorType>;

private:
    FixedList* _fixedList;
    bool _isSaveSequence;
    size_t _nodeID;
    SequenceType _sequence;
    IteratorSequence* _parent;

public:
    IteratorSequence(FixedList& fixedList, bool isSaveSeq, size_t nodeID) : 
        _fixedList(&fixedList), _isSaveSequence(isSaveSeq), _nodeID(nodeID) {}


    void initSequence() {
        auto fixedListIterator = _fixedList->begin();

        while (fixedListIterator != _fixedList->end()) {
            if (_isSaveSequence) _fixedList->referencePosition(fixedListIterator);
            _sequence.push_back(fixedListIterator);
            ++fixedListIterator;
        }
    }

    void generateSequence(const BlockList &blocklist, IteratorSequence &parentSeq) {
        size_t position;
        size_t length;
        size_t insertion;
        _parent = &parentSeq;


        for (auto it = blocklist.begin(); it != blocklist.end(); ++it) {
            position = (*it)[static_cast<int>(BLOCK::POSITION)];
            length = (*it)[static_cast<int>(BLOCK::LENGTH)];
            insertion = (*it)[static_cast<int>(BLOCK::INSERTION)];



            if (position==0 && length==1 && insertion==0) {
                continue;
            }


            // Copy length elements from parent
            for (size_t i = 0; i < length; i++) {
                if (_isSaveSequence) {
                    _fixedList->referencePosition(_parent->_sequence[position+i]);
                } 
                _sequence.push_back(_parent->_sequence[position+i]);
            }

            while (_parent->_sequence.size() == 0) _parent = _parent->_parent;

            auto insertAfterIt = _parent->_sequence[position+length];
            // if (!_sequence.empty()) {
            //     insertAfterIt = _parent->_sequence[position+length-1];
            //     ++insertAfterIt;
            // }

            // Insert new positions
            for (size_t i = 0; i < insertion; i++) {
                insertAfterIt = _fixedList->insertAfter(insertAfterIt, _isSaveSequence);
                _sequence.push_back(insertAfterIt);
                // ++insertAfterIt;
            }
        }
    }

    FixedList* getFixedList() {
        return _fixedList;
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

    iteratorType getPos(size_t pos) {
        return _sequence[pos];
    }

    void printSequence() {
        for(auto &it: _sequence) {
            std::cout << *it << " ";
        }
        std::cout << "\n";
    }

    bool checkSequenceValidity() {
        size_t maxSeqSize = _fixedList->size();
        for (size_t i = 1; i < maxSeqSize; i++) {
            size_t numberOfAppearances = 0;
            for (auto j: _sequence) {
                if (i==*j) numberOfAppearances++;
            }
            if (numberOfAppearances > 1) {
                std::cout << "position " << i << " appears " << numberOfAppearances << " times\n";
                return false;
            }
        }
        return true;
    }

    size_t getSequenceNodeID() {
        return _nodeID;
    }

    void clear() {
        _sequence.clear();
    }

    ~IteratorSequence() {}
};

#endif