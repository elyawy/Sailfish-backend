#ifndef _SEQUENCE
#define _SEQUENCE


#include <stddef.h>
#include <iostream>
#include <vector>
#include <list>
#include <numeric>
#include <iterator>

#include "SuperSequence.h"
#include "BlockTree.h"

class Sequence
{
    using iteratorType = std::list<SuperSequence::columnContainer>::iterator;
    using SequenceType = std::vector<iteratorType>;

private:
    SuperSequence* _superSequence;
    bool _isSaveSequence;
    size_t _nodeID;
    SequenceType _sequence;
    Sequence* _parent;

    // size_t _numLeaf;
public:

    Sequence(SuperSequence& superSeq, bool isSaveSeq, size_t nodeID) : 
        _superSequence(&superSeq), _isSaveSequence(isSaveSeq), _nodeID(nodeID) {}

    Sequence(const Sequence &seq) {
        for (size_t i = 0; i < seq._sequence.size(); i++) {
            _sequence.push_back(seq._sequence[i]);
        }
        _superSequence = seq._superSequence;
        _isSaveSequence = seq._isSaveSequence;
        _nodeID = seq._nodeID;
    }

    void initSequence() {
        auto superSeqIterator = _superSequence->begin();

        while (superSeqIterator != _superSequence->end()) {
            if (_isSaveSequence) _superSequence->referencePosition(superSeqIterator);
            _sequence.push_back(superSeqIterator);
            superSeqIterator++;
        }
    }

    void generateSequence (const BlockList &blocklist, Sequence &parentSeq) {

        size_t position;
        size_t length;
        size_t insertion;
        size_t randomPos = _superSequence->getRandomSequencePosition();
        _parent = &parentSeq;
        // std::cout << parentSeq.getSequenceNodeID() << "\n";
        for (auto it = blocklist.begin(); it != blocklist.end(); ++it) {
            position = (*it)[static_cast<int>(BLOCK::POSITION)];//(&it)->key();
            length = (*it)[static_cast<int>(BLOCK::LENGTH)];//(*it).length;
            insertion = (*it)[static_cast<int>(BLOCK::INSERTION)];//(*it).insertion;

            // std::cout << "current Block is: " << position <<"|" << length << "|" << insertion << "\n";

            if (position==0 && length==1 && insertion==0) {
                // _sequence.push_back(parentSeq._sequence[0]);
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
                superSeqIterator++;
                randomPos = _superSequence->incrementRandomSequencePosition();
            }
        }

        if (_isSaveSequence) _superSequence->incrementLeafNum();
    }

    SuperSequence* getSuperSequence() {
        return _superSequence;
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

    size_t getSequenceNodeID() {
        return _nodeID;
    }


    void clear() {
        _sequence.clear();
    }



    ~Sequence() {
    }
};


#endif 
