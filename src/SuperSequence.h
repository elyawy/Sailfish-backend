#include <cstdio>
#include <list>
#include <iterator> 
#include <cstddef> 
#include <limits>

#include "BlockTree.h"

#include "BlockTreeWithRates.h"
#include "SimulationContext.h"

// TODO: add the rateCategory assignment that should be passed from the Sequence object.

template<typename RngType = std::mt19937_64, typename BlockTreeType = BlockTree>
class SuperSequence {
public:
    struct columnContainer {
    // Phase 1 (build):  raw sequence position (1…N), never 0, never reaches 2^56
    // Phase 2 (output): overwritten by setAbsolutePositions() to hold absolute MSA index
    // High bits encode flags — safe because biological positions never approach 2^56
    //
    //  bit 63     : COLUMN_BIT  — this site appears in the MSA
    //  bits 55-62 : RATE_BITS   — rate category (0-255)
    //  bits  0-55 : position / absolutePosition

        uint64_t _data = 0;

        static constexpr uint64_t COLUMN_BIT = uint64_t(1) << 63;
        static constexpr uint64_t RATE_SHIFT = 55;
        static constexpr uint64_t RATE_MASK  = uint64_t(0xFF) << 55;
        static constexpr uint64_t POS_MASK   = ~(COLUMN_BIT | RATE_MASK);

        // Constructors
        columnContainer() = default;
        explicit columnContainer(uint64_t pos)
            : _data(pos & POS_MASK) {}
        columnContainer(uint64_t pos, uint8_t rate, bool isCol)
            : _data((pos & POS_MASK)
                | (uint64_t(rate) << RATE_SHIFT)
                | (isCol ? COLUMN_BIT : 0)) {}

        // Position accessors — valid meaning changes after setAbsolutePositions()
        uint64_t position()         const { return _data & POS_MASK; }
        void   setPosition(uint64_t pos)  { _data = (_data & ~POS_MASK) | (pos & POS_MASK); }

        // absolutePosition is an alias for position() post-setAbsolutePositions()
        uint64_t absolutePosition() const { return _data & POS_MASK; }

        // isColumn flag
        bool isColumn()           const { return _data & COLUMN_BIT; }
        void setColumn()                { _data |= COLUMN_BIT; }

        // Rate category (if enabled)
        uint8_t rateCategory()    const { return (_data & RATE_MASK) >> RATE_SHIFT; }
        void setRateCategory(uint8_t c) {
            _data = (_data & ~RATE_MASK) | (uint64_t(c) << RATE_SHIFT);
        }
    };
    
    using SequenceType = std::list<columnContainer>;

private:
    SequenceType _sequence;
    std::vector<typename SequenceType::iterator> _positionToIterator;
    size_t _nextSiteCounter;
    size_t _leafNum;
    size_t _numSequences;
    size_t _msaSeqLength;
    size_t _originalSequenceSize;

    BlockTreeType _blocks;
    RngType & _rng;

    CategorySampler* _rateCategorySampler;
    std::shared_ptr<std::vector<uint8_t>> _msaRateCategories;

public:
    SuperSequence<RngType, BlockTreeType>(size_t sequenceSize, SimulationContext<RngType> &simContext):
         _rng(simContext.getRng()),
         _rateCategorySampler(simContext.getCategorySampler())
         {
        _originalSequenceSize = sequenceSize;
        _msaSeqLength = 0;
        _leafNum = 0;
        _numSequences = simContext.getNumberOfNodesToSave();
        _positionToIterator.resize(sequenceSize + 1);

        for (size_t i = 1; i <= sequenceSize; ++i) {
            columnContainer column(i); // = {i, std::numeric_limits<size_t>::max(), false};
            _sequence.push_back(column);
            _positionToIterator[i] = std::prev(_sequence.end());
        }
        _nextSiteCounter = sequenceSize + 1;
    }

    void referencePosition(typename SequenceType::iterator position) {
        // if (position->position == 0) return;
        if (!(*position).isColumn()) {
            (*position).setColumn();
            ++_msaSeqLength;
        }
    }

    std::vector<size_t> setAbsolutePositions() {
        if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
            _msaRateCategories = std::make_shared<std::vector<uint8_t>>(_msaSeqLength);
        }
        std::vector<size_t> rootPositions(_originalSequenceSize, SIZE_MAX);
        size_t i = 0;
        for (auto &column: _sequence) {
            if (!column.isColumn()) continue;
            // column.absolutePosition = i;
            size_t currentPosition = column.position();
            column.setPosition(i);
            if (currentPosition <= _originalSequenceSize) {
                rootPositions[currentPosition - 1] = i;
            }
            if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
               (*_msaRateCategories)[i] = column.rateCategory();
            }
            ++i;

        }
        return rootPositions;
    }

    typename SequenceType::iterator insertItemAtPosition(typename SequenceType::iterator position, size_t item, bool isToSave) {
        // std::cout << "INSERT POS: " << *position << " " << item << "\n";
        // printSequence();
        // columnContainer newColumn = {item, std::numeric_limits<size_t>::max(), false};
        columnContainer newColumn(item);

        if (isToSave) {
            newColumn.setColumn();
            ++_msaSeqLength;
        }
        auto inserted_iterator = _sequence.insert(position ,newColumn);
        _positionToIterator.push_back(inserted_iterator);

        return inserted_iterator;
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


    typename SequenceType::iterator getIteratorByPosition(size_t position) {
        return _positionToIterator[position];
    }


    void printSequence() {
        for (auto &item: _sequence) {
            std::cout << item.position()  << " ";
        }
        std::cout << "\n";
    }

    bool checkSequenceValidity() {
        
        for (size_t i = 1; i < _nextSiteCounter; i++)
        {
            size_t numberOfAppearances = 0;
            for (auto j: _sequence) {
                if (i==j.position()) numberOfAppearances++;

            }
            if (numberOfAppearances!=1) {
                std::cout << "position " << i << " appears " << numberOfAppearances << " times\n";
                return false;
            }
        }
        return true;
    }

    void initBlockTree(size_t seqLength){ 
            _blocks.initTree(seqLength);
    }


    void initBlockTree(size_t seqLength, const std::vector<uint8_t> &rootRates){ 
        _blocks.initTree(seqLength, rootRates);
    }


    void logEventInBlockTree(Event ev) {
        // Conditional handling based on BlockTreeType
        if constexpr (std::is_same_v<BlockTreeType, BlockTreeWithRates>) {
            _blocks.handleEvent(ev, *_rateCategorySampler, _rng);
        } else {
            // Simple BlockTree doesn't need sampler/rng
            _blocks.handleEvent(ev.type, ev.position, ev.length);
        }
    }

    BlockTreeType& getBlockTree(){ return _blocks;}


    size_t sampleRootCategory() {
        if (_rateCategorySampler == nullptr) {
            errorMsg::reportError("Category sampler is nullptr, Please assign a valid sampler to the simulation context\n");
        }
        return _rateCategorySampler->drawSample(_rng);
    }

    std::shared_ptr<std::vector<uint8_t>> getMsaRateCategories() const {
        return _msaRateCategories;  // shared_ptr copy, increments ref count
    }



    ~SuperSequence() {
        _sequence.clear();
    };
};