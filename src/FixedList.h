#ifndef _FIXEDLIST
#define _FIXEDLIST


#include <vector>

class FixedList {
private:
    // Structure of Arrays - using vectors for convenience
    std::vector<size_t> _nextIndices;        // Next pointer in the linked structure
    std::vector<size_t> _traversalPositions; // Traversal position data
    std::vector<bool> _isColumns;            // Boolean flag - NOTE: specialized template!
    
    // List management
    size_t _count = 0;                       // Number of elements currently stored
    size_t _headIndex = INVALID;             // Index of first element in logical order
    size_t _tailIndex = INVALID;             // Index of last element in logical order
    size_t _msaSeqLength = INVALID;


    static constexpr size_t INVALID = SIZE_MAX;

public:
    // Constructor - pre-allocates vectors and creates initial node
    explicit FixedList(size_t max_size)
        : _nextIndices(max_size)
        , _traversalPositions(max_size)
        , _isColumns(max_size)
        , _count(1)      // Start with one element
        , _headIndex(0)
        , _tailIndex(0)
    {
        // Initialize the first element
        _nextIndices[0] = INVALID;
        _traversalPositions[0] = INVALID;
        _isColumns[0] = false;  // Default value for first element
    }
    
    // Copy and move are handled automatically by vector
    size_t size() const { return _count; }
    size_t max_size() const { return _nextIndices.size(); }
    bool empty() const { return _count == 0; }
    bool full() const { return _count >= _nextIndices.size(); }

    void initialize(size_t sequenceSize) {
        if (sequenceSize <= 0) return;
        batchInsertAfter(_headIndex, false, sequenceSize);
    }
    

    // Insert after element at index k
    // Returns index of newly inserted element, or INVALID if failed
    inline size_t insertAfter(size_t nodeK, bool isColumn) {


        // Check if we have space
        if (full()) {
            return INVALID;
        }
        
        // Check if k is valid
        if (nodeK >= _count) {
            return INVALID;
        }
        
        // Get the new element's index (next available slot)
        size_t new_index = _count++;
        
        // Set the new element's data with default values
        _isColumns[new_index] = isColumn;

        // Update links: insert after element k
        size_t old_next = _nextIndices[nodeK];
        _nextIndices[nodeK] = new_index;        // k now points to new element
        _nextIndices[new_index] = old_next;     // new element points to what k used to point to
        
        // Update tail if we inserted after the current tail
        if (nodeK == _tailIndex) {
            _tailIndex = new_index;
        }

        
        return new_index;
    }

    // Insert before element at index nodeK
    // Returns index of newly inserted element, or INVALID if failed
    // NOTE: Cannot insert before head (node 0)
    inline size_t insertBefore(size_t nodeK, bool isColumn) {
        // Check if we have space
        if (full()) {
            return INVALID;
        }
        
        // Check if nodeK is valid
        if (nodeK >= _count) {
            return INVALID;
        }
        
        // Cannot insert before head (anchor node)
        if (nodeK == _headIndex || nodeK == 0) {
            return INVALID;
        }
        
        // Find the node that points to nodeK
        size_t prev = _headIndex;
        while (prev != INVALID && _nextIndices[prev] != nodeK) {
            prev = _nextIndices[prev];
        }
        
        // If we didn't find a predecessor, nodeK is not in the list
        if (prev == INVALID) {
            return INVALID;
        }
        
        // Now insert after prev (which is before nodeK)
        return insertAfter(prev, isColumn);
    }

    size_t batchInsertAfter(size_t nodeK, bool isColumn, int numNodes) {
        size_t currentIndex = nodeK;
        for (int i = 0; i < numNodes; ++i) {
            currentIndex = insertAfter(currentIndex, isColumn);
            if (currentIndex == INVALID) return INVALID; // Failed
        }
        return currentIndex;
    }


    bool referencePosition(size_t nodeK) {
        // Check if k is valid
        // if (nodeK >= _count) {
        //     return false;
        // }
        
        _isColumns[nodeK] = true;
        return true;
    }


    void setAbsolutePositions() {
        size_t i = 0;
        size_t current = _headIndex;
        
        while (current != INVALID) {
            if (!_isColumns[current]) {
                current = _nextIndices[current];
                continue;
            }
            _traversalPositions[current] = i;
            ++i;
            current = _nextIndices[current];
        }
        _msaSeqLength = i; // Cache the computed length
    }

    size_t getMsaSequenceLength() {
        return _msaSeqLength;
    }


    size_t getAbsolutePosition(size_t index) const {
        return _traversalPositions[index];
    }

    

    size_t getIsColumn(size_t index) const {
        return _isColumns[index];
    }


    void printSequence() {
        size_t i = 0;
        size_t current = _headIndex;

        while (current != INVALID) {
            std::cout << current  << " ";

            ++i;
            current = _nextIndices[current];
        }
        std::cout << "\n";
    }

    void printTraversalVec() {
        for (size_t i = 0; i < this->size(); i++) {
            std::cout << _traversalPositions[i] << " ";
        }
        std::cout << "\n";
    }

    void printIndicesVector() {
        for (size_t i = 0; i < this->size(); i++) {
            std::cout << _nextIndices[i] << " ";
        }
        std::cout << "\n";
    }

    bool checkSequenceValidity() {
        
        for (size_t i = 1; i < this->size(); i++)
        {
            size_t numberOfAppearances = 0;
            for (auto &index: _nextIndices) {
                if (i==index) numberOfAppearances++;
            }
            if (numberOfAppearances!=1) {
                std::cout << "position " << i << " appears " << numberOfAppearances << " times\n";
                return false;
            }
        }
        return true;
    }

    class iterator {
    private:
        FixedList* _list;
        size_t _currentIndex;
        
    public:
        iterator(FixedList* list, size_t index) 
            : _list(list), _currentIndex(index) {}
        
        size_t operator*() const { return _currentIndex; }  // Return the index itself
        
        iterator& operator++() {
            _currentIndex = _list->_nextIndices[_currentIndex];
            return *this;
        }

        iterator operator++(int) {
            iterator temp = *this;
            _currentIndex = _list->_nextIndices[_currentIndex];
            return temp;
        }
        
        bool operator!=(const iterator& other) const {
            return _currentIndex != other._currentIndex;
        }
    };


    iterator begin() { return iterator(this, _headIndex); }
    iterator end() { return iterator(this, INVALID); }

    iterator insertAfter(const iterator& it, bool isColumn) {
        size_t newIndex = insertAfter(*it, isColumn);  // Call the size_t version
        return iterator(this, newIndex);  // Wrap the index in an iterator
    }

    iterator insertBefore(const iterator& it, bool isColumn) {
        size_t newIndex = insertBefore(*it, isColumn);
        return iterator(this, newIndex);
    }


    // Accept iterator
    bool referencePosition(const iterator& it) {
        return referencePosition(*it);
    }

    size_t getAbsolutePosition(const iterator& it) const {
        return _traversalPositions[*it];
    }
};

#endif
