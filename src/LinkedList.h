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
    size_t _leafNum;
    size_t _msaSeqLength;

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
        _msaSeqLength = 0;
        _leafNum = 0;
        if (sequenceSize <= 0) return;
        
        // Update count to include the anchor (0) plus sequenceSize positions
        _count = sequenceSize + 1;
        
        _nextIndices[_headIndex] = 1;  // anchor points to position 1
        for (size_t i = 1; i < sequenceSize; ++i) {
            if (i+1 == INVALID) break;
            _nextIndices[i] = i+1;
            _isColumns[i] = false;  // Initialize as not columns
        }
        _nextIndices[sequenceSize] = INVALID;
        _tailIndex = sequenceSize;
        
    }
    
    // Insert after element at index k
    // Returns index of newly inserted element, or INVALID if failed
    size_t insert_after(size_t nodeK, bool isColumn) {
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

    bool referencePosition(size_t k, size_t position) {
        // Check if k is valid
        if (k >= _count) {
            return false;
        }
        
        _isColumns[k] = true;
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
    }
};