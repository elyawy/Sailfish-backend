///////////////////////////////////////////////////////////////////////////////
// \author (c) Marco Paland (info@paland.com)
//             2017-2020, paland consult, Hannover, Germany
//
// \license The MIT License (MIT)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// \brief avl_array class
// This is an AVL tree implementation using an array as data structure.
// avl_array combines the insert/delete and find advantages (log n) of an AVL
// tree with a static allocated arrays and minimal storage overhead. If memory
// is critical the 'Fast' template parameter can be set to false which removes
// the parent member of every node. This saves sizeof(size_type) * Size bytes,
// but slowes down the insert and delete operation by factor 10 due to 'parent
// search'. The find opeartion is not affected cause finding doesn't need a
// parent.
//
// usage:
// #include "avl_array.h"
// avl_array<int, int, int, 1024> avl;
// avl.insert(1, 1);
//
///////////////////////////////////////////////////////////////////////////////

// This file was originally part of the 'avl_array' library by Marco Paland.
// It was modified to fit the Sailfish project requirements, specifically to
// handle the block tree data structure needed for fast indel simulations.
// The core AVL tree implementation remains intact, ensuring efficient
// performance for block tree operations.
// The changes were introduces by Elya Wygoda according to the algorithm described in:
// Elya Wygoda, Asher Moshe, Nimrod Serok, Edo Dotan, Noa Ecker, Naiel Jabareen, 
// Omer Israeli, Itsik Pe’er, Tal Pupko, 
// Efficient algorithms for simulating sequences along a phylogenetic tree, 
// Bioinformatics, Volume 42, Issue 1, January 2026, btaf686, 
// https://doi.org/10.1093/bioinformatics/btaf686

#ifndef _AVL_ARRAY_WITH_RATES_H_
#define _AVL_ARRAY_WITH_RATES_H_

#include <cstdint>
#include <sstream>
#include <vector>
#include <array>
#include <tuple>

#include "Event.h"
#include "CategorySampler.h"

struct Block {
    size_t length;
    size_t insertion;
    std::vector<size_t> rateCategories;
    
    // Handle rate category insertions within this block.
    template<typename RngType>
    void handleInsertion(size_t position, size_t insertLength, CategorySampler* sampler, RngType &rng) {
      if (rateCategories.size() == 0) {
        rateCategories.push_back(SIZE_MAX); // anchor position at rateCategories[0], all insertions happen after it.
        for (size_t i = 0; i < insertLength; i++) {
          rateCategories.push_back(sampler->drawSample(rng));
        }
        return;
      }
      size_t leftFlankCategory = rateCategories[position];
      size_t rightFlankCategory = rateCategories[position + 1];
      std::vector<size_t> newRates;
      // handle insertion at the beginning case.
      if (position == 0) {
        newRates = sampler.sampleRightSidedBridge(rightFlankCategory, insertLength);
      }
      // handle insertion at the end case.
      if (position == (rateCategories.size() - 1)) {
        newRates = sampler.sampleLeftSidedBridge(leftFlankCategory, insertLength);
      }

      // handle insertion in the middle case.
      if (position > 0 && position < (rateCategories.size() - 1)) {
          newRates = sampler.sampleBridge(leftFlankCategory, rightFlankCategory, insertLength);
      }
      // insert new positions in rateCategories vector.
      rateCategories.insert(rateCategories.begin()+(position+1), newRates.begin(), newRates.end());

    }
    
    // Handle rate category deletions within this block
    void handleDeletion(size_t position, size_t deleteLength) {
      auto itStart = rateCategories.begin()+ position;
      auto itEnd = (rateCategories.begin() + position) + deleteLength;
      rateCategories.erase(itStart, itEnd);
    } 
};

typedef std::vector<std::array<size_t, 3>> BlockList;


/**
 * \param Key The key type. The type (class) must provide a 'less than' and
 * 'equal to' operator \param T The Data type \param size_type Container size
 * type \param Size Container size \param Fast If true every node stores an
 * extra parent index. This increases memory but speed up insert/erase by
 * factor 10
 */
template <typename Key, typename size_type, const size_type Size,
          const bool Fast = true, typename RngType>
class avl_array_with_rates {
  // child index pointer class
  typedef struct tag_child_type {
    size_type left;
    size_type right;
  } child_type;


  // node storage, due to possible structure packing effects, single arrays are
  // used instead of a 'node' structure
  Key key_[Size];             // node key
  Block val_[Size];               // node value
  std::size_t length_[Size]; // subtree length

  std::int8_t balance_[Size]; // subtree balance
  child_type child_[Size];    // node childs
  size_type size_;            // actual size
  size_type root_;            // root node
  size_type
      parent_[Fast ? Size : 1]; // node parent, use one element if not needed
                                // (zero sized array is not allowed)

  // invalid index (like 'nullptr' in a pointer implementation)
  static const size_type INVALID_IDX = Size;

  // iterator class
  typedef class tag_avl_array_iterator {
    avl_array *instance_; // array instance
    size_type idx_;       // actual node

    friend avl_array; // avl_array may access index pointer

  public:
    // ctor
    tag_avl_array_iterator(avl_array *instance = nullptr, size_type idx = 0U)
        : instance_(instance), idx_(idx) {}

    inline tag_avl_array_iterator &
    operator=(const tag_avl_array_iterator &other) {
      instance_ = other.instance_;
      idx_ = other.idx_;
      return *this;
    }

    inline bool operator==(const tag_avl_array_iterator &rhs) const {
      return idx_ == rhs.idx_;
    }

    inline bool operator!=(const tag_avl_array_iterator &rhs) const {
      return !(*this == rhs);
    }

    // dereference - access value
    inline Block &operator*() const { return val(); }

    // access value
    inline Block &val() const { return instance_->val_[idx_]; }

    // access key
    inline Key &key() const { return instance_->key_[idx_]; }

    // access length
    inline size_t &length() const { return instance_->length_[idx_]; }


    // preincrement
    tag_avl_array_iterator &operator++() {
      // end reached?
      if (idx_ >= Size) {
        return *this;
      }
      // take left most child of right child, if not existent, take parent
      size_type i = instance_->child_[idx_].right;
      if (i != instance_->INVALID_IDX) {
        // successor is the furthest left node of right subtree
        for (; i != instance_->INVALID_IDX; i = instance_->child_[i].left) {
          idx_ = i;
        }
      } else {
        // have already processed the left subtree, and
        // there is no right subtree. move up the tree,
        // looking for a parent for which nodePtr is a left child,
        // stopping if the parent becomes NULL. a non-NULL parent
        // is the successor. if parent is NULL, the original node
        // was the last node inorder, and its successor
        // is the end of the list
        i = instance_->get_parent(idx_);
        while ((i != instance_->INVALID_IDX) &&
               (idx_ == instance_->child_[i].right)) {
          idx_ = i;
          i = instance_->get_parent(idx_);
        }
        idx_ = i;
      }
      return *this;
    }

    // postincrement
    inline tag_avl_array_iterator operator++(int) {
      tag_avl_array_iterator _copy = *this;
      ++(*this);
      return _copy;
    }
  } avl_array_iterator;

public:
  typedef Block value_type;
  typedef Block *pointer;
  typedef const Block *const_pointer;
  typedef Block &reference;
  typedef const Block &const_reference;
  typedef Key key_type;
  typedef avl_array_iterator iterator;

  // ctor
  avl_array() : size_(0U), root_(Size) {}

  // iterators
  inline iterator begin() {
    size_type i = INVALID_IDX;
    if (root_ != INVALID_IDX) {
      // find smallest element, it's the farthest node left from root
      for (i = root_; child_[i].left != INVALID_IDX; i = child_[i].left)
        ;
    }
    return iterator(this, i);
  }

  inline iterator end() { return iterator(this, INVALID_IDX); }

  // capacity
  inline size_type size() const { return size_; }

  inline bool empty() const { return size_ == static_cast<size_type>(0); }

  inline size_type max_size() const { return Size; }

  /**
   * Clear the container
   */
  inline void clear() {
    size_ = 0U;
    root_ = INVALID_IDX;
  }

  /**
   * Insert or update an element
   * \param key The key to insert. If the key already exists, it is updated
   * \param val Value to insert or update
   * \return True if the key was successfully inserted or updated, false if
   * container is full
   */
  bool insert(const key_type &key, const value_type &val, int added_length) {
    if (root_ == INVALID_IDX) {
      key_[size_] = key;
      val_[size_] = val;
      balance_[size_] = 0;
      length_[size_] = added_length;//val.insertion + val.length;
      child_[size_] = {INVALID_IDX, INVALID_IDX};
      set_parent(size_, INVALID_IDX);
      root_ = size_++;
      return true;
    }

    for (size_type i = root_; i != INVALID_IDX;
         i = (key < key_[i]) ? child_[i].left : child_[i].right) {
      length_[i] += added_length;//(val.insertion + val.length) - (old_val.insertion + old_val.length);
      if (key < key_[i]) {
        if (child_[i].left == INVALID_IDX) {
          if (size_ >= max_size()) {
            // container is full
            return false;
          }
          key_[size_] = key;
          val_[size_] = val;
          balance_[size_] = 0;
          length_[size_] = added_length;//val.insertion + val.length; 
          child_[size_] = {INVALID_IDX, INVALID_IDX};
          set_parent(size_, i);
          child_[i].left = size_++;
          insert_balance(i, 1);
          return true;
        }
      } else if (key_[i] == key) {
        // found same key, update node
        // length_[i] -= val_[i].length + val_[i].insertion;
        val_[i] = val;
        return true;
      } else {
        if (child_[i].right == INVALID_IDX) {
          if (size_ >= max_size()) {
            // container is full
            return false;
          }
          key_[size_] = key;
          val_[size_] = val;
          balance_[size_] = 0;
          length_[size_] = added_length;//val.insertion + val.length;
          child_[size_] = {INVALID_IDX, INVALID_IDX};
          set_parent(size_, i);
          child_[i].right = size_++;
          insert_balance(i, -1);
          return true;
        }
      }
    }
    // node doesn't fit (should not happen) - discard it anyway
    return false;
  }

  /**
   * get the block index
   * \param pos position of the event
   * \return True if the key was successfully inserted or updated, false if
   * container is full
   */
  // this function assumes that the length of the subtree is valid at each node, this is not the case yet.
  // currently this function handles the length of the subtree, this is wrong and should be handled in
  // the rotation functions as well.
  size_type get_block_index(size_t &pos) {
    if (root_ == INVALID_IDX) {
      return INVALID_IDX;
    }

    size_type left_node;
    size_type right_node;
    size_type root_node;

    size_type i = root_;

    // BUG : the left branch case is not working!
    while (i != INVALID_IDX){

      left_node = child_[i].left;
      root_node = i;
      right_node = child_[i].right;

      Block root_val = val_[root_node];
    //   std::cout << "left " << "subtree = " << length_[left_node] << "\n";
    //   std::cout << pos << "\n";
      // val_[root_node].subtree_length += size;
      if (left_node != INVALID_IDX) {
        if (length_[left_node] < pos) {
          pos = pos - length_[left_node];
        } else {
          i = left_node;
          continue;
        }
      }

    //   std::cout << "center " << "subtree = " << root_val.length + root_val.insertion << "\n";;
    //   std::cout << pos << "\n";

      if (root_val.length + root_val.insertion < pos ){
        pos = pos - (root_val.length + root_val.insertion);
      } else {
        return i;
      }

    //   std::cout << "right " << "subtree = " << length_[right_node] << "\n";
    //   std::cout << pos << "\n";

      i = right_node;

    }
    // val_[i].subtree_length += size;

    return i;
  }

  size_type get_next_block(size_type block_index) {
    size_type current_node = block_index;
    current_node = child_[current_node].right;
    while (current_node != INVALID_IDX) {
      if (child_[current_node].left == INVALID_IDX) return current_node;
      current_node = child_[current_node].left;
    }
    current_node = get_parent(block_index);
    if (child_[current_node].left == block_index) return current_node;
    
    size_type previous_node = current_node;
    current_node = get_parent(previous_node);
    while (current_node != INVALID_IDX) {
      if (child_[current_node].left == previous_node) return current_node;
      previous_node = current_node;
      current_node = get_parent(previous_node);
    }

    return INVALID_IDX;
  }

  size_type get_previous_block(size_type block_index) {
    size_type current_node = block_index;
    current_node = child_[current_node].left;
    while (current_node != INVALID_IDX) {
      if (child_[current_node].right == INVALID_IDX) return current_node;
      current_node = child_[current_node].right;
    }
    current_node = get_parent(block_index);
    size_type right_child_node = child_[current_node].right;
    // if (child_[current_node].right == block_index) return current_node;
    if (right_child_node != INVALID_IDX &&
        key_[right_child_node] == key_[block_index]) return current_node;

    size_type previous_node = current_node;
    current_node = get_parent(previous_node);
    while (current_node != INVALID_IDX) {
      if (child_[current_node].right == previous_node) return current_node;
      previous_node = current_node;
      current_node = get_parent(previous_node);
    }

    return INVALID_IDX;
  }


  
  bool split_block(const size_type block_index, size_t pos, size_t event_size) {
    Block event_block = val_[block_index];

    int original_size = event_block.length + event_block.insertion;
    // bool is_anchor_block = (key_[block_index] == 0);
    pos = pos + 1;
    // insertion in added part
    if (pos >= event_block.length ) {
        event_block.insertion = event_block.insertion + event_size;
        int new_size = event_block.insertion + event_block.length;
        int difference_in_length = new_size - original_size;

        return this->insert(key_[block_index], event_block, difference_in_length);
    } else if (pos < event_block.length) { // insertion in origianl part
        // if (pos == 0) pos = 1;

        Block potential_block = { event_block.length - pos, event_block.insertion};
        Block updated_block = {pos, event_size};

        int new_size = updated_block.insertion + updated_block.length;
        int difference_in_length = new_size - original_size;
        size_t potential_block_size = potential_block.insertion + potential_block.length;

        bool event_a = this->insert(key_[block_index], updated_block, difference_in_length);
        bool event_b = this->insert(key_[block_index] + pos, potential_block, potential_block_size);

        return event_a && event_b;
    }
    return false;
  }

  int block_difference(Block &a, Block &b) {
    return (a.length + a.insertion) - (b.length + b.insertion);
  }


  //   xxxxxxxxxxx 
  //  [------OP------|---AP---]
  bool remove_case_a(const size_type block_index, size_t position, size_t event_size, size_t length, size_t insertion) {

    key_type event_key = key_[block_index];
    bool is_valid = true;
    Block new_block = {length - event_size, insertion};
    if (key_[block_index] == 0) { // checking if this is the first block in the blocklist
      Block first_block = {1, 0};
      is_valid = this->insert(0 , first_block, 1 - (length + insertion));
    } else {
      is_valid = this->erase(key_[block_index], length + insertion);
    }
    return this->insert(event_key + event_size, new_block, (length + insertion) - event_size) && is_valid;
  }

  //   xxxxxxxxxxxxxx xxxxxxxx
  //  [------OP------|---AP---]
  bool remove_case_b(const size_type block_index, size_t position, size_t event_size, size_t length, size_t insertion) {

    if (key_[block_index] == 0) { // checking if this is the first block in the blocklist
      Block first_block = {1, 0};
      return this->insert(0, first_block, 1 - (length + insertion));
    } else {
      return this->erase(key_[block_index], length + insertion);
    }
  }

  //   xxxxxxxxxxxxxx xxxx
  //  [------OP------|---AP---]
  bool remove_case_c(const size_type block_index, size_t position, size_t event_size, size_t length, size_t insertion) {
    bool is_valid = true;
    size_t insertion_leftover = (length + insertion) - event_size;
    if (key_[block_index] == 0) { // checking if this is the first block in the blocklist
      Block first_block = {1, insertion_leftover};
      int new_size = first_block.insertion + first_block.length;
      int difference_in_length = new_size - (length + insertion);
      return this->insert(0, first_block, difference_in_length);
    } else {
      size_type previous_block_index = this->get_previous_block(block_index);
      Block previous_block = val_[previous_block_index];
      Block updated_block = {previous_block.length, previous_block.insertion + insertion_leftover};
      is_valid = this->erase(key_[block_index], length + insertion);
      return this->insert(key_[previous_block_index], updated_block, insertion_leftover) && is_valid;
    }
  }

  //      xxxxxxxx
  //  [------OP------|---AP---]
  bool remove_case_d(const size_type block_index, size_t position, size_t event_size, size_t length, size_t insertion) {
    bool is_valid = true;
    Block first_block = {position, 0};
    int new_size = first_block.insertion + first_block.length;
    int difference_in_length = new_size - (length + insertion);
    is_valid = this->insert(key_[block_index], first_block, difference_in_length);
    Block new_block = {length - (position + event_size), insertion};
    new_size = new_block.insertion + new_block.length;
    return this->insert(key_[block_index] + position + event_size, new_block, new_size) && is_valid;
  }


  //         xxxxxxxx
  //  [------OP------|---AP---]
  bool remove_case_e(const size_type block_index, size_t position, size_t event_size, size_t length, size_t insertion) {

    Block first_block = {position, insertion};
    int new_size = first_block.insertion + first_block.length;
    int difference_in_length = new_size - (length + insertion);
    return this->insert(key_[block_index], first_block, difference_in_length);
  }


  //             xxxx xxxx
  //  [------OP------|---AP---]
  bool remove_case_f(const size_type block_index, size_t position, size_t event_size, size_t length, size_t insertion) {
    position = position <= length ? position : length;

    Block first_block = {position, (length + insertion) - (position + event_size)};
    int new_size = first_block.insertion + first_block.length;
    int difference_in_length = new_size - (length + insertion);
    return this->insert(key_[block_index], first_block, difference_in_length);
  }



  bool remove_block(const size_type block_index, size_t position, size_t event_size) {
    Block event_block = val_[block_index];
    size_t length  = event_block.length;
    size_t insertion = event_block.insertion;
    size_t original_size = length + insertion;
    bool is_valid = true;
    // removing the first block still keep a single position, thus we remove 1 from the event size. 
    // this will be removed during the decoding of the blocklist.
    // this handles non spanning deletions.
    if (position + event_size <= original_size) {
      if (position == 0) {
        if (event_size == length + insertion) is_valid = remove_case_b(block_index, position, event_size, length, insertion);
        else if (event_size < length) is_valid = remove_case_a(block_index, position, event_size, length, insertion);
        else if (event_size >= length) is_valid = remove_case_c(block_index, position, event_size, length, insertion);
      } else {
        if ((position + event_size) < length) is_valid = remove_case_d(block_index, position, event_size, length, insertion);
        else if ((position + event_size) == length) is_valid = remove_case_e(block_index, position, event_size, length, insertion);
        else if ((position + event_size) > length) is_valid = remove_case_f(block_index, position, event_size, length, insertion);
      }
    } else {
        // this means the deletion spans multiple blocks, thus we split the deletion and apply
        // it recursively

        size_type next_block_index = this->get_next_block(block_index);
        // std::cout << "Affected blocks:\n";
        // std::cout << key_[block_index] << " " << key_[next_block_index] << "\n";
        size_t updated_size = original_size - position;
        is_valid = this->remove_block(block_index, position, updated_size);

        if (next_block_index != INVALID_IDX) {
          is_valid = this->remove_block(next_block_index, 0, event_size - updated_size) && is_valid; // it this correct?
        }
    }
      
  return is_valid;
  }


  /**
   * Find an element and return an iterator as result
   * \param key The key to find
   * \return Iterator if key was found, else end() is returned
   */
  inline iterator find(const key_type &key) {
    for (size_type i = root_; i != INVALID_IDX;) {
      if (key < key_[i]) {
        i = child_[i].left;
      } else if (key == key_[i]) {
        // found key
        return iterator(this, i);
      } else {
        i = child_[i].right;
      }
    }
    // key not found, return end() iterator
    return end();
  }

  /**
   * Count elements with a specific key
   * Searches the container for elements with a key equivalent to key and
   * returns the number of matches. Because all elements are unique, the
   * function can only return 1 (if the element is found) or zero (otherwise).
   * \param key The key to find/count
   * \return 0 if key was not found, 1 if key was found
   */
  inline size_type count(const key_type &key) {
    return find(key) != end() ? 1U : 0U;
  }

  /**
   * Remove element by key
   * \param key The key of the element to remove
   * \param added_length value that is added to the subtree length from root to target node
   * \return True if the element ws removed, false if key was not found
   */
  inline bool erase(const key_type &key, size_t added_length) { return erase(find(key), added_length); }

  /**
   * Remove element by iterator position
   * THIS ERASE OPERATION INVALIDATES ALL ITERATORS!
   * \param position The iterator position of the element to remove
   * \param added_length value that is added to the subtree length from root to target node
   * \return True if the element was successfully removed, false if error
   */
  bool erase(iterator position, size_t added_length) {
    if (empty() || (position == end())) {
      return false;
    }

    const size_type node = position.idx_;
    const size_type left = child_[node].left;
    const size_type right = child_[node].right;


    size_type predecessor = this->get_parent(node);
    length_[node] -= added_length;
    // if (node == root_) length_[root_] -= added_length;
    while (predecessor != INVALID_IDX) {
      length_[predecessor] -= added_length;
      predecessor = this->get_parent(predecessor);
    }



    if (left == INVALID_IDX) {
      if (right == INVALID_IDX) {
        const size_type parent = get_parent(node);
        if (parent != INVALID_IDX) {
          if (child_[parent].left == node) {
            child_[parent].left = INVALID_IDX;
            delete_balance(parent, -1);
          } else {
            child_[parent].right = INVALID_IDX;
            delete_balance(parent, 1);
          }
        } else {
          root_ = INVALID_IDX; //  this will not be possible since case (a) always keeps the 0 key block.
        }
      } else {
        const size_type parent = get_parent(node);
        if (parent != INVALID_IDX) {
          child_[parent].left == node ? child_[parent].left = right
                                      : child_[parent].right = right;
        } else {
          root_ = right;
        }
        set_parent(right, parent);
        delete_balance(right, 0);
      }
    } else if (right == INVALID_IDX) {
      const size_type parent = get_parent(node);
      if (parent != INVALID_IDX) {
        child_[parent].left == node ? child_[parent].left = left
                                    : child_[parent].right = left;
      } else {
        root_ = left;
      }
      set_parent(left, parent);
      delete_balance(left, 0);
    } else {
      size_type successor = right;
      if (child_[successor].left == INVALID_IDX) {
        const size_type parent = get_parent(node);
        child_[successor].left = left;
        balance_[successor] = balance_[node];
        length_[successor] = length_[node];

        set_parent(successor, parent);
        set_parent(left, successor);

        if (node == root_) {
          root_ = successor;
        } else {
          if (child_[parent].left == node) {
            child_[parent].left = successor;
          } else {
            child_[parent].right = successor;
          }
        }
        delete_balance(successor, 1);
      } else {
        while (child_[successor].left != INVALID_IDX) {
          successor = child_[successor].left;
        }

        size_t length_successor_only = val_[successor].length + val_[successor].insertion;
        size_type successor_length_removal = get_parent(successor);


        while (successor_length_removal != node) {
          length_[successor_length_removal] -= length_successor_only;
          successor_length_removal = get_parent(successor_length_removal);
        }


        const size_type parent = get_parent(node);
        const size_type successor_parent = get_parent(successor);
        const size_type successor_right = child_[successor].right;




        // size_t length_successor_parent_only =  val_[successor_parent].length + val_[successor_parent].insertion;

        // size_t sum_node_parent = parent != INVALID_IDX ? length_[parent] : 0;
        // size_t sum_successor_right = successor_right != INVALID_IDX ? length_[successor_right] : 0;



        if (child_[successor_parent].left == successor) {
          child_[successor_parent].left = successor_right;
        } else {
          child_[successor_parent].right = successor_right;
        }

        set_parent(successor_right, successor_parent);
        set_parent(successor, parent);
        set_parent(right, successor);
        set_parent(left, successor);

        child_[successor].left = left;
        child_[successor].right = right;
        balance_[successor] = balance_[node];
        length_[successor] = length_[right] + length_[left] + length_successor_only;




        if (node == root_) {
          root_ = successor;
        } else {
          if (child_[parent].left == node) {
            child_[parent].left = successor;
          } else {
            child_[parent].right = successor;
          }
        }
        delete_balance(successor_parent, -1);
      }
    }
    size_--;

    // relocate the node at the end to the deleted node, if it's not the
    // deleted one
    if (node != size_) {
      size_type parent = INVALID_IDX;
      if (root_ == size_) {
        root_ = node;
      } else {
        parent = get_parent(size_);
        child_[parent].left == size_ ? child_[parent].left = node
                                     : child_[parent].right = node;
      }

      // correct childs parent
      set_parent(child_[size_].left, node);
      set_parent(child_[size_].right, node);

      // move content
      key_[node] = key_[size_];
      val_[node] = val_[size_];
      balance_[node] = balance_[size_];
      child_[node] = child_[size_];
      length_[node] = length_[size_];
      set_parent(node, parent);
    }

    return true;
  }

  /**
   * Integrity (self) check
   * \return True if the tree intergity is correct, false if error (should not
   * happen normally)
   */
  bool check() const {
    // check root
    if (empty() && (root_ != INVALID_IDX)) {
      // invalid root
      return false;
    }
    if (size() && root_ >= size()) {
      // root out of bounds
      return false;
    }

    // check tree
    for (size_type i = 0U; i < size(); ++i) {
      if ((child_[i].left != INVALID_IDX) &&
          (!(key_[child_[i].left] < key_[i]) ||
           (key_[child_[i].left] == key_[i]))) {
        // wrong key order to the left
        return false;
      }
      if ((child_[i].right != INVALID_IDX) &&
          ((key_[child_[i].right] < key_[i]) ||
           (key_[child_[i].right] == key_[i]))) {
        // wrong key order to the right
        return false;
      }
      const size_type parent = get_parent(i);
      if ((i != root_) && (parent == INVALID_IDX)) {
        // no parent
        return false;
      }
      if ((i == root_) && (parent != INVALID_IDX)) {
        // invalid root parent
        return false;
      }
    }
    // check passed
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Helper functions
private:
  // find parent element
  inline size_type get_parent(size_type node) const {
    if (Fast) {
      return parent_[node];
    } else {
      const Key key_node = key_[node];
      for (size_type i = root_; i != INVALID_IDX;
           i = (key_node < key_[i]) ? child_[i].left : child_[i].right) {
        if ((child_[i].left == node) || (child_[i].right == node)) {
          // found parent
          return i;
        }
      }
      // parent not found
      return INVALID_IDX;
    }
  }

  // set parent element (only in Fast version)
  inline void set_parent(size_type node, size_type parent) {
    if (Fast) {
      if (node != INVALID_IDX) {
        parent_[node] = parent;
      }
    }
  }


  void insert_balance(size_type node, std::int8_t balance) {
    while (node != INVALID_IDX) {
      balance = (balance_[node] += balance);

      if (balance == 0) {
        return;
      } else if (balance == 2) {

        if (balance_[child_[node].left] == 1) {
          rotate_right(node);
        } else {
          rotate_left_right(node);
        }
        return;
      } else if (balance == -2) {

        if (balance_[child_[node].right] == -1) {
          rotate_left(node);
        } else {
          rotate_right_left(node);
        }
        return;
      }

      const size_type parent = get_parent(node);
      if (parent != INVALID_IDX) {
        balance = child_[parent].left == node ? 1 : -1;
      }
      node = parent;
    }
  }

  void delete_balance(size_type node, std::int8_t balance) {
    while (node != INVALID_IDX) {
      balance = (balance_[node] += balance);

      if (balance == -2) {
        if (balance_[child_[node].right] <= 0) {
          node = rotate_left(node);
          if (balance_[node] == 1) {
            return;
          }
        } else {
          node = rotate_right_left(node);
        }
      } else if (balance == 2) {
        if (balance_[child_[node].left] >= 0) {
          node = rotate_right(node);
          if (balance_[node] == -1) {
            return;
          }
        } else {
          node = rotate_left_right(node);
        }
      } else if (balance != 0) {
        return;
      }

      if (node != INVALID_IDX) {
        const size_type parent = get_parent(node);
        if (parent != INVALID_IDX) {
          balance = child_[parent].left == node ? -1 : 1;
        }
        node = parent;
      }
    }
  }

  size_type rotate_left(size_type node) {
    const size_type right = child_[node].right;
    const size_type right_left = child_[right].left;
    const size_type parent = get_parent(node);

    size_t sum_left = 0;
    size_t sum_right_right = 0;
    size_t sum_right_left = 0;
    if (child_[right].right != INVALID_IDX) {
      sum_right_right = length_[child_[right].right];
    } 
    if (child_[node].left != INVALID_IDX) {
      sum_left = length_[child_[node].left];
    }
    if (child_[right].left != INVALID_IDX) {
      sum_right_left = length_[child_[right].left];
    }


    set_parent(right, parent);
    set_parent(node, right);
    set_parent(right_left, node);
    child_[right].left = node;
    child_[node].right = right_left;

    if (node == root_) {
      root_ = right;
    } else if (child_[parent].right == node) {
      child_[parent].right = right;
    } else {
      child_[parent].left = right;
    }

    balance_[right]++;
    balance_[node] = -balance_[right];

    size_t length_node_only = val_[node].length + val_[node].insertion;
    size_t length_right_only = val_[right].length + val_[right].insertion;

    length_[node] = sum_left + sum_right_left + length_node_only;
    length_[right] = length_[node] + sum_right_right + length_right_only;



    return right;
  }

  size_type rotate_right(size_type node) {
    const size_type left = child_[node].left;
    const size_type left_right = child_[left].right;
    const size_type parent = get_parent(node);


    size_t sum_right = 0;
    size_t sum_left_left = 0;
    size_t sum_left_right = 0;
    if (child_[left].left != INVALID_IDX) {
      sum_left_left = length_[child_[left].left];
    } 
    if (child_[node].right != INVALID_IDX) {
      sum_right = length_[child_[node].right];
    }
    if (child_[left].right != INVALID_IDX) {
      sum_left_right = length_[child_[left].right];
    }


    set_parent(left, parent);
    set_parent(node, left);
    set_parent(left_right, node);
    child_[left].right = node;
    child_[node].left = left_right;

    if (node == root_) {
      root_ = left;
    } else if (child_[parent].left == node) {
      child_[parent].left = left;
    } else {
      child_[parent].right = left;
    }

    balance_[left]--;
    balance_[node] = -balance_[left];



    size_t length_node_only = val_[node].length + val_[node].insertion;
    size_t length_left_only = val_[left].length + val_[left].insertion;
    length_[node] = sum_right + sum_left_right + length_node_only;
    length_[left] = length_[node] + sum_left_left + length_left_only;


    return left;
  }

  size_type rotate_left_right(size_type node) {
    const size_type left = child_[node].left;
    const size_type left_right = child_[left].right;
    const size_type left_right_right = child_[left_right].right;
    const size_type left_right_left = child_[left_right].left;
    const size_type parent = get_parent(node);

    size_t sum_right = 0;
    size_t sum_left_left = 0;
    size_t sum_left_right_right = 0;
    size_t sum_left_right_left= 0;
    if (child_[node].right != INVALID_IDX) {
      sum_right = length_[child_[node].right];
    }
    if (child_[left].left != INVALID_IDX) {
      sum_left_left = length_[child_[left].left];
    }
    if (child_[left_right].right != INVALID_IDX) {
      sum_left_right_right = length_[child_[left_right].right];
    } 
    if (child_[left_right].left != INVALID_IDX) {
      sum_left_right_left = length_[child_[left_right].left];
    }

    set_parent(left_right, parent);
    set_parent(left, left_right);
    set_parent(node, left_right);
    set_parent(left_right_right, node);
    set_parent(left_right_left, left);
    child_[node].left = left_right_right;
    child_[left].right = left_right_left;
    child_[left_right].left = left;
    child_[left_right].right = node;

    if (node == root_) {
      root_ = left_right;
    } else if (child_[parent].left == node) {
      child_[parent].left = left_right;
    } else {
      child_[parent].right = left_right;
    }

    if (balance_[left_right] == 0) {
      balance_[node] = 0;
      balance_[left] = 0;
    } else if (balance_[left_right] == -1) {
      balance_[node] = 0;
      balance_[left] = 1;
    } else {
      balance_[node] = -1;
      balance_[left] = 0;
    }
    balance_[left_right] = 0;

    size_t length_node_only = val_[node].length + val_[node].insertion;
    size_t length_left_only = val_[left].length + val_[left].insertion;
    size_t length_left_right_only = val_[left_right].length + val_[left_right].insertion;
    length_[node] = sum_right + sum_left_right_right + length_node_only;
    length_[left] = sum_left_right_left + sum_left_left + length_left_only;
    length_[left_right] = length_[node] + length_[left] + length_left_right_only;


    return left_right;
  }

  size_type rotate_right_left(size_type node) {
    const size_type right = child_[node].right;
    const size_type right_left = child_[right].left;
    // const size_type right_right = child_[right].right;
    const size_type right_left_left = child_[right_left].left;
    const size_type right_left_right = child_[right_left].right;
    const size_type parent = get_parent(node);

    size_t sum_left = 0;
    size_t sum_right_right = 0;
    size_t sum_right_left_left = 0;
    size_t sum_right_left_right= 0;
    if (child_[node].left != INVALID_IDX) {
      sum_left = length_[child_[node].left];
    }
    if (child_[right].right != INVALID_IDX) {
      sum_right_right = length_[child_[right].right];
    }
    if (child_[right_left].left != INVALID_IDX) {
      sum_right_left_left = length_[child_[right_left].left];
    } 
    if (child_[right_left].right != INVALID_IDX) {
      sum_right_left_right = length_[child_[right_left].right];
    }

    set_parent(right_left, parent);
    set_parent(right, right_left);
    set_parent(node, right_left);
    set_parent(right_left_left, node);
    set_parent(right_left_right, right);
    child_[node].right = right_left_left;
    child_[right].left = right_left_right;
    child_[right_left].right = right;
    child_[right_left].left = node;

    if (node == root_) {
      root_ = right_left;
    } else if (child_[parent].right == node) {
      child_[parent].right = right_left;
    } else {
      child_[parent].left = right_left;
    }

    if (balance_[right_left] == 0) {
      balance_[node] = 0;
      balance_[right] = 0;
    } else if (balance_[right_left] == 1) {
      balance_[node] = 0;
      balance_[right] = -1;
    } else {
      balance_[node] = 1;
      balance_[right] = 0;
    }
    balance_[right_left] = 0;

    size_t length_node_only = val_[node].length + val_[node].insertion;
    size_t length_right_only = val_[right].length + val_[right].insertion;
    size_t length_right_left_only = val_[right_left].length + val_[right_left].insertion;
    length_[node] = sum_left + sum_right_left_left + length_node_only;
    length_[right] = sum_right_left_right + sum_right_right + length_right_only;
    length_[right_left] = length_[node] + length_[right] + length_right_left_only;


    return right_left;
  }



void printAVL(std::stringstream& ss, const std::string& prefix, const size_type node, bool isLeft) {
    // prefix = prefix + (isLeft ? "│   " : "    ");

    if (node == INVALID_IDX) {
      return;
    }

    ss << prefix;

    ss << (isLeft ? "├──" : "└──" );

    // print the value of the node
    print_block(ss, node);

    // enter the next tree level - left and right branch
    printAVL(ss, prefix + (isLeft ? "│   " : "    "), child_[node].left, true);
    printAVL(ss, prefix + (isLeft ? "│   " : "    "), child_[node].right, false);

    }
public:

void print_block(std::stringstream &ss, const size_type node) {
    key_type key = key_[node];
    Block block = val_[node];
    size_t subtree_length = length_[node];

    // std::string block_str = std::string("[") + std::to_string(key) + "|" + 
    //     std::to_string(block.length) + "|" + std::to_string(block.insertion) +
    //     "]->" + std::to_string(subtree_length) + "\n";
    // return block_str;
    ss << "[" << key << "|" << 
        block.length << "|" << block.insertion << "]->" << subtree_length << std::endl;
  }

std::string print_avl() {
    std::stringstream tree_stream;

    printAVL(tree_stream, "", root_, false);
    return tree_stream.str();
  };


// this function assumes that the length of the subtree is valid at each node, this is not the case yet.
bool handle_event(Event &ev, RngType &rng) {

    size_type block_index = this->get_block_index(ev.position);
    if (block_index == INVALID_IDX) {
        std::cout << "could not get event key!\n";
        return false;
    }

    if (ev.type == INSERTION) {
      return split_block(block_index, ev.position, ev.size);
    }

    if (ev.type == DELETION) {
      if (event_position == 0) throw std::out_of_range("event_position exceeds sequence");
      return remove_block(block_index, ev.position, ev.size);
    }



    return false;    
    
}



BlockList get_blocklist() {
    BlockList blocklist;
    for (auto it = this->begin(); it != this->end(); ++it) {
        std::array<size_t,3> current_block =  {(&it)->key(), (*it).length, (*it).insertion};
        // std::tuple<int, int, int> current_block ((&it)->key(), (*it).length, (*it).insertion);
        blocklist.push_back(current_block);
    }
    return blocklist;
}

bool init_tree(size_t sequence_length) {
  this->clear();

  Block root_block = {sequence_length, 0};
  return this->insert(0, root_block, sequence_length);
}


bool checkLength(size_type node) {
    bool left_length = true;
    bool right_length = true;
    size_t node_length_only = val_[node].length + val_[node].insertion;

    size_type left = child_[node].left;
    size_type right = child_[node].right;



    if (left == INVALID_IDX && right == INVALID_IDX) {
        return length_[node] == node_length_only;
    } 
    if (left != INVALID_IDX && right != INVALID_IDX) {
        left_length = (length_[node] == length_[left] + length_[right] + node_length_only);

        return checkLength(left) && checkLength(right) && left_length;
    }
    if (left != INVALID_IDX) {
        left_length = (length_[node] == length_[left] + node_length_only);
        return checkLength(left) && left_length;
    }
    if (right != INVALID_IDX) {
        right_length = (length_[node] == length_[right] + node_length_only);
        return checkLength(right) && right_length;
    }

    return false;
}

bool checkLength() {
    return checkLength(root_);
}

size_t getTotalLength() {
    return length_[root_];
}

size_t memoryUsage() {
  return (3*sizeof(size_type) + sizeof(Key)
           + sizeof(Block) + sizeof(size_t)
           + sizeof(std::int8_t) )*size();
}

};

#endif // _AVL_ARRAY_H_
