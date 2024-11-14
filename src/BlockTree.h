#ifndef _BLOCK_TREE_H_
#define _BLOCK_TREE_H_

#include <memory>
#include <unordered_map>
#include <tuple>
#include "AvlTree.h"


typedef std::unordered_map<size_t, std::tuple<BlockList, size_t>> BlockMap;


enum class BLOCK {
  POSITION = 0,
  LENGTH = 1,
  INSERTION = 2
};

enum class BLOCKLIST {
    BLOCKS = 0,
    LENGTH = 1
};

/**
 * A class representing a block tree data structure, handling insertion and deletion events
 * efficiently.
 */
class BlockTree
{
private:
  /**
   * Type alias for the AVL array template used to implement the block 
tree.
    *
    * Note: The size of 1000000U ensures that the tree has enough 
capacity for all blocks.
  */
  using TreeType = avl_array<std::uint32_t, std::uint32_t, 1000000U, true>;
  std::shared_ptr<TreeType> _avlTree;
public:
  /**
   * Constructor for BlockTree class.
   *
   * Initializes a new block tree with the specified first block size,
   * creating an AVL array tree and setting up its internal state.
   *
   * @param first_block_size The initial size of the first block in the 
tree.
   */
  BlockTree(int first_block_size) {
    _avlTree = std::make_shared<TreeType>();//new TreeType;
    _avlTree->init_tree(first_block_size + 1);
  }

  BlockTree(const BlockTree& otherTree): _avlTree(otherTree._avlTree) {}
  BlockTree() {}

  /**
 * Handles an event by delegating the operation to the AVL tree.
 *
 * @param ev The event type being handled.
 * @param event_position The position of the event within the BlockTree(current sequence).
 * @param event_size The length of the event e.g. "deletion of length 3".
 */
  void handleEvent (event ev, size_t event_position, size_t event_size) {
    _avlTree->handle_event(ev, event_position, event_size);
  }

  std::string printTree () {
    return _avlTree->print_avl();
  }

  BlockList getBlockList () {
    return _avlTree->get_blocklist();
  }

  TreeType::iterator begin () {
    return _avlTree->begin();
  }

  TreeType::iterator end () {
    return _avlTree->end();
  }

  size_t length() {
    return _avlTree->getTotalLength();
  }

  size_t memoryUsage() {
    return _avlTree->memoryUsage();
  }

/**
 * Usefull for debugging:
 * Checks if the overall length of the BlockTree is valid throughout.
 * The length of the tree is the sum of the length of the root block + left childs length
 * + right childs length. If any of these is not correct, this function will return false.
 *
 * @return TRUE if the length is valid, FALSE otherwise.
 */
  bool checkLength() {
    return _avlTree->checkLength();
  }


  ~BlockTree() {
  }
};


#endif