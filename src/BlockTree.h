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

class BlockTree
{
private:
 // consider creating different sized avl_arrays to reduce ram usage.
  using TreeType = avl_array<std::uint32_t, std::uint32_t, 1000000U, true>;
  std::shared_ptr<TreeType> _avlTree;
public:
  BlockTree(int first_block_size) {
    _avlTree = std::make_shared<TreeType>();//new TreeType;
    _avlTree->init_tree(first_block_size + 1);
  }

  BlockTree(const BlockTree& otherTree): _avlTree(otherTree._avlTree) {}
  BlockTree() {}

  void handleEvent (event ev, size_t event_position, size_t event_size) {
    if (event_size == 0) return;
    if(!_avlTree->handle_event(ev, event_position, event_size)) {
      throw std::out_of_range("event_position exceeds sequence");
    }
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

  bool checkLength() {
    return _avlTree->checkLength();
  }


  ~BlockTree() {
  }
};

// class BlockContainer {
// private:
//   std::map<std::string, BlockTree*> nodeToBlock;
// public:
//   BlockContainer() {}

// }

#endif