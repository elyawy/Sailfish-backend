#ifndef _BLOCK_TREE_H_
#define _BLOCK_TREE_H_

#include <memory>

#include "AvlTree.h"

class BlockTree
{
private:
 // consider creating different sized avl_arrays to reduce ram usage.
  using TreeType = avl_array<int, std::uint32_t, 100000U, true>;
  std::shared_ptr<TreeType> _avlTree;
public:
  BlockTree(int first_block_size) {
    _avlTree = std::make_shared<TreeType>();//new TreeType;
    _avlTree->init_tree(first_block_size + 1);
  }

  BlockTree(const BlockTree& otherTree): _avlTree(otherTree._avlTree) {}
  BlockTree() {}

  void handleEvent (event ev, size_t event_position, size_t event_size) {
    _avlTree->handle_event(ev, event_position, event_size);
  }

  std::string printTree () {
    return _avlTree->print_avl();
  }

  std::vector<std::tuple<int, int, int>> blockList () {
    return _avlTree->get_blocklist();
  }

  TreeType::iterator begin () {
    return _avlTree->begin();
  }

  const TreeType::iterator begin () const {
    return _avlTree->begin();
  }

  TreeType::iterator end () {
    return _avlTree->end();
  }

  const TreeType::iterator end () const {
    return _avlTree->end();
  }

  size_t length() const {
    return _avlTree->getTotalLength();
  }

  size_t memoryUsage() {
    return _avlTree->memoryUsage();
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