#ifndef _BLOCK_TREE_WITH_RATES_H_
#define _BLOCK_TREE_WITH_RATES_H_

#include <memory>
#include <unordered_map>
#include <tuple>
#include "AvlTreeWithRates.h"
#include "Event.h"
#include "BlockCommon.h"

class BlockTreeWithRates
{
private:
 // consider creating different sized avl_arrays to reduce ram usage.
  using TreeType = avl_array_with_rates<std::uint32_t, std::uint32_t, 1000000U, true>;
  std::unique_ptr<TreeType> _avlTree;
public:
  BlockTreeWithRates() {
    _avlTree = std::make_unique<TreeType>();
  }

  template<typename RngType = std::mt19937_64>
  void handleEvent (Event &ev, CategorySampler& sampler, RngType &rng) {
    if (ev.length == 0) return;
    if(!_avlTree->handle_event(ev, sampler, rng)) {
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

  void clear(){
    _avlTree->clear();
  }

  void initTree(int first_block_size, const std::vector<size_t>& rateCategories){
    _avlTree->clear();
    _avlTree->init_tree(first_block_size + 1, rateCategories);
  }


  ~BlockTreeWithRates() {
  }
};


#endif