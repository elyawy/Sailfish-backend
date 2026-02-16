#ifndef _BLOCK_COMMON_H_
#define _BLOCK_COMMON_H_

#include <unordered_map>
#include <tuple>
#include <vector>
#include <array>

typedef std::vector<std::array<size_t, 3>> BlockList;

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

#endif