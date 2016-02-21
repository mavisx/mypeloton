//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// bwtree.cpp
//
// Identification: src/backend/index/bwtree.cpp
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
// LeiQi
//===----------------------------------------------------------------------===//

#include "backend/index/bwtree.h"
#include "backend/common/logger.h"

namespace peloton {
namespace index {

// Add your function definitions here
    template<typename KeyType>
    PidType BWTree::Search(PidType pid, KeyType key) {
        auto itr = mapping_table.find( pid );
        if(itr == mapping_table.end())
            return -1;

        Node* node = itr->second;
        return Search(node, key);
    }

  template<typename KeyType>
  PidType BWTree::Search(Node* node, KeyType key) {
      //should always kep track of right key range even in delta node
      if(key < node->low_key || key >= node->high_key) {
          LOG_ERROR("Search Range Err: key not in range")
          return -1;
      }
      switch(node->node_type) {
          case LEAF:
              return pid;
          case RECORD_DELTA:
              break;
          case INDEX_ENTRY_DELTA: case DELETE_INDEX_TERM_DELTA:
              return Search(node->next, key);
          case REMOVE_NODE_DELTA:
              break;
          case MERGE_DELTA:
              break;
          default:
              break;
      }



  }



}  // End index namespace
}  // End peloton namespace