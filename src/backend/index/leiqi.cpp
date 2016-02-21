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
  std::stack<PidType> BWTree::search(PidType pid, KeyType key) {
    auto itr = mapping_table.find( pid );
    if(itr == mapping_table.end())
      return -1;

  BWTree::Node* node = itr->second;
    return search(node, key);
  }

  template<typename KeyType>
  PidType BWTree::search(Node* node, KeyType key, std::stack<PidType>& path) {
    //should always kep track of right key range even in delta node
    if(key < node->low_key || key >= node->high_key) {
      LOG_ERROR("Search Range Err: key not in range");
      return -1;
    }
    switch(node->node_type) {
      case LEAF:
      case RECORD_DELTA:
        return node->pid;

      case INDEX_ENTRY_DELTA:
      case DELETE_INDEX_TERM_DELTA:
        return search(node->next, key);

      case REMOVE_NODE_DELTA:
        LOG_INFO("Search Range Info: meet a removed node");
        return -1;

      case MERGE_DELTA:
      case SPLIT_DELTA:
        LOG_INFO("Search Range Info: meet split/merge delta");
        if(key >= ((MergeDelta)node)->Kp) {
          node = mapping_table.get(((MergeDelta)node).pQ);
          if(node == NULL) {
            LOG_ERROR("pid in split/merge delta in not exist");
            return -1;
          }
          return search()
        }
      default:
        break;
    }



  }



}  // End index namespace
}  // End peloton namespace