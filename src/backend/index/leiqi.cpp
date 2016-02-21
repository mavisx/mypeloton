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
    auto node = mapping_table.get( pid );
    if(node == nullptr)
    {
      std::stack<PidType> empty_res;
      return empty_res;
    }
    std::stack<PidType> path;
    path.push(pid);
    PidType res = search(node, key, path);
    if(res==-1) {
      std::stack<PidType> empty_res;
      return empty_res;
    }

    return path;

  }

  template<typename KeyType>
  PidType BWTree::search(Node* node, KeyType key, std::stack<PidType>& path) {
    //should always kep track of right key range even in delta node
    PidType pid;
//    if(key < node->low_key || key >= node->high_key) {
//      LOG_ERROR("Search Range Err: key not in range");
//      return -1;
//    }
    switch(node->node_type) {
      case LEAF:
      case RECORD_DELTA:
        return node->pid;

      case INDEX_ENTRY_DELTA:
      case DELETE_INDEX_TERM_DELTA:
        if(key >= ((IndexEntryDelta*)node)->Kp
            && key < ((IndexEntryDelta*)node)->Kq) {
          pid= ((MergeDelta *)node)->pQ;
          node = mapping_table.get(pid);
          if(node == nullptr) {
            LOG_ERROR("pid in split/merge delta not exist");
            return -1;
          }
          path.push(pid);
          return search(node, key, path);
        }

        return search(node->next, key, path);

      case REMOVE_NODE_DELTA:
        LOG_INFO("Search Range Info: meet a removed node");

        path.pop();
        if(path.empty()){
          LOG_INFO("Search Path empty");
          return -1;
        }
        pid = path.top();

        node = mapping_table.get(pid);
        if(node == nullptr) {
          LOG_ERROR("Pid in split/merge delta not exist");
          return -1;
        }
        else {
          return search(node, key, path);
        }
      case MERGE_DELTA:
        LOG_INFO("Search Range Info: meet merge delta");
        pid= ((MergeDelta *)node)->pQ;
        if(key >= ((MergeDelta *)node)->Kp) {
          node = mapping_table.get(pid);
          if(node == nullptr) {
            LOG_ERROR("pid in split delta not exist");
            return -1;
          }
          return search(node, key, path);
        }
        return search(node->next, key, path);
      case SPLIT_DELTA:
        LOG_INFO("Search Range Info: meet split/merge delta");
        pid= ((MergeDelta *)node)->pQ;
        if(key >= ((MergeDelta *)node)->Kp) {
          node = mapping_table.get(pid);
          if(node == nullptr) {
            LOG_ERROR("pid in split/merge delta not exist");
            return -1;
          }
          path.push(pid);
          return search(node, key, path);
        }
        return search(node->next, key, path);

      case INNER:
        if(node->slotuse == 0) {
          LOG_ERROR("empty inner node");
          return -1;
        } else {
          int i = 0;
          for(i = 0; i < node->slotuse; i++) {
            if(key >= ((InnerNode*)node)->slotkey[i]) continue;
            else break;
          }
          pid= ((InnerNode *)node)->childid[i];
          node = mapping_table.get(pid);
          if(node == nullptr) {
            LOG_ERROR("pid in inner node not exist");
            return -1;
          }
          path.push(pid);
          return search(node, key, path);
        }
      default:
        return -1;

    }
  }


  template<typename KeyType>
  bool BWTree::DeleteEntry( KeyType key){
    std::stack<PidType > path = search(BWTree::root, key);
    if(path.empty()) {
      LOG_ERROR("InsertEntry get empty tree");
    }
    PidType basic_pid = path.top();
    path.pop();

    Node* basic_node = mapping_table.get( basic_pid );
    RecordDelta* new_delta = new RecordDelta(basic_pid, RecordDelta::DELETE, key);


  };

  template<typename KeyType>
  bool BWTree::is_in( KeyType key, Node** nptr) {

  };




}  // End index namespace
}  // End peloton namespace