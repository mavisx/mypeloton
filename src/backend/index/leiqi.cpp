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
    if(node == nullptr) {
      std::stack<PidType> empty_res;
      return empty_res;
    }

    std::stack<PidType> path;
    path.push(pid);
    PidType res = search(node, key, path);
    if(res == -1) {
      std::stack<PidType> empty_res;
      return empty_res;
    }

    return path;

  }

  template<typename KeyType>
  PidType BWTree::search(Node* node, KeyType key, std::stack<PidType>& path) {
    //should always keep track of right key range even in delta node
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
          pid = ((IndexEntryDelta *)node)->pQ;
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

        if(key >= ((MergeDelta *)node)->Kp) {
          node = ((MergeDelta *)node)->orignal_node;
          if(node == nullptr) {
            LOG_ERROR("pid in split delta not exist");
            return -1;
          }
          return search(node, key, path);
        }
        return search(node->next, key, path);
      case SPLIT_DELTA:
        LOG_INFO("Search Range Info: meet split/merge delta");
        pid= ((SplitDelta*)node)->pQ;
        if(key >= ((SplitDelta *)node)->Kp) {
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
  bool BWTree::delete_entry( KeyType key){
    std::stack<PidType > path = search(BWTree::root, key);
    if(path.empty()) {
      LOG_ERROR("InsertEntry get empty tree");
    }
    PidType basic_pid = path.top();
    path.pop();

    Node* basic_node = mapping_table.get( basic_pid );
    if(!is_in(key, basic_node)) {
      LOG_INFO("DeleteEntry Not Exist");
      return false;
    }

    RecordDelta* new_delta = new RecordDelta(basic_pid, RecordDelta::DELETE, key);



  };

  template<typename KeyType>
  bool BWTree::apend_delete(KeyType key, Node* basic_pid){
    RecordDelta* new_delta = new RecordDelta(basic_pid, RecordDelta::INSERT, key);
    new_delta->high_key = basic_node->high_key;
    new_delta->low_key = basic_node->low_key;
  }

  template<typename KeyType>
  bool BWTree::is_in( KeyType key, Node* listhead) {
    if(listhead == nullptr) return false;

    Node* node = listhead;
    switch (node->node_type){
      case RECORD_DELTA:
        RecordDelta* rcd_node = (RecordDelta*)node;
        if (rcd_node->op_type == RecordDelta::INSERT
            && rcd_node->key == key) {
          return true;
        }
        else if(rcd_node->op_type == RecordDelta::DELETE
            && rcd_node->key == key) {
          return false;
        }
        return is_in(key, node->next);
      case LEAF:
        LeafNode* lf_node = (LeafNode*)node;
        for(int i = 0; i< (lf_node->slotuse); i++){
          if(lf_node->slotkey[i] == key) {
            return true;
          }
        }
        return false;
      case MERGE_DELTA:
      case SPLIT_DELTA:
        if(key >= ((MergeDelta *)node)->Kp) {
          PidType pid= ((MergeDelta *)node)->pQ;
          return is_in(key,mapping_table.get(pid));
        }
        return is_in(key, node->next);
      default:
        return false;
    }

  };



  template<typename KeyType>
  BWTree::RecordDelta::RecordDelta(PidType next, BWTree::RecordDelta::RecordType op, KeyType k)
      : Node(NodeType::RECORD_DELTA){
    op_type = op;
    key = k;
    // Get node* of original node form mapping_table
    Node* orig_node = mapping_table.get(next);
    prepend(this, orig_node);

    // update the slotuse of the new delta node
    if(op_type == INSERT) {
      slotuse = (unsigned short) (orig_node->slotuse + 1);
    }
    else if(op_type == DELETE) {
      slotuse = (unsigned short) (orig_node->slotuse - 1);
    }

  }


  template<typename KeyType>
  BWTree::RecordDelta::RecordDelta(PidType next, BWTree::RecordDelta::RecordType op, KeyType k, ValueType v)
      : Node(NodeType::RECORD_DELTA){
    RecordDelta(next, op,k);
    this->value = v;
  }

}  // End index namespace
}  // End peloton namespace