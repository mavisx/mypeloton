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
//===----------------------------------------------------------------------===//

#include "backend/index/bwtree.h"

namespace peloton {
namespace index {

MappingTable::MappingTable() {
  for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
    mappingtable_1[i] = nullptr;
  }
  nextPid = 0;
}

MappingTable::~MappingTable() {
  for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
    if (mappingtable_1[i] != nullptr) {
      delete[] mappingtable_1[i];
      mappingtable_1[i] = nullptr;
    }
  }
}

BWTree::Node* MappingTable::get(PidType pid) {
  long tier1_idx = GET_TIER1_INDEX(pid);
  long tier2_idx = GET_TIER2_INDEX(pid);

  if (mappingtable_1[tier1_idx] == nullptr) return nullptr;

  return mappingtable_1[tier1_idx][tier2_idx];
}

// cas set the item in mappintable[pid] as addr
bool MappingTable::set(PidType pid, void* addr) {
  long tier1_idx = GET_TIER1_INDEX(pid);
  long tier2_idx = GET_TIER2_INDEX(pid);

  if (mappingtable_1[tier1_idx] == nullptr) {
    return false;
  }

  void* expected = mappingtable_1[tier1_idx][tier2_idx];

  // using cas to set the value in mapping table
  return std::atomic_compare_exchange_strong(
      (std::atomic<void*>*)&mappingtable_1[tier1_idx][tier2_idx], &expected,
      addr);
}

// if correctly add this new addr, return the new pid,
// else return -1
long MappingTable::add(void* addr) {
  unsigned long new_pid = nextPid++;
  long tier1_idx = GET_TIER1_INDEX(new_pid);
  long tier2_idx = GET_TIER2_INDEX(new_pid);

  void* expectded = (void*)mappingtable_1[tier1_idx];

  // atomically add new secondary index
  if (mappingtable_1[tier1_idx] == nullptr) {
    void* desired = new void* [MAPPING_TABLE_SIZE];
    for (int i = 0; i < MAPPING_TABLE_SIZE; i++ ) {
      desired[i] = nullptr;
    }
    if (!std::atomic_compare_exchange_strong(
            (std::atomic<void*>*)&mappingtable_1[tier1_idx], &expectded,
            desired)) {
      delete[] desired;
    }
  }

  expectded = (void*)mappingtable_1[tier1_idx][tier2_idx];
  // atomically add addr to desired secondary index
  if (std::atomic_compare_exchange_strong(
          (std::atomic<void*>*)&mappingtable_1[tier1_idx][tier2_idx],
          &expectded, addr)) {
    return new_pid;
  } else {
    return -1;
  }
}

void MappingTable::remove(PidType pid) {
  long tier1_idx = GET_TIER1_INDEX(pid);
  long tier2_idx = GET_TIER2_INDEX(pid);

  mappingtable_1[tier1_idx][tier2_idx] = nullptr;
}

template <class KeyComparator, typename KeyEqualityChecker>
BWTree::BWTree(KeyComparator &kc, KeyEqualityChecker &ke):
    m_key_less(kc), m_key_equal(ke) {

  LeafNode* addr = new LeafNode();
  long newpid = mapping_table.add(addr);
  if (newpid >= 0) {
    // initialize the root, head and tail pid.
    root = newpid;
    headleaf = tailleaf = newpid;
  } else {
    LOG_ERROR("Can't create the initial leafNode!");
  }

}

BWTree::~BWTree() {}

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
      pid = ((SplitDelta*)node)->pQ;
      if(key >= ((SplitDelta *)node)->Kp) {
        node = mapping_table.get(pid);
        if(node == nullptr) {
          LOG_ERROR("pid in split/merge delta not exist");
          return -1;
        }
        // replace the top with our split node
        path.pop();
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
  while(!apend_delete(key)){
    LOG_INFO("delete_entry fail, retry...");
  }

  //TODO:apend merge_delta

  //TODO:apend delete_index_term_delta

  return false;
};

template<typename KeyType>
bool BWTree::apend_delete(KeyType key){
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
  return mapping_table.set(basic_node->pid,new_delta);
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
      if(key >= ((MergeDelta *)node)->Kp) {
        node= ((MergeDelta *)node)->orignal_node;
        return is_in(key,node);
      }
      return is_in(key, node->next);

    case SPLIT_DELTA:
      if(key >= ((SplitDelta *)node)->Kp) {
        PidType pid= ((SplitDelta *)node)->pQ;
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
BWTree::RecordDelta::RecordDelta(PidType next, BWTree::RecordDelta::RecordType op,
                                 KeyType k, ValueType v)
    : Node(NodeType::RECORD_DELTA){
  RecordDelta(next, op,k);
  this->value = v;
}

// Add your function definitions here
template <typename KeyType, typename ValueType, class KeyComparator>
bool BWTree::insert_entry(KeyType key, ValueType value) {

  std::stack<PidType > path = search(BWTree::root, key);
  if(path.empty()) {
    LOG_ERROR("InsertEntry get empty tree");
  }
  PidType basic_pid = path.top();
  path.pop();

  Node* basic_node = mapping_table.get( basic_pid );
  RecordDelta* new_delta = new RecordDelta(basic_pid, RecordDelta::INSERT, key, value);

  new_delta->high_key = basic_node->high_key;
  new_delta->low_key = basic_node->low_key;

  // check if the leaf node need to be split
  if(  new_delta->need_split() ) {

    KeyType pivotal;

    // create a slibling leaf node
    PidType  new_leaf_pid = create_leaf( new_delta, &pivotal );

    // ceate and prepend a split node
    SplitDelta* new_split = new SplitDelta(new_delta, pivotal, new_leaf_pid);

    // create index-entry-


  }



  return false;
}


template <typename KeyType, typename ValueType, class KeyComparator>
PidType BWTree::create_leaf<typename KeyType, typename ValueType>( Node* new_delta, KeyType* pivotal ){


  // put deltas into stack
  std::stack<Node *> delta_chain;
  Node* tmp_cur_node = new_delta;
  while (tmp_cur_node) {
    delta_chain.push(tmp_cur_node);
    tmp_cur_node = tmp_cur_node -> next;
  }

  // prepare two array to store what logical k-v pairs we have
  KeyType tmpkeys[leafslotmax+1];
  ValueType tmpvals[leafslotmax+1];

  // the first node must be the original leaf node itself
  LeafNode* orig_leaf_node = (LeafNode*)delta_chain.top();
  delta_chain.pop();
  for(int i=0; i< orig_leaf_node->slotuse; i++){
    tmpkeys[i] = orig_leaf_node->slotkey[i];
    tmpvals[i] = orig_leaf_node->slotdata[i];
  }

  // traverse the delta chain
  while( !delta_chain.empty() ){
    Node* cur_delta = delta_chain.top();
    delta_chain.pop();
    switch (cur_delta -> node_type ){
      case RECORD_DELTA:
        if( ((RecordDelta*)cur_delta) ->op_type == RecordDelta::INSERT){

          if( cur_delta->slotuse == 0 ){
            tmpkeys[0] = ((RecordDelta*)cur_delta)->key;
            tmpvals[0] = ((RecordDelta*)cur_delta)->value;
          }
          else {
            int target_pos = 0;
            for (int x = cur_delta->slotuse; x > 0; x--) {
              if ( ((RecordDelta *) cur_delta)->key < tmpkeys[x-1] ) {
                tmpkeys[x] = tmpkeys[x - 1];
                tmpvals[x] = tmpvals[x - 1];
              }
              else{
                target_pos = x;
                break;
              }
            }
            tmpkeys[target_pos] = ((RecordDelta *) cur_delta)->key;
            tmpvals[target_pos] = ((RecordDelta *) cur_delta)->value;
          }

        } else if ( ((RecordDelta*)cur_delta) ->op_type == RecordDelta::DELETE ){

          int target_pos = cur_delta->slotuse - 1;
          for( int x = 0; x < cur_delta->slotuse; x++ ){
            if( tmpkeys[x] == ((RecordDelta *) cur_delta)->key ){
              target_pos = x;
              break;
            }
          }
          for( int x = target_pos; x<cur_delta->slotuse - 1; x++ ){
            tmpkeys[x] = tmpkeys[x+1];
            tmpvals[x] = tmpvals[x+1];
          }

        } else if ( ((RecordDelta*)cur_delta) ->op_type == RecordDelta::UPDATE ){

          for( int x = 0; x < cur_delta->slotuse; x++ ){
            if( tmpkeys[x] == ((RecordDelta *) cur_delta)->key ){
              tmpvals[x] = ((RecordDelta *) cur_delta)-> value;
              break;
            }
          }
        }
        break;
      case SPLIT_DELTA:
        break;
      case MERGE_DELTA:
        break;
      case REMOVE_NODE_DELTA:
        break;
    }
  }

  LeafNode *new_leaf = new LeafNode();
  new_leaf->delta_list_len = 0;
  new_leaf->high_key = orig_leaf_node->high_key;
  for(int i=leafslotmax/2; i<leafslotmax; i++){
    new_leaf->slotkey[i-leafslotmax/2] = tmpkeys[i];
    new_leaf->slotdata[i-leafslotmax/2] = tmpvals[i];
  }
  new_leaf->high_key = orig_leaf_node->high_key;
  new_leaf->low_key = new_leaf->slotdata[0];
  new_leaf->slotuse = (unsigned short) (leafslotmax / 2);

  *pivotal = new_leaf->slotdata[0];
  PidType new_leaf_pid = mapping_table.add(new_leaf);
  return new_leaf_pid;
  // TODO: updated the leaf chain
  // TODO: register in the mapping table

}


}  // End index namespace
}  // End peloton namespace
