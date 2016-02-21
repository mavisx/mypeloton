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

BWTree::BWTree() {
  LeafNode* addr = new LeafNode();
  long newpid = mapping_table.add(addr);
  if (newpid > 0) {
    root = newpid;
  } else {
    LOG_ERROR("Can't create the first leafNode");
  }
}

BWTree::~BWTree() {}

}  // End index namespace
}  // End peloton namespace
