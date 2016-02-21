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
      mappingtable_1[i] = NULL;
    }
    nextPid = 0;
  }

  void* MappingTable::get(PidType pid) {
    long tier1_idx = GET_TIER1_INDEX(pid);
    long tier2_idx = GET_TIER2_INDEX(pid);

    if (mappingtable_1[tier1_idx] == NULL)
      return NULL;

    return mappingtable_1[tier1_idx][tier2_idx];
  }

  bool MappingTable::set(PidType pid, void *addr) {
    long tier1_idx = GET_TIER1_INDEX(pid);
    long tier2_idx = GET_TIER2_INDEX(pid);

    if (mappingtable_1[tier1_idx] == NULL) {
     return false;
    }

    void* expected = mappingtable_1[tier1_idx][tier2_idx];
    return true;
  }


  bool MappingTable::add(void * addr) {
    unsigned long new_pid = nextPid++;

  }

  bool MappingTable::remove(PidType pid) {
    long tier1_idx = GET_TIER1_INDEX(pid);
    long tier2_idx = GET_TIER2_INDEX(pid);

  }


  BWTree::BWTree() {

  }

}  // End index namespace
}  // End peloton namespace
