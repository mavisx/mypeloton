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
// Mavis
//===----------------------------------------------------------------------===//


#include "bwtree.h"
#include <atomic>

namespace peloton {
namespace index {

    // Add your function definitions here
    template <typename KeyType, typename ValueType, class KeyComparator>
    bool BwTree::InsertEntry(KeyType key, ValueType value) {


    }

    template <typename KeyType, typename ValueType, class KeyComparator>
    bool BwTree::Prepend( Node* delta_node, PidType orig_pid) {

        // Get node* of original node form mapping_table
        auto itr = mapping_table.find( orig_pid );
        if( itr == mapping_table.end() )
            return false;
        Node* orig_node = itr->second;

        // update the delta_list_len of the delta node
        delta_node->delta_list_len = orig_node->delta_list_len + 1;

        // update the slotuse of the new delta node
        if( delta_node->node_type == RECORD_DELTA ){
            switch ( ((struct RecordDelta*)delta_node)->op_type ){
                case RecordDelta::INSERT :
                    delta_node->slotuse = orig_node->slotuse + 1;
                    break;
                case RecordDelta::DELETE :
                    delta_node->slotuse = orig_node->slotuse - 1;
                    break;
                case RecordDelta::UPDATE :
                    delta_node->slotuse = orig_node->slotuse;
                    break;
            }
        }

        // CAS operation
        std::atomic<Node *> myCAS;
        while(!myCAS.compare_exchange_weak(orig_node, delta_node));

        return true;
    }


}  // End index namespace
}  // End peloton namespace
