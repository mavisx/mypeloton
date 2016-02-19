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

        PidType basic_node_pid = Search(BwTree::root, key);
        auto itr = mapping_table.find( basic_node_pid );
        if( itr == mapping_table.end() )
            return false;
        Node* basic_node = itr->second;
        RecordDelta* new_delta = new RecordDelta();
        new_delta->op_type = RecordDelta::INSERT;
        new_delta->value = value;
        new_delta->high_key = basic_node->high_key;
        new_delta->low_key = basic_node->low_key;
        if( Prepend(new_delta, basic_node_pid) )
            return false;

        Node* tmp_cur_node = new_delta;
        while( tmp_cur_node -> slotuse > innerslotmax){


        }




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

        // maintain next, prev pointer
        delta_node -> next = orig_node;
        orig_node -> prev = delta_node;

        // CAS operation
        while(std::atomic_compare_exchange_weak<Node*>
                ( &(itr->second), itr->second, delta_node));


        return true;
    }



}  // End index namespace
}  // End peloton namespace
