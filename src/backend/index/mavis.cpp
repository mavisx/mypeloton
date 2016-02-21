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


namespace peloton {
namespace index {

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
        if( ((LeafNode*)new_delta)->isfull() ) {

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
            new_leaf->slotuse = leafslotmax/2;

            // TODO: updated the leaf chain
            // TODO: register in the mapping table



        }
    }

    template <typename KeyType, typename ValueType, class KeyComparator>
    static bool BWTree::prepend( Node* delta_node, PidType orig_pid) {

        // Get node* of original node form mapping_table
        Node* orig_node = mapping_table.get(orig_pid);

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
        mapping_table.set(orig_pid, delta_node);


        return true;
    }



}  // End index namespace
}  // End peloton namespace