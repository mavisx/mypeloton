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
    bool BWTree::InsertEntry(KeyType key, ValueType value) {

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


        if( ((LeafNode*)new_delta)->isfull() ) {

            std::stack<Node *> delta_chain;
            Node* tmp_cur_node = new_delta;
            while (tmp_cur_node) {
                delta_chain.push(tmp_cur_node);
                tmp_cur_node = tmp_cur_node -> next;
            }

            KeyType tmpkeys[leafslotmax+1];
            ValueType tmpvals[leafslotmax+1];
            LeafNode* orig_leaf_node = (LeafNode*)delta_chain.top();
            delta_chain.pop();
            for(int i=0; i< leafslotmax; i++){
                tmpkeys[i] = orig_leaf_node->slotkey[i];
                tmpvals[i] = orig_leaf_node->slotdata[i];
            }
            while( !delta_chain.empty() ){
                Node* cur_delta = delta_chain.top();
                delta_chain.pop();
                switch (cur_delta -> node_type ){
                    case RECORD_DELTA:
                        if( ((RecordDelta*)cur_delta) ->op_type == RecordDelta::INSERT){

                            for(int x=leafslotmax; x>=0; x--){
                                if( ((RecordDelta*)cur_delta)->key > )
                                tmpkeys[x+1] = tmpkeys[x];
                            }

                        } else if ( ((RecordDelta*)cur_delta) ->op_type == RecordDelta::DELETE ){

                        } else if ( ((RecordDelta*)cur_delta) ->op_type == RecordDelta::UPDATE ){

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