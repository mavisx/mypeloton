//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// BWTree.h
//
// Identification: src/backend/index/BWTree.h
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <stddef.h>
#include <assert.h>
#include <atomic>
#include <unordered_map>
#include <stack>
#include <vector>
#include <unordered_set>
#include "../common/types.h"
#include "../storage/tuple.h"

// in bytes
#define BWTREE_NODE_SIZE 64

#define BWTREE_MAX(a, b) ((a) < (b) ? (b) : (a))

#define MAX_DELTA_CHAIN_LEN 8

// in bits
#define MAPPING_TABLE_SIZE_BITNUM 10
#define MAPPING_TABLE_SIZE (1 << (MAPPING_TABLE_SIZE_BITNUM))

#define GET_TIER1_INDEX(pid) ((pid) >> 10)
#define GET_TIER2_INDEX(pid) ((pid)&0x3ff)

#define NULL_PID -1

namespace peloton {
namespace index {

typedef long long PidType;

enum NodeType {
  LEAF = 0,
  INNER = 1,
  RECORD_DELTA = 2,
  INDEX_ENTRY_DELTA = 3,
  REMOVE_NODE_DELTA = 4,
  MERGE_DELTA = 5,
  DELETE_INDEX_TERM_DELTA = 6,
  SPLIT_DELTA = 7
};


// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
template <typename KeyType, typename ValueType, class KeyComparator,
          typename KeyEqualityChecker, class ValueEqualityChecker>
class BWTree {
 public:
  // *** Constructed Types

  /// Typedef of our own type
  typedef BWTree<KeyType, ValueType, KeyComparator, KeyEqualityChecker>
      BwTreeSelf;

  /// Size type used to count keys
  typedef size_t SizeType;

  /// The pair of KeyType and data_type, this may be different from value_type.
  typedef std::pair<KeyType, ValueType> PairType;

 public:
  // *** Static Constant Options and Values of the Bw Tree
  /// Base B2 tree parameter: The number of key/data slots in each leaf
  static const unsigned short leafslotmax =
      BWTREE_MAX(8, BWTREE_NODE_SIZE / (sizeof(KeyType) + sizeof(ValueType)));

  /// Base B+ tree parameter: The number of key slots in each inner node,
  /// this can differ from slots in each leaf.
  static const unsigned short innerslotmax =
      BWTREE_MAX(8, BWTREE_NODE_SIZE / (sizeof(KeyType) + sizeof(PidType)));

  /// Computed B+ tree parameter: The minimum number of key/data slots used
  /// in a leaf. If fewer slots are used, the leaf will be merged or slots
  /// shifted from it's siblings.
  static const unsigned short minleafslots = (leafslotmax / 2);

  /// Computed B+ tree parameter: The minimum number of key slots used
  /// in an inner node. If fewer slots are used, the inner node will be
  /// merged or slots shifted from it's siblings.
  static const unsigned short mininnerslots = (innerslotmax / 2);

  struct Node;
  struct MappingTable {
   private:
    Node** mappingtable_1[MAPPING_TABLE_SIZE];
    std::atomic<unsigned long> nextPid;

   public:
    MappingTable() {
      for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
        mappingtable_1[i] = nullptr;
      }
      nextPid = 0;
    }

    ~MappingTable() {
      for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
        if (mappingtable_1[i] != nullptr) {
          delete[] mappingtable_1[i];
          mappingtable_1[i] = nullptr;
        }
      }
    }

    Node* get(PidType pid) {
      if (pid == NULL_PID)
        return nullptr;

      long tier1_idx = GET_TIER1_INDEX(pid);
      long tier2_idx = GET_TIER2_INDEX(pid);

      if (mappingtable_1[tier1_idx] == nullptr) return nullptr;

      return mappingtable_1[tier1_idx][tier2_idx];
    }

    bool set(PidType pid, Node* addr, Node* expected) {
      long tier1_idx = GET_TIER1_INDEX(pid);
      long tier2_idx = GET_TIER2_INDEX(pid);

      if (mappingtable_1[tier1_idx] == nullptr) {
        return false;
      }

      // using cas to set the value in mapping table
      return std::atomic_compare_exchange_strong(
          (std::atomic<Node*>*)&mappingtable_1[tier1_idx][tier2_idx], &expected,
          addr);
    }


    // if correctly add this new addr, return the new pid,
    // else return -1
    long add(Node* addr) {
      unsigned long new_pid = nextPid++;
      long tier1_idx = GET_TIER1_INDEX(new_pid);
      long tier2_idx = GET_TIER2_INDEX(new_pid);

      Node** expectded = mappingtable_1[tier1_idx];

      // atomically add new secondary index
      if (mappingtable_1[tier1_idx] == nullptr) {
        Node** desired = new Node* [MAPPING_TABLE_SIZE];
        for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
          desired[i] = nullptr;
        }
        if (!std::atomic_compare_exchange_strong(
            (std::atomic<Node**>*)&mappingtable_1[tier1_idx], &expectded,
            desired)) {
          delete[] desired;
        }
      }

      Node* expectded2 = mappingtable_1[tier1_idx][tier2_idx];
      // atomically add addr to desired secondary index
      if (std::atomic_compare_exchange_strong(
          (std::atomic<Node*>*)&mappingtable_1[tier1_idx][tier2_idx],
          &expectded2, addr)) {
        return new_pid;
      } else {
        return -1;
      }
    }

    void remove(PidType pid) {
      long tier1_idx = GET_TIER1_INDEX(pid);
      long tier2_idx = GET_TIER2_INDEX(pid);

      mappingtable_1[tier1_idx][tier2_idx] = nullptr;
    }
  };


 private:
  // We need a root node
  PidType root;

  MappingTable mapping_table;

  /// Pointer to first leaf in the double linked leaf chain
  PidType headleaf;

  /// Pointer to last leaf in the double linked leaf chain
  PidType tailleaf;

 public:
  /**
   * The Node inheritance hierachy
   * **/
  struct Node {
    // reference to outer mapping table
    MappingTable& mapping_table;
    /// Number of key slotuse use, so number of valid children or data
    /// pointers
    unsigned short slotuse;

    // flag to indicate whether this chain belongs to a leaf node
    bool is_leaf;

    // Delta chain next pointer
    Node* next;

    // Node pid
    PidType pid;

    // Length of current delta chain
    size_t delta_list_len;

    // type of this node
    NodeType node_type;

    // minimal and maximal key in this node
    KeyType low_key, high_key;

    // Double linked list pointers to traverse the leaves
    //Node* prev_node;

    // Double linked list pointers to traverse the leaves
    PidType next_leafnode;

      // constructor
    Node(NodeType n, size_t delta_l, MappingTable& mt, PidType next_leaf) : mapping_table(mt) {
      node_type = n;
      slotuse = 0;
      next = nullptr;
      delta_list_len = delta_l;
      next_leafnode = next_leaf;
    }

    // True if this is a leaf node
    inline bool is_leaf_node() const { return (node_type == NodeType::LEAF); }

    inline bool need_split() const {
      if (is_leaf) {
        return slotuse >= leafslotmax;
      } else {
        return slotuse >= innerslotmax;
      }
    }

    inline bool need_merge() const {
      if (is_leaf) {
        return slotuse <= minleafslots;
      } else {
        return slotuse <= mininnerslots;
      }
    }
  };

 private:
  /// Extended structure of a inner node in-memory. Contains only keys and no
  /// data items.
  struct InnerNode : public Node {
    /// Define an related allocator for the inner_node structs.
    // typedef typename _Alloc::template rebind<inner_node>::other alloc_type

    /// Keys of children or data pointers,
    //  we plus one so as to avoid overflow when consolidation
    KeyType slotkey[innerslotmax + 1];

    /// Pointers to children,
    //  we plus one so as to avoid overflow when consolidation
    PidType childid[innerslotmax + 1 + 1];

    /// Set variables to initial values
    InnerNode(MappingTable& mapping_table, PidType next_leafnode)
        : Node(NodeType::INNER, 0, mapping_table, next_leafnode) {}

    /// True if the node's slots are full
    inline bool isfull() const { return (Node::slotuse == innerslotmax); }

    /// True if few used entries, less than half full
    inline bool isfew() const { return (Node::slotuse <= mininnerslots); }

    /// True if node has too few entries
    inline bool isunderflow() const { return (Node::slotuse < mininnerslots); }
  };

  /// Extended structure of a leaf node in memory. Contains pairs of keys and
  /// data items. Key and data slots are kept in separate arrays, because the
  /// key array is traversed very often compared to accessing the data items.
  struct LeafNode : public Node {
    /// Define an related allocator for the leaf_node structs.
    //    typedef typename _Alloc::template rebind<LeafNode>::other alloc_type;


    /// Keys of children or data pointers
    //  we plus one so as to avoid overflow when consolidation
    KeyType slotkey[leafslotmax + 1];

    /// Array of data, each bucket points to an vector of actual data
    //  we plus one so as to avoid overflow when consolidation
    std::vector<ValueType>* slotdata[leafslotmax + 1];

    LeafNode(MappingTable& mapping_table, PidType next_leafnode )
        : Node(NodeType::LEAF, 0, mapping_table, next_leafnode ){
      // initialize the value bucket
      for (int i = 0; i < leafslotmax + 1; i++) {
        slotdata[i] = nullptr;
      }
    }

    /// True if the node's slots are full
    inline bool isfull() const { return (Node::slotuse == leafslotmax); }

    /// True if few used entries, less than half full
    inline bool isfew() const { return (Node::slotuse <= minleafslots); }

    /// True if node has too few entries
    inline bool isunderflow() const { return (Node::slotuse < minleafslots); }

    /// Set the (key,data) pair in slot. Overloaded function used by
    /// bulk_load().
    inline void set_slot(unsigned short slot, const PairType& value) {
      assert(slot < Node::slotuse);
      slotkey[slot] = value.first;
      if (slotdata[slot] == nullptr) {
        slotdata[slot] = new std::vector<ValueType>();
      }
      slotdata[slot]->push_back(value.second);
    }

    /// Set the key pair in slot. Overloaded function used by
    /// bulk_load().
    inline void set_slot(unsigned short slot, const KeyType& key) {
      assert(slot < Node::slotuse);
      slotkey[slot] = key;
    }
  };

  // Delta Node for record update operation
  struct RecordDelta : public Node {
    // construction added -mavis
    enum RecordType { INSERT = 0, DELETE = 1, UPDATE = 2 };


    RecordDelta(PidType next, RecordType op, KeyType k, ValueType v,
                MappingTable& mapping_table,  PidType next_leafnode)
        : Node(NodeType::RECORD_DELTA, 0, mapping_table, next_leafnode ) {
      op_type = op;
      key = k;
      // Get node* of original node form mapping_table
      Node* orig_node = mapping_table.get(next);

      prepend(this, orig_node);

      // update the slotuse of the new delta node
     /* if (op_type == INSERT) {
        if( !key_is_in(k, orig_node) )
          this->slotuse = (unsigned short)(orig_node->slotuse + 1);
      }*/

      this->value = v;
    }
    RecordType op_type;
    KeyType key;
    ValueType value;
  };

  // Delta Node for spliting operation
  struct SplitDelta : public Node {
    SplitDelta(Node* next, KeyType Kp, PidType pQ, MappingTable& mapping_table,  PidType next_leafnode )
        : Node(NodeType::SPLIT_DELTA, 0, mapping_table, next_leafnode ), Kp(Kp), pQ(pQ) {}
    KeyType Kp;
    PidType pQ;
  };

  struct IndexEntryDelta : public Node {
    IndexEntryDelta(Node* next, KeyType Kp, KeyType Kq, PidType pQ,
                    MappingTable& mapping_table, PidType next_leafnode)
        : Node(NodeType::INDEX_ENTRY_DELTA, 0, mapping_table, next_leafnode),
          Kp(Kp),
          Kq(Kq),
          pQ(pQ) {}
    KeyType Kp, Kq;
    PidType pQ;
  };

  // Delta Node for merging operation
  struct RemoveDelta : public Node {
    RemoveDelta(Node* next, MappingTable& mapping_table, PidType next_leafnode)
        : Node(NodeType::REMOVE_NODE_DELTA, 0, mapping_table, next_leafnode ) {}
  };

  struct MergeDelta : public Node {
    MergeDelta(Node* next, KeyType Kp, Node* orignal_node,
               MappingTable& mapping_table,  PidType next_leafnode)
        : Node(NodeType::MERGE_DELTA, 0, mapping_table, next_leafnode ),
          Kp(Kp),
          orignal_node(orignal_node) {}
    KeyType Kp;
    Node* orignal_node;
  };

  struct DeleteIndexDelta : public Node {
    DeleteIndexDelta(Node* next, KeyType Kp, KeyType Kq, PidType pQ,
                     MappingTable& mapping_table, PidType next_leafnode)
        : Node(NodeType::INDEX_ENTRY_DELTA, 0, mapping_table, next_leafnode),
          Kp(Kp),
          Kq(Kq),
          pQ(pQ) {}
    KeyType Kp, Kq;
    PidType pQ;
  };

 public:
  // constructor

  BWTree(const KeyComparator& kc, const KeyEqualityChecker& ke, const ValueEqualityChecker& ve)
      : m_key_less(kc), m_key_equal(ke), m_value_equal(ve) {
    LeafNode* addr = new LeafNode(mapping_table, NULL_PID );
   long newpid = mapping_table.add(addr);
    if (newpid >= 0) {
      // initialize the root, head and tail pid.
      root = newpid;
      headleaf = tailleaf = newpid;
    } else {
      LOG_ERROR("Can't create the initial leafNode!");
    }
  }

  // destructor
  ~BWTree(){};

  /*
   ************************************************
   *    public method exposed to users -leiqi     *
   ************************************************
   */
 public:
  // True if a < b ? "constructed" from m_key_less()
  inline bool key_less(const KeyType& a, const KeyType b) const {
    return m_key_less(a, b);
  }

  // True if a <= b ? constructed from key_less()
  inline bool key_lessequal(const KeyType& a, const KeyType b) const {
    return !m_key_less(b, a);
  }

  // True if a > b ? constructed from key_less()
  inline bool key_greater(const KeyType& a, const KeyType& b) const {
    return m_key_less(b, a);
  }

  // True if a >= b ? constructed from key_less()
  inline bool key_greaterequal(const KeyType& a, const KeyType b) const {
    return !m_key_less(a, b);
  }

  // True if a == b ? constructed from key_less().
  inline bool key_equal(const KeyType& a, const KeyType& b) const {
    return m_key_equal(a, b);
  }

  inline bool value_equal(const ValueType& a, const ValueType& b) const {
    return m_value_equal(a, b);
  }
  /*
    ************************************************
    * private functions, invisible to users -leiqi *
    ************************************************
    */

 private:
  KeyComparator m_key_less;
  KeyEqualityChecker m_key_equal;
  ValueEqualityChecker m_value_equal;

  std::stack<PidType> search(PidType pid, KeyType key) {
    auto node = mapping_table.get(pid);
    if (node == nullptr) {
      std::stack<PidType> empty_res;
      return empty_res;
    }

    std::stack<PidType> path;
    path.push(pid);
    PidType res = search(node, key, path);
    if (res == -1) {
      std::stack<PidType> empty_res;
      return empty_res;
    }

    return path;
  }

  PidType search(Node* node, KeyType key, std::stack<PidType>& path) {
    // should always keep track of right key range even in delta node
    PidType pid;
    //    if(key < node->low_key || key >= node->high_key) {
    //      LOG_ERROR("Search Range Err: key not in range");
    //      return -1;
    //    }
    switch (node->node_type) {
      case LEAF:
      case RECORD_DELTA:
        return node->pid;

      case INDEX_ENTRY_DELTA:
      case DELETE_INDEX_TERM_DELTA:
        if (key_greaterequal(key, ((IndexEntryDelta*)node)->Kp) &&
            key_less(key, ((IndexEntryDelta*)node)->Kq)) {
          pid = ((IndexEntryDelta*)node)->pQ;
          node = mapping_table.get(pid);
          if (node == nullptr) {
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
        if (path.empty()) {
          LOG_INFO("Search Path empty");
          return -1;
        }
        pid = path.top();

        node = mapping_table.get(pid);
        if (node == nullptr) {
          LOG_ERROR("Pid in split/merge delta not exist");
          return -1;
        } else {
          return search(node, key, path);
        }
      case MERGE_DELTA:
        LOG_INFO("Search Range Info: meet merge delta");

        if (key_greaterequal(key, ((MergeDelta*)node)->Kp)) {
          node = ((MergeDelta*)node)->orignal_node;
          if (node == nullptr) {
            LOG_ERROR("pid in split delta not exist");
            return -1;
          }

          return search(node, key, path);
        }
        return search(node->next, key, path);
      case SPLIT_DELTA:
        LOG_INFO("Search Range Info: meet split/merge delta");
        pid = ((SplitDelta*)node)->pQ;
        if (key_greaterequal(key, ((SplitDelta*)node)->Kp)) {
          node = mapping_table.get(pid);
          if (node == nullptr) {
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
        if (node->slotuse == 0) {
          LOG_ERROR("empty inner node");
          return -1;
        } else {
          int i = 0;
          for (i = 0; i < node->slotuse; i++) {
            if (key_greaterequal(key, ((InnerNode*)node)->slotkey[i]))
              continue;
            else
              break;
          }
          pid = ((InnerNode*)node)->childid[i];
          node = mapping_table.get(pid);
          if (node == nullptr) {
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


  typedef std::unordered_set<ValueType,std::hash<ValueType>, m_value_equal> DelSet;

  inline bool key_is_in(KeyType key, Node* listhead) {
    DelSet deleted_set;
    return key_is_in(key, listhead, deleted_set);
  };

  bool key_is_in(KeyType key, Node* listhead, DelSet& deleted) {
    if (listhead == nullptr) return false;

    Node* node = listhead;
    switch (node->node_type) {
      case RECORD_DELTA:
        RecordDelta* rcd_node = (RecordDelta*)node;
        if (rcd_node->op_type == RecordDelta::INSERT
            && key_equal(rcd_node->key, key)
            && (!deleted.count(rcd_node->value))) {
          return true;
        }
        else if (rcd_node->op_type == RecordDelta::DELETE
                 && key_equal(rcd_node->key, key)) {
          deleted.insert(rcd_node->value);
        }
        return key_is_in(key, node->next, deleted);
      case LEAF:
        LeafNode* lf_node = (LeafNode*)node;
        for (int i = 0; i < (lf_node->slotuse); i++) {
          if (key_equal(lf_node->slotkey[i], key)) {
            for(auto val:(*lf_node->slotdata[i])){
              if(!deleted.count(val)) {
                return true;
              }
            }
            return false;
          }
        }
        return false;
      case MERGE_DELTA:
        if (key_greaterequal(key, ((MergeDelta*)node)->Kp)) {
          node = ((MergeDelta*)node)->orignal_node;
          return key_is_in(key, node, deleted);
        }
        return key_is_in(key, node->next, deleted);

      case SPLIT_DELTA:
        if (key_greaterequal(key, ((SplitDelta*)node)->Kp)) {
          PidType pid = ((SplitDelta*)node)->pQ;
          return key_is_in(key, mapping_table.get(pid), deleted);
        }
        return key_is_in(key, node->next, deleted);
      default:
        return false;
    }

  }

  std::pair<int, int> count_pair(KeyType key, Value value, Node* listhead) {
    int total_count = 0;
    int pair_count = 0;
    DelSet deleted_set;
    count_pair(key, value, listhead, total_count, pair_count, deleted_set);
    return std::pair<int, int>(total_count, pair_count);
  };

  //return pair nums, also calculates total val nums of the key in count
  int count_pair(KeyType key, Value value, Node* listhead,
                 int& total_count, int& pair_count, DelSet& deleted) {
    if (listhead == nullptr) return pair_count;

    Node* node = listhead;
    switch (node->node_type) {
      case RECORD_DELTA:
        RecordDelta* rcd_node = (RecordDelta*)node;
        if (rcd_node->op_type == RecordDelta::INSERT
            && key_equal(rcd_node->key, key)
            && (!deleted.count(rcd_node->value))) {
          total_count++;
          if(value_equal(rcd_node->value, value)) {
            pair_count++;
          }
        }
        else if (rcd_node->op_type == RecordDelta::DELETE &&
                   key_equal(rcd_node->key, key)) {
          deleted.insert(rcd_node->value);
        }
        break;
      case LEAF:
        LeafNode* lf_node = (LeafNode*)node;
        for (int i = 0; i < (lf_node->slotuse); i++) {
          if (key_equal(lf_node->slotkey[i], key)) {
            for (ValueType v: *(lf_node->slotdata[i])){
              if(deleted.count(v)) continue;
              total_count++;
              if(value_equal(v, value))
                pair_count++;
            }
            break;
          }
        }
        return pair_count;
      case MERGE_DELTA:
        if (key_greaterequal(key, ((MergeDelta*)node)->Kp)) {
          node = ((MergeDelta*)node)->Kp;
          return count_pair(key, value, node, total_count, pair_count, deleted);
        }
        break;

      case SPLIT_DELTA:
        if (key_greaterequal(key, ((SplitDelta*)node)->Kp)) {
          PidType pid = ((SplitDelta*)node)->pQ;
          return count_pair(key, value, mapping_table.get(pid),
                            total_count, pair_count, deleted);
        }
        break;
    }
    count_pair(key, value, node->next, total_count, pair_count, deleted);
  };

  bool append_delete( Node* basic_node, KeyType key, ValueType value, bool deletekey){

    RecordDelta* new_delta = new RecordDelta(basic_node->pid, RecordDelta::DELETE, key, value,
                        mapping_table, basic_node->next_leafnode);

    if(deletekey) {
      new_delta->slotuse -= 1;
    }

    return mapping_table.set(basic_node->pid, basic_node, new_delta);
  }

  bool apend_merge()
  {
    //TODO: A lot
    return false;
  }
  /*
   ************************************************
   *               end -leiqi                     *
   ************************************************
   */

  // private fuctions, invisible to users -mavis
  inline static bool prepend(Node* delta_node, Node* orig_node) {
    // update the delta_list_len of the delta node
    delta_node->delta_list_len = orig_node->delta_list_len + 1;

    // update the slotuse of the new delta node
    delta_node->slotuse = orig_node->slotuse;

    // maintain next, prev pointer
    delta_node->next = orig_node;

    delta_node->low_key = orig_node->low_key;
    delta_node->high_key = orig_node->high_key;

    return true;
  }
  // end -mavis

 public:


  void get_value(KeyType key, std::vector<ValueType>& result) {
    std::stack<PidType> path = search(root, key);
    if (path.empty()) {
      LOG_INFO("There is no result available in this tree");
      return;
    }

    PidType target_node = path.top();
    Node * next = mapping_table.get(target_node);

    if (!next->is_leaf) {
      LOG_ERROR("get_value's search result is not a leaf");
      return;
    }

    // traverse the delta chain from top down
    // to construct the correct result
    while (next) {
      switch (next->node_type) {
        case RECORD_DELTA : {
          RecordDelta *node = dynamic_cast<RecordDelta *>(next);
          if (key_equal(node->key, key)) {
            if (node->op_type == RecordDelta::RecordType::INSERT) {
              result.push_back(node->value);
            } else if (node->op_type == RecordDelta::RecordType::DELETE) {
              // if we meet a "delete" all later value associated with this key is invalid.
              return;
            }
          }
          next = node->next;
        }
        case LEAF: {
          LeafNode * leaf = dynamic_cast<LeafNode *>(next);
          for (int i=0; i < leaf->slotuse; i++) {
            if (key_equal(leaf->slotkey[i], key)) {
              result.insert(result.end(), leaf->slotdata[i]->begin(), leaf->slotdata[i]->end());
            }
          }
          next = nullptr;
        }
        case SPLIT_DELTA: {
          SplitDelta *split_delta = dynamic_cast<SplitDelta *>(next);
          // if key >= Kp, we go to the new split node
          if (key_greaterequal(key, split_delta->Kp)) {
            LOG_ERROR("search result direct to another node");
            next = mapping_table.get(split_delta->pQ);
          }
          // else we go alone this delta chain
          else {
            next = split_delta->next;
          }
        }
        case MERGE_DELTA: {
          MergeDelta *merge_delta = dynamic_cast<MergeDelta *>(next);
          // if key >= Kp, we go to the original node
          if (key_greaterequal(key, merge_delta->Kp)) {
            next = merge_delta->orignal_node;
          }
          // else we go alone this delta chain
          else {
            next = merge_delta->next;
          }
        }
        case REMOVE_NODE_DELTA:
          // if we meet remove node delta, we can search from the root again.
          result.clear();
          get_value(key, result);
        default:
          LOG_ERROR("meet wrong delta: %d during get_value", next->node_type);
      }
    }
  }

  // public method exposed to users -mavis
  bool insert_entry(KeyType key, ValueType value) {
    std::stack<PidType> path = search(BWTree::root, key);
    if (path.empty()) {
      LOG_ERROR("InsertEntry get empty tree");
      return false;
    }
    PidType basic_pid = path.top();
    path.pop();

    Node* basic_node = mapping_table.get(basic_pid);


    SplitDelta* new_split;
    // check if the leaf node need to be split before we add record delta
    if ( basic_node -> need_split() ) {
      KeyType pivotal;

      // create a slibling leaf node
      PidType new_leaf_pid = create_leaf(basic_pid, &pivotal);

      // ceate and prepend a split node
      new_split = new SplitDelta(basic_node, pivotal,
                                             new_leaf_pid, mapping_table, new_leaf_pid);

      if ( !mapping_table.set(basic_pid, new_split, basic_node)) {
        //TODO: if the CAS fails, release the new_split obj
        delete new_split;
      }

    }

    RecordDelta* new_delta = new RecordDelta(basic_pid, RecordDelta::INSERT,
                                             key, value, mapping_table, basic_node->next_leafnode);
    if( !key_is_in(key, new_split) )
      new_delta->slotuse = new_delta->next->slotuse+1;

    std::stack<PidType> now_path = search(BWTree::root, key);
    if ( now_path.empty() ) {
      LOG_ERROR("InsertEntry get empty tree");
      return false;
    }
    basic_pid = now_path.top();
    basic_node = mapping_table.get(basic_pid);
    new_delta->high_key = basic_node->high_key;
    new_delta->low_key = basic_node->low_key;

    //TODO: use CAS concatenate this new_delta to the delta chain
    if ( !mapping_table.set(basic_pid, new_delta, basic_node) ){
      return false;
    }

    return true;
  }

  bool delete_entry(KeyType key, ValueType value) {
    bool redo =true;

    //check and insert delete delta
    while(redo) {
      std::stack<PidType> path = search(BWTree::root, key);
      if (path.empty()) {
        LOG_ERROR("InsertEntry get empty tree");
        return false;
      }
      PidType basic_pid = path.top();
      path.pop();

      Node* basic_node = mapping_table.get(basic_pid);
      auto tv_count_pair = count_pair(key,value, basic_node);
      if (!tv_count_pair.second) {
        LOG_INFO("DeleteEntry Not Exist");
        return false;
      }

      bool deletekey = (tv_count_pair.second>=tv_count_pair.first);

      redo = append_delete(basic_node, key, value, deletekey);
    }
    // TODO:apend merge_delta

    // TODO:apend delete_index_term_delta

    return !redo;
  };

  bool update_entry(KeyType key, ValueType value);

  std::pair<std::vector<KeyType>*, std::vector<std::vector<ValueType>*>*>
      fake_consolidate( Node* new_delta) {

    std::stack<Node*> delta_chain;
    Node* tmp_cur_node = new_delta;
    while (tmp_cur_node) {
      delta_chain.push(tmp_cur_node);
      tmp_cur_node = tmp_cur_node->next;
    }

    // prepare two array to store what logical k-v pairs we have
    std::vector<KeyType>* tmpkeys = new std::vector<KeyType>();
    std::vector<std::vector<ValueType>*>* tmpvals = new std::vector<std::vector<ValueType>*>();

    // the first node must be the original leaf node itself
    LeafNode* orig_leaf_node = dynamic_cast<LeafNode*>(delta_chain.top());
    delta_chain.pop();

    // copy the data in the base node
    for (int i = 0; i < orig_leaf_node->slotuse; i++) {
      tmpkeys[i] = orig_leaf_node->slotkey[i];
      tmpvals = new std::vector<ValueType>(*(orig_leaf_node->slotdata[i]) );
    }

    // traverse the delta chain
    while (!delta_chain.empty()) {
      // get top delta node
      Node* cur_delta = delta_chain.top();
      delta_chain.pop();
      switch (cur_delta->node_type) {

        case RECORD_DELTA: {
          // first see the key has already existed
          bool no_dedup = true;
          RecordDelta *recordDelta = dynamic_cast<RecordDelta *>(cur_delta);

          if (recordDelta->op_type == RecordDelta::INSERT) {
            for (int x = 0; x < recordDelta->slotuse; x++) {
              if (key_equal(tmpkeys[x], recordDelta->key)) {
                (*tmpvals)[x]->push_back(recordDelta->value);
                no_dedup = false;
                break;
              }
            }
            // key not exists, need to insert somewhere
            if (no_dedup) {
              if (recordDelta->slotuse == 0) {
                tmpkeys[0] = recordDelta->key;
                (*tmpvals)[0] = new std::vector<ValueType>();
                (*tmpvals)[0]->push_back(recordDelta->value);
              } else {
                int target_pos = 0;
                for (int x = recordDelta->slotuse; x > 0; x--) {
                  if (key_less(recordDelta->key, tmpkeys[x - 1])) {
                    tmpkeys[x] = tmpkeys[x - 1];
                    tmpvals[x] = tmpvals[x - 1];
                  } else {
                    target_pos = x;
                    break;
                  }
                }
                tmpkeys[target_pos] = recordDelta->key;
                (*tmpvals)[target_pos] = new std::vector<ValueType>();
                (*tmpvals)[target_pos]->push_back(recordDelta->value);
              }
            } // end of RecordDelta::INSERT


          } else if (recordDelta->op_type == RecordDelta::DELETE) {

            int target_pos = recordDelta->slotuse - 1;
            bool need_remove = false;
            for (int x = 0; x < recordDelta->slotuse; x++) {
              if (key_equal(tmpkeys[x], recordDelta->key)) {
                // remove value in the vector
                (*tmpvals)[x]->erase( recordDelta->);

                // if vector is empty, needed to be removed
                if( (*tmpvals)[x]->size() == 0 ){
                  need_remove = true;
                }
                target_pos = x;
                break;
              }
            }

            if( need_remove ) {
              delete tmpvals[target_pos];
              for (int x = target_pos; x < cur_delta->slotuse - 1; x++) {
                tmpkeys[x] = tmpkeys[x + 1];
                tmpvals[x] = tmpvals[x + 1];
              }
            }

          }
          break;
        }
        case SPLIT_DELTA:
          break;
        case MERGE_DELTA:
          break;
        case REMOVE_NODE_DELTA:
          break;
      }
    }

    return std::make_pair(tmpkeys, tmpvals);

  }

  PidType create_leaf(PidType new_delta_pid, KeyType* pivotal) {

    Node* new_delta = mapping_table.get(new_delta_pid);
    LeafNode* new_leaf = new LeafNode(mapping_table, new_delta->next_leafnode);
    new_leaf->delta_list_len = 0;
    new_leaf->high_key = new_delta->high_key;

    std::pair< KeyType*, std::vector<ValueType>** > arrays = fake_consolidate(new_delta);

    for (int i = leafslotmax / 2; i < leafslotmax; i++) {
      new_leaf->slotkey[i - leafslotmax / 2] = arrays.first[i];
      new_leaf->slotdata[i - leafslotmax / 2] = arrays.second[i];
    }
    new_leaf->low_key = new_leaf->slotkey[0];
    new_leaf->slotuse = (unsigned short)(leafslotmax / 2);

    *pivotal = new_leaf->slotkey[0];
    PidType new_leaf_pid = mapping_table.add(new_leaf);


    // TODO: left right pointer
    return new_leaf_pid;

  }

  void scan_all(std::vector<ValueType>& v) {
    Node * node = mapping_table.get(headleaf);

    // scan the leaf nodes list from begin to the end
    while (node != nullptr) {
      auto all_key_value_pair =  fake_consolidate(node);

      v.insert(v.end(), all_key_value_pair.second->begin(), all_key_value_pair.second->end());

      node = mapping_table.get(node->next_leafnode);
    }
  }

  void scan(std::vector<KeyType> &keys_result,
      std::vector<ItemPointer> &values_result) {
    Node * node = mapping_table.get(headleaf);

    // scan the leaf nodes list from begin to the end
    while (node != nullptr) {
      auto all_key_value_pair =  fake_consolidate(node);
      auto keys = *all_key_value_pair.first;
      auto values = *all_key_value_pair.second;

      keys_result.insert(keys_result.end(), keys.begin(), keys.end());
      values_result.insert(values_result.end(), values.begin(), values.end());

      node = mapping_table.get(node->next_leafnode);
    }
  }

};

}  // End index namespace
}  // End peloton namespace
