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
#include "index.h"

//#define MY_PRINT_DEBUG

// in bytes
#define BWTREE_NODE_SIZE 4096

#define BWTREE_MAX(a, b) ((a) < (b) ? (b) : (a))

#define MAX_DELTA_CHAIN_LEN 8

// in bits
#define MAPPING_TABLE_SIZE_BITNUM 10
#define MAPPING_TABLE_SIZE (1 << (MAPPING_TABLE_SIZE_BITNUM))

#define GET_TIER1_INDEX(pid) ((pid) >> 10)
#define GET_TIER2_INDEX(pid) ((pid)&0x3ff)

#define NULL_PID -1

namespace std {
struct HashPair {
  size_t operator()(peloton::ItemPointer const& ptr) const {
    using std::hash;
    return hash<peloton::oid_t>()(ptr.block) ^
           hash<peloton::oid_t>()(ptr.offset);
  }
};
}

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

class ItemPointerEqualityChecker {
 public:
  inline bool operator()(const ItemPointer& lhs, const ItemPointer& rhs) const {
    return (lhs.block == rhs.block) && (lhs.offset == rhs.offset);
  }
};

class ItemPointerComparator {
 public:
  inline bool operator()(const ItemPointer& lhs, const ItemPointer& rhs) const {
    return (lhs.block < rhs.block) ||
           ((lhs.block == rhs.block) && (lhs.offset < rhs.offset));
  }
};

// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
template <typename KeyType, typename ValueType, class KeyComparator,
          typename KeyEqualityChecker>
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
    // The first level mapping table
    Node** mappingtable_1[MAPPING_TABLE_SIZE];
    std::atomic<unsigned long> nextPid;

   public:
    MappingTable() {
      for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
        mappingtable_1[i] = nullptr;
      }
      nextPid = 0;
    }

    bool delete_chain(Node* node) {
      Node* next = node;
      while (next) {
        node = next;
        if (node->node_type == MERGE_DELTA) {
          PidType temp = ((MergeDelta*)node)->orignal_node->pid;
          delete_chain(this->get(temp));
        }
        next = node->next;
        delete node;
      }
      return true;
    }

    ~MappingTable() {
      for (int i = 0; i < MAPPING_TABLE_SIZE; i++) {
        if (mappingtable_1[i] != nullptr) {
          for (int j = 0; j < MAPPING_TABLE_SIZE; j++) {
            if (mappingtable_1[i][j] != nullptr)
              delete_chain(mappingtable_1[i][j]);
          }
          delete[] mappingtable_1[i];
          mappingtable_1[i] = nullptr;
        }
      }
    }

    Node* get(PidType pid) {
      if (pid == NULL_PID) return nullptr;

      long tier1_idx = GET_TIER1_INDEX(pid);
      long tier2_idx = GET_TIER2_INDEX(pid);

      if (mappingtable_1[tier1_idx] == nullptr) return nullptr;

      return mappingtable_1[tier1_idx][tier2_idx];
    }

    bool set(PidType pid, Node* expected, Node* addr) {
      if (addr != nullptr) addr->pid = pid;

      long tier1_idx = GET_TIER1_INDEX(pid);
      long tier2_idx = GET_TIER2_INDEX(pid);

      if (mappingtable_1[tier1_idx] == nullptr) {
        LOG_ERROR("meet a null level1 idx");
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

      if (addr != nullptr) addr->pid = new_pid;

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
        return NULL_PID;
      }
    }

    void remove(PidType pid) {
      long tier1_idx = GET_TIER1_INDEX(pid);
      long tier2_idx = GET_TIER2_INDEX(pid);

      mappingtable_1[tier1_idx][tier2_idx] = nullptr;
    }
  };

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

    // check if lowkey == -inf, highkey == +inf
    bool inf_lowkey, inf_highkey;

    // Double linked list pointers to traverse the leaves
    // Node* prev_node;

    // Double linked list pointers to traverse the leaves
    PidType next_leafnode;

    // constructor
    Node(Node* ne, NodeType ntype, size_t delta_l, MappingTable& mt,
         PidType next_leaf, bool is_lf, KeyType lowkey, KeyType highkey,
         bool inf_low, bool inf_high)
        : mapping_table(mt),
          low_key(lowkey),
          high_key(highkey),
          inf_lowkey(inf_low),
          inf_highkey(inf_high) {
      node_type = ntype;
      slotuse = 0;
      next = ne;
      delta_list_len = delta_l;
      next_leafnode = next_leaf;
      is_leaf = is_lf;
    }

    virtual ~Node() {}

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
  // We need a root node
  PidType root;

  // the mapping table used in our BWTree
  MappingTable mapping_table;

  // another mapping table used for garbage collection
  MappingTable garbage_table;

  /// Pointer to first leaf in the double linked leaf chain
  PidType headleaf;

  /// Pointer to last leaf in the double linked leaf chain
  PidType tailleaf;

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
    InnerNode(MappingTable& mapping_table, KeyType lowkey, KeyType highkey,
              bool inf_low, bool inf_high)
        : Node(nullptr, NodeType::INNER, 0, mapping_table, NULL_PID, false,
               lowkey, highkey, inf_low, inf_high) {}

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

    LeafNode(MappingTable& mapping_table, PidType next_leafnode, KeyType lowkey,
             KeyType highkey, bool inf_low, bool inf_high)
        : Node(nullptr, NodeType::LEAF, 0, mapping_table, next_leafnode, true,
               lowkey, highkey, inf_low, inf_high) {
      // initialize the value bucket
      for (int i = 0; i < leafslotmax + 1; i++) {
        slotdata[i] = nullptr;
      }
    }

    ~LeafNode() {
      for (int i = 0; i < leafslotmax + 1; i++) {
        if (slotdata[i] != nullptr) delete slotdata[i];
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
    //    inline void set_slot(unsigned short slot, const PairType& value) {
    //      assert(slot < Node::slotuse);
    //      slotkey[slot] = value.first;
    //      if (slotdata[slot] == nullptr) {
    //        slotdata[slot] = new std::vector<ValueType>();
    //      }
    //      slotdata[slot]->push_back(value.second);
    //    }

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

    RecordDelta(Node* next, RecordType op, KeyType k, ValueType v,
                MappingTable& mapping_table, PidType next_leafnode,
                KeyType lowkey, KeyType highkey, bool inf_low, bool inf_high)
        : Node(nullptr, NodeType::RECORD_DELTA, 0, mapping_table, next_leafnode,
               next->is_leaf, lowkey, highkey, inf_low, inf_high),
          op_type(op),
          key(k),
          value(v) {
      prepend(this, next);
      // temporarily update the slotuse of the new delta node
      this->slotuse = next->slotuse;
    }
    RecordType op_type;
    KeyType key;
    ValueType value;
  };

  // Delta Node for spliting operation
  struct SplitDelta : public Node {
    SplitDelta(Node* next, KeyType Kp, PidType pQ, MappingTable& mapping_table,
               PidType next_leafnode, KeyType lowkey, KeyType highkey,
               bool inf_low, bool inf_high)
        : Node(next, NodeType::SPLIT_DELTA, 0, mapping_table, next_leafnode,
               next->is_leaf, lowkey, highkey, inf_low, inf_high),
          Kp(Kp),
          pQ(pQ) {
      prepend(this, next);
      // temporarily update the slotuse of the new delta node
      this->slotuse = next->slotuse / 2;
    }
    KeyType Kp;
    PidType pQ;
  };

  struct IndexEntryDelta : public Node {
    IndexEntryDelta(Node* next, KeyType Kp, KeyType Kq, bool Kq_is_inf,
                    PidType pQ, MappingTable& mapping_table,
                    PidType next_leafnode, KeyType lowkey, KeyType highkey,
                    bool inf_low, bool inf_high)
        : Node(next, NodeType::INDEX_ENTRY_DELTA, 0, mapping_table,
               next_leafnode, next->is_leaf, lowkey, highkey, inf_low,
               inf_high),
          Kp(Kp),
          Kq(Kq),
          pQ(pQ),
          inf_Kq(Kq_is_inf) {
      prepend(this, next);
      // update the slotuse of the new delta node
      this->slotuse = next->slotuse + 1;
    }
    KeyType Kp, Kq;
    PidType pQ;
    bool inf_Kq;
  };

  // Delta Node for merging operation
  struct RemoveDelta : public Node {
    RemoveDelta(Node* next, MappingTable& mapping_table, PidType next_leafnode,
                KeyType lowkey, KeyType highkey, bool inf_low, bool inf_high)
        : Node(next, NodeType::REMOVE_NODE_DELTA, 0, mapping_table,
               next_leafnode, next->is_leaf, lowkey, highkey, inf_low,
               inf_high) {}
  };

  struct MergeDelta : public Node {
    MergeDelta(Node* next, KeyType Kp, Node* orignal_node,
               MappingTable& mapping_table, PidType next_leafnode,
               KeyType lowkey, KeyType highkey, bool inf_low, bool inf_high)
        : Node(next, NodeType::MERGE_DELTA, 0, mapping_table, next_leafnode,
               next->is_leaf, lowkey, highkey, inf_low, inf_high),
          Kp(Kp),
          orignal_node(orignal_node) {}
    KeyType Kp;
    Node* orignal_node;
  };

  struct DeleteIndexDelta : public Node {
    DeleteIndexDelta(Node* next, KeyType Kp, KeyType Kq, bool Kq_is_inf,
                     PidType pQ, MappingTable& mapping_table,
                     PidType next_leafnode, KeyType lowkey, KeyType highkey,
                     bool inf_low, bool inf_high)
        : Node(next, NodeType::INDEX_ENTRY_DELTA, 0, mapping_table,
               next_leafnode, next->is_leaf, lowkey, highkey, inf_low,
               inf_high),
          Kp(Kp),
          Kq(Kq),
          pQ(pQ),
          inf_Kq(Kq_is_inf) {}
    KeyType Kp, Kq;
    PidType pQ;
    bool inf_Kq;
  };

 public:
  // constructor

  BWTree(const KeyComparator& kc, const KeyEqualityChecker& ke,
         peloton::index::IndexMetadata* metadata)
      : m_key_less(kc),
        m_key_equal(ke),
        m_value_equal(ItemPointerEqualityChecker()),
        m_metadata(metadata) {
    KeyType waste;
    LeafNode* addr =
        new LeafNode(mapping_table, NULL_PID, waste, waste, true, true);
    long newpid = mapping_table.add(addr);
    if (newpid >= 0) {
      // initialize the root, head and tail pid.
      root = newpid;
      headleaf = tailleaf = newpid;
    } else {
      delete addr;
      LOG_ERROR("Can't create the initial leafNode!");
    }

    LOG_INFO("leaf_max = %d, inner_max = %d", leafslotmax, innerslotmax);
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
  inline bool key_less(const KeyType& a, const KeyType b,
                       bool b_max_inf) const {
    if (b_max_inf) return true;
    return m_key_less(a, b);
  }

  // True if a <= b ? constructed from key_less()
  inline bool key_lessequal(const KeyType& a, const KeyType b,
                            bool b_max_inf) const {
    if (b_max_inf) return true;
    return !m_key_less(b, a);
  }

  // True if a > b ? constructed from key_less()
  inline bool key_greater(const KeyType& a, const KeyType& b,
                          bool b_min_inf) const {
    if (b_min_inf) return true;
    return m_key_less(b, a);
  }

  // True if a >= b ? constructed from key_less()
  inline bool key_greaterequal(const KeyType& a, const KeyType b,
                               bool b_min_inf) const {
    if (b_min_inf) return true;
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
  ItemPointerEqualityChecker m_value_equal;
  peloton::index::IndexMetadata* m_metadata;

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

    if (path.empty()) {
      LOG_ERROR("Search get empty tree");
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
      case RECORD_DELTA: {
        if (node->node_type == LEAF) {
          LOG_INFO("Search Range Info: meet a leaf, pid: %lld, slotuse %d",
                   node->pid, node->slotuse);
        } else {
          LOG_INFO("Search Range Info: meet a record, pid: %lld, slotuse %d",
                   node->pid, node->slotuse);
        }

        return node->pid;
      }
      case INDEX_ENTRY_DELTA:
      case DELETE_INDEX_TERM_DELTA: {
        IndexEntryDelta* indexEntryDelta = static_cast<IndexEntryDelta*>(node);
        if (key_greaterequal(key, indexEntryDelta->Kp, false) &&
            key_less(key, indexEntryDelta->Kq, indexEntryDelta->inf_Kq)) {
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
      }
      case REMOVE_NODE_DELTA: {
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
      }
      case MERGE_DELTA: {
        LOG_INFO("Search Range Info: meet merge delta");

        if (key_greaterequal(key, ((MergeDelta*)node)->Kp, false)) {
          node = ((MergeDelta*)node)->orignal_node;
          if (node == nullptr) {
            LOG_ERROR("pid in split delta not exist");
            return -1;
          }

          return search(node, key, path);
        }
        return search(node->next, key, path);
      }
      case SPLIT_DELTA: {
        LOG_INFO("Search Range Info: meet split/merge delta");
        pid = ((SplitDelta*)node)->pQ;
        if (key_greaterequal(key, ((SplitDelta*)node)->Kp, false)) {
          // TODO: should
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
      }
      case INNER: {
        if (node->slotuse == 0) {
          LOG_INFO("empty inner node");
          pid = ((InnerNode*)node)->childid[0];
          if (pid == NULL_PID) {
            LOG_ERROR("error leftmost pid -- NULL_PID");
            return -1;
          }

          path.push(pid);
          node = mapping_table.get(pid);
          if (node == nullptr) {
            LOG_ERROR("pid in inner node not exist");
            return -1;
          }

          return search(node, key, path);
        } else {
          int i = 0;
          for (i = 0; i < node->slotuse; i++) {
            if (key_less(key, ((InnerNode*)node)->slotkey[i], false))
              break;  // 0, 1 ,2
          }           // 0  1  2  3
          pid = ((InnerNode*)node)->childid[i];
          node = mapping_table.get(pid);
          if (node == nullptr) {
            LOG_ERROR("pid in inner node not exist");
            return -1;
          }
          path.push(pid);
          return search(node, key, path);
        }
      }
      default:
        return -1;
    }
    return -1;
  }

  typedef std::unordered_set<ValueType, std::HashPair,
                             ItemPointerEqualityChecker> DelSet;

  bool key_is_in(KeyType key, Node* listhead, DelSet& deleted) {
    if (listhead == nullptr) return false;

    Node* node = listhead;
    switch (node->node_type) {
      case RECORD_DELTA: {
        RecordDelta* rcd_node = (RecordDelta*)node;
        if (rcd_node->op_type == RecordDelta::INSERT &&
            key_equal(rcd_node->key, key) &&
            (!deleted.count(rcd_node->value))) {
          return true;
        } else if (rcd_node->op_type == RecordDelta::DELETE &&
                   key_equal(rcd_node->key, key)) {
          deleted.insert(rcd_node->value);
        }
        return key_is_in(key, node->next, deleted);
      }
      case LEAF: {
        LeafNode* lf_node = (LeafNode*)node;
        for (int i = 0; i < (lf_node->slotuse); i++) {
          if (key_equal(lf_node->slotkey[i], key)) {
            for (auto val : (*lf_node->slotdata[i])) {
              if (!deleted.count(val)) {
                return true;
              }
            }
            return false;
          }
        }
        return false;
      }
      case MERGE_DELTA: {
        if (key_greaterequal(key, ((MergeDelta*)node)->Kp, false)) {
          node = ((MergeDelta*)node)->orignal_node;
          return key_is_in(key, node, deleted);
        }
        return key_is_in(key, node->next, deleted);
      }
      case SPLIT_DELTA: {
        if (key_greaterequal(key, ((SplitDelta*)node)->Kp, false)) {
          PidType pid = ((SplitDelta*)node)->pQ;
          return key_is_in(key, mapping_table.get(pid), deleted);
        }
        return key_is_in(key, node->next, deleted);
      }
      default:
        return false;
    }
  }

  inline bool key_is_in(KeyType key, Node* listhead) {
    DelSet deleted_set;
    return key_is_in(key, listhead, deleted_set);
  };

  // return pair nums, also calculates total val nums of the key in count
  std::pair<int, int> count_pair(KeyType key, ValueType value, Node* listhead) {
    int total_count = 0;
    int pair_count = 0;
    DelSet deleted;
    Node* node = listhead;

    int len = 0;
    while (node != nullptr) {
      switch (node->node_type) {
        case RECORD_DELTA: {
          // TODO : this part is buggy!
          len++;
          RecordDelta* rcd_node = static_cast<RecordDelta*>(node);
          if (rcd_node->op_type == RecordDelta::INSERT) {
            if (key_equal(rcd_node->key, key) &&
                (!deleted.count(rcd_node->value))) {
              total_count++;
              if (value_equal(rcd_node->value, value)) {
                pair_count++;
              }
              //              printf("count_pair after inserting!\n");
            }
          } else if (rcd_node->op_type == RecordDelta::DELETE &&
                     key_equal(rcd_node->key, key)) {
            //            printf("count_pair before delete!\n");
            deleted.insert(rcd_node->value);
            //            printf("count_pair after delete!\n");
          }
          node = node->next;
          // TODO : this part is buggy!
          break;
        }
        case LEAF: {
          LeafNode* lf_node = static_cast<LeafNode*>(node);
          for (int i = 0; i < (lf_node->slotuse); i++) {
            if (key_equal(lf_node->slotkey[i], key)) {
              for (ValueType v : *(lf_node->slotdata[i])) {
                if (deleted.count(v)) continue;
                total_count++;
                if (value_equal(v, value)) pair_count++;
              }
              break;
            }
          }
          assert(node->next == nullptr);
          node = nullptr;
          break;
        }
        case MERGE_DELTA: {
          if (key_greaterequal(key, ((MergeDelta*)node)->Kp, false)) {
            node = ((MergeDelta*)node)->orignal_node;
          } else {
            node = node->next;
          }
          break;
        }
        case SPLIT_DELTA: {
          if (key_greaterequal(key, ((SplitDelta*)node)->Kp, false)) {
            PidType pid = ((SplitDelta*)node)->pQ;
            node = mapping_table.get(pid);
          } else {
            node = node->next;
          }
          break;
        }
        default:
          LOG_ERROR("count pair should not be here");
          break;
      }
    }

    return std::pair<int, int>(total_count, pair_count);
  }

  bool append_delete(Node* basic_node, KeyType key, ValueType value,
                     bool deletekey) {
    RecordDelta* new_delta = new RecordDelta(
        basic_node, RecordDelta::DELETE, key, value, mapping_table,
        basic_node->next_leafnode, basic_node->low_key, basic_node->high_key,
        basic_node->inf_lowkey, basic_node->inf_highkey);

    if (deletekey) {
      new_delta->slotuse -= 1;
    }

    if (mapping_table.set(basic_node->pid, basic_node, new_delta)) {
      return true;
    } else {
      LOG_INFO("CAS FAIL: redo delete recordDelta");
      delete new_delta;
      return false;
    };
  }

  bool apend_merge() {
    // TODO: A lot
    return false;
  }

  // private fuctions, invisible to users
  inline static bool prepend(Node* delta_node, Node* orig_node) {
    // update the delta_list_len of the delta node
    delta_node->delta_list_len = orig_node->delta_list_len + 1;

    // maintain next, prev pointer
    delta_node->next = orig_node;

    return true;
  }

 public:
  void get_value(KeyType key, std::vector<ValueType>& result) {
    std::stack<PidType> path = search(root, key);

    PidType target_node = path.top();
    Node* next = mapping_table.get(target_node);

    LOG_INFO("Search result: pid - %lld, slotuse: %d", next->pid,
             next->slotuse);

    if (!next->is_leaf) {
      LOG_ERROR("get_value's search result is not a leaf");
      return;
    }

    DelSet delset;

    // traverse the delta chain from top down
    // to construct the correct result
    while (next) {
      switch (next->node_type) {
        case RECORD_DELTA: {
          RecordDelta* node = static_cast<RecordDelta*>(next);
          if (key_equal(node->key, key)) {
            if (node->op_type == RecordDelta::RecordType::INSERT) {
              if (delset.find(node->value) == delset.end())
                result.push_back(node->value);
            } else if (node->op_type == RecordDelta::RecordType::DELETE) {
              // if we meet a "delete" all later value associated with this key
              // is invalid.
              delset.insert(node->value);
            }
          }
          next = node->next;
        } break;
        case LEAF: {
          LeafNode* leaf = static_cast<LeafNode*>(next);
          for (int i = 0; i < leaf->slotuse; i++) {
            if (key_equal(leaf->slotkey[i], key)) {
              unsigned long vsize = leaf->slotdata[i]->size();
              for (int j = 0; j < vsize; j++) {
                if (delset.find(leaf->slotdata[i]->at(j)) == delset.end())
                  result.push_back(leaf->slotdata[i]->at(j));
              }
            }
          }
          next = nullptr;
        } break;
        case SPLIT_DELTA: {
          SplitDelta* split_delta = static_cast<SplitDelta*>(next);
          // if key >= Kp, we go to the new split node
          if (key_greaterequal(key, split_delta->Kp, false)) {
            LOG_ERROR("search result direct to another node");
            next = mapping_table.get(split_delta->pQ);
          }
          // else we go alone this delta chain
          else {
            next = split_delta->next;
          }
        } break;
        case MERGE_DELTA: {
          MergeDelta* merge_delta = static_cast<MergeDelta*>(next);
          // if key >= Kp, we go to the original node
          if (key_greaterequal(key, merge_delta->Kp, false)) {
            next = merge_delta->orignal_node;
          }
          // else we go alone this delta chain
          else {
            next = merge_delta->next;
          }
        } break;
        case REMOVE_NODE_DELTA:
          // if we meet remove node delta, we can search from the root again.
          result.clear();
          get_value(key, result);
          break;
        default:
          LOG_ERROR("meet wrong delta: %d during get_value", next->node_type);
      }
    }
  }

  // perform split for the node containing (key, value)
  void split(KeyType key) {
    std::stack<PidType> path = search(BWTree::root, key);

    PidType check_split_pid = path.top();
    path.pop();
    Node* check_split_node = mapping_table.get(check_split_pid);

    // Step1: Check if we need to split
    while (check_split_node->need_split()) {
      LOG_INFO("pid = %llu, begin Split", check_split_pid);
      // Step 1.1 add splitDelta to current check_split_node
      SplitDelta* new_split;
      KeyType pivotal;

      //      std::vector<catalog::Column> columns;
      //      catalog::Column column1(VALUE_TYPE_INTEGER,
      //      GetTypeSize(VALUE_TYPE_INTEGER),
      //                              "A", true);
      //      catalog::Column column2(VALUE_TYPE_VARCHAR, 1024, "B", true);
      //      columns.push_back(column1);
      //      columns.push_back(column2);
      //      // INDEX KEY SCHEMA -- {column1, column2}
      //      catalog::Schema *key_schema = new catalog::Schema(columns);
      //      key_schema->SetIndexedColumns({0, 1});
      //      storage::Tuple * tp = (new storage::Tuple(key_schema, true));
      //      tp->SetValue(0, ValueFactory::GetIntegerValue(2500), nullptr);
      //      tp->SetValue(1, ValueFactory::GetStringValue("e"), nullptr);
      //      KeyType key4;
      //      key4.SetFromKey(tp);
      //
      //      ItemPointer item1(121, 7);
      //
      //
      //      Node * basic_node = check_split_node;
      //      {
      //        printf("check if 2500, e - item1  in original node:\n");
      //        auto res = leaf_fake_consolidate(basic_node);
      //        auto keys = res.first;
      //        auto vals = res.second;
      //        int count = 0;
      //        int ck = 0;
      //        for (int i=0; i < keys.size(); i++) {
      //          if (key_equal(keys[i], key4)) {
      //            ck = vals[i].size();
      //            for (int j=0; j<ck; j++) {
      //              if (value_equal(vals[i][j], item1)) {
      //                count ++;
      //              }
      //            }
      //            break;
      //          }
      //        }
      //        printf("k-v pair found by my sanity: %d, %d\n", ck, count);
      //      }

      // create a new right sibling node
      PidType new_node_pid;

      if (check_split_node->is_leaf)
        new_node_pid = create_leaf(check_split_pid, &pivotal);
      else
        new_node_pid = create_inner(check_split_pid, &pivotal);

      Node* new_node = mapping_table.get(new_node_pid);

      //      basic_node = new_node;
      //      {
      //        std::cout << "pivaotal: " <<
      //            pivotal.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
      //            << "," <<
      //            pivotal.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
      //            << std::endl;
      //
      //        printf("check if 2500, e - item1 in new leaf node: \n");
      //        std::cout << "new leaf lowkey: " <<
      //            basic_node->low_key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
      //            << "," <<
      //            basic_node->low_key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
      //            << std::endl;
      //
      //        auto res = leaf_fake_consolidate(basic_node);
      //        auto keys = res.first;
      //
      //        std::cout << "new leaf key[0]: " <<
      //            keys[0].GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
      //            << "," <<
      //            keys[0].GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
      //            << std::endl;
      //
      //        auto vals = res.second;
      //        int count = 0;
      //        int ck = 0;
      //        for (int i=0; i < keys.size(); i++) {
      //          if (key_equal(keys[i], key4)) {
      //            ck = vals[i].size();
      //            for (int j=0; j<ck; j++) {
      //              if (value_equal(vals[i][j], item1)) {
      //                count ++;
      //              }
      //            }
      //            break;
      //          }
      //        }
      //        printf("k-v pair found by my sanity: %d, %d\n", ck, count);
      //      }

      // create and prepend a split node
      new_split = new SplitDelta(check_split_node, pivotal, new_node_pid,
                                 mapping_table, new_node_pid,
                                 check_split_node->low_key, new_node->low_key,
                                 check_split_node->inf_lowkey, false);

      if (!mapping_table.set(check_split_pid, check_split_node, new_split)) {
        LOG_INFO("CAS FAIL: Split CAS fails");
        // clean created waste
        Node* old_ptr = mapping_table.get(new_node_pid);
        delete old_ptr;
        mapping_table.set(new_node_pid, old_ptr, nullptr);
        delete new_split;

        check_split_node = mapping_table.get(check_split_pid);
        continue;
      }

      if (check_split_node->is_leaf) {
        LOG_INFO("pid = %llu, Split finished, new leaf node %llu created",
                 check_split_pid, new_node_pid);
#ifdef MY_PRINT_DEBUG
        print_node_info(check_split_pid);
        print_node_info(new_node_pid);
#endif
      } else {
        LOG_INFO("pid = %llu, Split finished, new inner node %llu created",
                 check_split_pid, new_node_pid);
#ifdef MY_PRINT_DEBUG
        print_node_info(check_split_pid);
        print_node_info(new_node_pid);
#endif
      }

      check_split_node = mapping_table.get(check_split_pid);
      // check if we need to consolidate
      if (check_split_node->delta_list_len > MERGE_DELTA) {
        consolidate(check_split_pid);
      }

      //      basic_node = check_split_node;
      //      {
      //        printf("check if 2500, e - item1  in new original node:\n");
      //        std::cout << "new original highkey: " <<
      //            basic_node->high_key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
      //            << "," <<
      //            basic_node->high_key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
      //            << std::endl;
      //
      //        auto res = leaf_fake_consolidate(basic_node);
      //        auto keys = res.first;
      //
      //        std::cout << "new original key[slotuse-1]: " <<
      //            keys[keys.size()-1].GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
      //            << "," <<
      //            keys[keys.size()-1].GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
      //            << std::endl;
      //
      //        auto vals = res.second;
      //        int count = 0;
      //        int ck = 0;
      //        for (int i=0; i < keys.size(); i++) {
      //          if (key_equal(keys[i], key4)) {
      //            ck = vals[i].size();
      //            for (int j=0; j<ck; j++) {
      //              if (value_equal(vals[i][j], item1)) {
      //                count ++;
      //              }
      //            }
      //            break;
      //          }
      //        }
      //        printf("k-v pair found by my sanity: %d, %d\n", ck, count);
      //      }
      //      delete key_schema;
      //      delete tp;

      // Step 1.2 update our check_split_node as its parent (or create new root)
      if (path.empty()) {
        // create new root
        KeyType waste;
        InnerNode* new_root =
            new InnerNode(mapping_table, waste, waste, true, true);
        new_root->childid[0] = check_split_pid;
        check_split_pid = mapping_table.add(new_root);

        root = check_split_pid;
        LOG_INFO("new root = %llu, created", root);
      } else {
        // get the father node
        check_split_pid = path.top();
        path.pop();
        LOG_INFO("parent = %llu, created", check_split_pid);
      }

      check_split_node = mapping_table.get(check_split_pid);
      // check if we need to consolidate
      if (check_split_node->delta_list_len > MERGE_DELTA) {
        consolidate(check_split_pid);
      }

      // Step 1.3 add indexEntryDelta to current check_split_node
      bool redo = true;
      while (redo) {
        IndexEntryDelta* new_indexEntryDelta = new IndexEntryDelta(
            check_split_node, new_split->Kp, new_node->high_key,
            new_node->inf_highkey, new_split->pQ, mapping_table,
            check_split_node->next_leafnode, check_split_node->low_key,
            check_split_node->high_key, check_split_node->inf_lowkey,
            check_split_node->inf_highkey);

        //        std::cout << "indexEntryDelta Kp: " <<
        //            new_indexEntryDelta->Kp.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
        //            << "," <<
        //            new_indexEntryDelta->Kp.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
        //            << " pQ = "<< new_indexEntryDelta->pQ << std::endl;

        if (mapping_table.set(check_split_pid, check_split_node,
                              new_indexEntryDelta)) {
          LOG_INFO("new indexEntryDelta added to pid = %llu", check_split_pid);
          redo = false;
        } else {
          LOG_INFO("CAS FAIL: redo add indexEntryDelta");
          delete new_indexEntryDelta;
          check_split_node = mapping_table.get(check_split_pid);
        }
      }
      check_split_node = mapping_table.get(check_split_pid);
      // check if we need to consolidate
      if (check_split_node->delta_list_len > MERGE_DELTA) {
        consolidate(check_split_pid);
      }

#ifdef MY_PRINT_DEBUG
      print_node_info(check_split_pid);
#endif

      // if now the path is empty and current node
      // doesn't need splitting we stop split
      if (path.empty() && !check_split_node->need_split()) {
        break;
      }
    }
  }

  // public method exposed to users -mavis
  bool insert_entry(KeyType key, ValueType value) {
    // Step1: perform split if necessary
    split(key);

    // search again to make sure we are inserting the right node
    std::stack<PidType> path = search(BWTree::root, key);

    PidType basic_pid = path.top();
    path.pop();
    Node* basic_node = mapping_table.get(basic_pid);

    // Step2: Check whether we need to consolidate
    if (basic_node->delta_list_len > MAX_DELTA_CHAIN_LEN) {
      consolidate(basic_pid);
      //      print_node_info(basic_pid);
    }

    // Step3: Add the insert record delta to current delta chain
    // update the basic_node before adding record delta
    basic_node = mapping_table.get(basic_pid);

    bool redo = true;
    while (redo) {
      // check whether we can insert duplicate key
      if (m_metadata->HasUniqueKeys()) {
        LOG_INFO("unique key required!");
        auto count_res = count_pair(key, value, basic_node);
        if (count_res.first > 0) return false;
      }

      RecordDelta* new_delta = new RecordDelta(
          basic_node, RecordDelta::INSERT, key, value, mapping_table,
          basic_node->next_leafnode, basic_node->low_key, basic_node->high_key,
          basic_node->inf_lowkey, basic_node->inf_highkey);

      new_delta->pid = basic_pid;
      if (!key_is_in(key, new_delta->next))
        new_delta->slotuse = new_delta->next->slotuse + 1;

      basic_node = mapping_table.get(basic_pid);
      new_delta->high_key = basic_node->high_key;
      new_delta->low_key = basic_node->low_key;

      // TODO: use CAS concatenate this new_delta to the delta chain
      redo = !mapping_table.set(basic_pid, basic_node, new_delta);
      if (redo) {
        LOG_INFO("CAS FAIL: redo add insert record");
        delete new_delta;
        path = search(BWTree::root, key);

        basic_pid = path.top();
        path.pop();
        basic_node = mapping_table.get(basic_pid);
      }
      LOG_INFO("success add a new insert record delta, current delta len = %lu",
               mapping_table.get(basic_pid)->delta_list_len);
    }
    //    print_node_info(basic_pid);

    //    basic_node = mapping_table.get(basic_pid);
    //    if (basic_node->slotuse >= 24 && basic_node->slotuse <= 26 ) {
    //      std::cout << "inserted: " <<
    //      key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
    //          << "," <<
    //          key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
    //          << " " << value.block << "," << value.offset << " " <<
    //          std::endl;
    //      auto res = leaf_fake_consolidate(basic_node);
    //      auto keys = res.first;
    //      auto vals = res.second;
    //      int count = 0;
    //      int ck = 0;
    //      for (int i=0; i < keys.size(); i++) {
    //        if (key_equal(keys[i], key)) {
    //          ck = vals[i].size();
    //          for (int j=0; j<vals.size(); j++) {
    //            if (value_equal(vals[i][j], value)) {
    //              count ++;
    //            }
    //          }
    //          break;
    //        }
    //      }
    //      printf("k-v pair found by my sanity: %d, %d\n", ck, count);
    //    }
    return true;
  }

  bool delete_entry(KeyType key, ValueType value) {
    // Step1: perform split if necessary
    split(key);

    // search again to make sure we are deleting from the correct node
    std::stack<PidType> path = search(BWTree::root, key);
    PidType basic_pid = path.top();
    path.pop();
    Node* basic_node = mapping_table.get(basic_pid);

    // Step2: Check whether we need to consolidate
    if (basic_node->delta_list_len > MAX_DELTA_CHAIN_LEN) {
      consolidate(basic_pid);
      //      print_node_info(basic_pid);
    }

    // Step3: Add the insert record delta to current delta chain
    // update the basic_node before adding record delta
    bool redo = true;
    // check and insert delete delta
    while (redo) {
      std::stack<PidType> path = search(root, key);

      PidType basic_pid = path.top();
      path.pop();

      Node* basic_node = mapping_table.get(basic_pid);
      //      printf("before count pair in pid = %llu\n", basic_pid);
      //      if (basic_node->slotuse == 59) {
      //          std::cout << "want to delete: " <<
      //          key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
      //          << "," <<
      //          key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
      //              << " " << value.block << "," << value.offset << " " <<
      //              std::endl;
      //          auto res = leaf_fake_consolidate(basic_node);
      //          auto keys = res.first;
      //          auto vals = res.second;
      //          int count = 0;
      //          int ck = 0;
      //          for (int i=0; i < keys.size(); i++) {
      //            if (key_equal(keys[i], key)) {
      //              ck = vals[i].size();
      //              for (int j=0; j<vals.size(); j++) {
      //                if (value_equal(vals[i][j], value)) {
      //                  count ++;
      //                }
      //              }
      //              break;
      //            }
      //          }
      //          printf("k-v pair found by my sanity: %d, %d\n", ck, count);
      //      }
      auto tv_count_pair = count_pair(key, value, basic_node);

      //      if (basic_node->slotuse == 59)
      //        printf("k-v pair found by count_pair: %d %d\n",
      //        tv_count_pair.first, tv_count_pair.second);

      if (!tv_count_pair.second) {
        LOG_INFO("DeleteEntry Not Exist");
        return false;
      }

      if (tv_count_pair.second > tv_count_pair.first) {
        LOG_ERROR("error!! count pair second > first");
      }

      bool deletekey = (tv_count_pair.second == tv_count_pair.first);

      //      printf("before append_delete\n");
      redo = !append_delete(basic_node, key, value, deletekey);
      //      printf("after append_delete\n");
    }
    // TODO:apend merge_delta
    // TODO:apend delete_index_term_delta
    //    print_node_info(basic_pid);
    return true;
  };

  bool consolidate(PidType pid) {
    Node* orinode = mapping_table.get(pid);

    LOG_INFO("begin consolidation!");

    if (orinode->is_leaf) {
      // We need to consolidate a leaf node
      LeafNode* new_leaf = new LeafNode(
          mapping_table, orinode->next_leafnode, orinode->low_key,
          orinode->high_key, orinode->inf_lowkey, orinode->inf_highkey);

      auto res = leaf_fake_consolidate(orinode);
      auto keys = res.first;
      auto vals = res.second;

      if (keys.size() != vals.size() || keys.size() > leafslotmax) {
        LOG_ERROR("wrong consolidated leaf key size!");
      }

      for (int i = 0; i < keys.size(); i++) {
        new_leaf->slotkey[i] = keys[i];
        new_leaf->slotdata[i] = new std::vector<ValueType>(vals[i]);
      }
      new_leaf->slotuse = keys.size();

      // CAS replace the item with new leaf and throw old delta chain away
      if (mapping_table.set(pid, orinode, new_leaf)) {
        while (garbage_table.add(orinode) == NULL_PID)
          ;
        LOG_INFO("leaf consolidation finished!");
        return true;
      } else {
        delete new_leaf;
        LOG_INFO(
            "CAS FAIL: unnecessary consolidate, remove just created "
            "consolidated node");
        return false;
      }
    } else {
      // We need to consolidate an inner node
      InnerNode* new_inner =
          new InnerNode(mapping_table, orinode->low_key, orinode->high_key,
                        orinode->inf_lowkey, orinode->inf_highkey);

      auto res = inner_fake_consolidate(orinode);
      auto keys = res.first;
      auto childs = res.second;

      if (keys.size() != childs.size() - 1 || keys.size() > innerslotmax) {
        LOG_ERROR("wrong consolidated inner key size!");
      }

      int i;
      for (i = 0; i < keys.size(); i++) {
        new_inner->slotkey[i] = keys[i];
        new_inner->childid[i] = childs[i];
      }
      new_inner->childid[i] = childs[i];
      new_inner->slotuse = keys.size();

      // CAS replace the item with new leaf and throw old delta chain away
      if (mapping_table.set(pid, orinode, new_inner)) {
        while (garbage_table.add(orinode) == NULL_PID)
          ;
        LOG_INFO("inner consolidation finished!");
        return true;
      } else {
        delete new_inner;
        LOG_INFO(
            "CAS FAIL: unnecessary consolidate, remove just created "
            "consolidated node");
        return false;
      }
    }

    return false;
  }

  std::pair<std::vector<KeyType>, std::vector<std::vector<ValueType>>>
  leaf_fake_consolidate(Node* new_delta) {
    std::stack<Node*> delta_chain;
    Node* tmp_cur_node = new_delta;
    while (tmp_cur_node) {
      delta_chain.push(tmp_cur_node);
      tmp_cur_node = tmp_cur_node->next;
    }

    // prepare two array to store what logical k-v pairs we have
    std::vector<KeyType> tmpkeys;

    std::vector<std::vector<ValueType>> tmpvals;

    // TODO: move the leaf-node branch inside switch and add inner-node branch.
    // the first node must be the original leaf node itself
    assert(delta_chain.top()->node_type == LEAF);

    LeafNode* orig_leaf_node = static_cast<LeafNode*>(delta_chain.top());
    delta_chain.pop();

    // copy the data in the base node
    for (int i = 0; i < orig_leaf_node->slotuse; i++) {
      tmpkeys.push_back(orig_leaf_node->slotkey[i]);
      tmpvals.push_back(std::vector<ValueType>(*(orig_leaf_node->slotdata[i])));
    }

    LOG_INFO("Consolidate: delta_chain.len = %lu", delta_chain.size());

    //    int c = 0;
    // traverse the delta chain
    while (!delta_chain.empty()) {
      // get top delta node
      Node* cur_delta = delta_chain.top();
      delta_chain.pop();
      switch (cur_delta->node_type) {
        case RECORD_DELTA: {
          // first see the key has already existed
          bool no_key = true;
          RecordDelta* recordDelta = static_cast<RecordDelta*>(cur_delta);

          if (recordDelta->op_type == RecordDelta::INSERT) {
            for (int x = 0; x < tmpkeys.size(); x++) {
              if (key_equal(tmpkeys[x], recordDelta->key)) {
                tmpvals[x].push_back(recordDelta->value);
                no_key = false;
                break;
              }
            }

            //            LOG_INFO("insert %d", ++c);
            //            std::cout <<
            //            recordDelta->key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(0)
            //            << " "
            //                <<
            //                recordDelta->key.GetTupleForComparison(m_metadata->GetKeySchema()).GetValue(1)
            //                << std::endl;

            // key not exists, need to insert somewhere
            if (no_key) {
              if (recordDelta->slotuse == 0) {
                tmpkeys.push_back(recordDelta->key);
                tmpvals.push_back(
                    std::vector<ValueType>(1, recordDelta->value));
              } else {
                int target_pos = 0;
                for (int x = ((int)tmpkeys.size()) - 1; x >= 0; x--) {
                  if (key_greaterequal(recordDelta->key, tmpkeys[x], false)) {
                    target_pos = x + 1;
                    break;
                  }
                }

                tmpkeys.insert(tmpkeys.begin() + target_pos, recordDelta->key);
                tmpvals.insert(tmpvals.begin() + target_pos,
                               std::vector<ValueType>(1, recordDelta->value));
              }
              assert(tmpvals.size() == recordDelta->slotuse);
            }  // end of RecordDelta::INSERT

          } else if (recordDelta->op_type == RecordDelta::DELETE) {
            for (int x = 0; x < tmpkeys.size(); x++) {
              if (key_equal(tmpkeys[x], recordDelta->key)) {
                // remove value in the vector
                for (int j = ((int)tmpvals[x].size() - 1); j >= 0; j--) {
                  if (m_value_equal(tmpvals[x][j], recordDelta->value)) {
                    tmpvals[x].erase(tmpvals[x].begin() + j);
                  }
                }

                // if vector is empty, needed to be removed
                if (tmpvals[x].size() == 0) {
                  tmpkeys.erase(tmpkeys.begin() + x);
                  tmpvals.erase(tmpvals.begin() + x);
                }
                break;
              }
            }
            assert(tmpvals.size() == recordDelta->slotuse);
          }
          break;
        }
        case SPLIT_DELTA: {
          SplitDelta* splitDelta = static_cast<SplitDelta*>(cur_delta);
          // truncate all the values whose key is >= Kp
          for (int i = 0; i < tmpkeys.size(); i++) {
            if (key_greaterequal(tmpkeys[i], splitDelta->Kp, false)) {
              tmpkeys.resize(i);
              tmpvals.resize(i);
              break;
            }
          }
        } break;
        case MERGE_DELTA:
          break;
        case REMOVE_NODE_DELTA:
          break;
        case LEAF:
          LOG_ERROR("Wrong leaf delta!");
        default:
          break;
      }
    }

    return std::make_pair(tmpkeys, tmpvals);
  }

  std::pair<std::vector<KeyType>, std::vector<PidType>> inner_fake_consolidate(
      Node* new_delta) {
    std::stack<Node*> delta_chain;
    Node* tmp_cur_node = new_delta;
    while (tmp_cur_node) {
      delta_chain.push(tmp_cur_node);
      tmp_cur_node = tmp_cur_node->next;
    }

    // prepare two array to store what logical k-v pairs we have
    std::vector<KeyType> tmpkeys;

    std::vector<PidType> tmpchilds;

    // the first node must be the original leaf node itself
    assert(delta_chain.top()->node_type == INNER);

    InnerNode* orig_inner_node = static_cast<InnerNode*>(delta_chain.top());
    delta_chain.pop();

    // copy the data in the base node
    for (int i = 0; i < orig_inner_node->slotuse; i++) {
      tmpkeys.push_back(orig_inner_node->slotkey[i]);
      tmpchilds.push_back(orig_inner_node->childid[i]);
    }
    tmpchilds.push_back(orig_inner_node->childid[orig_inner_node->slotuse]);

    LOG_INFO("Consolidate: delta_chain.len = %lu", delta_chain.size());

    // traverse the delta chain
    while (!delta_chain.empty()) {
      // get top delta node
      Node* cur_delta = delta_chain.top();
      delta_chain.pop();

      switch (cur_delta->node_type) {
        case INDEX_ENTRY_DELTA: {
          IndexEntryDelta* indexEntryDelta =
              static_cast<IndexEntryDelta*>(cur_delta);
          int pos = 0;
          for (pos = 0; pos < tmpkeys.size(); pos++) {
            if (key_less(indexEntryDelta->Kp, tmpkeys[pos], false)) {
              break;
            }
          }

          tmpkeys.insert(tmpkeys.begin() + pos, indexEntryDelta->Kp);
          tmpchilds.insert(tmpchilds.begin() + pos + 1, indexEntryDelta->pQ);
        } break;

        case SPLIT_DELTA: {
          SplitDelta* splitDelta = static_cast<SplitDelta*>(cur_delta);
          // truncate all the values whose key is >= Kp
          for (int i = 0; i < tmpkeys.size(); i++) {
            if (key_greaterequal(tmpkeys[i], splitDelta->Kp, false)) {
              tmpkeys.resize(i);
              tmpchilds.resize(i + 1);
              break;
            }
          }
        } break;

        default:
          LOG_ERROR("impossible type on inner node");
          break;
      }
    }

    return std::make_pair(tmpkeys, tmpchilds);
  }

  PidType create_leaf(PidType check_split_pid, KeyType* pivotal) {
    Node* check_split_node = mapping_table.get(check_split_pid);
    KeyType waste;
    LeafNode* new_leaf = new LeafNode(
        mapping_table, check_split_node->next_leafnode, waste,
        check_split_node->high_key, false, check_split_node->inf_highkey);

    auto res = leaf_fake_consolidate(check_split_node);
    int orisize = check_split_node->slotuse;

    for (int i = orisize / 2; i < orisize; i++) {
      new_leaf->slotkey[i - orisize / 2] = res.first[i];
      new_leaf->slotdata[i - orisize / 2] =
          new std::vector<ValueType>(res.second[i]);
    }
    new_leaf->low_key = new_leaf->slotkey[0];
    new_leaf->slotuse = (unsigned short)((orisize + 1) / 2);

    *pivotal = new_leaf->slotkey[0];
    PidType new_leaf_pid = mapping_table.add(new_leaf);

    return new_leaf_pid;
  }

  PidType create_inner(PidType check_split_pid, KeyType* pivotal) {
    Node* check_split_node = mapping_table.get(check_split_pid);
    KeyType waste;
    InnerNode* new_inner =
        new InnerNode(mapping_table, waste, check_split_node->high_key, false,
                      check_split_node->inf_highkey);

    auto res = inner_fake_consolidate(check_split_node);
    int orisize = check_split_node->slotuse;

    new_inner->childid[0] = NULL_PID;
    for (int i = orisize / 2; i < orisize; i++) {
      new_inner->slotkey[i - orisize / 2] = res.first[i];
      new_inner->childid[i - orisize / 2 + 1] = res.second[i + 1];
    }
    new_inner->low_key = new_inner->slotkey[0];
    new_inner->slotuse = (unsigned short)((orisize + 1) / 2);

    *pivotal = new_inner->slotkey[0];
    PidType new_inner_pid = mapping_table.add(new_inner);

    return new_inner_pid;
  }

  void scan_all(std::vector<ValueType>& v) {
    Node* node = mapping_table.get(headleaf);

    // scan the leaf nodes list from begin to the end
    while (node != nullptr) {
      auto all_key_value_pair = leaf_fake_consolidate(node);

      LOG_INFO("fake_consolidate size: %lu", all_key_value_pair.second.size());

      for (auto const& value : all_key_value_pair.second) {
        v.insert(v.end(), value.begin(), value.end());
      }

      node = mapping_table.get(node->next_leafnode);
    }
  }

  void scan(std::vector<KeyType>& keys_result,
            std::vector<std::vector<ItemPointer>>& values_result) {
    Node* node = mapping_table.get(headleaf);

    // scan the leaf nodes list from begin to the end
    while (node != nullptr) {
      auto all_key_value_pair = leaf_fake_consolidate(node);
      std::vector<KeyType>& keys = all_key_value_pair.first;
      std::vector<std::vector<ValueType>>& values = all_key_value_pair.second;

      keys_result.insert(keys_result.end(), keys.begin(), keys.end());
      values_result.insert(values_result.end(), values.begin(), values.end());

      node = mapping_table.get(node->next_leafnode);
    }
  }

  void print_node_info(PidType pid) {
    Node* node = mapping_table.get(pid);
    size_t total_len = node->delta_list_len;
    printf("pid - %lld, delta_chain_len: %ld, slotuse: %d  ", pid, total_len,
           node->slotuse);

    // print out deltas added to this node in order
    for (int i = 0; i <= total_len; i++) {
      if (node->delta_list_len != total_len - i) {
        LOG_ERROR("Wrong delta chain length!");
      }
      if (node->node_type == RECORD_DELTA) {
        if (((RecordDelta*)node)->op_type == RecordDelta::INSERT)
          printf("insert->");
        else if (((RecordDelta*)node)->op_type == RecordDelta::DELETE)
          printf("delete->");
      } else if (node->node_type == SPLIT_DELTA) {
        printf("split(pQ=%llu)->", ((SplitDelta*)node)->pQ);
      } else if (node->node_type == LEAF) {
        printf("leaf");
      } else if (node->node_type == INNER) {
        printf("inner");
      }
      node = node->next;
    }

    printf("\n");
  }
};

}  // End index namespace
}  // End peloton namespace
