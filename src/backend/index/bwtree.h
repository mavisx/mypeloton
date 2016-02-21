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

namespace peloton {
namespace index {

typedef long long PidType;

class MappingTable {
 private:
  BWTree::Node **mappingtable_1[MAPPING_TABLE_SIZE];
  static std::atomic<unsigned long> nextPid;

 public:
  MappingTable();
  ~MappingTable();
  BWTree::Node *get(PidType pid);
  bool set(PidType pid, void *addr);
  long add(void *addr);
  void remove(PidType pid);
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

 private:
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

  // We need a root node
  PidType root;

  static MappingTable mapping_table;

 public:
  /**
   * The Node inheritance hierachy
   * **/
  struct Node {
    /// Number of key slotuse use, so number of valid children or data
    /// pointers
    unsigned short slotuse;

    // flag to indicate whether this chain belongs to a leaf node
    bool is_leaf;

    // Delta chain next pointer
    Node *next;

    // Node pid
    PidType pid;

    // Length of current delta chain
    size_t delta_list_len;

    // type of this node
    NodeType node_type;

    // minimal and maximal key in this node
    KeyType low_key, high_key;

    // constructor
    Node(NodeType n) {
      node_type = n;
      slotuse = 0;
      next = nullptr;
      delta_list_len = -1;
    }

    // True if this is a leaf node
    inline bool isleafnode() const { return (node_type == NodeType::LEAF); }
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
    InnerNode() : Node(NodeType::INNER) {}

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

    /// Double linked list pointers to traverse the leaves
    LeafNode *prevleaf;

    /// Double linked list pointers to traverse the leaves
    LeafNode *nextleaf;

    /// Keys of children or data pointers
    //  we plus one so as to avoid overflow when consolidation
    KeyType slotkey[leafslotmax + 1];

    /// Array of data
    //  we plus one so as to avoid overflow when consolidation
    ValueType slotdata[leafslotmax + 1];

    LeafNode() : Node(NodeType::LEAF), prevleaf(nullptr), nextleaf(nullptr) {}

    /// True if the node's slots are full
    inline bool isfull() const { return (Node::slotuse == leafslotmax); }

    /// True if few used entries, less than half full
    inline bool isfew() const { return (Node::slotuse <= minleafslots); }

    /// True if node has too few entries
    inline bool isunderflow() const { return (Node::slotuse < minleafslots); }

    /// Set the (key,data) pair in slot. Overloaded function used by
    /// bulk_load().
    inline void set_slot(unsigned short slot, const PairType &value) {
      assert(slot < Node::slotuse);
      slotkey[slot] = value.first;
      slotdata[slot] = value.second;
    }

    /// Set the key pair in slot. Overloaded function used by
    /// bulk_load().
    inline void set_slot(unsigned short slot, const KeyType &key) {
      assert(slot < Node::slotuse);
      slotkey[slot] = key;
    }
  };

  // Delta Node for record update operation
  struct RecordDelta : public Node {
    // construction added -mavis
    RecordDelta(PidType next, RecordType op, KeyType k, ValueType v)
        : Node(NodeType::RECORD_DELTA) {
      this->op_type = op;
      this->key = k;
      this->value = v;
      prepend(this, next);
    }

    enum RecordType { INSERT = 0, DELETE = 1, UPDATE = 2 };

    RecordType op_type;
    KeyType key;
    ValueType value;
  };

  // Delta Node for spliting operation
  struct SplitDelta : public Node {
    SplitDelta(Node *next, KeyType Kp, PidType pQ)
        : Node(NodeType::SPLIT_DELTA), Kp(Kp), pQ(pQ) {}
    KeyType Kp;
    PidType pQ;
  };

  struct IndexEntryDelta : public Node {
    IndexEntryDelta(Node *next, KeyType Kp, KeyType Kq, PidType pQ)
        : Node(NodeType::INDEX_ENTRY_DELTA), Kp(Kp), Kq(Kq), pQ(pQ) {}
    KeyType Kp, Kq;
    PidType pQ;
  };

  // Delta Node for merging operation
  struct RemoveDelta : public Node {
    RemoveDelta(Node *next) : Node(NodeType::REMOVE_NODE_DELTA) {}
  };

  struct MergeDelta : public Node {
    MergeDelta(Node *next) : Node(NodeType::MERGE_DELTA), Kp(Kp), pQ(pQ) {}
    KeyType Kp;
    PidType pQ;
  };

  struct DeleteIndexDelta : public Node {
    DeleteIndexDelta(Node *next) : Node(NodeType::DELETE_INDEX_TERM_DELTA) {}
  };

 public:
  // constructor
  BWTree(){};

  // destructor
  ~BWTree(){};

  /*
   ************************************************
   * private functions, invisible to users -leiqi *
   ************************************************
   */

  KeyComparator m_key_less;

  /*
   ************************************************
   *    public method exposed to users -leiqi     *
   ************************************************
   */

  std::stack<PidType> search<typename KeyType>(PidType rootpid, KeyType key);
  PidType search<typename KeyType>(Node *node, KeyType key,
                                   std::stack<PidType> &path);

  bool is_in<typename KeyType>( KeyType key, Node** nptr);

  // True if a < b ? "constructed" from m_key_less()
  inline bool operator<(const KeyType &a, const KeyType b) const {
    return m_key_less(a, b);
  }

  // True if a <= b ? constructed from key_less()
  inline bool operator<=(const KeyType &a, const KeyType b) const {
    return !m_key_less(b, a);
  }

  // True if a > b ? constructed from key_less()
  inline bool operator>(const KeyType &a, const KeyType &b) const {
    return m_key_less(b, a);
  }

  // True if a >= b ? constructed from key_less()
  inline bool operator>=(const KeyType &a, const KeyType b) const {
    return !m_key_less(a, b);
  }

  // True if a == b ? constructed from key_less().
  inline bool operator==(const KeyType &a, const KeyType &b) const {
    return !m_key_less(a, b) && !m_key_less(b, a);
  }

  /*
   ************************************************
   *               end -leiqi                     *
   ************************************************
   */


  // public method exposed to users -mavis
  bool InsertEntry<typename KeyType, typename ValueType>(KeyType key,
                                                         ValueType value);
  bool DeleteEntry<typename KeyType, typename ValueType>(KeyType key);
  bool UpdateEntry<typename KeyType, typename ValueType>(KeyType key,
                                                         ValueType value);
  // interfaces of SCAN to be added -mavis

  // private fuctions, invisible to users -mavis
  static bool prepend(Node *delta_node, PidType orig_pid);
  // end -mavis
};

}  // End index namespace
}  // End peloton namespace
