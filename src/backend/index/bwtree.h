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


#include <bits/stl_pair.h>
#include <bits/unordered_map.h>
#include <stddef.h>
#include <assert.h>
#include <unordered_map>
#include "../common/types.h"
#include "../storage/tuple.h"

#define BWTREE_NODE_SIZE  64

#define BWTREE_MAX(a,b)       ((a) < (b) ? (b) : (a))

#define MAX_DELTA_CHAIN_LEN   8


namespace peloton {
namespace index {

// Look up the stx btree interface for background.
// peloton/third_party/stx/btree.h
template <typename KeyType, typename ValueType, class KeyComparator, typename KeyEqualityChecker>
class BwTree {


 public:
  // *** Constructed Types

  /// Typedef of our own type
  typedef BwTree<KeyType, ValueType, KeyComparator> BwTreeSelf;

  /// Size type used to count keys
  typedef size_t                              SizeType;

  /// The pair of key_type and data_type, this may be different from value_type.
  typedef std::pair<KeyType, ValueType>      PairType;


 public:
  // *** Static Constant Options and Values of the Bw Tree
  /// Base B2 tree parameter: The number of key/data slots in each leaf
  static const unsigned short  leafslotmax = BWTREE_MAX( 8, BWTREE_NODE_SIZE / (sizeof(KeyType) + sizeof(ValueType)) );

  /// Base B+ tree parameter: The number of key slots in each inner node,
  /// this can differ from slots in each leaf.
  static const unsigned short  innerslotmax = BWTREE_MAX( 8, BWTREE_NODE_SIZE / (sizeof(KeyType) + sizeof(void*)) );

  /// Computed B+ tree parameter: The minimum number of key/data slots used
  /// in a leaf. If fewer slots are used, the leaf will be merged or slots
  /// shifted from it's siblings.
  static const unsigned short minleafslots = (leafslotmax / 2);

  /// Computed B+ tree parameter: The minimum number of key slots used
  /// in an inner node. If fewer slots are used, the inner node will be
  /// merged or slots shifted from it's siblings.
  static const unsigned short mininnerslots = (innerslotmax / 2);

  // We need a root node
  PidType  root;

 private:
  struct Node;
typedef long long PidType;

  enum NodeType
  {
      LEAF = 0,
      INNER = 1,
      RECORD_DELTA = 2,
      INDEX_ENTRY_DELTA = 3,
      REMOVE_NODE_DELTA = 4,
      MERGE_DELTA = 5,
      DELETE_INDEX_TERM_DELTA = 6
  };

  std::unordered_map<PidType, Node*> mapping_table;


  /**
   * The Node inheritance hierachy
   * **/
  struct Node
  {
    /// Number of key slotuse use, so number of valid children or data
    /// pointers
    unsigned short  slotuse;

    // Delta chain next pointer
    Node *next;
    // Do we need a prev pointer? -mavis
    Node *prev;

    // Length of current delta chain
    size_t delta_list_len;

    //type of this node
    NodeType node_type;

    // minimal and maximal key in this node
    KeyType low_key, high_key;

    /// Delayed initialisation of constructed node
    inline void initialize(NodeType n)
    {
      node_type = n;
      slotuse = 0;
    }

    /// True if this is a leaf node
    inline bool isleafnode() const
    {
      return (node_type == NodeType:: LEAF);
    }
  };

  /// Extended structure of a inner node in-memory. Contains only keys and no
  /// data items.
  struct InnerNode : public Node
  {
    /// Define an related allocator for the inner_node structs.
    // typedef typename _Alloc::template rebind<inner_node>::other alloc_type

    /// Keys of children or data pointers,
    //  we plus one so as to avoid overflow when consolidation
    KeyType        slotkey[innerslotmax + 1];

    /// Pointers to children,
    //  we plus one so as to avoid overflow when consolidation
    PidType        childid[innerslotmax+1 + 1];

    /// Set variables to initial values
    inline void initialize()
    {
      Node::initialize(NodeType:: INNER);
    }

    /// True if the node's slots are full
    inline bool isfull() const
    {
      return (Node::slotuse == innerslotmax);
    }

    /// True if few used entries, less than half full
    inline bool isfew() const
    {
      return (Node::slotuse <= mininnerslots);
    }

    /// True if node has too few entries
    inline bool isunderflow() const
    {
      return (Node::slotuse < mininnerslots);
    }
  };

  /// Extended structure of a leaf node in memory. Contains pairs of keys and
  /// data items. Key and data slots are kept in separate arrays, because the
  /// key array is traversed very often compared to accessing the data items.
  struct LeafNode : public Node
  {
    /// Define an related allocator for the leaf_node structs.
//    typedef typename _Alloc::template rebind<LeafNode>::other alloc_type;

    /// Double linked list pointers to traverse the leaves
    LeafNode       *prevleaf;

    /// Double linked list pointers to traverse the leaves
    LeafNode       *nextleaf;

    /// Keys of children or data pointers
    //  we plus one so as to avoid overflow when consolidation
    LeafNode        slotkey[leafslotmax + 1];

    /// Array of data
    //  we plus one so as to avoid overflow when consolidation
    ValueType       slotdata[leafslotmax + 1];

    /// Set variables to initial values
    inline void initialize()
    {
      Node::initialize(NodeType::LEAF);
      prevleaf = nextleaf = NULL;
    }

    /// True if the node's slots are full
    inline bool isfull() const
    {
      return (Node::slotuse == leafslotmax);
    }

    /// True if few used entries, less than half full
    inline bool isfew() const
    {
      return (Node::slotuse <= minleafslots);
    }

    /// True if node has too few entries
    inline bool isunderflow() const
    {
      return (Node::slotuse < minleafslots);
    }

    /// Set the (key,data) pair in slot. Overloaded function used by
    /// bulk_load().
    inline void set_slot(unsigned short slot, const PairType& value)
    {
      assert(slot < Node::slotuse);
      slotkey[slot] = value.first;
      slotdata[slot] = value.second;
    }

    /// Set the key pair in slot. Overloaded function used by
    /// bulk_load().
    inline void set_slot(unsigned short slot, const KeyType& key)
    {
      assert(slot < Node::slotuse);
      slotkey[slot] = key;
    }
  };

  // Delta Node for record update operation
  struct RecordDelta : public Node {

    //construction added -mavis
    RecordDelta(){ this->node_type = RECORD_DELTA; }

    enum RecordType
    {
      INSERT = 0,
      DELETE = 1,
      UPDATE = 2
    };

    RecordType op_type;
    KeyType key;
    ValueType value;
  };

  // Delta Node for spliting operation
  struct SplitDelta : public Node {
    KeyType Kp;
    PidType pQ;
  };

  struct IndexEntryDelta : public Node {
    KeyType Kp, Kq;
    PidType pQ;
  };

  // Delta Node for merging operation
  struct RemoveDelta : public Node {

  };

  struct MergeDelta : public Node {

  };

  struct DeleteIndexDelta : public Node {

  };

  //private functions, invisible to users -leiqi
  PidType Search<typename KeyType>(PidType rootpid, KeyType key);

  //public method exposed to users -mavis
  bool InsertEntry<typename KeyType, typename ValueType>( KeyType key, ValueType value );
  bool DeleteEntry<typename KeyType, typename ValueType>( KeyType key);
  bool UpdateEntry<typename KeyType, typename ValueType>( KeyType key, ValueType value );
  //interfaces of SCAN to be added -mavis

  //private fuctions, invisible to users -mavis
  bool Prepend(Node* delta_node, PidType orig_pid);
  //end -mavis

};

}  // End index namespace
}  // End peloton namespace
