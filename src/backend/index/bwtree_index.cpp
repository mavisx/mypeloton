//===----------------------------------------------------------------------===//
//
//                         PelotonDB
//
// btree_index.cpp
//
// Identification: src/backend/index/btree_index.cpp
//
// Copyright (c) 2015, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "backend/common/logger.h"
#include "backend/index/bwtree_index.h"
#include "backend/index/index_key.h"
#include "backend/storage/tuple.h"

namespace peloton {
namespace index {

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
BWTreeIndex<KeyType, ValueType, KeyComparator, KeyEqualityChecker,
            ValueComparator, ValueEqualityChecker>::BWTreeIndex(
    IndexMetadata *metadata)
    : Index(metadata),
      container(KeyComparator(metadata), KeyEqualityChecker(metadata), ValueEqualityChecker()),
      equals(metadata),
      comparator(metadata) {}

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
BWTreeIndex<KeyType, ValueType, KeyComparator,
            KeyEqualityChecker, ValueComparator, ValueEqualityChecker>::~BWTreeIndex() {
  // Add your implementation here
}

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
bool BWTreeIndex<KeyType, ValueType, KeyComparator,
                 KeyEqualityChecker,
                 ValueComparator, ValueEqualityChecker>::InsertEntry(__attribute__((unused))
                                                  const storage::Tuple *key,
                                                  __attribute__((unused))
                                                  const ItemPointer location) {
  KeyType index_key;
  index_key.SetFromKey(key);

  auto key_pair = std::pair<KeyType, ValueType>(index_key, location);

  return container.insert_entry(key_pair.first, key_pair.second);
}

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
bool BWTreeIndex<KeyType, ValueType, KeyComparator,
                 KeyEqualityChecker,
                 ValueComparator, ValueEqualityChecker>::DeleteEntry(__attribute__((unused))
                                                  const storage::Tuple *key,
                                                  __attribute__((unused))
                                                  const ItemPointer location) {
  KeyType index_key;
  index_key.SetFromKey(key);

  return container.delete_entry(index_key, location);
}

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
std::vector<ItemPointer>
BWTreeIndex<KeyType, ValueType, KeyComparator, KeyEqualityChecker,
            ValueComparator, ValueEqualityChecker>::Scan(
    __attribute__((unused)) const std::vector<Value> &values,
    __attribute__((unused)) const std::vector<oid_t> &key_column_ids,
    __attribute__((unused)) const std::vector<ExpressionType> &expr_types,
    __attribute__((unused)) const ScanDirectionType &scan_direction) {
  std::vector<KeyType> keys_result;
  std::vector<ItemPointer> values_result;
  std::vector<ItemPointer> result;

  {
    switch(scan_direction){
      case SCAN_DIRECTION_TYPE_FORWARD:
      case SCAN_DIRECTION_TYPE_BACKWARD:
        container.scan(keys_result, values_result);

        unsigned long vector_size = values_result.size();
        for (int i=0; i<vector_size; i++){
          auto tuple = keys_result[i].GetTupleForComparison(metadata->GetKeySchema());

          // Compare the current key in the scan with "values" based on "expression types"
          // For instance, "5" EXPR_GREATER_THAN "2" is true
          if (Compare(tuple, key_column_ids, expr_types, values) == true) {
            result.push_back(values_result[i]);
          }
        }

        break;

      case SCAN_DIRECTION_TYPE_INVALID:
      default:
        throw Exception("Invalid scan direction \n");
        break;
    }
  }

  return result;
}

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
std::vector<ItemPointer> BWTreeIndex<KeyType, ValueType, KeyComparator,
                                     KeyEqualityChecker,
                                     ValueComparator, ValueEqualityChecker>::ScanAllKeys() {
  std::vector<ItemPointer> result;

  container.scan_all(result);

  return result;
}

/**
 * @brief Return all locations related to this key.
 */
template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
std::vector<ItemPointer>
BWTreeIndex<KeyType, ValueType, KeyComparator, KeyEqualityChecker,
            ValueComparator, ValueEqualityChecker>::ScanKey(
    __attribute__((unused)) const storage::Tuple *key) {

  std::vector<ItemPointer> result;
  KeyType index_key;
  index_key.SetFromKey(key);

  container.get_value(index_key, result);

  return result;
}

template <typename KeyType, typename ValueType, class KeyComparator,
          class KeyEqualityChecker, class ValueComparator, class ValueEqualityChecker>
std::string BWTreeIndex<KeyType, ValueType, KeyComparator,
                        KeyEqualityChecker,
                        ValueComparator, ValueEqualityChecker>::GetTypeName() const {
  return "BWTree";
}

// Explicit template instantiation
template class BWTreeIndex<IntsKey<1>, ItemPointer, IntsComparator<1>,
                           IntsEqualityChecker<1>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<IntsKey<2>, ItemPointer, IntsComparator<2>,
                           IntsEqualityChecker<2>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<IntsKey<3>, ItemPointer, IntsComparator<3>,
                           IntsEqualityChecker<3>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<IntsKey<4>, ItemPointer, IntsComparator<4>,
                           IntsEqualityChecker<4>, ItemPointerComparator, ItemPointerEqualityChecker>;

template class BWTreeIndex<GenericKey<4>, ItemPointer, GenericComparator<4>,
                           GenericEqualityChecker<4>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<8>, ItemPointer, GenericComparator<8>,
                           GenericEqualityChecker<8>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<12>, ItemPointer, GenericComparator<12>,
                           GenericEqualityChecker<12>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<16>, ItemPointer, GenericComparator<16>,
                           GenericEqualityChecker<16>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<24>, ItemPointer, GenericComparator<24>,
                           GenericEqualityChecker<24>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<32>, ItemPointer, GenericComparator<32>,
                           GenericEqualityChecker<32>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<48>, ItemPointer, GenericComparator<48>,
                           GenericEqualityChecker<48>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<64>, ItemPointer, GenericComparator<64>,
                           GenericEqualityChecker<64>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<96>, ItemPointer, GenericComparator<96>,
                           GenericEqualityChecker<96>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<128>, ItemPointer, GenericComparator<128>,
                           GenericEqualityChecker<128>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<256>, ItemPointer, GenericComparator<256>,
                           GenericEqualityChecker<256>, ItemPointerComparator, ItemPointerEqualityChecker>;
template class BWTreeIndex<GenericKey<512>, ItemPointer, GenericComparator<512>,
                           GenericEqualityChecker<512>, ItemPointerComparator, ItemPointerEqualityChecker>;

template class BWTreeIndex<TupleKey, ItemPointer, TupleKeyComparator,
                           TupleKeyEqualityChecker, ItemPointerComparator, ItemPointerEqualityChecker>;

}  // End index namespace
}  // End peloton namespace
