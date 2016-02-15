/*-------------------------------------------------------------------------
 *
 * hash_join.cpp
 * file description
 *
 * Copyright(c) 2015, CMU
 *
 * /peloton/src/executor/hash_join_executor.cpp
 *
 *-------------------------------------------------------------------------
 */

#include <vector>

#include "backend/common/types.h"
#include "backend/common/logger.h"
#include "backend/executor/logical_tile_factory.h"
#include "backend/executor/hash_join_executor.h"
#include "backend/expression/abstract_expression.h"
#include "backend/expression/container_tuple.h"

namespace peloton {
namespace executor {

/**
 * @brief Constructor for hash join executor.
 * @param node Hash join node corresponding to this executor.
 */
HashJoinExecutor::HashJoinExecutor(const planner::AbstractPlan *node,
                                   ExecutorContext *executor_context)
    : AbstractJoinExecutor(node, executor_context) {}

bool HashJoinExecutor::DInit() {
  assert(children_.size() == 2);

  auto status = AbstractJoinExecutor::DInit();
  if (status == false) return status;

  assert(children_[1]->GetRawNode()->GetPlanNodeType() == PLAN_NODE_TYPE_HASH);

  hash_executor_ = reinterpret_cast<HashExecutor *>(children_[1]);

  return true;
}

/**
 * @brief Creates logical tiles from the two input logical tiles after applying
 * join predicate.
 * @return true on success, false otherwise.
 */
bool HashJoinExecutor::DExecute() {
  LOG_INFO("********** Hash Join %s Join executor :: 2 children \n",
           GetJoinTypeString());

  // If there is still something in the buffer, we return the first tile in it
  if (buffered_output_tiles.size() > 0) {
    SetOutput(buffered_output_tiles[0]);
    buffered_output_tiles.pop_front();
    return true;
  }

  //===--------------------------------------------------------------------===//
  // Pick all right tiles
  //===--------------------------------------------------------------------===//

  if (right_child_done_ != true) {
    // build hash table on right tiles
    while (hash_executor_->Execute() == true) {
      BufferRightTile(hash_executor_->GetOutput());
    }
    right_child_done_ = true;
  }

  if (right_result_tiles_.empty()) {
    assert(left_result_tiles_.empty());
    LOG_TRACE("Right child returned nothing. Exit.");
    return false;
  }

  // Get the hash table from the hash executor
  HashExecutor::HashMapType &hash_table_ = hash_executor_->GetHashTable();
  // Get the address of column_ids
  auto col_ids_address = &(hash_executor_->GetHashKeyIds());

  // Loop until we have non-empty result join logical tile or exit
  for (;;) {
    // Build outer join output when done
    if (left_child_done_) {
      assert(right_child_done_ == true);
      return BuildOuterJoinOutput();
    }

    //===--------------------------------------------------------------------===//
    // Pick next left tiles
    //===--------------------------------------------------------------------===//
    LogicalTile *left_tile = nullptr;

    // Left child is finished, no more tiles
    if (children_[0]->Execute() == false) {
      LOG_TRACE("Left child is exhausted. Returning false.");

      // Left child exhausted.
      // Release cur left tile. Clear right child's result buffer and return.
      assert(right_result_tiles_.size() > 0);
      left_child_done_ = true;

      // hash_table_.clear();
      return BuildOuterJoinOutput();
    }
      // Buffer the left child's result
    else {
      LOG_TRACE("Advance the left child.");
      BufferLeftTile(children_[0]->GetOutput());
    }

    left_tile = left_result_tiles_.back().get();

    //===--------------------------------------------------------------------===//
    // Build Join Tile
    //===--------------------------------------------------------------------===//
    std::unordered_map<size_t, std::unique_ptr<LogicalTile>> output_tile_map;
    std::unordered_map<size_t, LogicalTile::PositionListsBuilder>
        pos_lists_builder_map;

    // Go over the left logical tile
    for (auto left_tile_row_itr : *left_tile) {
      auto target_key = HashExecutor::HashMapType::key_type(
          left_tile, left_tile_row_itr, col_ids_address);

      auto got = hash_table_.find(target_key);

      // There is no corresponding right tuples. Skip this left tuple and
      // continue.
      if (got == hash_table_.end()) {
        continue;
      }

      // Go over the matching right tuples
      // Key : container tuple with a subset of tuple attributes
      // Value : < child_tile offset, tuple offset >
      auto right_matched_tuples = got->second;
      for (auto value : right_matched_tuples) {
        // For Right and Full Outer Join
        RecordMatchedRightRow(value.first, value.second);

        // If this left-right pair has not been recorded
        if (output_tile_map.find(value.first) == output_tile_map.end()) {
          // we have to firstly update the right_tile,
          // because this right tuple can belong to a tile which is
          // different from the one we use to initialize our pos_lists_builder
          LogicalTile *current_right_tile =
              right_result_tiles_[value.first].get();

          // Build and record output join logical tile
          output_tile_map[value.first] =
              BuildOutputLogicalTile(left_tile, current_right_tile);

          // Build and record corresponding position lists
          LogicalTile::PositionListsBuilder pos_lists_builder(
              left_tile, current_right_tile);
          pos_lists_builder_map[value.first] = pos_lists_builder;
        }

        // Insert a tuple into the correct output logical tile in the hashmap
        pos_lists_builder_map[value.first].AddRow(left_tile_row_itr,
                                                  value.second);
      }

      // For Left and Full Outer Join
      // if we get some matched right tuples for this left tuple,
      // we remove it from the no_matching_left_row_sets_
      RecordMatchedLeftRow(left_result_tiles_.size() - 1, left_tile_row_itr);
    }

    // fill the buffered_output_tiles with values in hashmap
    for (auto iter = output_tile_map.begin(); iter != output_tile_map.end();
         ++iter) {
      iter->second->SetPositionListsAndVisibility(
          pos_lists_builder_map[iter->first].Release());
      buffered_output_tiles.push_back(iter->second.release());
    }

    // Check if we have any join tuples
    if (buffered_output_tiles.size() > 0) {
      SetOutput(buffered_output_tiles[0]);
      buffered_output_tiles.pop_front();
      return true;
    }

    LOG_TRACE("This left tile produces empty join result. Continue the loop.");
  }
}

}  // namespace executor
}  // namespace peloton
