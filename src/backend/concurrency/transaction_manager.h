//===----------------------------------------------------------------------===//
//
//                         Peloton
//
// transaction_manager.h
//
// Identification: src/backend/concurrency/transaction_manager.h
//
// Copyright (c) 2015-16, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include <atomic>
#include <unordered_map>

#include "backend/common/platform.h"
#include "backend/common/types.h"
#include "backend/concurrency/transaction.h"
#include "backend/storage/data_table.h"
#include "backend/storage/tile_group.h"
#include "backend/storage/tile_group_header.h"

#include "libcuckoo/cuckoohash_map.hh"

namespace peloton {
namespace concurrency {

extern thread_local Transaction *current_txn;

#define RUNNING_TXN_BUCKET_NUM 10

class TransactionManager {
 public:
  TransactionManager() {
    next_txn_id_ = ATOMIC_VAR_INIT(START_TXN_ID);
    next_cid_ = ATOMIC_VAR_INIT(START_CID);
  }

  virtual ~TransactionManager() {}

  txn_id_t GetNextTransactionId() { return next_txn_id_++; }

  cid_t GetNextCommitId() { return next_cid_++; }

  virtual bool IsVisible(
      const storage::TileGroupHeader *const tile_group_header,
      const oid_t &tuple_id) = 0;

  bool IsVisbleOrDirty(__attribute__((unused)) const storage::Tuple *key, const ItemPointer &position) {
    auto tile_group_header = catalog::Manager::GetInstance().GetTileGroup(position.block)->GetHeader();
    auto tuple_id = position.offset;

    txn_id_t tuple_txn_id = tile_group_header->GetTransactionId(tuple_id);
    cid_t tuple_begin_cid = tile_group_header->GetBeginCommitId(tuple_id);
    cid_t tuple_end_cid = tile_group_header->GetEndCommitId(tuple_id);
    if (tuple_txn_id == INVALID_TXN_ID) {
      // the tuple is not available.
      return false;
    }
    bool own = (current_txn->GetTransactionId() == tuple_txn_id);

    // there are exactly two versions that can be owned by a transaction.
    // unless it is an insertion.
    if (own == true) {
      if (tuple_begin_cid == MAX_CID && tuple_end_cid != INVALID_CID) {
        assert(tuple_end_cid == MAX_CID);
        // the only version that is visible is the newly inserted one.
        return true;
      } else {
        // the older version is not visible.
        return false;
      }
    } else {
//      // if the tuple is owned by other transactions or not owned by anyone
//      bool activated = (current_txn->GetBeginCommitId() >= tuple_begin_cid);
//      bool invalidated = (current_txn->GetBeginCommitId() >= tuple_end_cid);
//      if ((!invalidated) && (activated || tuple_begin_cid == MAX_CID)) {
//        // return true when:
//        //  1) valid and activated version owned by no one
//        //  2) valid and activated version owned by others
//        //  3) dirty insert and update owned by others, for an update of others we see two versions
//        return true;
//      }
//      // return false when:
//      //    1) Dirty delete
//      //    2) Old versions
//      return false;
      bool activated = (current_txn->GetBeginCommitId() >= tuple_begin_cid);
      bool invalidated = (current_txn->GetBeginCommitId() >= tuple_end_cid);
      if (tuple_txn_id != INITIAL_TXN_ID) {
        // if the tuple is owned by other transactions.
        if (tuple_begin_cid == MAX_CID) {
          // uncommitted version.
          if (tuple_end_cid == INVALID_CID) {
            // dirty delete is invisible
            return false;
          } else {
            // dirty update or insert is visible
            return true;
          }
        } else {
          // the older version may be visible.
          if (activated && !invalidated) {
            return true;
          } else {
            return false;
          }
        }
      } else {
        // if the tuple is not owned by any transaction.
        if (activated && !invalidated) {
          return true;
        } else {
          return false;
        }
      }
    }
  }

  virtual bool IsOwner(const storage::TileGroupHeader *const tile_group_header,
                       const oid_t &tuple_id) = 0;

  virtual bool IsOwnable(
      const storage::TileGroupHeader *const tile_group_header,
      const oid_t &tuple_id) = 0;

  virtual bool AcquireOwnership(
      const storage::TileGroupHeader *const tile_group_header,
      const oid_t &tile_group_id, const oid_t &tuple_id) = 0;

  virtual void SetOwnership(const oid_t &tile_group_id,
                            const oid_t &tuple_id) = 0;

  virtual bool PerformInsert(const oid_t &tile_group_id,
                             const oid_t &tuple_id) = 0;

  virtual bool PerformRead(const oid_t &tile_group_id,
                           const oid_t &tuple_id) = 0;

  virtual bool PerformUpdate(const oid_t &tile_group_id, const oid_t &tuple_id,
                             const ItemPointer &new_location) = 0;

  virtual bool PerformDelete(const oid_t &tile_group_id, const oid_t &tuple_id,
                             const ItemPointer &new_location) = 0;

  virtual void PerformUpdate(const oid_t &tile_group_id,
                             const oid_t &tuple_id) = 0;

  virtual void PerformDelete(const oid_t &tile_group_id,
                             const oid_t &tuple_id) = 0;

  void SetTransactionResult(const Result result) {
    current_txn->SetResult(result);
  }

  //for use by recovery
  void SetNextCid(cid_t cid) { next_cid_ = cid; }

  virtual Transaction *BeginTransaction() {
    txn_id_t txn_id = GetNextTransactionId();
    cid_t begin_cid = GetNextCommitId();
    Transaction *txn = new Transaction(txn_id, begin_cid);
    current_txn = txn;
    
    RegisterTransaction(txn_id, begin_cid);
    return txn;
  }

  virtual void EndTransaction() {
    txn_id_t txn_id = current_txn->GetTransactionId();
    DeregisterTransaction(txn_id);
    
    delete current_txn;
    current_txn = nullptr;
  }

  virtual Result CommitTransaction() = 0;

  virtual Result AbortTransaction() = 0;

  void ResetStates() {
    next_txn_id_ = START_TXN_ID;
    next_cid_ = START_CID;
  }


  void RegisterTransaction(const txn_id_t &txn_id, const cid_t &begin_cid) {
    running_txn_buckets_[txn_id % RUNNING_TXN_BUCKET_NUM][txn_id] = begin_cid;
  }

  void DeregisterTransaction(const txn_id_t &txn_id) {
    running_txn_buckets_[txn_id % RUNNING_TXN_BUCKET_NUM].erase(txn_id);
  }

  // this function generates the maximum commit id of committed transactions.
  // please note that this function only returns a "safe" value instead of a precise value.
  cid_t GetMaxCommittedCid() {
    cid_t min_running_cid = 0;
    for (size_t i = 0; i < RUNNING_TXN_BUCKET_NUM; ++i) {
      {
        auto iter = running_txn_buckets_[i].lock_table();
        for (auto &it : iter) {
          if (min_running_cid == 0 || it.second < min_running_cid) {
            min_running_cid = it.second;
          }
        }
      }
    }
    assert(min_running_cid > 0);
    return min_running_cid - 1;
  }

 private:
  std::atomic<txn_id_t> next_txn_id_;
  std::atomic<cid_t> next_cid_;

  cuckoohash_map<txn_id_t, cid_t> running_txn_buckets_[RUNNING_TXN_BUCKET_NUM];
  
};
}  // End storage namespace
}  // End peloton namespace
