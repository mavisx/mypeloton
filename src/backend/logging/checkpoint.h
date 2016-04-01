/*-------------------------------------------------------------------------
 *
 * checkpoint.h
 * file description
 *
 * Copyright(c) 2015, CMU
 *
 * /peloton/src/backend/logging/checkpoint.h
 *
 *-------------------------------------------------------------------------
 */

#pragma once
#include <cassert>
#include <string>
#include <sys/stat.h>

#include "backend/logging/log_manager.h"
#include "backend/logging/backend_logger.h"
#include "backend/common/pool.h"

namespace peloton {
namespace logging {

//===--------------------------------------------------------------------===//
// Checkpoint
//===--------------------------------------------------------------------===//
class Checkpoint {
 public:
  Checkpoint() { pool.reset(new VarlenPool(BACKEND_TYPE_MM)); }

  virtual ~Checkpoint(void) {}

  // Do checkpoint periodically
  virtual void DoCheckpoint() = 0;

  // Initialize resources for checkpoint
  virtual void Init() = 0;

  // Do recovery from most recent version of checkpoint
  virtual bool DoRecovery() = 0;

 protected:
  std::string ConcatFileName(std::string checkpoint_dir, int version);

  void InitDirectory();

  // variable length memory pool
  std::unique_ptr<VarlenPool> pool;

  std::string checkpoint_dir = "pl_checkpoint";

  // the version of next checkpoint
  int checkpoint_version = -1;

  // prefix for checkpoint file name
  const std::string FILE_PREFIX = "peloton_checkpoint_";

  // suffix for checkpoint file name
  const std::string FILE_SUFFIX = ".log";
};

}  // namespace logging
}  // namespace peloton
