//===----------------------------------------------------------------------===//
//
//                         Peloton
//
// wal_frontend_logger.h
//
// Identification: src/backend/logging/loggers/wal_frontend_logger.h
//
// Copyright (c) 2015-16, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include "backend/logging/frontend_logger.h"
#include "backend/logging/records/tuple_record.h"
#include "backend/logging/log_file.h"
#include <dirent.h>
#include <vector>

namespace peloton {

class VarlenPool;

namespace concurrency {
class Transaction;
}

namespace logging {

//===--------------------------------------------------------------------===//
// Write Ahead Frontend Logger
//===--------------------------------------------------------------------===//

class WriteAheadFrontendLogger : public FrontendLogger {
 public:
  WriteAheadFrontendLogger(void);

  ~WriteAheadFrontendLogger(void);

  void FlushLogRecords(void);

  //===--------------------------------------------------------------------===//
  // Recovery
  //===--------------------------------------------------------------------===//

  void DoRecovery(void);

  void StartTransactionRecovery(cid_t commit_id);

  void CommitTransactionRecovery(cid_t commit_id);

  void InsertTuple(TupleRecord *recovery_txn);

  void DeleteTuple(TupleRecord *recovery_txn);

  void UpdateTuple(TupleRecord *recovery_txn);

  void AbortActiveTransactions();

  void InitLogFilesList();

  void CreateNewLogFile(bool);

  bool FileSwitchCondIsTrue();

  void OpenNextLogFile();

  LogRecordType GetNextLogRecordTypeForRecovery(FILE *, size_t);

  void TruncateLog(int);

  void SetLogDirectory(char *);

  void InitLogDirectory();

  std::string GetFileNameFromVersion(int);

 private:
  std::string GetLogFileName(void);

  //===--------------------------------------------------------------------===//
  // Member Variables
  //===--------------------------------------------------------------------===//

  // File pointer and descriptor
  FILE *log_file;
  int log_file_fd;

  // Size of the log file
  size_t log_file_size;

  // Txn table during recovery
  std::map<txn_id_t, std::vector<TupleRecord *>> recovery_txn_table;

  // Keep tracking max oid for setting next_oid in manager
  // For active processing after recovery
  oid_t max_oid = 0;

  cid_t max_cid = 0;

  // pool for allocating non-inlined values
  VarlenPool *recovery_pool;

  // abj1 adding code here!
  std::vector<LogFile *> log_files_;

  int log_file_counter_;

  int log_file_cursor_;

  std::string peloton_log_directory = "peloton_log";

  std::string LOG_FILE_PREFIX = "peloton_log_";

  std::string LOG_FILE_SUFFIX = ".log";
};

}  // namespace logging
}  // namespace peloton
