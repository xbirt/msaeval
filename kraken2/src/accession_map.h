/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_ACCESSION_MAP_H_
#define KRAKEN2_ACCESSION_MAP_H_

#include "kraken2_headers.h"
#include "mmap_file.h"
#include <string>
#include <unordered_map>
#include <mutex>

namespace kraken2 {

class AccessionMap {
public:
  AccessionMap();
  AccessionMap(const std::string &filename, bool memory_mapping = false);
  ~AccessionMap();

  // Get accession ID for a given minimizer hash
  std::string GetAccession(uint64_t minimizer) const;
  
  // Add a mapping from minimizer to accession ID (thread-safe)
  void AddMapping(uint64_t minimizer, const std::string &accession_id);
  
  // Write the mapping to a file
  void WriteToFile(const std::string &filename) const;
  
  // Load mapping from file
  void LoadFromFile(const std::string &filename, bool memory_mapping = false);
  
  // Get the number of mappings
  size_t Size() const { 
    std::lock_guard<std::mutex> lock(mapping_mutex_);
    return mappings_.size(); 
  }
  
  // Check if empty
  bool Empty() const { 
    std::lock_guard<std::mutex> lock(mapping_mutex_);
    return mappings_.empty(); 
  }

  // Debug function to check object state
  bool IsValid() const { return initialized_; }
  
  // Enable/disable debug output
  static void SetDebugMode(bool debug) { debug_mode_ = debug; }

private:
  std::unordered_map<uint64_t, std::string> mappings_;
  bool file_backed_;
  MMapFile backing_file_;
  char* file_data_;
  size_t file_size_;
  bool initialized_;
  
  mutable std::mutex mapping_mutex_;  // Thread safety for mappings_
  static bool debug_mode_;
  
  void LoadFromMemoryMappedFile();
  void DebugPrint(const std::string &message) const;
  AccessionMap(const AccessionMap &rhs) = delete;
  AccessionMap& operator=(const AccessionMap &rhs) = delete;
};

}  // end namespace

#endif 