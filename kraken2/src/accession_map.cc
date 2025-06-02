/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "accession_map.h"
#include <fstream>
#include <iostream>
#include <cstring>

namespace kraken2 {

// Static member initialization
bool AccessionMap::debug_mode_ = false;

AccessionMap::AccessionMap() : file_backed_(false), file_data_(nullptr), file_size_(0), initialized_(true) {
  DebugPrint("AccessionMap: Default constructor called");
}

AccessionMap::AccessionMap(const std::string &filename, bool memory_mapping) 
    : file_backed_(false), file_data_(nullptr), file_size_(0), initialized_(false) {
  DebugPrint("AccessionMap: Constructor with filename: " + filename);
  try {
    LoadFromFile(filename, memory_mapping);
    initialized_ = true;
    DebugPrint("AccessionMap: Successfully loaded from file");
  } catch (const std::exception &e) {
    DebugPrint("AccessionMap: Exception during construction: " + std::string(e.what()));
    initialized_ = false;
  }
}

AccessionMap::~AccessionMap() {
  DebugPrint("AccessionMap: Destructor called");
  if (file_backed_ && file_data_) {
    try {
      backing_file_.CloseFile();
      DebugPrint("AccessionMap: Memory-mapped file closed");
    } catch (const std::exception &e) {
      DebugPrint("AccessionMap: Exception during destruction: " + std::string(e.what()));
    }
  }
}

void AccessionMap::DebugPrint(const std::string &message) const {
  if (debug_mode_) {
    std::cerr << "[AccessionMap DEBUG] " << message << std::endl;
  }
}

std::string AccessionMap::GetAccession(uint64_t minimizer) const {
  if (!initialized_) {
    DebugPrint("GetAccession called on uninitialized object");
    return "";
  }

  //DebugPrint("GetAccession called for minimizer: " + std::to_string(minimizer));
  
  if (file_backed_) {
    // For memory-mapped files, we need to search through the file data
    // Format: [count][minimizer1][accession_len1][accession1]...
    if (!file_data_ || file_size_ < sizeof(uint64_t)) {
      DebugPrint("Invalid file data in memory-mapped mode");
      return "";
    }
    
    try {
      char* ptr = file_data_;
      uint64_t count = *reinterpret_cast<uint64_t*>(ptr);
      ptr += sizeof(uint64_t);
      
      DebugPrint("Searching through " + std::to_string(count) + " entries in memory-mapped file");
      
      for (size_t i = 0; i < count && (ptr - file_data_) < (ptrdiff_t)file_size_; i++) {
        if ((ptr - file_data_) + sizeof(uint64_t) + sizeof(uint32_t) > (ptrdiff_t)file_size_) {
          DebugPrint("Memory bounds check failed at entry " + std::to_string(i));
          break;
        }
        
        uint64_t stored_minimizer = *reinterpret_cast<uint64_t*>(ptr);
        ptr += sizeof(uint64_t);
        
        uint32_t accession_len = *reinterpret_cast<uint32_t*>(ptr);
        ptr += sizeof(uint32_t);
        
        if ((ptr - file_data_) + accession_len > (ptrdiff_t)file_size_) {
          DebugPrint("Accession length would exceed file bounds");
          break;
        }
        
        if (stored_minimizer == minimizer) {
          std::string result(ptr, accession_len);
          DebugPrint("Found accession: " + result);
          return result;
        }
        
        ptr += accession_len;
      }
      DebugPrint("Minimizer not found in memory-mapped file");
      return "";
    } catch (const std::exception &e) {
      DebugPrint("Exception in GetAccession (memory-mapped): " + std::string(e.what()));
      return "";
    }
  } else {
    std::lock_guard<std::mutex> lock(mapping_mutex_);
    auto it = mappings_.find(minimizer);
    if (it != mappings_.end()) {
      //DebugPrint("Found accession in memory map: " + it->second);
      return it->second;
    } else {
      DebugPrint("Minimizer not found in memory map");
      return "";
    }
  }
}

void AccessionMap::AddMapping(uint64_t minimizer, const std::string &accession_id) {
  if (!initialized_) {
    DebugPrint("AddMapping called on uninitialized object");
    return;
  }

  if (file_backed_) {
    DebugPrint("Cannot add mappings to file-backed AccessionMap");
    errx(EX_SOFTWARE, "Cannot add mappings to file-backed AccessionMap");
  }
  
  //DebugPrint("Adding mapping: " + std::to_string(minimizer) + " -> " + accession_id);
  
  try {
    std::lock_guard<std::mutex> lock(mapping_mutex_);
    mappings_[minimizer] = accession_id;
    //DebugPrint("Successfully added mapping (total: " + std::to_string(mappings_.size()) + ")");
  } catch (const std::exception &e) {
    DebugPrint("Exception in AddMapping: " + std::string(e.what()));
  }
}

void AccessionMap::WriteToFile(const std::string &filename) const {
  if (!initialized_) {
    DebugPrint("WriteToFile called on uninitialized object");
    errx(EX_SOFTWARE, "Cannot write uninitialized AccessionMap to file");
  }

  if (file_backed_) {
    DebugPrint("Cannot write file-backed AccessionMap to file");
    errx(EX_SOFTWARE, "Cannot write file-backed AccessionMap to file");
  }
  
  DebugPrint("Writing AccessionMap to file: " + filename);
  
  try {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
      DebugPrint("Failed to create output file");
      errx(EX_CANTCREAT, "Cannot create accession map file: %s", filename.c_str());
    }
    
    std::lock_guard<std::mutex> lock(mapping_mutex_);
    
    // Write the number of entries
    uint64_t count = mappings_.size();
    DebugPrint("Writing " + std::to_string(count) + " entries");
    outfile.write(reinterpret_cast<const char*>(&count), sizeof(count));
    
    // Write each mapping
    for (const auto &pair : mappings_) {
      // Write minimizer
      outfile.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
      
      // Write accession length and accession string
      uint32_t accession_len = pair.second.length();
      outfile.write(reinterpret_cast<const char*>(&accession_len), sizeof(accession_len));
      outfile.write(pair.second.c_str(), accession_len);
    }
    
    outfile.close();
    if (!outfile) {
      DebugPrint("Error occurred while writing file");
      errx(EX_IOERR, "Error writing accession map file: %s", filename.c_str());
    }
    
    DebugPrint("Successfully wrote AccessionMap to file");
  } catch (const std::exception &e) {
    DebugPrint("Exception in WriteToFile: " + std::string(e.what()));
    throw;
  }
}

void AccessionMap::LoadFromFile(const std::string &filename, bool memory_mapping) {
  DebugPrint("Loading AccessionMap from file: " + filename + " (memory_mapping=" + (memory_mapping ? "true" : "false") + ")");
  
  file_backed_ = memory_mapping;
  
  try {
    if (memory_mapping) {
      DebugPrint("Using memory mapping");
      backing_file_.OpenFile(filename);
      file_data_ = backing_file_.fptr();
      file_size_ = backing_file_.filesize();
      DebugPrint("Memory-mapped file size: " + std::to_string(file_size_) + " bytes");
    } else {
      DebugPrint("Loading into memory");
      std::ifstream infile(filename, std::ios::binary);
      if (!infile) {
        DebugPrint("Cannot open input file");
        errx(EX_NOINPUT, "Cannot open accession map file: %s", filename.c_str());
      }
      
      // Read the number of entries
      uint64_t count;
      infile.read(reinterpret_cast<char*>(&count), sizeof(count));
      if (!infile) {
        DebugPrint("Error reading file header");
        errx(EX_DATAERR, "Error reading accession map file header: %s", filename.c_str());
      }
      
      DebugPrint("Reading " + std::to_string(count) + " entries from file");
      
      std::lock_guard<std::mutex> lock(mapping_mutex_);
      
      // Read each mapping
      for (uint64_t i = 0; i < count; i++) {
        uint64_t minimizer;
        infile.read(reinterpret_cast<char*>(&minimizer), sizeof(minimizer));
        
        uint32_t accession_len;
        infile.read(reinterpret_cast<char*>(&accession_len), sizeof(accession_len));
        
        if (accession_len > 10000) {  // Sanity check
          DebugPrint("Suspicious accession length: " + std::to_string(accession_len));
          errx(EX_DATAERR, "Suspicious accession length %u at entry %llu in %s", 
               accession_len, (unsigned long long)i, filename.c_str());
        }
        
        std::string accession_id(accession_len, '\0');
        infile.read(&accession_id[0], accession_len);
        
        if (!infile) {
          DebugPrint("Error reading entry " + std::to_string(i));
          errx(EX_DATAERR, "Error reading accession map file entry %llu: %s", 
               (unsigned long long)i, filename.c_str());
        }
        
        mappings_[minimizer] = accession_id;
      }
      
      infile.close();
      DebugPrint("Successfully loaded " + std::to_string(mappings_.size()) + " mappings");
    }
  } catch (const std::exception &e) {
    DebugPrint("Exception in LoadFromFile: " + std::string(e.what()));
    throw;
  }
}

void AccessionMap::LoadFromMemoryMappedFile() {
  // This method is called internally when using memory mapping
  // The actual loading is done in LoadFromFile
  DebugPrint("LoadFromMemoryMappedFile called (no-op)");
}

}  // end namespace 