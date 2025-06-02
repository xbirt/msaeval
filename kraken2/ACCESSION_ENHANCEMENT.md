# Kraken2 Accession ID Enhancement

This enhancement adds accession ID tracking and reporting to Kraken2, allowing users to see which database sequence accession was the best match for classified sequences.

## Changes Made

### 1. New AccessionMap Class
- **Files**: `src/accession_map.h`, `src/accession_map.cc`
- **Purpose**: Stores and manages mappings between minimizer hashes and sequence accession IDs
- **Features**:
  - Memory-mapped file support for efficient loading
  - Binary file format for compact storage
  - Thread-safe access during classification

### 2. Database Building Enhancement
- **File**: `src/build_db.cc`
- **Changes**:
  - Added `-a` option to specify accession map output filename
  - Tracks accession IDs (first sequence ID from FASTA headers) during database creation
  - Creates accession map file alongside the main database

### 3. Classification Enhancement
- **File**: `src/classify.cc`
- **Changes**:
  - Added `-A` option to specify accession map input filename
  - Loads accession map during database loading
  - Tracks best accession ID match during sequence classification
  - Includes accession ID in classified sequence headers

### 4. Output Format Enhancement
- **Classified sequences** now include accession ID in their headers when available:
  - Format: `kraken:taxid|TAXID|ACCESSION_ID`
  - Falls back to: `kraken:taxid|TAXID` when no accession available

## Usage

### Building Database with Accession Mapping
```bash
build_db [existing options] -a accession_map.k2map
```

### Classification with Accession Information
```bash
classify [existing options] -A accession_map.k2map
```

## File Format

The accession map file (`.k2map`) uses a binary format:
- Header: 8-byte count of entries
- For each entry:
  - 8-byte minimizer hash
  - 4-byte accession string length
  - Variable-length accession string

## Benefits

1. **Sequence Provenance**: Users can trace classified sequences back to specific database entries
2. **Quality Assessment**: Helps evaluate classification confidence by examining source sequences
3. **Debugging**: Facilitates troubleshooting classification results
4. **Downstream Analysis**: Enables more detailed analysis of classification results

## Performance Impact

- **Memory**: Minimal additional memory usage during classification
- **Speed**: Very small performance impact due to hash lookups
- **Storage**: Additional accession map file (size depends on database complexity)

## Backward Compatibility

- All changes are backward compatible
- Accession mapping is optional (controlled by command-line flags)
- Existing workflows continue to work unchanged 