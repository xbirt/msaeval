# rRNA 16S Sequence Alignment Benchmarking Suite

This repository contains a comprehensive suite of tools for benchmarking and evaluating RNA sequence alignment tools. The project includes installation scripts for various alignment tools, benchmarking scripts, and data processing utilities.

## Project Structure

### Installation Scripts (`/install`)
Contains shell scripts for installing various sequence alignment tools:
- Kalign
- ProbCons
- Clustal Omega
- FastSP
- MAFFT
- Mothur
- MUSCLE
- NAST
- PASTA
- TPMA
- UPP
- T-Coffee
- Synthetic Data Tools

### Benchmarking Scripts (`/benchmark`)
Contains scripts for running benchmarks on different alignment tools:
- Main benchmark script (`benchmark.sh`)
- Individual tool benchmark scripts:
  - Kalign
  - ProbCons
  - PASTA
  - T-Coffee
  - Clustal Omega
  - MAFFT
  - MUSCLE

### Data Processing (`/data`)
Contains Python scripts and shell scripts for data processing and evaluation:

#### Python Scripts
- `score.py`: Main scoring module for evaluating alignments
- `evaluate.py`: Evaluation utilities for alignment results
- `compute_coverage.py`: Coverage computation utilities
- `reference_dataset.py`: Reference dataset handling
- `reference_aligned_dataset.py`: Aligned reference dataset handling
- `abundance_model.py`: Abundance modeling utilities

#### Shell Scripts
- `generate_data.sh`: Data generation utilities
- `compute_coverage.sh`: Coverage computation script
- `build_references.sh`: Reference sequence building
- `align_references.sh`: Reference sequence alignment
- `filter_species.sh`: Species filtering utilities
- `porechop.sh`: Porechop processing script
- `get_silva_db.sh`: SILVA database retrieval

## Features

1. **Multiple Alignment Tool Support**
   - Comprehensive support for popular alignment tools
   - Easy installation scripts for each tool
   - Individual benchmarking capabilities

2. **Data Processing**
   - Reference dataset generation and management
   - Alignment evaluation and scoring
   - Coverage computation
   - Species filtering
   - Abundance modeling

3. **Benchmarking**
   - Standardized benchmarking across tools
   - Performance metrics collection
   - Comparative analysis capabilities

## Usage

### Installation

To install a specific alignment tool, use the corresponding installation script:

```bash
./install/install_[tool_name].sh
```

### Running Benchmarks

To run benchmarks for any tool:

```bash
./benchmark/benchmark.sh
```

For individual tool benchmarks:

```bash
./benchmark/benchmark_[tool_name].sh
```

### Data Processing

To generate synthetic data:

```bash
./data/generate_data.sh
```

To compute coverage:

```bash
./data/compute_coverage.sh
```

## Dependencies

- Python 3.x
- BioPython
- Various alignment tools (installed via provided scripts)
- Shell environment (bash)

## Notes

- The scoring system uses a custom RNA scoring matrix optimized for RNA sequence alignment
- Gap penalties are configurable (default: -10 for opening, -0.5 for extension)
- The system supports both standard scoring and alternative scoring with built-in gap penalties
