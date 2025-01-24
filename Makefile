# Variables
PYTHON := python
CONDA_ENV_NAME := smorfenv
CONDA_ENV_FILE := environment.yml
BASE_DIR := $(shell pwd)
DATA_DIR := data
FASTA_DIR := $(DATA_DIR)/large_fasta/fna_files
GFF_DIR := $(DATA_DIR)/large_gff/gff_files
OUTPUT_DIR := $(DATA_DIR)/large_processed
LOGS_DIR := logs
N_WORKERS := 8
BATCH_SIZE := 100
SCRIPT_ARGS := --fasta_dir $(FASTA_DIR) --gff_dir $(GFF_DIR) --output_dir $(OUTPUT_DIR) --workers $(N_WORKERS) --batch_size $(BATCH_SIZE)
PYTHONPATH := $(BASE_DIR)
.PHONY: setup install clean pipeline run-quick run-all help

# Setup project structure
setup:
	$(PYTHON) -c "from src.preprocessing.utils import setup_directories; setup_directories('$(BASE_DIR)')"
	@echo "Directory structure created"

setup-env:
	@echo "Setting up conda environment..."
	conda env create -f $(CONDA_ENV_FILE) || conda env update -f $(CONDA_ENV_FILE)
	@echo "Activating environment..."
	@echo "Run 'conda activate $(CONDA_ENV_NAME)' to activate the environment"

# Update the install target
install:
	python -m pip install --upgrade pip
	python -m pip install -r requirements.txt
	python -m pip install -e .

# Clean generated files
clean:
	rm -rf $(OUTPUT_DIR)/*
	rm -rf $(LOGS_DIR)/*
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.mmap" -delete
	find . -type f -name "*.log" -delete
	@echo "Cleaned generated files"

# Run main pipeline
pipeline:
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) src/pipeline.py $(SCRIPT_ARGS)

# Quick run (skip setup/install)
run-quick:
	@echo "Running pipeline..."
	$(MAKE) pipeline

# Complete run with setup
run-all:
	@echo "=== Starting papa-smORF Pipeline with Memory Mapping ==="
	@echo "1. Cleaning and setting up..."
	$(MAKE) clean setup
	@echo "\n2. Installing dependencies..."
	$(MAKE) install
	@echo "\n3. Running pipeline..."
	$(MAKE) pipeline
	@echo "\n=== Pipeline Complete ==="
	@echo "Check:"
	@echo "- logs/ for execution logs"
	@echo "- data/processed/ for memory-mapped data files (.mmap)"

# Help target
help:
	@echo "Available targets:"
	@echo "  setup      - Create necessary directories"
	@echo "  install    - Install dependencies"
	@echo "  clean      - Remove generated files"
	@echo "  pipeline   - Run main processing"
	@echo "  run-quick  - Run pipeline without setup"
	@echo "  run-all    - Run complete pipeline with fresh setup"
	@echo "  help       - Show this help message"
	@echo "\nEnvironment:"
	@echo "  Workers: $(N_WORKERS)"
	@echo "  Batch size: $(BATCH_SIZE)"
