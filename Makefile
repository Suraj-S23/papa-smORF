# Variables
PYTHON := python
PIP := pip3  # or python -m pip
BASE_DIR := $(shell pwd)
DATA_DIR := data
FASTA_DIR := $(DATA_DIR)/fasta
GFF_DIR := $(DATA_DIR)/gff
OUTPUT_DIR := $(DATA_DIR)/processed
LOGS_DIR := logs
N_WORKERS := 8
SCRIPT_ARGS := --fasta_dir $(FASTA_DIR) --gff_dir $(GFF_DIR) --output_dir $(OUTPUT_DIR) --workers $(N_WORKERS)

.PHONY: setup install clean pipeline run-quick run-all

# Setup project structure
setup:
	$(PYTHON) -c "from src.preprocessing.utils import setup_directories; setup_directories('$(BASE_DIR)')"

# Install dependencies and package
install:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e .

# Clean generated files
clean:
	rm -rf $(OUTPUT_DIR)/*
	rm -rf $(LOGS_DIR)/*
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.log" -delete

# Run main pipeline
pipeline:
	$(PYTHON) src/pipeline.py $(SCRIPT_ARGS)

# Quick run (skip setup/install)
run-quick:
	$(MAKE) pipeline

# Complete run with setup
run-all:
	@echo "=== Starting papa-smORF Pipeline ==="
	@echo "1. Cleaning and setting up..."
	$(MAKE) clean setup
	@echo "\n2. Installing dependencies..."
	$(MAKE) install
	@echo "\n3. Running pipeline..."
	$(MAKE) pipeline
	@echo "\n=== Pipeline Complete ==="
	@echo "Check:"
	@echo "- logs/ for execution logs"
	@echo "- data/processed/ for processed data"

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
