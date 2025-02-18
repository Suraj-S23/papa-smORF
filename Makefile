# Variables
PYTHON := python
CONDA_ENV_NAME := smorfenv
CONDA_ENV_FILE := environment.yml
BASE_DIR := $(shell pwd)
DATA_DIR := data
FASTA_DIR := $(DATA_DIR)/large_fasta/fna_files
GFF_DIR := $(DATA_DIR)/large_gff/gff_files
OUTPUT_DIR := $(DATA_DIR)/processed
LOGS_DIR := run_logs
N_WORKERS := 8
BATCH_SIZE := 100
SCRIPT_ARGS := --fasta_dir $(DATA_DIR)/10_fasta --gff_dir $(DATA_DIR)/10_gff --output_dir $(OUTPUT_DIR) --workers $(N_WORKERS) --batch_size $(BATCH_SIZE)
PYTHONPATH := $(BASE_DIR)
SMALL_FASTA_DIR := $(DATA_DIR)/10_fasta
SMALL_GFF_DIR := $(DATA_DIR)/10_gff
SMALL_PROCESSED_DIR := $(DATA_DIR)/processed
PROCESSED_DIR := $(DATA_DIR)/large_processed_2
ENCODED_DIR := $(DATA_DIR)/large_encoded

# Add HDF5 specific cleanup patterns
.PHONY: setup install pipeline run-quick run-all help clean-hdf5

# Setup project structure
setup:
	$(PYTHON) -c "from src.preprocessing.utils import setup_directories; setup_directories('$(BASE_DIR)')"
	@echo "Directory structure created"

setup-env:
	@echo "Setting up conda environment..."
	conda env create -f $(CONDA_ENV_FILE) || conda env update -f $(CONDA_ENV_FILE)
	@echo "Activating environment..."
	@echo "Run 'conda activate $(CONDA_ENV_NAME)' to activate the environment"

install:
	python -m pip install --upgrade pip
	python -m pip install -r requirements.txt
	python -m pip install -e .
	# Add h5py installation if not in requirements.txt
	python -m pip install h5py

# Clean generated files
clean:
	rm -rf $(OUTPUT_DIR)/*
	rm -rf $(LOGS_DIR)/*
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.h5" -delete
	find . -type f -name "*.log" -delete
	@echo "Cleaned generated files"

# Clean only HDF5 files
clean-hdf5:
	find . -type f -name "*.h5" -delete
	@echo "Cleaned HDF5 files"

# Run main pipeline
pipeline:
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) src/pipeline.py $(SCRIPT_ARGS)

# Quick run (skip setup/install)
run-quick:
	@echo "Running pipeline..."
	$(MAKE) pipeline

# Complete run with setup
run-all:
	@echo "=== Starting papa-smORF Pipeline with HDF5 Storage ==="
	@echo "1. Cleaning and setting up..."
	$(MAKE) clean setup
	@echo "\n2. Installing dependencies..."
	$(MAKE) install
	@echo "\n3. Running pipeline..."
	$(MAKE) pipeline
	@echo "\n=== Pipeline Complete ==="
	@echo "Check:"
	@echo "- logs/ for execution logs"
	@echo "- data/processed/ for HDF5 data files (.h5)"

# Process small dataset
process-small:
	@echo "Processing small dataset..."
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) src/pipeline.py \
	--fasta_dir $(SMALL_FASTA_DIR) \
	--gff_dir $(SMALL_GFF_DIR) \
	--output_dir $(SMALL_PROCESSED_DIR) \
	--workers $(N_WORKERS) \
	--batch_size $(BATCH_SIZE)

# Process large dataset
process-large:
	@echo "Processing large dataset..."
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) src/pipeline.py \
	--fasta_dir $(FASTA_DIR) \
	--gff_dir $(GFF_DIR) \
	--output_dir $(PROCESSED_DIR) \
	--workers $(N_WORKERS) \
	--batch_size $(BATCH_SIZE)

# Rest of the Makefile remains the same...

# === Encoding Variables ===
ENCODED_DIR := $(DATA_DIR)/large_encoded
K_VALUE := 3
VOCAB_SIZE := 50
MAX_SEQUENCES := 1000

# Encoding setup
encode-setup:
	@echo "Setting up encoding directories..."
	@mkdir -p $(ENCODED_DIR)/ohe
	@mkdir -p $(ENCODED_DIR)/kmer
	@mkdir -p $(ENCODED_DIR)/bpe
	@mkdir -p run_logs

# Main encoding target
# Main encoding target
encode: encode-setup
	@echo "Running encoding..."
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) -m src.encoding.encode \
		--encoding_type $(encoding) \
		--processed_dir $(PROCESSED_DIR) \
		--output_dir $(ENCODED_DIR) \
		--workers $(N_WORKERS) \
		--batch_size $(BATCH_SIZE) \
		--vocab_size $(VOCAB_SIZE) \
		--max_sequences $(MAX_SEQUENCES) \
		$(if $(K_VALUE),--k $(K_VALUE),)


# Clean only encoded data
clean-encoded:
	@echo "Cleaning encoded data..."
	@rm -rf $(ENCODED_DIR)/*
	@echo "Cleaned encoded data"


# Add these targets
process-small:
	@echo "Processing small dataset..."
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) src/pipeline.py \
		--fasta_dir $(SMALL_FASTA_DIR) \
		--gff_dir $(SMALL_GFF_DIR) \
		--output_dir $(SMALL_PROCESSED_DIR) \
		--workers $(N_WORKERS) \
		--batch_size $(BATCH_SIZE)

encode-small: encode-setup
	@echo "Encoding small processed dataset..."
	PYTHONPATH=$(PYTHONPATH) $(PYTHON) -m src.encoding.encode \
		--encoding_type $(encoding) \
		--processed_dir $(SMALL_PROCESSED_DIR) \
		--output_dir $(SMALL_ENCODED_DIR) \
		--vocab_size $(VOCAB_SIZE) \
		--workers $(N_WORKERS) \
		--batch_size $(BATCH_SIZE)

clean-small:
	@echo "Cleaning small dataset results..."
	@rm -rf $(SMALL_PROCESSED_DIR)
	@rm -rf $(SMALL_ENCODED_DIR)
