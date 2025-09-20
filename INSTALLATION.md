# PanGenomePlus Installation Guide

A comprehensive guide to installing PanGenomePlus and all its dependencies across different operating systems.

## Table of Contents

1. [System Requirements](#system-requirements)
2. [External Tool Installation](#external-tool-installation)
3. [PanGenomePlus Installation](#pangenomeplus-installation)
4. [Installation Verification](#installation-verification)
5. [Platform-Specific Instructions](#platform-specific-instructions)
6. [Advanced Installation Options](#advanced-installation-options)
7. [Troubleshooting](#troubleshooting)
8. [Alternative Installation Methods](#alternative-installation-methods)

## System Requirements

### Minimum Requirements
- **Python**: 3.9 or later
- **Memory**: 8 GB RAM (16 GB+ recommended)
- **Storage**: 10 GB available space
- **CPU**: Multi-core processor (4+ cores recommended)

### Recommended Specifications
- **Memory**: 32-128 GB RAM for large datasets (100+ genomes)
- **Storage**: SSD with 100+ GB for large analyses
- **CPU**: 8-16 cores for optimal performance
- **GPU**: NVIDIA GPU with CUDA support (optional, for acceleration)

### Operating System Support
- **Linux**: Ubuntu 18.04+, CentOS 7+, RHEL 7+
- **macOS**: 10.15 (Catalina) or later
- **Windows**: Windows 10/11 with WSL2

## External Tool Installation

PanGenomePlus requires several external bioinformatics tools. Here's how to install them:

### 1. MMseqs2 (Required)

MMseqs2 is the core clustering engine for PanGenomePlus.

#### Option A: Conda (Recommended)
```bash
# Install conda if not already installed
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Install MMseqs2
conda install -c conda-forge mmseqs2

# Verify installation
mmseqs version
```

#### Option B: Homebrew (macOS)
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install MMseqs2
brew install mmseqs2

# Verify installation
mmseqs version
```

#### Option C: From Source
```bash
# Download latest release
wget https://github.com/soedinglab/MMseqs2/releases/latest/download/mmseqs-linux-x86_64.tar.gz
tar xvzf mmseqs-linux-x86_64.tar.gz

# Add to PATH (add this to your ~/.bashrc or ~/.zshrc)
export PATH="$(pwd)/mmseqs/bin:$PATH"

# Verify installation
mmseqs version
```

### 2. Prodigal (Required)

Prodigal is used for gene prediction in bacterial genomes.

#### Option A: Conda
```bash
conda install -c bioconda prodigal
```

#### Option B: Homebrew (macOS)
```bash
brew install prodigal
```

#### Option C: From Source
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install prodigal

# CentOS/RHEL
sudo yum install prodigal

# From source
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal
make install
```

#### Verify Installation
```bash
prodigal -v
```

### 3. tRNAscan-SE (Optional)

For tRNA gene detection.

#### Option A: Conda
```bash
conda install -c bioconda trnascan-se
```

#### Option B: From Source
```bash
# Download from http://lowelab.ucsc.edu/tRNAscan-SE/
wget http://lowelab.ucsc.edu/software/tRNAscan-SE-2.0.9.tar.gz
tar -xzf tRNAscan-SE-2.0.9.tar.gz
cd tRNAscan-SE-2.0
./configure --prefix=/usr/local
make
sudo make install
```

#### Verify Installation
```bash
tRNAscan-SE -h
```

### 4. Barrnap (Optional)

For ribosomal RNA gene detection.

#### Option A: Conda
```bash
conda install -c bioconda barrnap
```

#### Option B: From Source
```bash
git clone https://github.com/tseemann/barrnap.git
cd barrnap
sudo cp bin/barrnap /usr/local/bin/
sudo cp -r db /usr/local/share/barrnap/
```

#### Verify Installation
```bash
barrnap --help
```

### 5. MINCED (Optional)

For CRISPR array detection.

#### Option A: Conda
```bash
conda install -c bioconda minced
```

#### Option B: From Source
```bash
git clone https://github.com/ctSkennerton/minced.git
cd minced
make
sudo cp minced /usr/local/bin/
```

#### Verify Installation
```bash
minced -h
```

## PanGenomePlus Installation

### Method 1: From GitHub (Recommended)

```bash
# Clone the repository
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus

# Create and activate virtual environment
python -m venv pangenome_env
source pangenome_env/bin/activate  # On Windows: pangenome_env\Scripts\activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install -r requirements.txt

# Install PanGenomePlus in development mode
pip install -e .
```

### Method 2: Direct Download

```bash
# Download and extract
wget https://github.com/mol-evol/pangenomeplus/archive/main.zip
unzip main.zip
cd pangenomeplus-main

# Create virtual environment and install
python -m venv pangenome_env
source pangenome_env/bin/activate
pip install -r requirements.txt
pip install -e .
```

### Method 3: Using Conda Environment

```bash
# Create conda environment with Python
conda create -n pangenomeplus python=3.9
conda activate pangenomeplus

# Clone and install
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus
pip install -r requirements.txt
pip install -e .
```

## Installation Verification

### Quick Test
```bash
# Activate your environment
source pangenome_env/bin/activate  # or conda activate pangenomeplus

# Test PanGenomePlus installation
python -m pangenomeplus --help

# Check version
python -c "import pangenomeplus; print('PanGenomePlus installed successfully')"
```

### External Tools Check
```bash
# Test all required tools
echo "Testing external tools..."

# MMseqs2 (required)
if command -v mmseqs &> /dev/null; then
    echo "[OK] MMseqs2 found: $(mmseqs version | head -1)"
else
    echo "[ERROR] MMseqs2 not found"
fi

# Prodigal (required)
if command -v prodigal &> /dev/null; then
    echo "[OK] Prodigal found: $(prodigal -v 2>&1 | head -1)"
else
    echo "[ERROR] Prodigal not found"
fi

# tRNAscan-SE (optional)
if command -v tRNAscan-SE &> /dev/null; then
    echo "[OK] tRNAscan-SE found"
else
    echo "[OPTIONAL] tRNAscan-SE not found (optional)"
fi

# Barrnap (optional)
if command -v barrnap &> /dev/null; then
    echo "[OK] Barrnap found"
else
    echo "[OPTIONAL] Barrnap not found (optional)"
fi

# MINCED (optional)
if command -v minced &> /dev/null; then
    echo "[OK] MINCED found"
else
    echo "[OPTIONAL] MINCED not found (optional)"
fi
```

### Full Functionality Test
```bash
# Create test directory
mkdir pangenomeplus_test
cd pangenomeplus_test

# Download test data (optional - you can use your own small genomes)
# For now, create a minimal test
python -m pangenomeplus --init-config test_config.yaml
echo "If this command works, installation is successful!"
```

## Platform-Specific Instructions

### Ubuntu/Debian Linux

```bash
# Update package manager
sudo apt-get update

# Install system dependencies
sudo apt-get install -y python3 python3-pip python3-venv git wget curl

# Install build tools (needed for some packages)
sudo apt-get install -y build-essential

# Install external tools via apt (alternative to conda)
sudo apt-get install -y prodigal

# Follow general installation instructions above
```

### CentOS/RHEL Linux

```bash
# Update package manager
sudo yum update

# Install system dependencies
sudo yum install -y python3 python3-pip git wget curl

# Install development tools
sudo yum groupinstall -y "Development Tools"

# Install EPEL for additional packages
sudo yum install -y epel-release

# Follow general installation instructions above
```

### macOS

```bash
# Install Xcode command line tools
xcode-select --install

# Install Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Python and tools via Homebrew
brew install python git

# Install external bioinformatics tools
brew install mmseqs2 prodigal

# Follow general installation instructions above
```

### Windows (WSL2)

```bash
# Install WSL2 with Ubuntu
# Follow Microsoft's WSL2 installation guide

# In WSL2 Ubuntu terminal:
sudo apt-get update
sudo apt-get install -y python3 python3-pip python3-venv git wget curl build-essential

# Follow Ubuntu instructions above
```

## Advanced Installation Options

### HPC/Cluster Environment

For high-performance computing environments:

```bash
# Load required modules (example for SLURM)
module load python/3.9
module load gcc/9.3.0

# Create environment in your home directory
python -m venv ~/pangenomeplus_env
source ~/pangenomeplus_env/bin/activate

# Install with specific optimizations
pip install --no-cache-dir -r requirements.txt
pip install -e .

# Install external tools in user space
mkdir -p ~/local/bin
export PATH="$HOME/local/bin:$PATH"

# Install MMseqs2 to user directory
wget https://github.com/soedinglab/MMseqs2/releases/latest/download/mmseqs-linux-x86_64.tar.gz
tar xvzf mmseqs-linux-x86_64.tar.gz
cp mmseqs/bin/* ~/local/bin/
```

### Docker Installation

```dockerfile
# Dockerfile for PanGenomePlus
FROM ubuntu:20.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git \
    wget \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install external tools
RUN wget -O mmseqs.tar.gz https://github.com/soedinglab/MMseqs2/releases/latest/download/mmseqs-linux-x86_64.tar.gz \
    && tar xzf mmseqs.tar.gz \
    && cp mmseqs/bin/* /usr/local/bin/ \
    && rm -rf mmseqs*

# Install Prodigal
RUN apt-get update && apt-get install -y prodigal

# Install PanGenomePlus
WORKDIR /app
COPY . .
RUN pip install -r requirements.txt && pip install -e .

# Set entrypoint
ENTRYPOINT ["python", "-m", "pangenomeplus"]
```

Build and run:
```bash
docker build -t pangenomeplus .
docker run -v /path/to/genomes:/data pangenomeplus /data output/
```

### Singularity/Apptainer

```bash
# Create Singularity definition file
cat > pangenomeplus.def << EOF
Bootstrap: docker
From: ubuntu:20.04

%post
    apt-get update
    apt-get install -y python3 python3-pip git wget curl build-essential prodigal

    # Install MMseqs2
    wget -O mmseqs.tar.gz https://github.com/soedinglab/MMseqs2/releases/latest/download/mmseqs-linux-x86_64.tar.gz
    tar xzf mmseqs.tar.gz
    cp mmseqs/bin/* /usr/local/bin/

    # Install PanGenomePlus
    git clone https://github.com/mol-evol/pangenomeplus.git
    cd pangenomeplus
    pip install -r requirements.txt
    pip install -e .

%runscript
    python -m pangenomeplus "$@"
EOF

# Build container
singularity build pangenomeplus.sif pangenomeplus.def

# Run
singularity run pangenomeplus.sif genomes/ output/
```

## Troubleshooting

### Common Issues

#### 1. Python Version Issues
```bash
# Check Python version
python --version
python3 --version

# If Python 3.9+ not available, install via pyenv
curl https://pyenv.run | bash
pyenv install 3.9.16
pyenv global 3.9.16
```

#### 2. Permission Errors
```bash
# If getting permission errors, use user installation
pip install --user -r requirements.txt
pip install --user -e .

# Or fix ownership
sudo chown -R $USER:$USER ~/.local/
```

#### 3. Tool Not Found Errors
```bash
# Check PATH
echo $PATH

# Find where tools are installed
which mmseqs
which prodigal

# Add to PATH if needed (add to ~/.bashrc)
export PATH="/path/to/tool/bin:$PATH"
source ~/.bashrc
```

#### 4. Memory Issues During Installation
```bash
# Install with no cache to reduce memory usage
pip install --no-cache-dir -r requirements.txt

# Install dependencies one by one
pip install numpy
pip install pandas
pip install biopython
# etc.
```

#### 5. Network/Firewall Issues
```bash
# If behind corporate firewall, configure pip
pip install --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host files.pythonhosted.org -r requirements.txt

# Or download packages manually
pip download -r requirements.txt
pip install --no-index --find-links . -r requirements.txt
```

### Platform-Specific Troubleshooting

#### macOS Issues

**Issue**: "Developer cannot be verified" error
```bash
# For downloaded binaries, remove quarantine
xattr -d com.apple.quarantine /path/to/binary
```

**Issue**: Command line tools not found
```bash
# Reinstall Xcode command line tools
sudo xcode-select --reset
xcode-select --install
```

#### Linux Issues

**Issue**: "No module named '_ssl'" error
```bash
# Install SSL development libraries
sudo apt-get install libssl-dev  # Ubuntu/Debian
sudo yum install openssl-devel   # CentOS/RHEL

# Rebuild Python if using pyenv
pyenv install 3.9.16
```

#### Windows WSL Issues

**Issue**: Performance is slow
```bash
# Move project to WSL filesystem (not Windows filesystem)
cp -r /mnt/c/pangenomeplus ~/pangenomeplus
cd ~/pangenomeplus
```

## Alternative Installation Methods

### Using Poetry
```bash
# Install Poetry
curl -sSL https://install.python-poetry.org | python3 -

# Clone and install
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus
poetry install
poetry shell
```

### Using Pipenv
```bash
# Install pipenv
pip install pipenv

# Clone and install
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus
pipenv install
pipenv shell
```

### System-wide Installation (Not Recommended)
```bash
# Only if you understand the implications
sudo pip install -e .
```

## Post-Installation Setup

### Environment Variables
Add these to your shell configuration file (~/.bashrc, ~/.zshrc):

```bash
# PanGenomePlus environment
export PANGENOMEPLUS_HOME="/path/to/pangenomeplus"
export PATH="$PANGENOMEPLUS_HOME/bin:$PATH"

# Tool paths (if needed)
export MMSEQS_PATH="/path/to/mmseqs"
export PRODIGAL_PATH="/path/to/prodigal"
```

### Shell Completion (Optional)
```bash
# Add tab completion for PanGenomePlus commands
python -m pangenomeplus --install-completion
```

### Performance Tuning
```bash
# Set optimal number of threads (usually = number of CPU cores)
export OMP_NUM_THREADS=8

# For large datasets, increase memory limits
ulimit -v unlimited
```

---

## Quick Installation Summary

For most users, this one-liner approach will work:

```bash
# Complete installation (choose one method)

# Method 1: Conda-based
conda create -n pangenomeplus python=3.9
conda activate pangenomeplus
conda install -c conda-forge mmseqs2
conda install -c bioconda prodigal trnascan-se barrnap minced
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus
pip install -r requirements.txt
pip install -e .

# Method 2: Homebrew-based (macOS)
brew install python mmseqs2 prodigal
git clone https://github.com/mol-evol/pangenomeplus.git
cd pangenomeplus
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

**Need help?** Check the [TUTORIAL.md](TUTORIAL.md) for step-by-step usage guide or open an issue on [GitHub](https://github.com/mol-evol/pangenomeplus/issues).