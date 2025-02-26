#!/bin/bash

# Add required Conda channels
echo "Adding Conda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install dependencies
echo "Installing dependencies..."
mamba install -y --file requirements.txt

# Verify installation
echo "Verifying installed packages..."
conda list
