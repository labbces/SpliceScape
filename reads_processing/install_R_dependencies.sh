echo "Installing dependencies..."
conda install -y --file requirements.txt

# Verify installation
echo "Verifying installed packages..."
conda list

echo "All dependencies installed successfully!"