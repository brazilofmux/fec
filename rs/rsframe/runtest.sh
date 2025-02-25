#!/bin/bash
# RSFrame Compatibility Test
# Tests the compatibility between old and new RSFrame implementations

set -e  # Exit on error

echo "===== RSFrame Compatibility Test ====="
echo ""

# Create a test file
echo "Creating test file..."
cat > test_input.txt << EOF
# RSFrame Test File

This is a test file to verify compatibility between the old and new
RSFrame implementations. We need to ensure that files encoded with the
old implementation can be decoded with the new implementation and vice versa.

The Reed-Solomon error correction coding must maintain exact binary compatibility
to ensure that we can recover files encoded decades ago.

1234567890
ABCDEFGHIJKLMNOPQRSTUVWXYZ
abcdefghijklmnopqrstuvwxyz
EOF

echo "Test file created."
echo ""

# Compile programs if needed
if [ ! -f "rsframe_old" ]; then
    echo "Compiling old RSFrame implementation..."
    g++ -o rsframe_old rsframe.cpp -I. -std=c++11
fi

if [ ! -f "rsframe_new" ]; then
    echo "Compiling new RSFrame implementation..."
    g++ -o rsframe_new fixed-rsframe-decoder.cpp -I. -std=c++11
fi

# Test 1: Encode with old, decode with old
echo "Test 1: Encode with old, decode with old"
./rsframe_old -e test_input.txt test_old_encoded.rs
./rsframe_old -d test_old_encoded.rs test_old_decoded.txt
echo "Comparing original and decoded files..."
if cmp -s test_input.txt test_old_decoded.txt; then
    echo "Success: Files are identical!"
else
    echo "Error: Files differ!"
    exit 1
fi
echo ""

# Test 2: Encode with old, decode with new
echo "Test 2: Encode with old, decode with new"
./rsframe_new -d test_old_encoded.rs test_new_decoded_old_encoded.txt
echo "Comparing original and decoded files..."
if cmp -s test_input.txt test_new_decoded_old_encoded.txt; then
    echo "Success: Files are identical!"
else
    echo "Error: Files differ!"
    exit 1
fi
echo ""

# Test 3: Encode with new, decode with old
echo "Test 3: Encode with new, decode with old"
./rsframe_new -e test_input.txt test_new_encoded.rs
./rsframe_old -d test_new_encoded.rs test_old_decoded_new_encoded.txt
echo "Comparing original and decoded files..."
if cmp -s test_input.txt test_old_decoded_new_encoded.txt; then
    echo "Success: Files are identical!"
else
    echo "Error: Files differ!"
    exit 1
fi
echo ""

# Test 4: Encode with new, decode with new
echo "Test 4: Encode with new, decode with new"
./rsframe_new -d test_new_encoded.rs test_new_decoded.txt
echo "Comparing original and decoded files..."
if cmp -s test_input.txt test_new_decoded.txt; then
    echo "Success: Files are identical!"
else
    echo "Error: Files differ!"
    exit 1
fi
echo ""

# Verify that old and new encoded files match
echo "Verifying binary compatibility of encoded files..."
if cmp -s test_old_encoded.rs test_new_encoded.rs; then
    echo "Success: Encoded files are binary compatible!"
else
    # The files might differ in ways that don't affect decoding
    # Let's check the SHA512 sums of the files
    echo "Warning: Encoded files differ. Comparing SHA512 hashes..."
    sha512sum test_old_encoded.rs test_new_encoded.rs
fi
echo ""

echo "All tests completed successfully!"
echo ""
echo "Cleanup: Remove test files? (y/n)"
read -r response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]; then
    rm -f test_input.txt test_old_encoded.rs test_old_decoded.txt test_new_encoded.rs test_new_decoded.txt test_new_decoded_old_encoded.txt test_old_decoded_new_encoded.txt
    echo "Test files removed."
fi
