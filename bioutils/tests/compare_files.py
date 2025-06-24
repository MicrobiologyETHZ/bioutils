#!/usr/bin/env python3

import os
import sys
import hashlib
import argparse
from collections import defaultdict


def calculate_md5(filepath):
    """Calculate MD5 hash of a file."""
    hash_md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        # Read in chunks to handle large files efficiently
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def find_identical_files(dir1, dir2, output_format='table'):
    """Find identical files between two directories."""
    # Create dictionaries to store file hashes
    dir1_hashes = defaultdict(list)
    dir2_hashes = defaultdict(list)

    # Calculate hashes for files in directory 1
    print(f"Calculating MD5 checksums for files in {dir1}...")
    for root, _, files in os.walk(dir1):
        for filename in files:
            filepath = os.path.join(root, filename)
            try:
                file_hash = calculate_md5(filepath)
                dir1_hashes[file_hash].append(filepath)
            except (IOError, PermissionError) as e:
                print(f"Error processing {filepath}: {e}", file=sys.stderr)

    # Calculate hashes for files in directory 2
    print(f"Calculating MD5 checksums for files in {dir2}...")
    for root, _, files in os.walk(dir2):
        for filename in files:
            filepath = os.path.join(root, filename)
            try:
                file_hash = calculate_md5(filepath)
                dir2_hashes[file_hash].append(filepath)
            except (IOError, PermissionError) as e:
                print(f"Error processing {filepath}: {e}", file=sys.stderr)

    # Find and output identical files
    print("Identifying identical files...")

    if output_format == 'table':
        print(f"{'File in ' + dir1:<50}\t{'Identical file in ' + dir2:<50}")
        print(f"{'-' * 50}\t{'-' * 50}")
    elif output_format == 'csv':
        print("file1,file2")

    # Count matches for summary
    total_matches = 0

    # For each hash in directory 1
    for file_hash in dir1_hashes:
        # If the same hash exists in directory 2
        if file_hash in dir2_hashes:
            # For each file with this hash in directory 1
            for file1 in dir1_hashes[file_hash]:
                # For each file with this hash in directory 2
                for file2 in dir2_hashes[file_hash]:
                    total_matches += 1
                    if output_format == 'table':
                        print(f"{file1:<50}\t{file2:<50}")
                    elif output_format == 'csv':
                        print(f"\"{file1}\",\"{file2}\"")

    print(f"\nFound {total_matches} identical file pairs.")
    return total_matches


def main():
    parser = argparse.ArgumentParser(
        description='Find identical files between two directories using MD5 checksums.')
    parser.add_argument('dir1', help='First directory to compare')
    parser.add_argument('dir2', help='Second directory to compare')
    parser.add_argument('--format', choices=['table', 'csv'], default='table',
                        help='Output format (default: table)')
    parser.add_argument('--output', '-o', help='Output file (default: stdout)')

    args = parser.parse_args()

    # Validate directories
    if not os.path.isdir(args.dir1):
        print(f"Error: {args.dir1} is not a directory", file=sys.stderr)
        return 1

    if not os.path.isdir(args.dir2):
        print(f"Error: {args.dir2} is not a directory", file=sys.stderr)
        return 1

    # Redirect output if specified
    if args.output:
        sys.stdout = open(args.output, 'w')

    try:
        find_identical_files(args.dir1, args.dir2, args.format)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.", file=sys.stderr)
        return 1
    finally:
        if args.output:
            sys.stdout.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
