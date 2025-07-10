#!/usr/bin/env python3
"""
Robust Bowtie2 log parser with proper error handling and extensibility.
Usage: python parse_bowtie_logs.py log1.log log2.log ... [--format json|csv|tsv]
"""

import sys
import re
import json
import csv
import argparse
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import List, Optional, Dict, Any
from enum import Enum


class OutputFormat(Enum):
    TSV = "tsv"
    CSV = "csv" 
    JSON = "json"


@dataclass
class AlignmentStats:
    """Type-safe container for alignment statistics."""
    sample: str
    total_reads: int
    concordant_0_count: int
    concordant_0_pct: float
    concordant_1_count: int
    concordant_1_pct: float
    concordant_multi_count: int
    concordant_multi_pct: float
    overall_alignment_rate: float
    
    def validate(self) -> None:
        """Validate that the statistics make sense."""
        if self.total_reads <= 0:
            raise ValueError(f"Invalid total reads: {self.total_reads}")
        
        # Check that counts sum to total
        total_concordant = self.concordant_0_count + self.concordant_1_count + self.concordant_multi_count
        if total_concordant != self.total_reads:
            raise ValueError(f"Concordant counts ({total_concordant}) don't sum to total reads ({self.total_reads})")
        
        # Check that percentages are reasonable
        total_pct = self.concordant_0_pct + self.concordant_1_pct + self.concordant_multi_pct
        if not (99.0 <= total_pct <= 101.0):  # Allow for rounding errors
            raise ValueError(f"Concordant percentages sum to {total_pct}%, expected ~100%")
        
        if not (0.0 <= self.overall_alignment_rate <= 100.0):
            raise ValueError(f"Invalid overall alignment rate: {self.overall_alignment_rate}%")


class BowtieLogParser:
    """Robust parser for Bowtie2 alignment logs."""
    
    # Compiled regex patterns for efficiency
    PATTERNS = {
        'total_reads': re.compile(r'(\d+) reads; of these:'),
        'concordant_0': re.compile(r'(\d+) \(([0-9.]+)%\) aligned concordantly 0 times'),
        'concordant_1': re.compile(r'(\d+) \(([0-9.]+)%\) aligned concordantly exactly 1 time'),
        'concordant_multi': re.compile(r'(\d+) \(([0-9.]+)%\) aligned concordantly >1 times'),
        'overall_rate': re.compile(r'([0-9.]+)% overall alignment rate')
    }
    
    @staticmethod
    def extract_sample_name(log_path: Path) -> str:
        """Extract sample name from log file path."""
        return log_path.stem.replace('.bowtie', '')
    
    @classmethod
    def parse_log_content(cls, content: str) -> Dict[str, Any]:
        """Extract statistics from log file content."""
        results = {}
        
        # Extract total reads
        if match := cls.PATTERNS['total_reads'].search(content):
            results['total_reads'] = int(match.group(1))
        else:
            raise ValueError("Could not find total reads count")
        
        # Extract concordant alignment stats
        if match := cls.PATTERNS['concordant_0'].search(content):
            results['concordant_0_count'] = int(match.group(1))
            results['concordant_0_pct'] = float(match.group(2))
        else:
            raise ValueError("Could not find concordant 0 times alignment stats")
        
        if match := cls.PATTERNS['concordant_1'].search(content):
            results['concordant_1_count'] = int(match.group(1))
            results['concordant_1_pct'] = float(match.group(2))
        else:
            raise ValueError("Could not find concordant 1 time alignment stats")
        
        if match := cls.PATTERNS['concordant_multi'].search(content):
            results['concordant_multi_count'] = int(match.group(1))
            results['concordant_multi_pct'] = float(match.group(2))
        else:
            raise ValueError("Could not find concordant >1 times alignment stats")
        
        if match := cls.PATTERNS['overall_rate'].search(content):
            results['overall_alignment_rate'] = float(match.group(1))
        else:
            raise ValueError("Could not find overall alignment rate")
        
        return results
    
    @classmethod
    def parse_log_file(cls, log_path: Path) -> AlignmentStats:
        """Parse a single Bowtie2 log file."""
        if not log_path.exists():
            raise FileNotFoundError(f"Log file not found: {log_path}")
        
        if not log_path.is_file():
            raise ValueError(f"Path is not a file: {log_path}")
        
        try:
            content = log_path.read_text()
        except UnicodeDecodeError as e:
            raise ValueError(f"Could not decode log file {log_path}: {e}")
        
        if not content.strip():
            raise ValueError(f"Log file is empty: {log_path}")
        
        sample_name = cls.extract_sample_name(log_path)
        extracted_data = cls.parse_log_content(content)
        
        stats = AlignmentStats(
            sample=sample_name,
            **extracted_data
        )
        
        # Validate the parsed data
        stats.validate()
        
        return stats


class OutputFormatter:
    """Handle different output formats for alignment statistics."""
    
    @staticmethod
    def to_tsv(stats_list: List[AlignmentStats]) -> str:
        """Format statistics as TSV."""
        if not stats_list:
            return ""
        
        # Get headers from first object
        headers = list(asdict(stats_list[0]).keys())
        lines = ['\t'.join(headers)]
        
        for stats in stats_list:
            values = [str(v) for v in asdict(stats).values()]
            lines.append('\t'.join(values))
        
        return '\n'.join(lines)
    
    @staticmethod
    def to_csv(stats_list: List[AlignmentStats]) -> str:
        """Format statistics as CSV."""
        if not stats_list:
            return ""
        
        import io
        output = io.StringIO()
        
        headers = list(asdict(stats_list[0]).keys())
        writer = csv.DictWriter(output, fieldnames=headers)
        writer.writeheader()
        
        for stats in stats_list:
            writer.writerow(asdict(stats))
        
        return output.getvalue()
    
    @staticmethod
    def to_json(stats_list: List[AlignmentStats]) -> str:
        """Format statistics as JSON."""
        return json.dumps([asdict(stats) for stats in stats_list], indent=2)


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse Bowtie2 log files and extract alignment statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        'log_files',
        nargs='+',
        type=Path,
        help="Bowtie2 log files to parse"
    )
    
    parser.add_argument(
        '--format',
        type=OutputFormat,
        choices=list(OutputFormat),
        default=OutputFormat.TSV,
        help="Output format (default: tsv)"
    )
    
    parser.add_argument(
        '--strict',
        action='store_true',
        help="Exit on first parsing error (default: skip bad files with warnings)"
    )
    
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_arguments()
    
    successful_stats = []
    failed_files = []
    
    for log_file in args.log_files:
        try:
            stats = BowtieLogParser.parse_log_file(log_file)
            successful_stats.append(stats)
        except Exception as e:
            error_msg = f"Failed to parse {log_file}: {e}"
            failed_files.append((log_file, str(e)))
            
            if args.strict:
                print(f"Error: {error_msg}", file=sys.stderr)
                sys.exit(1)
            else:
                print(f"Warning: {error_msg}", file=sys.stderr)
    
    if not successful_stats:
        print("Error: No log files were successfully parsed", file=sys.stderr)
        sys.exit(1)
    
    # Output results
    formatter = OutputFormatter()
    if args.format == OutputFormat.TSV:
        output = formatter.to_tsv(successful_stats)
    elif args.format == OutputFormat.CSV:
        output = formatter.to_csv(successful_stats)
    elif args.format == OutputFormat.JSON:
        output = formatter.to_json(successful_stats)
    
    print(output)
    
    # Report summary
    if failed_files:
        print(f"\nSummary: {len(successful_stats)} files parsed successfully, "
              f"{len(failed_files)} files failed", file=sys.stderr)


if __name__ == "__main__":
    main()
