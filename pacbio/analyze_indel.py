#!/usr/bin/env python3
import sys

def robust_find_deletions(coord_file, min_len=500, max_len=200000, merge_gap=500):
    """
    Core analysis function
    coord_file: Path to input .coords file
    min_len: Minimum length for candidate deletions (bp)
    max_len: Maximum length for candidate deletions (bp)
    merge_gap: Maximum gap between adjacent blocks to merge (bp). Gaps smaller than this will be merged.
    """
    # 1. Read and parse data
    with open(coord_file, 'r') as f:
        lines = f.readlines()
    
    # Find data start line (skip file paths, empty lines, and headers)
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() and line[0].isdigit():  # Data lines start with numbers
            data_start = i
            break
    
    if data_start == 0:
        # Alternative: try common header line counts
        data_start = 5
    
    # 2. Data structure initialization
    chromosome_data = {}  # {chr_name: {'length': int, 'blocks': [(start, end), ...]}}
    
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith('='):
            continue
        
        parts = line.split()
        if len(parts) < 18:
            continue
        
        try:
            # Parse key fields
            ref_start = int(parts[0])
            ref_end = int(parts[1])
            ref_chr = parts[-2]  # Chromosome identifier
            
            # Chromosome length ([LEN R] column, index needs verification based on your file)
            # Common position is parts[12], but let's find it dynamically
            chr_len = None
            for i in range(10, 15):
                if i < len(parts) and parts[i].isdigit():
                    chr_len = int(parts[i])
                    break
            
            if chr_len is None:
                continue  # Unable to get length, skip
            
            # Initialize or update chromosome data
            if ref_chr not in chromosome_data:
                chromosome_data[ref_chr] = {
                    'length': chr_len,
                    'blocks': []
                }
            else:
                # Ensure we use the correct length (should be consistent)
                if chromosome_data[ref_chr]['length'] != chr_len:
                    print(f"Warning: Inconsistent length for chromosome {ref_chr}: "
                          f"{chromosome_data[ref_chr]['length']} vs {chr_len}")
            
            # Add alignment block
            chromosome_data[ref_chr]['blocks'].append((min(ref_start, ref_end), 
                                                      max(ref_start, ref_end)))
            
        except (ValueError, IndexError) as e:
            # Silently skip parsing errors
            continue
    
    print(f"Successfully parsed {len(chromosome_data)} chromosomes")
    
    # 3. Intelligently merge alignment blocks for each chromosome and find gaps
    all_candidates = []
    
    for chr_name, data in chromosome_data.items():
        chr_len = data['length']
        blocks = data['blocks']
        
        if not blocks:
            # Entire chromosome is uncovered
            if chr_len >= min_len:
                all_candidates.append((chr_name, 1, chr_len, chr_len, "ENTIRE_CHROMOSOME"))
            continue
        
        # Sort and merge blocks
        blocks.sort(key=lambda x: x[0])
        merged_blocks = []
        
        for start, end in blocks:
            if not merged_blocks:
                merged_blocks.append([start, end])
                continue
            
            last_start, last_end = merged_blocks[-1]
            
            # If current block overlaps or is very close to previous block, merge them
            if start <= last_end + merge_gap:
                merged_blocks[-1][1] = max(last_end, end)
            else:
                merged_blocks.append([start, end])
        
        # 4. Find gaps (candidate deletions)
        # Check beginning of chromosome
        current_pos = 1
        for i, (block_start, block_end) in enumerate(merged_blocks):
            # Gap before block start
            if block_start > current_pos:
                gap_len = block_start - current_pos
                if min_len <= gap_len <= max_len:
                    # Determine position type
                    pos_type = classify_position(current_pos, block_start-1, chr_len)
                    all_candidates.append((chr_name, current_pos, block_start-1, 
                                         gap_len, pos_type))
            
            current_pos = max(current_pos, block_end + 1)
        
        # Check end of chromosome
        if current_pos <= chr_len:
            gap_len = chr_len - current_pos + 1
            if min_len <= gap_len <= max_len:
                pos_type = classify_position(current_pos, chr_len, chr_len)
                all_candidates.append((chr_name, current_pos, chr_len, 
                                     gap_len, pos_type))
    
    return all_candidates, chromosome_data

def classify_position(start, end, chr_len):
    """Simple classification of gap positions"""
    midpoint = (start + end) / 2
    relative_pos = midpoint / chr_len
    
    if relative_pos < 0.1:
        return "TELOMERE_PROXIMAL"
    elif relative_pos > 0.9:
        return "TELOMERE_PROXIMAL"
    elif 0.4 < relative_pos < 0.6:
        return "CENTROMERE_REGION"
    elif 0.2 < relative_pos < 0.8:
        return "CHROMOSOME_ARM"
    else:
        return "INTERMEDIATE"

def analyze_and_report(candidates, chromosome_data, output_file):
    """Analyze results and generate report"""
    
    if not candidates:
        print("No qualified candidate deletion fragments found.")
        return
    
    # Sort by length
    sorted_by_len = sorted(candidates, key=lambda x: x[3], reverse=True)
    
    # Sort by position (chromosome + start position)
    sorted_by_pos = sorted(candidates, key=lambda x: (x[0], x[1]))
    
    # Statistics
    total_candidates = len(candidates)
    total_bp = sum(c[3] for c in candidates)
    avg_len = total_bp / total_candidates
    
    # Statistics by type
    type_stats = {}
    for _, _, _, length, pos_type in candidates:
        type_stats[pos_type] = type_stats.get(pos_type, 0) + 1
    
    # Statistics by chromosome
    chr_stats = {}
    for chr_name, start, end, length, pos_type in candidates:
        if chr_name not in chr_stats:
            chr_stats[chr_name] = {'count': 0, 'total_len': 0, 'types': {}}
        chr_stats[chr_name]['count'] += 1
        chr_stats[chr_name]['total_len'] += length
        chr_stats[chr_name]['types'][pos_type] = chr_stats[chr_name]['types'].get(pos_type, 0) + 1
    
    # Output report
    with open(output_file, 'w') as f:
        f.write("# Designed deletion candidate fragment analysis report (based on sensitive parameter alignment)\n")
        f.write(f"# Input file: SY14_repeat_sensitive.coords\n")
        f.write(f"# Total: {total_candidates} fragments, {total_bp:,} bp\n")
        f.write(f"# Average length: {avg_len:,.1f} bp\n\n")
        
        f.write("#CHROM\tSTART\tEND\tLENGTH(bp)\tTYPE\tCHR_POSITION\n")
        for chr_name, start, end, length, pos_type in sorted_by_pos:
            chr_len = chromosome_data.get(chr_name, {}).get('length', 'unknown')
            rel_pos = f"{(start+end)/2/chr_len*100:.1f}%" if isinstance(chr_len, (int, float)) else "unknown"
            f.write(f"{chr_name}\t{start}\t{end}\t{length}\t{pos_type}\t{rel_pos}\n")
    
    # Console output statistics
    print(f"\n{'='*60}")
    print("Designed deletion candidate fragment analysis report (based on sensitive parameter alignment)")
    print(f"{'='*60}")
    print(f"Total candidate fragments: {total_candidates}")
    print(f"Total missing base pairs: {total_bp:,} bp")
    print(f"Average length: {avg_len:,.1f} bp")
    print(f"Length range: {min(c[3] for c in candidates):,} - {max(c[3] for c in candidates):,} bp")
    
    print(f"\nDistribution by position type:")
    for pos_type, count in sorted(type_stats.items()):
        percentage = count / total_candidates * 100
        print(f"  {pos_type:<20} {count:>4} fragments ({percentage:>5.1f}%)")
    
    print(f"\nDistribution by chromosome:")
    for chr_name in sorted(chr_stats.keys()):
        stats = chr_stats[chr_name]
        chr_len = chromosome_data.get(chr_name, {}).get('length', 'unknown')
        len_str = f"{chr_len:,}" if isinstance(chr_len, (int, float)) else "unknown"
        print(f"  {chr_name:<15} {stats['count']:>3} fragments, "
              f"total {stats['total_len']:>9,} bp (chromosome length: {len_str} bp)")
    
    print(f"\nTop 15 longest candidate fragments:")
    for i, (chr_name, start, end, length, pos_type) in enumerate(sorted_by_len[:15]):
        print(f"  {i+1:2d}. {chr_name:<15} {start:>8,}-{end:>8,} "
              f"({length:>6,} bp) [{pos_type}]")
    
    print(f"\nDetailed results saved to: {output_file}")
    
    # Provide adjustment suggestions
    print(f"\n{'='*60}")
    print("Parameter adjustment suggestions:")
    if total_candidates < 50:
        print("  • Few candidate fragments, you can try:")
        print(f"    - Lower the minimum length threshold (current: {MIN_LEN})")
        print(f"    - Increase block merge gap (current: {MERGE_GAP})")
    elif total_candidates > 80:
        print("  • Too many candidate fragments, you can try:")
        print(f"    - Increase minimum length threshold to filter small fragments (current: {MIN_LEN})")
        print("    - Check for false positives (e.g., gaps caused by alignment errors)")
    
    # Biological plausibility check prompts
    print(f"\nBiological validation suggestions:")
    print("  1. Check if CENTROMERE_REGION type fragments are indeed near the centromere of corresponding chromosome")
    print("  2. Check if TELOMERE_PROXIMAL type fragments are located at chromosome ends")
    print("  3. Compare with previous standard parameter results to see if sensitive parameters captured more/larger fragments")

def main():
    # ============== Adjustable parameter section ==============
    COORD_FILE = "SY14_repeat_sensitive.coords"  # Input file
    OUTPUT_FILE = "design_deletion_candidates_sensitive.txt"  # Output file
    
    MIN_LEN = 500      # Minimum length (bp)
    MAX_LEN = 200000   # Maximum length (bp) - Significantly increased to capture large fragment deletions
    MERGE_GAP = 500    # Maximum gap to merge blocks (bp)
    # ============================================
    
    print("Starting analysis of designed deletion candidate fragments (sensitive parameter version)...")
    print(f"Input file: {COORD_FILE}")
    print(f"Parameter settings: min_length={MIN_LEN}bp, max_length={MAX_LEN:,}bp, merge_gap={MERGE_GAP}bp")
    print("-" * 60)
    
    try:
        # Run core analysis
        candidates, chromosome_data = robust_find_deletions(
            COORD_FILE, MIN_LEN, MAX_LEN, MERGE_GAP
        )
        
        # Generate report
        analyze_and_report(candidates, chromosome_data, OUTPUT_FILE)
        
        # Additional analysis: length distribution histogram data
        if candidates:
            print(f"\nLength distribution (for histogram creation):")
            lengths = [c[3] for c in candidates]
            bins = [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000]
            hist = {}
            for length in lengths:
                for i in range(len(bins)-1):
                    if bins[i] <= length < bins[i+1]:
                        bin_label = f"{bins[i]:,}-{bins[i+1]:,}"
                        hist[bin_label] = hist.get(bin_label, 0) + 1
                        break
                else:
                    if length >= bins[-1]:
                        bin_label = f">={bins[-1]:,}"
                        hist[bin_label] = hist.get(bin_label, 0) + 1
            
            for bin_label, count in sorted(hist.items()):
                percentage = count / len(lengths) * 100
                print(f"  {bin_label:<15} bp: {count:>3} fragments ({percentage:>5.1f}%)")
        
    except FileNotFoundError:
        print(f"Error: File not found {COORD_FILE}")
        print("Please ensure the file exists, or modify the COORD_FILE variable to the correct path.")
        print("Recommended command to generate this file:")
        print("  nucmer --mum -l 100 -c 500 ./01_reference_genome/BY4742_S288C_R64-2-1.fasta ./04_results/flye_assembly_full/assembly.fasta -p SY14_repeat_sensitive")
        print("  delta-filter -1 SY14_repeat_sensitive.delta > SY14_repeat_sensitive.filtered.delta")
        print("  show-coords -r -c -l SY14_repeat_sensitive.filtered.delta > SY14_repeat_sensitive.coords")
    except Exception as e:
        print(f"Error occurred during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
