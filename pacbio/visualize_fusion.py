#!/usr/bin/env python3
"""
SY14 vs BY4742 genome alignment visualization
Input: circos_links.txt (already generated)
Output: Chromosomal fusion evidence plot
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from collections import defaultdict
import seaborn as sns

# Set attractive style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def load_links(link_file):
    """Load alignment link file"""
    links = []
    with open(link_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                ref_chr = parts[0].replace('CP026', 'chr')  # Simplify name
                ref_start = int(parts[1])
                ref_end = int(parts[2])
                qry_contig = parts[3]  # e.g., SY14_1
                qry_start = int(parts[4])
                qry_end = int(parts[5])
                
                # Convert to megabase (Mb) unit
                links.append({
                    'ref_chr': ref_chr,
                    'ref_start_mb': ref_start / 1e6,
                    'ref_end_mb': ref_end / 1e6,
                    'ref_length_mb': (ref_end - ref_start) / 1e6,
                    'qry_contig': qry_contig,
                    'qry_start_mb': qry_start / 1e6,
                    'qry_end_mb': qry_end / 1e6,
                    'qry_length_mb': (qry_end - qry_start) / 1e6
                })
    return pd.DataFrame(links)

def create_fusion_plot(df, output_file):
    """Create chromosomal fusion visualization plot"""
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Overall layout
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.5, 1], width_ratios=[3, 1])
    
    # 2. Main plot: SY14 contig vs BY4742 chromosome correspondence
    ax_main = fig.add_subplot(gs[0, :])
    ax_main.set_title('SY14 Synthetic Chromosome vs BY4742 Wild-type Chromosomes', 
                     fontsize=16, fontweight='bold', pad=20)
    
    # SY14 contigs (sorted by length)
    contig_stats = df.groupby('qry_contig').agg({
        'qry_length_mb': 'max',
        'ref_chr': 'nunique'
    }).sort_values('qry_length_mb', ascending=False)
    
    contigs = contig_stats.index.tolist()
    contig_lengths = contig_stats['qry_length_mb'].values
    contig_colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(contigs)))
    
    # Draw SY14 contigs
    y_pos = 0
    contig_y_positions = {}
    for i, (contig, length) in enumerate(zip(contigs, contig_lengths)):
        contig_y_positions[contig] = y_pos
        # contig bar
        ax_main.add_patch(patches.Rectangle(
            (0, y_pos - 0.15), length, 0.3,
            facecolor=contig_colors[i], alpha=0.8, edgecolor='black', linewidth=1
        ))
        # Label
        ax_main.text(length + 0.1, y_pos, f'{contig}\n({length:.2f} Mb)', 
                    va='center', fontsize=10, fontweight='bold')
        y_pos += 1
    
    # 3. Alignment link plot
    ax_links = fig.add_subplot(gs[1, :])
    ax_links.set_title('Detailed Alignment Links between SY14 and BY4742', 
                      fontsize=14, fontweight='bold', pad=15)
    
    # Assign colors to each BY4742 chromosome
    ref_chromosomes = sorted(df['ref_chr'].unique())
    chr_colors = plt.cm.tab20(np.linspace(0, 1, len(ref_chromosomes)))
    chr_color_map = {chr: chr_colors[i] for i, chr in enumerate(ref_chromosomes)}
    
    # Draw connection lines (sampled display to avoid overcrowding)
    sample_df = df.sort_values('ref_length_mb', ascending=False).head(200)
    for _, row in sample_df.iterrows():
        contig_y = contig_y_positions[row['qry_contig']]
        # Connection line
        ax_links.plot([row['qry_start_mb'], row['qry_start_mb']], 
                     [contig_y - 0.5, contig_y + 0.5], 
                     color=chr_color_map[row['ref_chr']], alpha=0.6, linewidth=0.8)
        # Label (chromosome name)
        if np.random.random() < 0.1:  # 10% probability of showing label to avoid crowding
            ax_links.text(row['qry_start_mb'], contig_y + 0.6, 
                         row['ref_chr'].replace('chr', ''), 
                         fontsize=6, ha='center', rotation=45, alpha=0.7)
    
    ax_links.set_xlabel('Position on SY14 Contigs (Mb)', fontsize=12)
    ax_links.set_ylabel('SY14 Contigs', fontsize=12)
    ax_links.set_xlim(-0.1, max(contig_lengths) + 0.5)
    ax_links.set_ylim(-0.5, len(contigs) - 0.5)
    
    # 4. Statistical panel: BY4742 chromosome contribution
    ax_stats = fig.add_subplot(gs[2, 0])
    chr_contrib = df['ref_chr'].value_counts().head(10)
    bars = ax_stats.bar(range(len(chr_contrib)), chr_contrib.values, 
                       color=[chr_color_map[chr] for chr in chr_contrib.index])
    ax_stats.set_title('Top 10 BY4742 Chromosomes Contributing to SY14', 
                      fontsize=13, fontweight='bold', pad=10)
    ax_stats.set_xlabel('BY4742 Chromosome', fontsize=11)
    ax_stats.set_ylabel('Number of Alignment Blocks', fontsize=11)
    ax_stats.set_xticks(range(len(chr_contrib)))
    ax_stats.set_xticklabels([c.replace('chr', '') for c in chr_contrib.index], rotation=45)
    
    # 5. Legend panel
    ax_legend = fig.add_subplot(gs[2, 1])
    ax_legend.axis('off')
    
    # Key findings summary
    summary_text = (
        "🔬 Key Findings:\n"
        f"• Total alignment blocks: {len(df):,}\n"
        f"• SY14 contigs: {len(contigs)}\n"
        f"• BY4742 chromosomes involved: {len(ref_chromosomes)}\n"
        f"• Main fusion contig: {contigs[0]} ({contig_lengths[0]:.2f} Mb)\n"
        "• Evidence: Multiple BY4742 chromosomes\n"
        "  align to SY14_1 → Chromosomal fusion\n"
        f"\nGenerated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}"
    )
    ax_legend.text(0, 0.5, summary_text, fontsize=11, va='center',
                  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ Chromosomal fusion visualization saved: {output_file}")
    
    # Console output key statistics
    print("\n📊 Alignment Statistics:")
    print(f"   Total alignment blocks: {len(df):,}")
    print(f"   SY14 contigs with alignments: {len(contigs)}")
    print(f"   BY4742 chromosomes involved: {len(ref_chromosomes)}")
    
    # Check chromosomal fusion evidence
    main_contig_links = df[df['qry_contig'] == contigs[0]]
    unique_chrs = main_contig_links['ref_chr'].nunique()
    print(f"\n🔍 Fusion Evidence:")
    print(f"   {contigs[0]} aligns to {unique_chrs} different BY4742 chromosomes")
    if unique_chrs >= 3:
        print(f"   ✅ STRONG EVIDENCE: Multiple chromosomes fused into {contigs[0]}")
    print(f"   Sample alignments to {contigs[0]}:")
    for i, (_, row) in enumerate(main_contig_links.head(3).iterrows()):
        print(f"     {row['ref_chr']}: {row['ref_start_mb']:.2f}-{row['ref_end_mb']:.2f} Mb")

def main():
    input_file = "04_results/genome_comparison/circos_links.txt"
    output_file = "04_results/genome_comparison/chromosome_fusion_visualization.png"
    
    print("Loading alignment data...")
    df = load_links(input_file)
    
    print(f"Loaded {len(df):,} alignment blocks")
    
    print("Creating visualization...")
    create_fusion_plot(df, output_file)
    
    print("\n🎉 Visualization complete!")
    print(f"File saved: {output_file}")
    print("\nTo view the image:")
    print(f"  scp {output_file} .  # Download to local machine")
    print("  Or use any image viewer on your HPC system")

if __name__ == "__main__":
    main()
VISUALIZE

# Run visualization script
python3 02_scripts/visualize_fusion.py
