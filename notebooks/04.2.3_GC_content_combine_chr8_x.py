from Bio import SeqIO
import pandas as pd
import numpy as np
import logging
import sys
import os
from datetime import datetime

from genominterv.remapping import remap
from genominterv.remapping import interval_distance, genomic
from genominterv.remapping import remap_interval_data
import seaborn as sns

import pyBigWig
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm  # optional, for progress bar

# Set up logging
def setup_logging(log_file='gc_analysis.log'):
    """Set up logging to both file and console"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# Initialize logging
logger = setup_logging()
logger.info("Starting GC content analysis script")

# Set matplotlib backend for headless operation
plt.switch_backend('Agg')
logger.info("Set matplotlib backend to 'Agg' for headless operation")

try:
    file_path = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/gc5Base_chrX.txt"
    
    # Load the file — tab-delimited
    logger.info(f"Loading GC data from: {file_path}")
    df_gc = pd.read_csv(file_path, sep="\t", header=None)

    # Name the columns (standard BedGraph format)
    df_gc.columns = ["chrom", "start", "end", "value"]
    logger.info(f"GC data loaded successfully. Shape: {df_gc.shape}")

    def parse_compartment_data(file_name):
        logger.info(f"Parsing compartment data from: {file_name}")
        e1_100kb = pd.read_csv(file_name)

        # remove na (experiment)
        e1_100kb.dropna(inplace=True)

        e1_100kb['start'] = [i*100_000 for i in range(e1_100kb.index.size)]
        e1_100kb['end'] = e1_100kb.start + 100_000
        e1_100kb['sign'] = np.sign(e1_100kb.e1)
        e1_100kb['segment_id'] = ((e1_100kb.sign.shift() != e1_100kb.sign)).cumsum()
        
        comp = e1_100kb.groupby('segment_id', as_index=False).agg(dict(
             e1=['mean', 'sum'], 
             start='min', 
             end='max', 
             segment_id='mean', 
             sign='mean'
        ))
        comp.columns = ['_'.join(col).strip() for col in comp.columns.values]
        comp = comp.rename(
            columns={'start_min':'start',
                     'end_max':'end', 
                     'segment_id_mean':'segment_id', 
                     'sign_mean':'sign'}
        )
        comp['comp'] = ['A' if x > 0 else 'B' for x in comp.sign]
        comp = comp.reset_index()
        comp['chrom'] = 'chrX'
        
        _comp = comp.copy()
        for i in range(1, _comp.index.size-1):
            if np.isnan(_comp.loc[i-1, 'e1_mean']):
                _comp.loc[i, 'start'] = np.nan
            if np.isnan(_comp.loc[i+1, 'e1_mean']):
                _comp.loc[i, 'end'] = np.nan
        _comp = _comp.loc[~_comp.e1_mean.isnull(), :]
        _comp = _comp.reset_index()
        compartment_edges = pd.concat([_comp.start, _comp.end]).sort_values().unique()
        
        compartments = comp.loc[~comp.e1_mean.isnull()].copy()
        compartments['start'] = compartments.start.astype(int)
        compartments['end'] = compartments.end.astype(int)

        return compartments, compartment_edges

    def edge_segments(compartment_edges, flank):
        compartment_edge_segm = pd.DataFrame(np.column_stack((compartment_edges, compartment_edges+flank)), columns=['start', 'end'])
        compartment_edge_segm['chrom'] = 'chrX'
        return compartment_edge_segm

    # Load data
    eigentrack_dir = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/eigentracks"
    logger.info(f"Loading eigentrack data from: {eigentrack_dir}")
    
    eigentrack_files = [
        f for f in os.listdir(eigentrack_dir) if f.endswith("_10Mb.csv")
    ]
    logger.info(f"Found {len(eigentrack_files)} eigentrack files: {eigentrack_files}")

    comps_dict = {}
    edges_dict = {}
    generated_comps=[]
    generated_edges=[]
    a_and_b_comps = []

    for filename in eigentrack_files:
        filepath = os.path.join(eigentrack_dir, filename)
        base = os.path.splitext(filename)[0]
        comps_var = f"{base}_comps"
        edges_var = f"{base}_edges"
        comps, edges = parse_compartment_data(filepath)
        comps_dict[comps_var] = comps
        edges_dict[edges_var] = edges
        generated_comps.append(comps_var)
        generated_edges.append(edges_var)

    logger.info(f"Generated compartment variables: {generated_comps}")
    logger.info(f"Generated edge variables: {generated_edges}")

    for edge_var in generated_edges:
        # Get the numpy array of edges from the dictionary
        edges = edges_dict[edge_var]
        
        # Create the DataFrame using edge_segments with flank=1
        seg_name = f"{edge_var}_interval"
        seg_df = edge_segments(edges, 1)

        # Merge compartment assignment
        comps_var = edge_var.replace("_edges", "_comps")
        if comps_var in comps_dict:
            comps_df = comps_dict[comps_var]
            
            comp_df = pd.DataFrame({
                'comp': comps_df['comp'].reset_index(drop=True),
                'start': seg_df['start'].reset_index(drop=True),
                'end': seg_df['end'].reset_index(drop=True),
                'chrom': seg_df['chrom'].reset_index(drop=True)
            })

            # Save full merged comp_df
            comp_full_name = f"{edge_var}_interval_comp"
                # Save combined A and B compartments
            comp_AB_name = f"{edge_var}_AB"
            globals()[comp_AB_name] = comp_df
            a_and_b_comps.append(comp_AB_name)

            # Split into compartments A and B
            comp_A = comp_df[comp_df['comp'] == 'A'].reset_index(drop=True)
            comp_B = comp_df[comp_df['comp'] == 'B'].reset_index(drop=True)

            # Save A and B splits as new variables
            comp_A_name = f"{edge_var}_A"
            comp_B_name = f"{edge_var}_B"
            globals()[comp_A_name] = comp_A
            globals()[comp_B_name] = comp_B
            a_and_b_comps.append(comp_A_name)
            a_and_b_comps.append(comp_B_name)

    logger.info(f"Generated compartment A and B variables: {a_and_b_comps}")

    # Load remapped data
    logger.info("Loading remapped data files...")
    
    file_path1 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/fibroblast_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    fibroblast = pd.read_csv(file_path1, sep="\t")
    logger.info(f"Loaded fibroblast data: {fibroblast.shape}")

    file_path2 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/sperm_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    sperm = pd.read_csv(file_path2, sep="\t")
    logger.info(f"Loaded sperm data: {sperm.shape}")

    file_path3 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/round_spermatid_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    round_spermatid = pd.read_csv(file_path3, sep="\t")
    logger.info(f"Loaded round_spermatid data: {round_spermatid.shape}")

    file_path4 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/pachytene_spermatocyte_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    pachytene_spermatocyte = pd.read_csv(file_path4, sep="\t")
    logger.info(f"Loaded pachytene_spermatocyte data: {pachytene_spermatocyte.shape}")

    file_path5 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/spermatogonia_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    spermatogonia = pd.read_csv(file_path5, sep="\t")
    logger.info(f"Loaded spermatogonia data: {spermatogonia.shape}")

    # Load chr8 data
    logger.info("Loading chr8 data files...")
    
    file_path1 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/chr8/fibroblast_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    fibroblast_chr8 = pd.read_csv(file_path1, sep="\t")
    logger.info(f"Loaded fibroblast_chr8 data: {fibroblast_chr8.shape}")

    file_path2 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/chr8/sperm_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    sperm_chr8 = pd.read_csv(file_path2, sep="\t")
    logger.info(f"Loaded sperm_chr8 data: {sperm_chr8.shape}")

    file_path3 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/chr8/round_spermatid_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    round_spermatid_chr8 = pd.read_csv(file_path3, sep="\t")
    logger.info(f"Loaded round_spermatid_chr8 data: {round_spermatid_chr8.shape}")

    file_path4 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/chr8/pachytene_spermatocyte_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    pachytene_spermatocyte_chr8 = pd.read_csv(file_path4, sep="\t")
    logger.info(f"Loaded pachytene_spermatocyte_chr8 data: {pachytene_spermatocyte_chr8.shape}")

    file_path5 = "/home/johanulstrup/johan_gpn/people/johanulsrup/johan_gpn/data/macaca/remapping/chr8/spermatogonia_e1_100kb_10Mb_edges_AB_gc_remap_drop_na.txt"
    spermatogonia_chr8 = pd.read_csv(file_path5, sep="\t")
    logger.info(f"Loaded spermatogonia_chr8 data: {spermatogonia_chr8.shape}")

    def filter_into_A_B_compartments(result):
        """
        Splits the input DataFrame into A and B compartments based on sign of 'start' and compartment identity.
        
        Parameters:
            result (pd.DataFrame): Must include columns 'comp' and 'start'.
        
        Returns:
            A_val (pd.DataFrame): Entries where position is consistent with Compartment A.
            B_val (pd.DataFrame): Entries where position is consistent with Compartment B.
        """
        # Entries within A compartment (left of origin if A, right if B)
        A_val = result[
            ((result['comp'] == 'A') & (result['start'] < 0)) |     # is greater than 0
            ((result['comp'] == 'B') & (result['start'] > 0))       # is less than 0
        ].copy()

        # Entries within B compartment (right of origin if A, left if B)
        B_val = result[
            ((result['comp'] == 'A') & (result['start'] > 0)) |     # is less than 0
            ((result['comp'] == 'B') & (result['start'] < 0))       # is greater than 0
        ].copy()

        return A_val, B_val

    fibroblast_A, fibroblast_B = filter_into_A_B_compartments(fibroblast)
    logger.info(f"Filtered fibroblast data: A={fibroblast_A.shape}, B={fibroblast_B.shape}")

    # --- helpers ---------------------------------------------------------------
    # Colors for chrX
    # Colors for main set
    A_color = "#117733"  # deep teal green
    B_color = "#E69F00"  # rich amber/orange

    # Colors for chr8 set
    A_color_chr8 = "#56B4E9"  # bright sky blue
    B_color_chr8 = "#D55E00"  # bold vermilion red


    def _bin_edges(max_dist=2_000_000, bins=20):
        edges = np.linspace(0, max_dist, bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        return edges, centers

    def _edge_availability(edge_df, max_dist=2_000_000, bins=20):
        """Return reverse cumulative fraction of available edges per bin for comp A and B."""
        edges, _ = _bin_edges(max_dist, bins)
        out = {}
        for comp in ("A", "B"):
            starts = edge_df.loc[edge_df["comp"] == comp, "start"].values
            if starts.size < 2:
                counts = np.zeros(bins, dtype=int)
            else:
                half_d = (starts[1:] - starts[:-1]) / 2.0
                counts, _ = np.histogram(half_d, bins=edges)

            # Reverse cumulative sum
            rev_cum = counts[::-1].cumsum()[::-1]

            # Normalize so bin 0 starts at 1.0 (100% availability)
            frac = rev_cum / rev_cum[0] if rev_cum[0] else np.zeros_like(rev_cum, dtype=float)

            out[comp] = frac
        return out  # {"A": frac, "B": frac}

    def _comp_sides(df):
        """Split into A-side and B-side using your sign logic (no changes)."""
        A_side = df[((df['comp'] == 'A') & (df['start'] < 0)) |
                    ((df['comp'] == 'B') & (df['start'] > 0))]
        B_side = df[((df['comp'] == 'A') & (df['start'] > 0)) |
                    ((df['comp'] == 'B') & (df['start'] < 0))]
        return A_side, B_side

    def _bin_stats(gc_df, edges, edge_frac):
        """Mean and 95% CI using SEM, scaled by edge availability as effective N."""
        means = []
        errs  = []
        counts= []
        for lo, hi, f in zip(edges[:-1], edges[1:], edge_frac):
            sel = gc_df[(gc_df["absmid"] >= lo) & (gc_df["absmid"] < hi)]["value"].dropna().to_numpy()
            n = sel.size
            counts.append(n)
            if n == 0:
                means.append(np.nan); errs.append(np.nan); continue
            m = float(sel.mean())
            if n == 1:
                # single-point guardrail: plot point, small default error
                e = 0.05
            else:
                # effective N down-weights sparse edge availability
                eff_n = max(n * max(f, 0.01), 1.0)
                sd = sel.std(ddof=1)
                sem = sd / np.sqrt(eff_n)
                e = 1.96 * sem  # approx 95% CI
            means.append(m); errs.append(e)
        return np.array(means), np.array(errs), np.array(counts)

    # --- plotting --------------------------------------------------------------

    def plot_gc_by_compartment_simple(ax, df, edge_df, title,
                                      max_dist=2_000_000, bins=20):
        """One clean plot: A/B means with edge-weighted 95% CI error bars."""
        # prep
        df = df.loc[df["absmid"] <= max_dist].copy()
        edges, centers = _bin_edges(max_dist, bins)
        avail = _edge_availability(edge_df, max_dist, bins)
        A_side, B_side = _comp_sides(df)

        # stats
        A_m, A_e, A_n = _bin_stats(A_side, edges, avail["A"])
        B_m, B_e, B_n = _bin_stats(B_side, edges, avail["B"])

        # draw (minimal styling)
        ax.errorbar(centers, A_m, yerr=A_e, marker='o', linestyle='-',
                    label='Compartment A', capsize=3, color =A_color)
        ax.errorbar(centers, B_m, yerr=B_e, marker='s', linestyle='--',
                    label='Compartment B', capsize=3, color =B_color)

        ax.set_title(title)
        ax.set_xlabel("Distance to compartment edge (bp)")
        ax.set_ylabel("Average GC content (%)")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best")

        return ax

    def plot_gc_by_compartment_combined_chr(ax, df_chrx, df_chr8, edge_df, title,
                                           max_dist=2_000_000, bins=20):
        """Plot both chrX and chr8 data on the same plot with different line styles."""
        # prep chrX data
        df_chrx = df_chrx.loc[df_chrx["absmid"] <= max_dist].copy()
        edges, centers = _bin_edges(max_dist, bins)
        avail = _edge_availability(edge_df, max_dist, bins)
        A_side_x, B_side_x = _comp_sides(df_chrx)

        # prep chr8 data
        df_chr8 = df_chr8.loc[df_chr8["absmid"] <= max_dist].copy()
        A_side_8, B_side_8 = _comp_sides(df_chr8)

        # stats for chrX
        A_m_x, A_e_x, A_n_x = _bin_stats(A_side_x, edges, avail["A"])
        B_m_x, B_e_x, B_n_x = _bin_stats(B_side_x, edges, avail["B"])
        
        # stats for chr8
        A_m_8, A_e_8, A_n_8 = _bin_stats(A_side_8, edges, avail["A"])
        B_m_8, B_e_8, B_n_8 = _bin_stats(B_side_8, edges, avail["B"])

        # draw chrX data (solid lines)
        ax.errorbar(centers, A_m_x, yerr=A_e_x, marker='o', linestyle='-',
                    label='Compartment A (chrX)', capsize=3, color=A_color, linewidth=2)
        ax.errorbar(centers, B_m_x, yerr=B_e_x, marker='s', linestyle='-',
                    label='Compartment B (chrX)', capsize=3, color=B_color, linewidth=2)
        
        # draw chr8 data (dashed lines with different colors)
        ax.errorbar(centers, A_m_8, yerr=A_e_8, marker='o', linestyle='--',
                    label='Compartment A (chr8)', capsize=3, color=A_color_chr8, linewidth=2)
        ax.errorbar(centers, B_m_8, yerr=B_e_8, marker='s', linestyle='--',
                    label='Compartment B (chr8)', capsize=3, color=B_color_chr8, linewidth=2)

        ax.set_title(title)
        ax.set_xlabel("Distance to compartment edge (bp)")
        ax.set_ylabel("Average GC content (%)")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize='small')

        return ax

    def plot_many_datasets_combined_chr(datasets_chrx, datasets_chr8, edge_datasets,
                                        max_dist=2_000_000, bins=20, ncols=2, figsize=(15, 12)):
        """
        Plot chrX and chr8 data together for each tissue type.
        datasets_chrx: list of (title, df_chrx)
        datasets_chr8: list of (title, df_chr8)
        edge_datasets: list of edge_df (same order)
        """
        n = len(datasets_chrx)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.atleast_1d(axes).ravel()

        for i, ((title_x, df_x), (title_8, df_8), edge_df) in enumerate(zip(datasets_chrx, datasets_chr8, edge_datasets)):
            plot_gc_by_compartment_combined_chr(axes[i], df_x, df_8, edge_df, title_x.replace(" (chrX)", ""),
                                               max_dist=max_dist, bins=bins)
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

        fig.tight_layout()
        return fig

    logger.info("Creating combined chrX and chr8 plot...")
    fig = plot_many_datasets_combined_chr(
        datasets_chrx=[
            ("Fibroblast (chrX)", fibroblast),
            ("Sperm (chrX)", sperm),
            ("Round Spermatid (chrX)", round_spermatid),
            ("Pachytene Spermatocyte (chrX)", pachytene_spermatocyte),
            ("Spermatogonia (chrX)", spermatogonia),
        ],
        datasets_chr8=[
            ("Fibroblast (chr8)", fibroblast_chr8),
            ("Sperm (chr8)", sperm_chr8),
            ("Round Spermatid (chr8)", round_spermatid_chr8),
            ("Pachytene Spermatocyte (chr8)", pachytene_spermatocyte_chr8),
            ("Spermatogonia (chr8)", spermatogonia_chr8),
        ],
        edge_datasets=[
            fibroblast_e1_100kb_10Mb_edges_AB,
            sperm_e1_100kb_10Mb_edges_AB,
            round_spermatid_e1_100kb_10Mb_edges_AB,
            pachytene_spermatocyte_e1_100kb_10Mb_edges_AB,
            spermatogonia_e1_100kb_10Mb_edges_AB,
        ],
        bins=20,
        max_dist=2_000_000,
    )
    
    # Save the plot with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    plot_filename = f"gc_compartment_analysis_chrX_chr8_combined_{timestamp}.png"
    fig.savefig(plot_filename, dpi=300, bbox_inches='tight')
    logger.info(f"Combined plot saved as: {plot_filename}")
    
    # Also save as PDF
    pdf_filename = f"gc_compartment_analysis_chrX_chr8_combined_{timestamp}.pdf"
    fig.savefig(pdf_filename, bbox_inches='tight')
    logger.info(f"Combined plot also saved as: {pdf_filename}")
    
    plt.close(fig)  # Free up memory
    
    logger.info("Analysis completed successfully!")

except Exception as e:
    logger.error(f"An error occurred: {str(e)}", exc_info=True)
    sys.exit(1)