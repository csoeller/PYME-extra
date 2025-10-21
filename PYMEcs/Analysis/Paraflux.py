import matplotlib.pyplot as plt
import numpy as np
import logging
logger = logging.getLogger(__file__)
import pandas as pd
import os

from PYMEcs.pyme_warnings import warn

# note that currently everything is hard coded for the standard 3D sequence
# TODO: should really check the metadata we have for (recent) zarr and check the CCRLimit array etc
#       and figure things out from that!!


# --- start Alex B contributed functions to (possibly save) and plot ITR stats (Paraflux like) ---
def paraflux_mk_df_fm(mfxdata):
    # Convert structured array (original mfx data) to DataFrame (containing all data from MFX experiment)
    df_mfx = pd.DataFrame()
    for name in mfxdata.dtype.names:
        col = mfxdata[name]
        if col.ndim == 1:
            df_mfx[name] = col
        else:
            n_dim = col.shape[1]
            expanded = pd.DataFrame(col.tolist(), index=df_mfx.index if not df_mfx.empty else None)
            if name in ['loc', 'lnc']:
                labels = ["x", "y", "z"]
                expanded.columns = [f"{name}_{labels[i]}" for i in range(n_dim)]
            elif name == 'dcr':
                expanded.columns = [f"{name}_{i+1}" for i in range(n_dim)]
            else:
                expanded.columns = [f"{name}_{i}" for i in range(n_dim)]
            df_mfx = pd.concat([df_mfx, expanded], axis=1)

    # Create failure_map (used to interpret failure reasons in analyze_failures)
    failure_map = {
        1: "Valid final", 2: "Valid not final",
        4: "Derived iteration", 5: "Reserved",
        6: "CFR failure", 8: "No signal",
        9: "DAC out of range", 11: "Background measurement"
    }

    return (df_mfx, failure_map)

# ==================================================
# --- Analysis functions ( Paraflux-like) ---
# ==================================================

# Create a df with list of vld tids per iteration + additional basic stats
def build_valid_df(df_mfx):
    if not isinstance(df_mfx, pd.DataFrame): # Convert to DataFrame if input is structured array
        df = pd.DataFrame(df_mfx)
    else:
        df = df_mfx.copy()
    df_valid = df[df['vld']] # Select only valid localizations
    vld_itr = df_valid.groupby('itr')['tid'].apply(lambda x: list(set(x))).reset_index() # Get list of unique tids per iteration
    vld_itr['Axis'] = np.where(vld_itr['itr'] % 2 == 0, 'x,y', 'z') # Add a col with axis of each iteration 
    vld_itr['vld loc count'] = vld_itr['tid'].apply(len) # Count valid locs per iteration
    vld_itr['failed loc count'] = vld_itr['vld loc count'].shift(1, fill_value=vld_itr['vld loc count'].iloc[0]) - vld_itr['vld loc count'] # Calculate failed loc count per iteration
    vld_itr.loc[0, 'failed loc count'] = 0 # Set failed loc count of first iteration to 0 (instead of NaN)
    vld_itr['failed loc cum sum'] = vld_itr['failed loc count'].cumsum() # Cumulative sum of failed locs
    return vld_itr

# Compute percentages of passed and failed localizations (from build_valid_df)
def compute_percentages(vld_itr):
    initial_count = vld_itr['vld loc count'].iloc[0] # Percentage calculations are based on initial count of valid locs
    vld_itr['passed itr %'] = (vld_itr['vld loc count'] * 100 / initial_count).round(1)
    vld_itr['failed % per itr'] = (vld_itr['failed loc count'] * 100 / initial_count).round(1)
    pair_sums = {} # This is done to mimic results from Paraflux
    for i in range(1, len(vld_itr), 2):
        pair_sums[i] = vld_itr.loc[i-1:i, 'failed % per itr'].sum().round(1)
    vld_itr['failed % per itr pairs'] = vld_itr.index.map(pair_sums)
    vld_itr['failed cum sum %'] = vld_itr['failed % per itr'].cumsum().round(1)
    return vld_itr

# Analyze failures between consecutive iterations and categorize them based on failure_map (found on wiki from Abberior)
def analyze_failures(vld_itr, df_mfx, failure_map):
    def analyze_failures_single_steps(vld_itr, df, itr_from, itr_to, failure_map):
        tids_from = set(vld_itr.loc[vld_itr['itr'] == itr_from, 'tid'].iloc[0]) # Select valid tids of the previous iteration
        tids_to   = set(vld_itr.loc[vld_itr['itr'] == itr_to, 'tid'].iloc[0]) # Select valid tids of the current iteration
        failed_tids = tids_from - tids_to # Determine tids that failed in the current iteration
        failed_df = df[df['tid'].isin(failed_tids) & (df['itr'] == itr_to)] # Create a df with only failed tids in the current iteration
        counts = failed_df['sta'].value_counts().rename_axis("sta").reset_index(name="count") # Count failure reasons
        counts["reason"] = counts["sta"].map(failure_map).fillna("Other") # Map failure reasons using failure_map
        counts.insert(0, "itr", itr_to) # Add iteration column
        return counts

    pairs = [(i, i+1) for i in range(vld_itr['itr'].max())] # Create pairs of consecutive iterations
    # Analyze failures for each pair and concatenate results
    failure_results = pd.concat(
        [analyze_failures_single_steps(vld_itr, df_mfx, i_from, i_to, failure_map) for i_from, i_to in pairs],
        ignore_index=True
    ) 
    # Pivot the results to have failure reasons as columns
    failure_pivot = failure_results.pivot_table( 
        index="itr", columns="reason", values="count", fill_value=0
    ).reset_index() 
    return vld_itr.merge(failure_pivot, on="itr", how="left")

# Compute percentages for failure reasons
def add_failure_metrics(vld_itr, initial_count):
    cfr_map = {5: 4, 7: 6} #map ITR where CFR failures occurs
    vld_itr['CFR failure %'] = np.nan # Initialize column with NaNs
    # Calculate CFR failure percentages based on cfr_map
    for target_itr, source_itr in cfr_map.items():
        if not vld_itr.loc[vld_itr['itr'] == source_itr, 'CFR failure'].empty: # Check if CFR failure data exists for the source iteration
            val = vld_itr.loc[vld_itr['itr'] == source_itr, 'CFR failure'].values[0] # Get the CFR failure count
            vld_itr.loc[vld_itr['itr'] == target_itr, 'CFR failure %'] = (val / initial_count * 100).round(1) # Calculate percentage and assign to target iteration
    # Calculate No signal percentage for each iteration       
    vld_itr['No signal %'] = (vld_itr['No signal'] * 100 / initial_count).round(1)
    # Define groups of iterations for No signal percentage calculation
    no_signal_groups = {1: [0, 1],3: [2, 3], 5: [4, 5], 7: [6, 7], 9: [8, 9]}
    # Calculate No signal percentage for each group and map to iterations
    no_signal_pct = {
        target_itr: (vld_itr.loc[vld_itr['itr'].isin(group), 'No signal'].sum() / initial_count * 100).round(1)
        for target_itr, group in no_signal_groups.items()
    }
    vld_itr['No signal % per itr pairs'] = vld_itr['itr'].map(no_signal_pct)
    return vld_itr

# Plot like in Paraflux
def paraflux_itr_plot(vld_paraflux):
    # Mapping rules
    label_map = {
        "passed": "Passed",
        "CFR": "CFR-filtered",
        "No signal": "Dark",
        "DAC": "Out of range",
        "Other": "Other",
    }

    # Function to get pretty label based on the label map (substring matching from vld_paraflux col names)
    def pretty_label(colname):
        """Map colname to user-friendly label based on substring rules."""
        for key, label in label_map.items():
            if key.lower() in colname.lower():
                return label
        return colname  # fallback: keep original name

    # Keep only odd iterations (Paraflux style, i.e., 1, 3, 5, 7, 9)
    vld_paraflux = vld_paraflux[vld_paraflux['itr'] % 2 == 1]

    # Set figure size
    plt.figure(figsize=(8, 6))

    # Base positions
    r1 = np.arange(len(vld_paraflux)) # Define the positions for each bar
    names = vld_paraflux['itr'] # Names of group
    barWidth = 0.85 # Bar width

    # Colors (extendable if more cols are added)
    colors = ["#0072B2", "#009E73", "#D55E00", "#E69F00", "#CC79A7"]

    # Define the bottom position for stacking
    bottompos = np.zeros(len(vld_paraflux))

    # Plot each column as a stacked bar
    for i, col in enumerate(vld_paraflux.columns[1:]): # Enumerate over all columns except 'itr'
        vals = vld_paraflux[col].fillna(0) # Get the values for the current column, filling NaNs with 0
        labels = pretty_label(col) # Get the pretty label for the legend

        plt.bar(
            r1, vals, bottom=bottompos,
            color=colors[i % len(colors)],
            edgecolor="white", width=barWidth, label=labels
        ) # Create the bar

        # Add labels inside each bar
        for j, v in enumerate(vals):
            if v > 0:
                plt.text(r1[j], bottompos[j] + v / 2, f"{v:.1f}%", # Only add text if value > 0
                        ha="center", va="center",
                        color="black",
                        fontsize=9)

        bottompos += vals.values

    # X/Y labels
    plt.xticks(r1, names)
    plt.xlabel("Iteration")
    plt.ylabel("Events (%)")
    plt.axhline(y=100, color="gray", linestyle="--", linewidth=1)

    plt.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    plt.show()

# Main function to compute stats of failed and valid localizations and save results
def compute_vld_stats(df_mfx, failure_map, store_path=None, ts=None):
    # Run the analysis steps
    vld_itr = build_valid_df(df_mfx)
    vld_itr = compute_percentages(vld_itr)
    vld_itr = analyze_failures(vld_itr, df_mfx, failure_map)
    initial_count = vld_itr['vld loc count'].iloc[0]
    vld_itr = add_failure_metrics(vld_itr, initial_count)
    vld_paraflux = paraflux_itr_plot(vld_itr[['itr', 'passed itr %', 'CFR failure %', 'No signal % per itr pairs']])

    if store_path is not None:
        if ts is None:
            warn("No timestamp found in metadata, saving as default.")
            timestamp = "no_ts"

        vld_itr = vld_itr.drop(columns='tid', errors='ignore')
        default_dir = str(store_path.parent)
        full_path = os.path.join(default_dir, f"{timestamp}_iteration_stats_full.csv")
        paraflux_path = os.path.join(default_dir, f"{timestamp}_iteration_stats_Paraflux_only.csv")

        vld_itr.to_csv(full_path, index=False)
        keep_cols = ["itr", 'passed itr %', 'CFR failure %', 'No signal % per itr pairs']
        vld_paraflux = vld_itr[keep_cols]
        vld_paraflux.to_csv(paraflux_path, index=False)

        logger.debug(f"✔ Saved full results to: {full_path}")
        logger.debug(f"✔ Saved cleaned results to: {paraflux_path}")

    return vld_itr

### --- End of Alex B contributed functions ---
