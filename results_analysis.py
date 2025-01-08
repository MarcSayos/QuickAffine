import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

color_swg = "#d62728"
color_windowed = "#4fc94f"
color_banded = "#2ca02c"
color_parascan = "#ff7f0e"
color_paradiag = "#1f77b4"

def read_input_file(filename):
    """
    Read the input file with comma-separated values and return a pandas DataFrame
    """
    columns = [
        'window_size', 'overlap_size', 'penalty_set', 'avg_query_length', 'dataset_type',
        'elapsed_SWG', 'memory_SWG', 'cells_SWG', 'score_SWG',
        'elapsed_windowed', 'memory_windowed', 'cells_windowed', 'score_windowed',
        'elapsed_banded', 'memory_banded', 'cells_banded', 'score_banded',
        'elapsed_parasail_scan', 'memory_parasail_scan', 'cells_parasail_scan', 'score_parasail_scan',
        'elapsed_parasail_diag', 'memory_parasail_diag', 'cells_parasail_diag', 'score_parasail_diag'
    ]
    
    # Read CSV file with comma separator
    df = pd.read_csv(filename, header=None, names=columns, skipinitialspace=True)
    
    # Convert numeric columns 
    numeric_columns = [
        'window_size', 'overlap_size', 'avg_query_length',
        'elapsed_SWG', 'memory_SWG', 'cells_SWG', 'score_SWG',
        'elapsed_windowed', 'memory_windowed', 'cells_windowed', 'score_windowed',
        'elapsed_banded', 'memory_banded', 'cells_banded', 'score_banded',
        'elapsed_parasail_scan', 'memory_parasail_scan', 'cells_parasail_scan', 'score_parasail_scan',
        'elapsed_parasail_diag', 'memory_parasail_diag', 'cells_parasail_diag', 'score_parasail_diag'
    ]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')
    
    return df

def create_figure_2_plot(df, output_dir):
    """
    Create a plot with Elapsed time for SWG and Windowed.
    Windowed includes 6 configurations: [32,8], [32,16], [64,16], [64,32], [128,32], [128,64].
    """
    os.makedirs(output_dir, exist_ok=True)

    windowed_configs = [(32, 8), (32, 16), (64, 16), (64, 32), (128, 32), (128, 64)]

    swg_data = df.groupby('avg_query_length')['elapsed_SWG'].mean()

    windowed_data = []
    for config in windowed_configs:
        ws, oss = config
        config_data = df[(df['window_size'] == ws) & (df['overlap_size'] == oss)].groupby('avg_query_length')['elapsed_windowed'].mean()
        windowed_data.append(config_data)

    combined_data = pd.DataFrame({'SWG': swg_data})
    for idx, config_data in enumerate(windowed_data):
        combined_data[f'Windowed {windowed_configs[idx]}'] = config_data

    bar_width = 0.1
    x = np.arange(len(combined_data.index))

    plt.figure(figsize=(14, 8))

    plt.bar(x - 3 * bar_width, combined_data['SWG'], width=bar_width, label='SWG', color=color_swg)
    colors = ['#53afed', '#3d97d4', '#297fba', '#18689e', '#0d5485', '#064570']
    for idx, config in enumerate(windowed_configs):
        plt.bar(x - 2 * bar_width + idx * bar_width, combined_data[f'Windowed {config}'], width=bar_width, label=f'Windowed {config}', color=colors[idx])

    plt.xlabel('Sequence Length', fontsize=12)
    plt.ylabel('Elapsed Time (seconds)', fontsize=12)
    plt.title('Elapsed Time: SWG vs Windowed Configurations', fontsize=14)
    plt.xticks(x, ['Simulated 100b', 'Illumina ~250b', 'Simulated 1000b', 'Nanopore ~5000b', 'Simulated 10000b'], fontsize=10)
    plt.legend(title='Algorithm and Configurations')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'figure_2.png'))
    plt.close()
    return 1

def create_figure_3_table(df, output_dir):
    """
    Create a table summarizing computed cells for specified configurations and datasets.
    Generates both a CSV and a LaTeX table.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Define configurations and dataset types
    configs = [(32, 8), (32, 16), (64, 16), (64, 32), (128, 32), (128, 64)]
    dataset_types = ['Sim 100b', 'Sim 1000b', 'Sim 10Kb', 'Illumina', 'Nanopore']

    # Map dataset types for filtering and final naming
    dataset_mapping = {
        'Sim 100b': {'dataset_type': 'Simulated', 'avg_query_length': 100},
        'Sim 1000b': {'dataset_type': 'Simulated', 'avg_query_length': 1000},
        'Sim 10Kb': {'dataset_type': 'Simulated', 'avg_query_length': 10000},
        'Illumina': {'dataset_type': 'Illumina'},
        'Nanopore': {'dataset_type': 'Nanopore'}
    }

    # Prepare the table data
    table_data = []
    for ws, oss in configs:
        row = []
        total_cells = 0
        for dataset, filters in dataset_mapping.items():
            subset = df[(df['window_size'] == ws) & (df['overlap_size'] == oss)]
            if 'avg_query_length' in filters:
                subset = subset[(subset['dataset_type'] == filters['dataset_type']) & (subset['avg_query_length'] == filters['avg_query_length'])]
            else:
                subset = subset[subset['dataset_type'] == filters['dataset_type']]
            total_cells = subset['cells_windowed'].sum() if not subset.empty else 0
            row.append(total_cells)
        table_data.append(row)

    # Create DataFrame for the table
    table_df = pd.DataFrame(table_data, columns=dataset_types, index=[f'{ws}, {oss}' for ws, oss in configs])

    # Save as CSV
    csv_path = os.path.join(output_dir, 'table_3.csv')
    table_df.to_csv(csv_path, float_format='%.0f')

    # Save as LaTeX table
    latex_path = os.path.join(output_dir, 'table_3.txt')
    with open(latex_path, 'w') as f:
        f.write("\\begin{table}[ht]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{l " + " ".join(["c"] * len(dataset_types)) + "}\\hline\n")
        f.write("\\rowcolor[HTML]{AFAFAF} \\textbf{Configuration} & \\textbf{" + "} & \\textbf{".join(dataset_types) + "} \\\\ \\hline\n")
        for idx, row in table_df.iterrows():
            f.write(f"{idx} & " + " & ".join([f"{val:.0f}" for val in row]) + " \\\\\n")
        f.write("\\hline \n") 
        f.write("\\textbf{Real n. of cells} & \\textbf{2000} & \\textbf{2000} & \\textbf{2000} & \\textbf{1546} & \\textbf{715} \\\\ \\hline \\end{tabular}\n")
        f.write("\\caption{Computed Cells (millions) in bounding by Configuration and Dataset}\n")
        f.write("\\label{tab:computed_cells_bound}\n")
        f.write("\\end{table}\n")

    return 1

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

def create_figure_4_plot(output_dir):
    """
    Create a plot with % accuracy of bound (Score deviation) for multiple datasets.
    Generates 6 subplots with 6 lines each, representing different window/overlap configurations.
    """
    # Read the input file
    columns = ['window_size', 'overlap_size', 'number_of_bases', 'dataset_type', 'real_score', 'windowed_score']
    data = pd.read_csv('scores.out', sep=',', header=None, names=columns)

    # Calculate error
    data['error'] = (data['windowed_score'] - data['real_score']) / data['windowed_score']

    # Normalize error to ensure it's in [0, 1]
    data['normalized_error'] = data['error'].clip(lower=0, upper=1)

    # Define the parts
    parts = [
        ("Simulated", 100),
        ("Simulated", 1000),
        ("Simulated", 10000),
        ("Illumina", None),
        ("Nanopore", None),
        ("PacBio", None),
    ]

    window_overlap_pairs = [(32, 8), (32, 16), (64, 16), (64, 32), (128, 32), (128, 64)]

    # Create the figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, (dataset_type, number_of_bases) in enumerate(parts):
        ax = axes[i]
        
        if number_of_bases is not None:
            part_data = data[(data['dataset_type'] == dataset_type) & (data['number_of_bases'] == number_of_bases)]
            title = f"Simulated {number_of_bases} Bases"
        else:
            part_data = data[data['dataset_type'] == dataset_type]
            title = dataset_type

        for window_size, overlap_size in window_overlap_pairs:
            pair_data = part_data[(part_data['window_size'] == window_size) & (part_data['overlap_size'] == overlap_size)]
            
            if pair_data.empty:
                # Ensure an empty line still goes from (0,0) to (1,1)
                sorted_errors = np.array([0, 1])
                cumulative_count = np.array([0, 1])
            else:
                # Sort by normalized error
                sorted_errors = np.sort(pair_data['normalized_error'])
                
                # Generate cumulative distribution
                cumulative_count = np.arange(1, len(sorted_errors) + 1) / len(sorted_errors)
                
                # Scale cumulative_count to ensure it reaches 1.0
                cumulative_count = cumulative_count / cumulative_count[-1]
                
                # Add endpoints
                sorted_errors = np.concatenate(([0], sorted_errors, [1]))
                cumulative_count = np.concatenate(([0], cumulative_count, [1]))

            # Plot
            ax.plot(sorted_errors, cumulative_count, label=f"({window_size}, {overlap_size})")



        # Customize the plot
        ax.set_title(title)
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True)

    # Set common labels
    fig.text(0.5, 0.02, 'Normalized Score', ha='center', fontsize=14)
    fig.text(0.02, 0.5, 'Cumulative Count', va='center', rotation='vertical', fontsize=14)

    # Add a single legend for all plots
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    # Remove duplicate entries from the legend
    unique_handles_labels = dict(zip(labels, handles))
    fig.legend(unique_handles_labels.values(), unique_handles_labels.keys(), loc='upper center', ncol=3, fontsize=12, bbox_to_anchor=(0.5, 0.999))


    # Adjust layout
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])

    # Save plot
    plt.savefig(os.path.join(output_dir, 'figure_4.png'))
    plt.close()
    return 1

    
def create_figures_4_5_6_plot(output_dir, bases):
    """
    Create a plot with % accuracy of bound (Score deviation) for {bases} bases.
    Generate 6 subplots for each combination of window_size and dataset_type.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the scores.out file
    scores_file = 'scores.out'
    columns = ['window_size', 'number_of_bases', 'dataset_type', 'real_score', 'windowed_score']
    df = pd.read_csv(scores_file, header=None, names=columns)

    # Filter rows with number_of_bases = 100 and real_score <= windowed_score
    filtered_df = df.loc[(df['number_of_bases'] == bases) & (df['real_score'] <= df['windowed_score'])].copy()

    # Calculate the error ((windowed_score - real_score) / windowed_score)
    filtered_df['error'] = (filtered_df['windowed_score'] - filtered_df['real_score']) / filtered_df['windowed_score']

    # Create the subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True)
    fig.suptitle(f'Cumulative Accuracy of Bound ({bases} Bases)', fontsize=16)

    # Define window_size and dataset_type combinations
    window_sizes = [32, 64, 128]
    dataset_types = ['Simulated', 'Real']

    # Iterate over combinations and plot each subplot
    for i, dataset_type in enumerate(dataset_types):
        for j, window_size in enumerate(window_sizes):
            ax = axes[i, j]
            subset = filtered_df[(filtered_df['dataset_type'] == dataset_type) & (filtered_df['window_size'] == window_size)]

            if not subset.empty:
                # Sort by error and calculate cumulative proportion
                subset = subset.sort_values('error')
                subset['cumulative'] = np.arange(1, len(subset) + 1) / len(subset)

                # Add (0, 0) and (1, 1) to ensure the line spans the entire range
                subset = pd.concat([
                    pd.DataFrame({'error': [0.0], 'cumulative': [0.0]}),
                    subset,
                    pd.DataFrame({'error': [1.0], 'cumulative': [1.0]})
                ]).reset_index(drop=True)

                # Plot the line
                ax.plot(subset['error'], subset['cumulative'], label=f'{dataset_type}, Window Size {window_size}', color=color_paradiag)

            ax.set_xlim(-0.05, 1.05)
            ax.set_ylim(-0.05, 1.05)
            ax.set_title(f'{dataset_type}, Window Size {window_size}', fontsize=12)
            ax.grid(True, linestyle='--', alpha=0.6)

    # Set common labels
    fig.text(0.5, 0.02, 'Normalized Score', ha='center', fontsize=14)
    fig.text(0.02, 0.5, 'Cumulative Count', va='center', rotation='vertical', fontsize=14)

    # Adjust layout
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])

    # Save plot
    if bases == 100: plt.savefig(os.path.join(output_dir, 'figure_4.png'))
    if bases == 1000: plt.savefig(os.path.join(output_dir, 'figure_5.png'))
    if bases == 10000: plt.savefig(os.path.join(output_dir, 'figure_6.png'))
    plt.close()
    return 1

def create_figure_7_plot(df, output_dir):
    """
    Create a plot with Elapsed time Windowed vs Banded.
    Separated by datasets and window-overlap configurations, normalized within each group.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Define x-axis labels and configurations
    datasets = ['Simulated 100b', 'Illumina ~250b', 'Simulated 1000b', 'Nanopore ~5000b', 'Simulated 10000b']
    configs = [(32, 8), (32, 16), (64, 16), (64, 32), (128, 32), (128, 64)]

    # Prepare data
    grouped_data = []
    for dataset, length in zip(datasets, [100, 248, 1000, 4847, 10000]):
        dataset_data = []
        for ws, oss in configs:
            subset = df[(df['avg_query_length'] == length) & 
                        (df['window_size'] == ws) & 
                        (df['overlap_size'] == oss)]
            if not subset.empty:
                elapsed_windowed = subset['elapsed_windowed'].mean()
                elapsed_banded = subset['elapsed_banded'].mean()
                dataset_data.append((elapsed_banded, elapsed_windowed))
            else:
                dataset_data.append((0, 0))
        # Normalize by the largest value in the group
        max_value = max((max(pair) for pair in dataset_data if max(pair) > 0), default=1)
        dataset_data = [(banded / max_value if max_value > 0 else 0, 
                         windowed / max_value if max_value > 0 else 0) for banded, windowed in dataset_data]
        grouped_data.append(dataset_data)

    # Plot
    plt.figure(figsize=(16, 8))
    bar_width = 0.1
    x = np.arange(len(datasets))

    for i, config in enumerate(configs):
        banded_values = [grouped_data[j][i][0] for j in range(len(datasets))]
        windowed_values = [grouped_data[j][i][1] for j in range(len(datasets))]
        bottom_values = banded_values  # Stack banded at the bottom

        bars = plt.bar(x + i * bar_width - 3 * bar_width, bottom_values, 
                       width=bar_width, color=color_banded, edgecolor='white', label=f'QuickAffine (Align)' if i == 0 else "")
        plt.bar(x + i * bar_width - 3 * bar_width, windowed_values, 
                width=bar_width, bottom=bottom_values, color=color_windowed, edgecolor='white', label=f'QuickAffine (Bound)' if i == 0 else "")

        # Add text on top of the bars
        for j, bar in enumerate(bars):
            plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f'  {config}', 
                     ha='center', va='bottom', fontsize=8, rotation=90)

    # Customize plot
    plt.xlabel('Dataset and Sequence Length', fontsize=12)
    plt.ylabel('Normalized Elapsed Time (seconds)', fontsize=12)
    plt.title('Normalized Elapsed Time: Windowed + Banded by Dataset and Configuration', fontsize=14)
    plt.xticks(x, datasets, fontsize=10)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'figure_7.png'))
    plt.close()
    return 1

def create_figure_9_plot(df, output_dir):
    """
    Create a plot with Elapsed time SWG vs Windowed vs Banded vs ParaScan vs ParaDiag with Simulated Datasets.
    """
    # Filter for simulated datasets
    filtered_df = df.loc[(df['dataset_type'] == 'Simulated') & (df['window_size'] == 32) & (df['overlap_size'] == 8)].copy()
    # filtered_df = df[df['dataset_type'] == 'Simulated']

    # Group by sequence length and compute mean elapsed times
    grouped = filtered_df.groupby('avg_query_length')[['elapsed_SWG', 'elapsed_windowed', 'elapsed_banded', 'elapsed_parasail_scan', 'elapsed_parasail_diag']].mean().reset_index()

    # Plot
    plt.figure(figsize=(12, 6))
    x = np.arange(len(grouped['avg_query_length']))
    width = 0.15

    plt.bar(x - 1.5*width, grouped['elapsed_SWG'], width, label='SWG', color=color_swg, edgecolor='white')
    plt.bar(x - 0.5*width, grouped['elapsed_parasail_diag'], width, label='ParaDiag', color=color_paradiag, edgecolor='white')
    plt.bar(x + 0.5*width, grouped['elapsed_parasail_scan'], width, label='ParaScan', color=color_parascan, edgecolor='white')
    plt.bar(x + 1.5*width, grouped['elapsed_windowed'] + grouped['elapsed_banded'], width, label='QuickAffine (Bound)', color=color_windowed, edgecolor='white')
    plt.bar(x + 1.5*width, grouped['elapsed_banded'], width, label='QuickAffine (Align)', color=color_banded, edgecolor='white')

    plt.xlabel('Dataset and Sequence Length', fontsize=12)
    plt.ylabel('Elapsed Time (seconds)', fontsize=12)
    plt.title('Elapsed Time: Algorithms with Simulated Datasets', fontsize=14)
    plt.xticks(x, ['Simulated 100b', 'Simulated 1000b', 'Simulated 10000b'], fontsize=10)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'figure_9.png'))
    plt.close()
    return 1

def create_figure_10_plot(df, output_dir):
    """
    Create a plot with Elapsed time SWG vs Windowed vs Banded vs ParaScan vs ParaDiag with Real Datasets.
    """
    # Filter for real datasets
    filtered_df = df.loc[(df['dataset_type'].isin(['Nanopore', 'Illumina', 'PacBio'])) & (df['window_size'] == 32) & (df['overlap_size'] == 8)].copy()
    # filtered_df = df[df['dataset_type'].isin(['Nanopore', 'Illumina', 'PacBio'])]

    # Group by sequence length and compute mean elapsed times
    grouped = filtered_df.groupby('avg_query_length')[['elapsed_SWG', 'elapsed_windowed', 'elapsed_banded', 'elapsed_parasail_scan', 'elapsed_parasail_diag']].mean().reset_index()

    # Plot
    plt.figure(figsize=(12, 6))
    x = np.arange(len(grouped['avg_query_length']))
    width = 0.15

    plt.bar(x - 1.5*width, grouped['elapsed_SWG'], width, label='SWG', color=color_swg, edgecolor='white')
    plt.bar(x - 0.5*width, grouped['elapsed_parasail_diag'], width, label='ParaDiag', color=color_paradiag, edgecolor='white')
    plt.bar(x + 0.5*width, grouped['elapsed_parasail_scan'], width, label='ParaScan', color=color_parascan, edgecolor='white')
    plt.bar(x + 1.5*width, grouped['elapsed_windowed'] + grouped['elapsed_banded'], width, label='QuickAffine (Bound)', color=color_windowed, edgecolor='white')
    plt.bar(x + 1.5*width, grouped['elapsed_banded'], width, label='QuickAffine (Align)', color=color_banded, edgecolor='white')

    plt.xlabel('Dataset and Sequence Length', fontsize=12)
    plt.ylabel('Elapsed Time (seconds)', fontsize=12)
    plt.title('Elapsed Time: Algorithms with Real Datasets', fontsize=14)
    plt.xticks(x, ['Illumina ~250b', 'Nanopore ~5000b'], fontsize=10)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'figure_10.png'))
    plt.close()
    return 1

def create_figure_12_plot(df, output_dir):
    """
    Create a plot with Memory consumed SWG vs Windowed vs Banded, with logarithmic scaling.
    The memory values are scaled logarithmically with base 10.
    """
    # Group by sequence length and compute mean memory consumed
    grouped = df.groupby('avg_query_length')[['memory_SWG', 'memory_windowed', 'memory_banded']].mean().reset_index()

    # Apply log10 transformation to memory consumed values
    grouped[['memory_SWG', 'memory_windowed', 'memory_banded']] = grouped[['memory_SWG', 'memory_windowed', 'memory_banded']].apply(lambda x: np.log10(x), axis=1)

    # Plot
    plt.figure(figsize=(10, 6))
    x = np.arange(len(grouped['avg_query_length']))
    width = 0.3

    plt.bar(x - width, grouped['memory_SWG'], width, label='SWG', color=color_swg, edgecolor='white')
    plt.bar(x, grouped['memory_windowed'], width, label='QuickAffine (Bound)', color=color_windowed, edgecolor='white')
    plt.bar(x + width, grouped['memory_banded'], width, label='QuickAffine (Align)', color=color_banded, edgecolor='white')

    plt.xlabel('Dataset and Sequence Length', fontsize=12)
    plt.ylabel('Log10 Memory Consumed (MB)', fontsize=12)
    plt.title('Log10 Memory Consumed: SWG vs Windowed vs Banded', fontsize=14)
    plt.xticks(x, ['Simulated 100b', 'Illumina ~250b', 'Simulated 1000b', 'Nanopore ~5000b', 'Simulated 10000b'], fontsize=10)
    plt.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'figure_12.png'))
    plt.close()
    return 1

def create_figure_13_table(df, output_dir):
    """
    Create a table summarizing computed cells for specified algorithms and datasets.
    Generates both a CSV and a LaTeX table.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Define algorithms and dataset types
    algorithms = ["SWG (real)", "QuickAffine (Bound)", "QuickAffine (Align)", "QuickAffine (Total)"]
    dataset_types = ['Sim 100b', 'Sim 1000b', 'Sim 10Kb', 'Illumina', 'Nanopore']

    # Map dataset types for filtering and final naming
    dataset_mapping = {
        'Sim 100b': {'dataset_type': 'Simulated', 'avg_query_length': 100},
        'Sim 1000b': {'dataset_type': 'Simulated', 'avg_query_length': 1000},
        'Sim 10Kb': {'dataset_type': 'Simulated', 'avg_query_length': 10000},
        'Illumina': {'dataset_type': 'Illumina'},
        'Nanopore': {'dataset_type': 'Nanopore'}
    }

    # Prepare the table data
    table_data = []
    for algorithm in algorithms:
        row = []
        for dataset, filters in dataset_mapping.items():
            subset = df[(df['window_size'] == 32) & (df['overlap_size'] == 8)]

            if 'avg_query_length' in filters:
                subset = subset[(subset['dataset_type'] == filters['dataset_type']) & (subset['avg_query_length'] == filters['avg_query_length'])]
            else:
                subset = subset[subset['dataset_type'] == filters['dataset_type']]

            if algorithm == "SWG (real)":
                total_cells = subset['cells_SWG'].sum() if not subset.empty else 0
            elif algorithm == "QuickAffine (Bound)":
                total_cells = subset['cells_windowed'].sum() if not subset.empty else 0
            elif algorithm == "QuickAffine (Align)":
                total_cells = subset['cells_banded'].sum() if not subset.empty else 0
            row.append(total_cells)

        if algorithm == "QuickAffine (Total)":
            row = [table_data[-2][i] + table_data[-1][i] for i in range(len(row))]
        table_data.append(row)

    # Create DataFrame for the table
    table_df = pd.DataFrame(table_data, columns=dataset_types, index=algorithms)

    # Save as CSV
    csv_path = os.path.join(output_dir, 'table_13.csv')
    table_df.to_csv(csv_path, float_format='%.0f')

    # Save as LaTeX table
    latex_path = os.path.join(output_dir, 'table_13.txt')
    with open(latex_path, 'w') as f:
        f.write("\\begin{table}[ht]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{l " + " ".join(["c"] * len(dataset_types)) + "}\\hline\n")
        f.write("\\rowcolor[HTML]{AFAFAF} \\textbf{Algorithm} & \\textbf{" + "} & \\textbf{".join(dataset_types) + "} \\\\ \\hline\n")
        for idx, row in table_df.iterrows():
            if idx == "QuickAffine (Total)":
                f.write(f"\\textbf{{{idx}}} & " + " & ".join([f"\\textbf{{{val:.0f}}}" for val in row]) + " \\\\ \\hline\n")
            else:
                f.write(f"{idx} & " + " & ".join([f"{val:.0f}" for val in row]) + " \\\\ \\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\caption{Computed Cells by Algorithm and Dataset}\n")
        f.write("\\label{tab:computed_cells_algorithm}\n")
        f.write("\\end{table}\n")

    return 1



def generate_plots(df, figures_created, output_dir, selected_figures):
    if 2 in selected_figures:
        figures_created[1] = create_figure_2_plot(df, output_dir)
    if 3 in selected_figures:
        figures_created[2] = create_figure_3_table(df, output_dir)
    if 4 in selected_figures:
        figures_created[3] = create_figure_4_plot(output_dir)
    if 7 in selected_figures:
        figures_created[6] = create_figure_7_plot(df, output_dir)
    if 9 in selected_figures:
        figures_created[8] = create_figure_9_plot(df, output_dir)
    if 10 in selected_figures:
        figures_created[9] = create_figure_10_plot(df, output_dir)
    if 12 in selected_figures:
        figures_created[11] = create_figure_12_plot(df, output_dir)
    if 13 in selected_figures:
        figures_created[12] = create_figure_13_table(df, output_dir)

def main():
    # Input and output file names
    input_filename = os.path.expanduser('~/Documents/GitHub/QuickAffine/res.out')
    
    # Read input file
    df = read_input_file(input_filename)
    
    num_of_figures = 13
    figures_created = [0] * num_of_figures

    # Ensure the output directory exists
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)

    # Request the figures to generate
    print("Enter the figure numbers to generate (e.g., 1,3,5-7):")
    user_input = input("> ")
    
    # Parse input into a list of selected figures
    selected_figures = []
    for part in user_input.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            selected_figures.extend(range(start, end + 1))
        else:
            selected_figures.append(int(part))
    
    selected_figures = list(set(selected_figures))  # Remove duplicates
    selected_figures.sort()

    # Generate plots
    generate_plots(df, figures_created, output_dir, selected_figures)
        
    print("Analysis complete. Output files generated:")
    for i in selected_figures:
        if figures_created[i - 1]:
            print(f"- figure_{i}.png")

if __name__ == "__main__":
    main()