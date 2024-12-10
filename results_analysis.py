import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def read_input_file(filename):
    """
    Read the input file with comma-separated values and return a pandas DataFrame
    """
    columns = [
        'window_size', 'overlap_size', 'penalty_set', 'avg_query_length', 'dataset_type',
        'elapsed_SWG', 'memory_SWG', 'cells_SWG', 'score_SWG',
        'elapsed_windowed', 'memory_windowed', 'cells_windowed', 'score_windowed',
        'elapsed_banded', 'memory_banded', 'cells_banded', 'score_banded'
    ]
    
    # Read CSV file with comma separator
    df = pd.read_csv(filename, header=None, names=columns, skipinitialspace=True)
    
    # Convert numeric columns 
    numeric_columns = [
        'window_size', 'overlap_size', 'avg_query_length',
        'elapsed_SWG', 'memory_SWG', 'cells_SWG', 'score_SWG',
        'elapsed_windowed', 'memory_windowed', 'cells_windowed', 'score_windowed',
        'elapsed_banded', 'memory_banded', 'cells_banded', 'score_banded'
    ]
    df[numeric_columns] = df[numeric_columns].apply(pd.to_numeric, errors='coerce')
    
    return df

def group_by_query_length(df, tolerance=0.1):
    """
    Group data by query length with a tolerance for similarity
    """
    def find_group(length, existing_groups):
        for group in existing_groups:
            if abs(length - group) / group <= tolerance:
                return group
        return length
    
    # Create a copy of the dataframe to avoid modifying the original
    df_grouped = df.copy()
    
    # Group query lengths
    df_grouped['query_length_group'] = df_grouped['avg_query_length'].apply(
        lambda x: find_group(x, df_grouped['avg_query_length'].unique())
    )
    
    return df_grouped

def create_comprehensive_plot(df, metric_columns, metric_name, output_filename):
    """
    Create a comprehensive plot with grouped bars for different algorithms
    
    Args:
    df (pandas.DataFrame): Input dataframe
    metric_columns (list): List of column names for different algorithms
    metric_name (str): Name of the metric for plot title
    output_filename (str): Output filename for the plot
    """
    # Group by query length
    df_grouped = group_by_query_length(df)
    
    # Aggregate data
    aggregated_data = []
    for column in metric_columns:
        agg = df_grouped.groupby(['query_length_group'])[column].mean()
        aggregated_data.append(agg)
    
    # Combine aggregated data
    plot_df = pd.DataFrame({col: data for col, data in zip(metric_columns, aggregated_data)})
    plot_df.reset_index(inplace=True)
    
    # Plotting
    plt.figure(figsize=(14, 7))
    
    # Number of algorithms
    num_algorithms = len(metric_columns)
    
    # Width of a bar 
    bar_width = 0.25
    
    # Position of bars on X axis
    r = np.arange(len(plot_df))
    
    # Plot bars for each algorithm
    plt.bar(r - bar_width, plot_df[metric_columns[0]], 
            width=bar_width, label='SWG', 
            edgecolor='white')
    plt.bar(r, plot_df[metric_columns[1]], 
            width=bar_width, label='Windowed', 
            edgecolor='white')
    plt.bar(r, plot_df[metric_columns[2]], 
            width=bar_width, label='Banded', 
            edgecolor='white')
    
    # Customize the plot
    plt.xlabel('Query Length', fontsize=12)
    plt.ylabel(metric_name, fontsize=12)
    plt.title(f'{metric_name} Comparison Across Algorithms', fontsize=15)
    
    # X-axis ticks
    plt.xticks(r, [f'{int(x)}' for x in plot_df['query_length_group']])
    
    # Legend
    plt.legend(title='Algorithm')
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_filename)
    plt.close()

def compare_scores(df):
    """
    Compare SWG scores with Windowed scores and calculate percentage difference
    
    Args:
    df (pandas.DataFrame): Input dataframe
    
    Returns:
    pandas.DataFrame: Dataframe with score comparison metrics
    """
    # Create a copy of the dataframe
    comparison_df = df.copy()
    
    # Calculate percentage difference for each row
    # Percentage difference = (Windowed Score - SWG Score) / SWG Score * 100
    comparison_df['score_percentage_difference'] = (
        (comparison_df['score_windowed'] - comparison_df['score_SWG']) / 
        comparison_df['score_SWG'] * 100
    )
    
    # Group by query length with tolerance
    def find_group(length, existing_groups, tolerance=0.1):
        for group in existing_groups:
            if abs(length - group) / group <= tolerance:
                return group
        return length
    
    comparison_df['query_length_group'] = comparison_df['avg_query_length'].apply(
        lambda x: find_group(x, comparison_df['avg_query_length'].unique())
    )
    
    # Aggregate percentage differences by query length group
    error_summary = comparison_df.groupby('query_length_group').agg({
        'score_percentage_difference': ['mean', 'std', 'min', 'max'],
        'score_SWG': 'mean',
        'score_windowed': 'mean',
        'score_banded': 'mean'
    }).reset_index()
    
    # Flatten multi-level column names
    error_summary.columns = [
        'query_length_group', 
        'percentage_diff_mean', 
        'percentage_diff_std', 
        'percentage_diff_min', 
        'percentage_diff_max',
        'swg_score_mean',
        'windowed_score_mean',
        'banded_score_mean'
    ]
    
    return comparison_df, error_summary

def plot_score_comparison(df):
    """
    Create a comprehensive plot comparing SWG and Windowed scores
    
    Args:
    df (pandas.DataFrame): Input dataframe
    """
    # Get detailed comparison and summary data
    detailed_df, summary_df = compare_scores(df)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # Plot percentage differences
    plot_data = detailed_df.groupby('query_length_group')['score_percentage_difference']
    
    # Boxplot of percentage differences
    ax1.boxplot([
        detailed_df[detailed_df['query_length_group'] == group]['score_percentage_difference'] 
        for group in sorted(detailed_df['query_length_group'].unique())
    ])
    
    ax1.set_xlabel('Query Length Group')
    ax1.set_ylabel('Percentage Difference (%)')
    ax1.set_title('Percentage Difference in Windowed Score Compared to SWG')
    ax1.set_xticklabels([f'{int(group)}' for group in sorted(detailed_df['query_length_group'].unique())])
    
    # Bar plot of average scores
    x = np.arange(len(summary_df))
    bar_width = 0.35
    
    ax2.bar(x - bar_width/2, summary_df['swg_score_mean'], 
            width=bar_width, label='SWG Score', color='blue', edgecolor='white')
    ax2.bar(x + bar_width/2, summary_df['windowed_score_mean'], 
            width=bar_width, label='Windowed Score', color='orange', edgecolor='white')
    ax2.bar(x + bar_width/2, summary_df['banded_score_mean'], 
            width=bar_width, label='Banded Score', color='green', edgecolor='white')
    
    ax2.set_xlabel('Query Length')
    ax2.set_ylabel('Average Score')
    ax2.set_title('Comparison of Average SWG and Windowed Scores')
    ax2.set_xticks(x)
    ax2.set_xticklabels([f'{int(ql)}' for ql in summary_df['query_length_group']])
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig('score_comparison.png')
    plt.close()
    
    # Print detailed summary
    print("Score Comparison Summary:")
    print(summary_df.to_string(index=False))
    
    return detailed_df, summary_df

def generate_plots(df):
    """
    Generate comprehensive plots
    """
    # Elapsed Time Plot
    create_comprehensive_plot(
        df, 
        ['elapsed_SWG', 'elapsed_windowed', 'elapsed_banded'], 
        'Elapsed Time (milli-seconds)', 
        'elapsed_time_comparison.png'
    )
    
    # Cells Computation Plot
    create_comprehensive_plot(
        df, 
        ['cells_SWG', 'cells_windowed', 'cells_banded'], 
        'Cells Computation (millions)', 
        'cells_computation_comparison.png'
    )
    
    # Memory Usage Plot
    create_comprehensive_plot(
        df, 
        ['memory_SWG', 'memory_windowed', 'memory_banded'], 
        'Memory Usage (MB)', 
        'memory_usage_comparison.png'
    )
    plot_score_comparison(df)

def main():
    # Input and output file names
    input_filename = os.path.expanduser('~/Documents/GitHub/QuickAffine/res.out')
    
    # Read input file
    df = read_input_file(input_filename)
        
    # Generate plots
    generate_plots(df)
        
    print("Analysis complete. Output files generated:")
    print("- elapsed_time_comparison.png")
    print("- cells_computation_comparison.png")
    print("- memory_usage_comparison.png")

if __name__ == "__main__":
    main()