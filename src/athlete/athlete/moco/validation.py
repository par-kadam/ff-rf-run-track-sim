import os
import glob
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from collections import defaultdict
import re
import xml.etree.ElementTree as ET
import opensim as osim
from datetime import datetime
from collections import Counter
from typing import Optional, Dict, List, Tuple


def analyze_gait_cycles_comprehensive(data_directory, coordinate_name=None, output_dir=None, 
                                    compute_errors=True):
    """
    Top-level function that performs comprehensive gait cycle analysis across all trial types.
    
    Parameters:
    -----------
    data_directory : str
        Path to directory containing ik_move.sto files
    coordinate_name : str, optional
        Specific coordinate to analyze and plot. If None, will analyze all coordinates
    output_dir : str, optional
        Directory to save plots. If None, plots will be shown but not saved
    compute_errors : bool
        Whether to load error metrics from _ik_errors.log files
        
    Returns:
    --------
    dict
        Dictionary containing all processed data and statistics for each trial type
    """
    print(f"Starting comprehensive gait cycle analysis in: {data_directory}")
    
    # Step 1: Categorize all files and load error metrics from log files
    categorized_files = categorize_ik_files(data_directory, compute_errors=compute_errors)
    
    # Step 2: Process all gait cycles for each category
    all_trial_data = {}
    
    for trial_type, file_data in categorized_files.items():
        if not file_data['files']:
            print(f"No files found for trial type: {trial_type}")
            continue
            
        print(f"Processing {len(file_data['files'])} files for trial type: {trial_type}")
        trial_data = process_trial_category(file_data['files'], trial_type, 
                                          error_metrics=file_data['error_metrics'])
        all_trial_data[trial_type] = trial_data
    
    # Step 3: Generate plots for specific coordinate or all coordinates
    if coordinate_name:
        plot_coordinate_across_trials(all_trial_data, coordinate_name, output_dir)
    else:
        # Find all available coordinates
        all_coordinates = set()
        for trial_data in all_trial_data.values():
            if trial_data['rms_data']:
                all_coordinates.update(trial_data['rms_data'].keys())
        
        print(f"Found {len(all_coordinates)} coordinates to analyze")
        for coord in sorted(all_coordinates):
            plot_coordinate_across_trials(all_trial_data, coord, output_dir)
    
    return all_trial_data

def categorize_ik_files(data_directory, compute_errors=True):
    """
    Level 1: Categorizes all ik.sto files based on naming string criteria and loads error XML files.
    
    Parameters:
    -----------
    data_directory : str
        Path to directory containing ik.sto files and ik_errors.xml files
        
    Returns:
    --------
    dict
        Dictionary with trial types as keys and lists of file paths as values
        Trial types: 'ft_02', 'ft_03', 'rt_02', 'rt_03', 'fc_02', 'fc_03', 'rc_02', 'rc_03'
    """
    # Find all ik.sto files and XML error files recursively
    ik_files = []
    xml_files = []
    
    # Search recursively for files
    for root, dirs, files in os.walk(data_directory):
        for file in files:
            if file.endswith('_ik_move.sto'):
                ik_files.append(os.path.join(root, file))
            if compute_errors and file.endswith('_ik_errors.xml'):  # FIXED: was .log, now .xml
                xml_files.append(os.path.join(root, file))
    
    print(f"Found {len(ik_files)} _ik_move.sto files")
    if compute_errors:
        print(f"Found {len(xml_files)} _ik_errors.xml files")
    
    # Initialize categories with both files and error_metrics lists
    categories = {
        'ft_02': {'files': [], 'error_metrics': []},
        'ft_03': {'files': [], 'error_metrics': []},
        'rt_02': {'files': [], 'error_metrics': []},
        'rt_03': {'files': [], 'error_metrics': []},
        'fc_02': {'files': [], 'error_metrics': []},
        'fc_03': {'files': [], 'error_metrics': []},
        'rc_02': {'files': [], 'error_metrics': []},
        'rc_03': {'files': [], 'error_metrics': []}
    }
    
    # Match STO files with XML files if error computation is enabled
    file_to_xml_map = {}
    if compute_errors:
        file_to_xml_map = match_sto_files_with_error_xml(ik_files, xml_files)
        print(f"Successfully matched {len(file_to_xml_map)} .sto files with .xml error files")
    
    for file_path in ik_files:
        filename = os.path.basename(file_path)
        
        # Extract stride type (ft, rt, fc, rc)
        stride_type = None
        if 'stri' in filename:
            # Find the two characters immediately after 'stri'
            stri_index = filename.find('stri')
            if stri_index != -1 and stri_index + 6 <= len(filename):
                stride_chars = filename[stri_index + 4:stri_index + 6]
                if stride_chars in ['ft', 'rt', 'fc', 'rc']:
                    stride_type = stride_chars
        
        # Extract speed (02 or 03)
        speed_type = None
        if 'speed' in filename:
            # Find the two characters after 'speed'
            speed_index = filename.find('speed')
            if speed_index != -1:
                # Look for speed pattern (could be speedp0200, speed02000, etc.)
                speed_section = filename[speed_index + 5:]
                # Extract first two digits after 'speed'
                speed_match = re.search(r'(\d{2})', speed_section)
                if speed_match:
                    speed_digits = speed_match.group(1)
                    if speed_digits in ['02', '03']:
                        speed_type = speed_digits
        
        # Categorize file
        if stride_type and speed_type:
            category_key = f"{stride_type}_{speed_type}"
            categories[category_key]['files'].append(file_path)
            
            # Add corresponding error metrics if available
            if compute_errors and file_path in file_to_xml_map:
                xml_file_path = file_to_xml_map[file_path]
                error_metrics = parse_ik_error_xml(xml_file_path)  # NEW: Parse XML file
                categories[category_key]['error_metrics'].append(error_metrics)
                
                if error_metrics:
                    print(f"Categorized {filename} as {category_key} with error data (RMS: {error_metrics['mean_rms_error']:.4f}, Max: {error_metrics['mean_max_error']:.4f})")
                else:
                    print(f"Categorized {filename} as {category_key} (XML parsing failed)")
            else:
                if compute_errors:
                    categories[category_key]['error_metrics'].append(None)
                    print(f"Categorized {filename} as {category_key} (no matching XML file found)")
                else:
                    print(f"Categorized {filename} as {category_key}")
        else:
            print(f"Could not categorize {filename} (stride: {stride_type}, speed: {speed_type})")
    
    # Print summary
    for category, data in categories.items():
        files_count = len(data['files'])
        if compute_errors:
            valid_errors = sum(1 for e in data['error_metrics'] if e is not None)
            print(f"Category {category}: {files_count} files, {valid_errors} with valid error data")
        else:
            print(f"Category {category}: {files_count} files")
    
    return categories


def parse_ik_error_xml(xml_file_path):
    """Parse IK error XML file to extract summary statistics."""
    try:
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
        
        metadata = root.find('Metadata')
        if metadata is None:
            return None
            
        summary = metadata.find('Summary')
        if summary is None:
            return None
        
        summary_attrs = summary.attrib
        trial_info = metadata.find('TrialInfo')
        trial_attrs = trial_info.attrib if trial_info is not None else {}
        
        result = {
            'mean_rms_error': float(summary_attrs.get('rms_error_mean', 'nan')),
            'mean_max_error': float(summary_attrs.get('max_error_mean', 'nan')),
            'rms_error_std': float(summary_attrs.get('rms_error_std', 'nan')),
            'max_error_std': float(summary_attrs.get('max_error_std', 'nan')),
            'total_frames': int(summary_attrs.get('total_frames', '0')),
            'trial_name': trial_attrs.get('name', os.path.basename(xml_file_path)),
            'xml_file_path': xml_file_path
        }
        
        return result
        
    except Exception as e:
        print(f"Error parsing XML file {xml_file_path}: {e}")
        return None

def match_sto_files_with_error_xml(sto_files, xml_files):
    """Match .sto files with their corresponding _ik_errors.xml files."""
    matches = {}
    
    for sto_file in sto_files:
        sto_basename = os.path.basename(sto_file)
        
        if sto_basename.endswith('_ik_move.sto'):
            base_name = sto_basename[:-12]  # Remove '_ik_move.sto'
        else:
            continue
            
        for xml_file in xml_files:
            xml_basename = os.path.basename(xml_file)
            
            if xml_basename.startswith(base_name) and xml_basename.endswith('_ik_errors.xml'):
                matches[sto_file] = xml_file
                break
    
    return matches

def normalize_gait_cycle(file_path):
    """
    Level 2: Normalizes time series data from a single gait cycle file to 100 points.
    
    Parameters:
    -----------
    file_path : str
        Path to the ik_move.sto file
        
    Returns:
    --------
    dict
        Dictionary with coordinate names as keys and normalized arrays (100 points) as values
    """
    try:
        # Load the STO file
        df = load_sto_file(file_path)
        
        if df is None or df.empty:
            print(f"Warning: Could not load data from {file_path}")
            return {}
        
        # Get time column
        time_col = None
        for col in df.columns:
            if col.lower() == 'time':
                time_col = col
                break
        
        if time_col is None:
            print(f"Warning: No time column found in {file_path}")
            return {}
        
        time_data = df[time_col].values
        
        # Create normalized time points (0 to 100% of gait cycle)
        normalized_time = np.linspace(0, 100, 100)
        original_time_percent = np.linspace(0, 100, len(time_data))
        
        # Normalize all coordinate columns
        normalized_data = {}
        
        for column in df.columns:
            if column != time_col:
                try:
                    # Get coordinate data
                    coord_data = df[column].values
                    
                    # Check for valid data
                    if len(coord_data) < 2 or np.all(np.isnan(coord_data)):
                        continue
                    
                    # Interpolate to 100 points
                    # Remove NaN values for interpolation
                    valid_indices = ~np.isnan(coord_data)
                    if np.sum(valid_indices) < 2:
                        continue
                    
                    valid_time = original_time_percent[valid_indices]
                    valid_data = coord_data[valid_indices]
                    
                    # Create interpolation function
                    interp_func = interp1d(valid_time, valid_data, kind='linear', 
                                         bounds_error=False, fill_value='extrapolate')
                    
                    # Interpolate to normalized time points
                    normalized_coord = interp_func(normalized_time)
                    
                    # Clean coordinate name
                    clean_coord_name = extract_coordinate_name(column)
                    normalized_data[clean_coord_name] = normalized_coord
                    
                except Exception as e:
                    print(f"Warning: Could not normalize coordinate {column} in {file_path}: {e}")
                    continue
        
        return normalized_data
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return {}

def process_trial_category(file_list, trial_type, error_metrics=None):
    """
    Level 3: Processes all files in a trial category to compute RMS and standard deviation.
    
    Parameters:
    -----------
    file_list : list
        List of file paths for this trial category
    trial_type : str
        Name of the trial type (e.g., 'ft_02')
    error_metrics : list, optional
        List of error metric dictionaries corresponding to files
        
    Returns:
    --------
    dict
        Dictionary containing RMS data, standard deviation data, raw data, and error metrics
    """
    print(f"Processing trial category: {trial_type} with {len(file_list)} files")
    
    # Store all normalized data
    all_cycles_data = defaultdict(list)  # coordinate_name -> list of 100-point arrays
    
    # Process each file
    for file_path in file_list:
        normalized_data = normalize_gait_cycle(file_path)
        
        for coord_name, coord_data in normalized_data.items():
            all_cycles_data[coord_name].append(coord_data)
    
    # Compute statistics for each coordinate
    rms_data = {}
    std_data = {}
    mean_data = {}
    raw_data = {}
    
    for coord_name, cycles_list in all_cycles_data.items():
        if not cycles_list:
            continue
            
        # Convert to numpy array (cycles x time_points)
        cycles_array = np.array(cycles_list)
        raw_data[coord_name] = cycles_array
        
        # Compute statistics at each time point
        mean_values = np.mean(cycles_array, axis=0)
        std_values = np.std(cycles_array, axis=0)
        rms_values = np.sqrt(np.mean(cycles_array**2, axis=0))
        
        mean_data[coord_name] = mean_values
        std_data[coord_name] = std_values
        rms_data[coord_name] = rms_values
        
        print(f"  Processed {coord_name}: {len(cycles_list)} cycles")
    
    # Process error metrics if provided
    processed_error_metrics = None
    if error_metrics:
        valid_metrics = [m for m in error_metrics if m is not None]
        if valid_metrics:
            mean_rms_errors = [m['mean_rms_error'] for m in valid_metrics if not np.isnan(m['mean_rms_error'])]
            mean_max_errors = [m['mean_max_error'] for m in valid_metrics if not np.isnan(m['mean_max_error'])]
            
            processed_error_metrics = {
                'individual_metrics': valid_metrics,
                'mean_rms_errors': mean_rms_errors,
                'mean_max_errors': mean_max_errors,
                'avg_mean_rms': np.mean(mean_rms_errors) if mean_rms_errors else np.nan,
                'avg_mean_max': np.mean(mean_max_errors) if mean_max_errors else np.nan
            }
    
    return {
        'trial_type': trial_type,
        'num_files': len(file_list),
        'file_list': file_list,
        'rms_data': rms_data,
        'std_data': std_data,
        'mean_data': mean_data,
        'raw_data': raw_data,
        'error_metrics': processed_error_metrics,
        'time_points': np.linspace(0, 100, 100)
    }

def plot_error_box_plots(all_trial_data, ax_rms, ax_max):
    """
    Create box and whisker plots for RMS and max errors across trial types.
    
    Parameters:
    -----------
    all_trial_data : dict
        Dictionary containing processed data for all trial types
    ax_rms : matplotlib.axes.Axes
        Axis for RMS error box plot
    ax_max : matplotlib.axes.Axes
        Axis for max error box plot
    """
    # Collect error data for each trial type
    trial_types = ['ft_02', 'ft_03', 'rt_02', 'rt_03', 'fc_02', 'fc_03', 'rc_02', 'rc_03']
    rms_data = []
    max_data = []
    labels = []
    
    for trial_type in trial_types:
        if trial_type in all_trial_data and all_trial_data[trial_type]['error_metrics']:
            error_metrics = all_trial_data[trial_type]['error_metrics']
            
            if error_metrics['mean_rms_errors']:
                rms_data.append(error_metrics['mean_rms_errors'])
                labels.append(trial_type)
            
            if error_metrics['mean_max_errors']:
                max_data.append(error_metrics['mean_max_errors'])
    
    # Plot RMS error box plot
    if rms_data:
        ax_rms.boxplot(rms_data, labels=labels)
        ax_rms.set_title('Mean RMS Errors by Trial Type')
        ax_rms.set_xlabel('Trial Type')
        ax_rms.set_ylabel('Mean RMS Error')
        ax_rms.tick_params(axis='x', rotation=45)
        ax_rms.grid(True, alpha=0.3)
    else:
        ax_rms.text(0.5, 0.5, 'No RMS Error Data Available', 
                   ha='center', va='center', transform=ax_rms.transAxes)
        ax_rms.set_title('Mean RMS Errors by Trial Type')
    
    # Plot max error box plot
    if max_data:
        ax_max.boxplot(max_data, labels=labels)
        ax_max.set_title('Mean Max Errors by Trial Type')
        ax_max.set_xlabel('Trial Type')
        ax_max.set_ylabel('Mean Max Error')
        ax_max.tick_params(axis='x', rotation=45)
        ax_max.grid(True, alpha=0.3)
    else:
        ax_max.text(0.5, 0.5, 'No Max Error Data Available', 
                   ha='center', va='center', transform=ax_max.transAxes)
        ax_max.set_title('Mean Max Errors by Trial Type')



def plot_coordinate_across_trials(all_trial_data, coordinate_name, output_dir=None):
    """
    Level 4: Plots mean time series with shaded standard deviation regions and error box plots.
    Creates 4 subplots grouped by speed and first character of stride type, plus error analysis.
    Each subplot compares trajectories differing only by second character (t vs c).
    
    Parameters:
    -----------
    all_trial_data : dict
        Dictionary containing processed data for all trial types
    coordinate_name : str
        Name of the coordinate to plot
    output_dir : str, optional
        Directory to save plots
    """
    print(f"Plotting coordinate: {coordinate_name}")
    
    # Create figure with subplots (2x3 layout to include error plots)
    fig = plt.figure(figsize=(18, 12))
    
    # Create subplot layout: 2x3 grid
    ax1 = plt.subplot(2, 3, 1)  # f-type, Speed 02
    ax2 = plt.subplot(2, 3, 2)  # f-type, Speed 03
    ax3 = plt.subplot(2, 3, 3)  # RMS Error Box Plot
    ax4 = plt.subplot(2, 3, 4)  # r-type, Speed 02
    ax5 = plt.subplot(2, 3, 5)  # r-type, Speed 03
    ax6 = plt.subplot(2, 3, 6)  # Max Error Box Plot
    
    fig.suptitle(f'Gait Cycle Analysis: {coordinate_name}', fontsize=16)
    
    time_points = np.linspace(0, 100, 100)
    
    # Define subplot groupings and colors
    # Each subplot compares 't' vs 'c' conditions
    subplot_configs = [
        {'ax': ax1, 'title': 'f-type, Speed 02', 'trials': ['ft_02', 'fc_02'], 'colors': ['blue', 'green']},
        {'ax': ax2, 'title': 'f-type, Speed 03', 'trials': ['ft_03', 'fc_03'], 'colors': ['lightblue', 'lightgreen']},
        {'ax': ax4, 'title': 'r-type, Speed 02', 'trials': ['rt_02', 'rc_02'], 'colors': ['red', 'purple']},
        {'ax': ax5, 'title': 'r-type, Speed 03', 'trials': ['rt_03', 'rc_03'], 'colors': ['lightcoral', 'plum']}
    ]
    
    plotted_any = False
    
    # Plot time series data
    for config in subplot_configs:
        ax = config['ax']
        ax.set_title(config['title'])
        ax.set_xlabel('Gait Cycle (%)')
        ax.set_ylabel('Coordinate Value')
        ax.grid(True, alpha=0.3)
        
        # Plot each trial type in this subplot
        for i, trial_type in enumerate(config['trials']):
            if trial_type in all_trial_data and coordinate_name in all_trial_data[trial_type]['mean_data']:
                trial_data = all_trial_data[trial_type]
                color = config['colors'][i]
                
                mean_values = trial_data['mean_data'][coordinate_name]
                std_values = trial_data['std_data'][coordinate_name]
                
                # Determine label based on second character (t vs c)
                second_char = trial_type[1]  # 't' or 'c'
                label = f'{second_char}-type (n={trial_data["num_files"]})'
                
                # Plot mean line
                ax.plot(time_points, mean_values, color=color, linewidth=2, label=label)
                
                # Add shaded standard deviation region
                ax.fill_between(time_points, 
                               mean_values - std_values,
                               mean_values + std_values,
                               color=color, alpha=0.3)
                
                plotted_any = True
        
        # Add legend to each subplot
        ax.legend(loc='best')
    
    # Create error box plots
    plot_error_box_plots(all_trial_data, ax3, ax6)
    
    if not plotted_any:
        print(f"Warning: No data found for coordinate '{coordinate_name}'")
        plt.close(fig)
        return
    
    plt.tight_layout()
    
    # Save or show plot
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        filename = f"gait_analysis_{coordinate_name.replace('/', '_')}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Plot saved: {filepath}")
        plt.close(fig)
    else:
        plt.show()


def load_sto_file(filename):
    """
    Load STO file and return pandas DataFrame.
    Reuses the existing function from retrieve_results.py
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the line with 'endheader'
        endheader_idx = None
        for i, line in enumerate(lines):
            if 'endheader' in line.lower():
                endheader_idx = i
                break
        
        if endheader_idx is None:
            return None
        
        # Headers are in the line after 'endheader'
        header_line = lines[endheader_idx + 1]
        headers = [h.strip() for h in header_line.split('\t') if h.strip()]
        
        # Data starts two lines after 'endheader'
        data_start_idx = endheader_idx + 2
        
        # Parse the data
        data_dict = {header: [] for header in headers}
        
        for i in range(data_start_idx, len(lines)):
            line = lines[i].strip()
            if not line:  # Skip empty lines
                continue
            
            values = line.split('\t')
            if len(values) >= len(headers):
                for j, header in enumerate(headers):
                    try:
                        data_dict[header].append(float(values[j]))
                    except (ValueError, IndexError):
                        data_dict[header].append(np.nan)
        
        return pd.DataFrame(data_dict)
        
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def extract_coordinate_name(column_name):
    """
    Extract coordinate name from different column naming conventions.
    Reuses the existing function from retrieve_results.py
    """
    # Remove common prefixes and suffixes
    name = column_name.replace('/jointset/', '').replace('/value', '').replace('/speed', '')
    
    # Handle path-style naming (e.g., '/jointset/ground_pelvis/pelvis_tx/value')
    if '/' in name:
        parts = name.split('/')
        # Usually the coordinate name is the last part or second to last
        if len(parts) >= 2:
            return parts[-1]  # Take the last part as coordinate name
    
    return name

# Example usage functions:
def analyze_single_coordinate(data_directory, coordinate_name, output_dir=None, 
                            compute_errors=True):
    """
    Convenience function to analyze a single coordinate across all trial types.
    """
    return analyze_gait_cycles_comprehensive(data_directory, coordinate_name, output_dir, 
                                           compute_errors)

def analyze_all_coordinates(data_directory, output_dir=None, compute_errors=True):
    """
    Convenience function to analyze all available coordinates.
    """
    return analyze_gait_cycles_comprehensive(data_directory, None, output_dir, 
                                           compute_errors)

def get_available_coordinates(data_directory):
    """
    Returns a list of all available coordinates in the dataset.
    """
    categorized_files = categorize_ik_files(data_directory, compute_errors=False)
    all_coordinates = set()
    
    for trial_type, file_data in categorized_files.items():
        if file_data['files']:
            # Sample one file to get coordinate names
            sample_data = normalize_gait_cycle(file_data['files'][0])
            all_coordinates.update(sample_data.keys())
    
    return sorted(list(all_coordinates))

def get_error_summary(data_directory):
    """
    Get a summary of error metrics across all trial types from existing log files.
    
    Parameters:
    -----------
    data_directory : str
        Path to directory containing ik_move.sto files and corresponding _ik_errors.log files
        
    Returns:
    --------
    dict
        Dictionary containing error summary statistics for each trial type
    """
    categorized_files = categorize_ik_files(data_directory, compute_errors=True)
    
    error_summary = {}
    
    for trial_type, file_data in categorized_files.items():
        if file_data['error_metrics']:
            valid_metrics = [m for m in file_data['error_metrics'] if m is not None]
            
            if valid_metrics:
                mean_rms_errors = [m['mean_rms_error'] for m in valid_metrics 
                                 if not np.isnan(m['mean_rms_error'])]
                mean_max_errors = [m['mean_max_error'] for m in valid_metrics 
                                 if not np.isnan(m['mean_max_error'])]
                
                error_summary[trial_type] = {
                    'num_files': len(file_data['files']),
                    'num_valid_analyses': len(valid_metrics),
                    'mean_rms_errors': mean_rms_errors,
                    'mean_max_errors': mean_max_errors,
                    'avg_rms_error': np.mean(mean_rms_errors) if mean_rms_errors else np.nan,
                    'std_rms_error': np.std(mean_rms_errors) if mean_rms_errors else np.nan,
                    'avg_max_error': np.mean(mean_max_errors) if mean_max_errors else np.nan,
                    'std_max_error': np.std(mean_max_errors) if mean_max_errors else np.nan
                }
            else:
                error_summary[trial_type] = {
                    'num_files': len(file_data['files']),
                    'num_valid_analyses': 0,
                    'error': 'No valid error log files found'
                }
    
    return error_summary

def get_negative_muscle_forces_setup(model, solution, output_path=None):
    """
    Get the setup for negative muscle forces from the model and solution.
    
    Parameters:
    -----------
    model : Model
        The OpenSim model object
    solution : Solution
        The solution object containing results
    
    Returns:
    --------
    dict
        Dictionary with negative muscle forces setup
    """
    model.initSystem()
    outputs = osim.analyze(model, solution.exportToStatesTable(), 
                           solution.exportToControlsTable(), ['.*\|tendon_force'])
    def simtkmin(simtkvec):
        lowest = np.inf
        for i in range(simtkvec.size()):
            if simtkvec[i] < lowest:
                lowest = simtkvec[i]
        return lowest
    
    neg_forces = list()
    muscle_names = list()
    neg_forces_data = []

    with open(output_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Muscle Name', 'Negative Force (F_iso)']) 
        print(f'Negative force report for {model.getName()}:')

        for imusc in range(model.getMuscles().getSize()):
            muscle = model.getMuscles().get(imusc)
            max_iso = muscle.get_max_isometric_force()
            tendon_force = outputs.getDependentColumn(
                muscle.getAbsolutePathString() + '|tendon_force')
            neg = simtkmin(tendon_force) / max_iso
            if neg < 0:
                neg_forces.append(neg)
                muscle_names.append(muscle.getName())
                # for the console printout
                print(f' {muscle.getName()}: {neg} F_iso')
                # for the CSV file
                neg_forces_data.append((muscle.getName(), neg))
            #write the whole array into the CSV file
            writer.writerows(neg_forces_data)

        if len(neg_forces) == 0:
            print('No negative muscle forces found.')
        else:
            imin = np.argmin(neg_forces)
            print(f'Largest negative force: {muscle_names[imin]} '
                    f'({neg_forces[imin]} F_iso)')
        return min([0] + neg_forces)
    
def get_negative_muscle_forces(model, solution, output_path=None):
    print(f'Negative force report for ')
    model.initSystem()
    return get_negative_muscle_forces_setup(model, solution, output_path=output_path)

def calc_muscle_mechanics(model, solution):
    """
    Calculate muscle mechanics from the model and solution.
    
    Parameters:
    -----------
    model : Model
        The OpenSim model object
    solution : Solution
        The solution moco trajcetory containing results
    
    Returns:
    --------
    dict
        Dictionary with muscle mechanics data
    """
    output_list = list()
    for output in ['normalized_fiber_length', 'normalized_fiber_velocity', 
                       'active_fiber_force', 'passive_fiber_force', 
                       'tendon_force', 'activation', 'cos_pennation_angle', 
                       'active_force_length_multiplier', 
                       'force_velocity_multiplier', 'passive_force_multiplier',
                       'tendon_length', 'tendon_strain']:
        for imusc in range(model.getMuscles().getSize()):
            musc = model.updMuscles().get(imusc)
            output_list.append(f'.*{musc.getName()}.*\|{output}')
    for output in ['length', 'lengthening_speed']:
        output_list.append(f'.*\|{output}')

    outputs = osim.analyze(model, solution.exportToStatesTable(),
                           solution.exportToControlsTable(), output_list)
    print(f'Analyzing muscle mechanics for {model.getName()}')
    print(f'Found {outputs.getNumRows()} rows and {outputs.getNumColumns()} columns in output table')
    return outputs        
    
    # Process outputs to extract muscle mechanics data
    # This is a placeholder for actual implementation
    # You would typically extract forces, lengths, velocities, etc. here

def parse_ik_errors_batch(directory: str, 
                    subject_id: Optional[str] = None,
                    output_filename: Optional[str] = None) -> Dict:
    """
    Parse multiple IK errors XML files and compute batch statistics.
    
    Parameters:
    -----------
    directory : str
        Directory containing IK errors XML files
    subject_id : str, optional
        Subject ID (2 characters) to filter files. If specified, only processes
        files with matching subject ID in filename (subxx format)
    output_filename : str, optional
        Name of output report file. If None, auto-generates based on subject_id
        
    Returns:
    --------
    dict
        Dictionary containing all computed statistics and metadata
    """
    
    def parse_single_xml(xml_file: str) -> Tuple[pd.DataFrame, Dict]:
        """Parse a single IK errors XML file"""
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            
            # Extract metadata
            metadata = {'filename': os.path.basename(xml_file)}
            trial_info = root.find('.//TrialInfo')
            if trial_info is not None:
                metadata.update(trial_info.attrib)
            
            summary = root.find('.//Summary')
            if summary is not None:
                metadata.update(summary.attrib)
            
            # Extract frame data
            frames = []
            for frame in root.findall('.//Frame'):
                frame_data = {
                    'frame_number': int(frame.get('number')),
                    'time': float(frame.get('time')),
                    'rms_error': float(frame.get('rms_error')),
                    'max_error': float(frame.get('max_error')),
                    'max_error_marker': frame.get('max_error_marker')
                }
                frames.append(frame_data)
            
            df = pd.DataFrame(frames)
            return df, metadata
            
        except Exception as e:
            print(f"Error parsing {xml_file}: {e}")
            return None, None
    
    def extract_subject_from_filename(filename: str) -> Optional[str]:
        """Extract subject ID from filename (subxx format)"""
        match = re.search(r'sub(\w{2})', filename)
        return match.group(1) if match else None
    
    # Find all XML files in directory
    if not os.path.exists(directory):
        raise ValueError(f"Directory does not exist: {directory}")
    
    xml_files = []
    for filename in os.listdir(directory):
        if filename.endswith('_ik_errors.xml'):
            # Check subject ID filter if specified
            if subject_id is not None:
                file_subject = extract_subject_from_filename(filename)
                if file_subject != subject_id:
                    continue
            
            filepath = os.path.join(directory, filename)
            if os.path.isfile(filepath):
                xml_files.append(filepath)
    
    if not xml_files:
        error_msg = f"No IK errors XML files found in {directory}"
        if subject_id:
            error_msg += f" for subject {subject_id}"
        raise ValueError(error_msg)
    
    print(f"Found {len(xml_files)} IK errors XML files to process:")
    for f in xml_files:
        print(f"  {os.path.basename(f)}")
    
    # Parse all files and collect data
    all_frames = []
    all_metadata = []
    file_statistics = []
    marker_error_counts = Counter()
    
    for xml_file in xml_files:
        df, metadata = parse_single_xml(xml_file)
        
        if df is not None and not df.empty:
            # Add filename to each frame for tracking
            df['source_file'] = os.path.basename(xml_file)
            all_frames.append(df)
            all_metadata.append(metadata)
            
            # Compute per-file statistics
            file_stats = {
                'filename': os.path.basename(xml_file),
                'total_frames': len(df),
                'mean_rms_error': df['rms_error'].mean(),
                'std_rms_error': df['rms_error'].std(),
                'mean_max_error': df['max_error'].mean(),
                'std_max_error': df['max_error'].std(),
                'time_range': f"{df['time'].min():.3f} - {df['time'].max():.3f}",
                'most_common_marker': df['max_error_marker'].mode().iloc[0] if not df['max_error_marker'].mode().empty else 'N/A'
            }
            file_statistics.append(file_stats)
            
            # Count marker error occurrences
            marker_counts = df['max_error_marker'].value_counts()
            marker_error_counts.update(marker_counts.to_dict())
    
    if not all_frames:
        raise ValueError("No valid data found in any XML files")
    
    # Combine all frame data
    combined_df = pd.concat(all_frames, ignore_index=True)
    
    # Compute overall statistics
    overall_stats = {
        'total_files_processed': len(xml_files),
        'total_frames_analyzed': len(combined_df),
        'mean_rms_error': combined_df['rms_error'].mean(),
        'std_rms_error': combined_df['rms_error'].std(),
        'mean_max_error': combined_df['max_error'].mean(),
        'std_max_error': combined_df['max_error'].std(),
        'min_rms_error': combined_df['rms_error'].min(),
        'max_rms_error': combined_df['rms_error'].max(),
        'min_max_error': combined_df['max_error'].min(),
        'max_max_error': combined_df['max_error'].max()
    }
    
    # Get top 3 most common error markers
    top_3_markers = marker_error_counts.most_common(3)
    
    # Compute per-file means and their statistics
    file_means_rms = [stats['mean_rms_error'] for stats in file_statistics]
    file_means_max = [stats['mean_max_error'] for stats in file_statistics]
    
    file_level_stats = {
        'mean_of_file_rms_means': np.mean(file_means_rms),
        'std_of_file_rms_means': np.std(file_means_rms, ddof=1) if len(file_means_rms) > 1 else 0.0,
        'mean_of_file_max_means': np.mean(file_means_max),
        'std_of_file_max_means': np.std(file_means_max, ddof=1) if len(file_means_max) > 1 else 0.0
    }
    
    # Generate output filename if not provided
    if output_filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if subject_id:
            output_filename = f"ik_errors_batch_analysis_sub{subject_id}_{timestamp}.txt"
        else:
            output_filename = f"ik_errors_batch_analysis_all_subjects_{timestamp}.txt"
    
    output_path = os.path.join(directory, output_filename)
    
    # Write comprehensive report
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("IK ERRORS BATCH ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n")
        f.write(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Directory: {directory}\n")
        f.write(f"Subject Filter: {'sub' + subject_id if subject_id else 'All subjects'}\n")
        f.write(f"Files Processed: {overall_stats['total_files_processed']}\n")
        f.write(f"Total Frames Analyzed: {overall_stats['total_frames_analyzed']}\n")
        f.write("\n")
        
        # Overall Frame-Level Statistics
        f.write("OVERALL FRAME-LEVEL STATISTICS\n")
        f.write("=" * 50 + "\n")
        f.write(f"Mean RMS Error:         {overall_stats['mean_rms_error']:.6f} ± {overall_stats['std_rms_error']:.6f}\n")
        f.write(f"Mean Max Error:         {overall_stats['mean_max_error']:.6f} ± {overall_stats['std_max_error']:.6f}\n")
        f.write(f"RMS Error Range:        {overall_stats['min_rms_error']:.6f} - {overall_stats['max_rms_error']:.6f}\n")
        f.write(f"Max Error Range:        {overall_stats['min_max_error']:.6f} - {overall_stats['max_max_error']:.6f}\n")
        f.write("\n")
        
        # File-Level Statistics (Means of Means)
        f.write("FILE-LEVEL STATISTICS (MEANS OF FILE MEANS)\n")
        f.write("=" * 50 + "\n")
        f.write(f"Mean of RMS Error Means: {file_level_stats['mean_of_file_rms_means']:.6f} ± {file_level_stats['std_of_file_rms_means']:.6f}\n")
        f.write(f"Mean of Max Error Means: {file_level_stats['mean_of_file_max_means']:.6f} ± {file_level_stats['std_of_file_max_means']:.6f}\n")
        f.write("\n")
        
        # Top 3 Most Common Error Markers
        f.write("TOP 3 MOST COMMON ERROR MARKERS\n")
        f.write("=" * 50 + "\n")
        for i, (marker, count) in enumerate(top_3_markers, 1):
            percentage = (count / overall_stats['total_frames_analyzed']) * 100
            f.write(f"{i}. {marker}: {count} occurrences ({percentage:.1f}%)\n")
        f.write("\n")
        
        # Per-File Statistics Summary
        f.write("PER-FILE STATISTICS SUMMARY\n")
        f.write("=" * 80 + "\n")
        f.write(f"{'Filename':<50} {'Frames':<8} {'RMS Mean':<12} {'RMS Std':<12} {'Max Mean':<12} {'Max Std':<12} {'Top Marker':<15}\n")
        f.write("-" * 130 + "\n")
        
        for stats in file_statistics:
            f.write(f"{stats['filename']:<50} "
                   f"{stats['total_frames']:<8} "
                   f"{stats['mean_rms_error']:<12.6f} "
                   f"{stats['std_rms_error']:<12.6f} "
                   f"{stats['mean_max_error']:<12.6f} "
                   f"{stats['std_max_error']:<12.6f} "
                   f"{stats['most_common_marker']:<15}\n")
        
        f.write("\n")
        
        # Detailed File Information
        f.write("DETAILED FILE INFORMATION\n")
        f.write("=" * 50 + "\n")
        for i, stats in enumerate(file_statistics, 1):
            f.write(f"{i}. {stats['filename']}\n")
            f.write(f"   Frames: {stats['total_frames']}\n")
            f.write(f"   Time Range: {stats['time_range']} seconds\n")
            f.write(f"   RMS Error: {stats['mean_rms_error']:.6f} ± {stats['std_rms_error']:.6f}\n")
            f.write(f"   Max Error: {stats['mean_max_error']:.6f} ± {stats['std_max_error']:.6f}\n")
            f.write(f"   Most Common Error Marker: {stats['most_common_marker']}\n")
            f.write("\n")
        
        # All Marker Error Counts
        f.write("ALL MARKER ERROR COUNTS\n")
        f.write("=" * 50 + "\n")
        sorted_markers = sorted(marker_error_counts.items(), key=lambda x: x[1], reverse=True)
        for marker, count in sorted_markers:
            percentage = (count / overall_stats['total_frames_analyzed']) * 100
            f.write(f"{marker}: {count} occurrences ({percentage:.1f}%)\n")
        
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")
    
    print(f"\nBatch analysis completed!")
    print(f"Report saved to: {output_path}")
    
    # Prepare return dictionary
    results = {
        'overall_statistics': overall_stats,
        'file_level_statistics': file_level_stats,
        'top_3_markers': top_3_markers,
        'all_marker_counts': dict(marker_error_counts),
        'per_file_statistics': file_statistics,
        'combined_dataframe': combined_df,
        'output_file': output_path,
        'files_processed': [os.path.basename(f) for f in xml_files]
    }
    
    return results


def plot_ik_errors_boxplot(txt_files: List[str], 
                          output_directory: Optional[str] = None,
                          output_filename: Optional[str] = None,
                          figure_size: Tuple[int, int] = (12, 8),
                          show_plot: bool = True,
                          custom_labels: Optional[List[str]] = None) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create box and whisker plots from IK errors batch analysis report files.
    
    Parameters:
    -----------
    txt_files : list of str
        List of paths to .txt report files generated by parse_ik_errors_batch
    output_filename : str, optional
        Filename to save the plot. If None, auto-generates filename
    figure_size : tuple, optional
        Figure size as (width, height) in inches
    show_plot : bool, optional
        Whether to display the plot
    custom_labels : list of str, optional
        Custom labels for each dataset. If None, uses filenames
        
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects
    """
    
    def parse_report_file(txt_file: str) -> Dict:
        """Parse a single IK errors report file to extract statistics"""
        
        if not os.path.exists(txt_file):
            raise FileNotFoundError(f"Report file not found: {txt_file}")
        
        results = {
            'filename': os.path.basename(txt_file),
            'subject_filter': None,
            'files_processed': 0,
            'total_frames': 0,
            'frame_level_rms_mean': None,
            'frame_level_rms_std': None,
            'frame_level_max_mean': None,
            'frame_level_max_std': None,
            'file_level_rms_mean': None,
            'file_level_rms_std': None,
            'file_level_max_mean': None,
            'file_level_max_std': None,
            'per_file_data': []
        }
        
        with open(txt_file, 'r') as f:
            content = f.read()
        
        # Extract header information
        subject_match = re.search(r'Subject Filter:\s*(.+)', content)
        if subject_match:
            results['subject_filter'] = subject_match.group(1).strip()
        
        files_processed_match = re.search(r'Files Processed:\s*(\d+)', content)
        if files_processed_match:
            results['files_processed'] = int(files_processed_match.group(1))
        
        total_frames_match = re.search(r'Total Frames Analyzed:\s*(\d+)', content)
        if total_frames_match:
            results['total_frames'] = int(total_frames_match.group(1))
        
        # Extract overall frame-level statistics
        frame_rms_match = re.search(r'Mean RMS Error:\s*([\d.]+)\s*±\s*([\d.]+)', content)
        if frame_rms_match:
            results['frame_level_rms_mean'] = float(frame_rms_match.group(1))
            results['frame_level_rms_std'] = float(frame_rms_match.group(2))
        
        frame_max_match = re.search(r'Mean Max Error:\s*([\d.]+)\s*±\s*([\d.]+)', content)
        if frame_max_match:
            results['frame_level_max_mean'] = float(frame_max_match.group(1))
            results['frame_level_max_std'] = float(frame_max_match.group(2))
        
        # Extract file-level statistics (means of means)
        file_rms_match = re.search(r'Mean of RMS Error Means:\s*([\d.]+)\s*±\s*([\d.]+)', content)
        if file_rms_match:
            results['file_level_rms_mean'] = float(file_rms_match.group(1))
            results['file_level_rms_std'] = float(file_rms_match.group(2))
        
        file_max_match = re.search(r'Mean of Max Error Means:\s*([\d.]+)\s*±\s*([\d.]+)', content)
        if file_max_match:
            results['file_level_max_mean'] = float(file_max_match.group(1))
            results['file_level_max_std'] = float(file_max_match.group(2))
        
        # Extract per-file data from the summary table
        # Look for the table section
        table_start = content.find('PER-FILE STATISTICS SUMMARY')
        if table_start != -1:
            # Find the data lines (after the header line with dashes)
            table_section = content[table_start:]
            lines = table_section.split('\n')
            
            # Find where data starts (after the dash line)
            data_start_idx = None
            for i, line in enumerate(lines):
                if '---' in line:
                    data_start_idx = i + 1
                    break
            
            if data_start_idx is not None:
                # Parse each data line
                for line in lines[data_start_idx:]:
                    line = line.strip()
                    if not line or line.startswith('DETAILED'):
                        break
                    
                    # Split by whitespace and extract numeric values
                    parts = line.split()
                    if len(parts) >= 7:  # Ensure we have enough columns
                        try:
                            file_data = {
                                'filename': parts[0],
                                'frames': int(parts[1]),
                                'rms_mean': float(parts[2]),
                                'rms_std': float(parts[3]),
                                'max_mean': float(parts[4]),
                                'max_std': float(parts[5])
                            }
                            results['per_file_data'].append(file_data)
                        except (ValueError, IndexError):
                            continue  # Skip malformed lines
        
        return results
    
    def create_box_data(mean: float, std: float) -> Dict:
        """Create box plot data structure from mean and standard deviation"""
        # For box plot with std as whiskers:
        # - Box: mean ± 0.5*std (represents IQR-like range)
        # - Whiskers: mean ± std
        # - Median: mean
        
        return {
            'med': mean,           # Median line
            'q1': mean - 0.5*std,  # First quartile (bottom of box)
            'q3': mean + 0.5*std,  # Third quartile (top of box)
            'whislo': mean - std,  # Lower whisker
            'whishi': mean + std,  # Upper whisker
            'fliers': []           # No outliers for this synthetic data
        }
    
    # Validate inputs
    if not txt_files:
        raise ValueError("At least one txt file must be provided")
    
    for txt_file in txt_files:
        if not os.path.exists(txt_file):
            raise FileNotFoundError(f"File not found: {txt_file}")
    
    # Parse all report files
    all_results = []
    for txt_file in txt_files:
        try:
            results = parse_report_file(txt_file)
            all_results.append(results)
            print(f"Parsed: {os.path.basename(txt_file)}")
        except Exception as e:
            print(f"Error parsing {txt_file}: {e}")
            continue
    
    if not all_results:
        raise ValueError("No valid report files could be parsed")
    
    # Prepare data for plotting
    rms_box_data = []
    max_box_data = []
    labels = []
    
    for i, results in enumerate(all_results):
        # Use custom labels if provided, otherwise use subject filter or filename
        # if custom_labels and i < len(custom_labels):
        #     label = custom_labels[i]
        # elif results['subject_filter'] and results['subject_filter'] != 'All subjects':
        #     label = results['subject_filter']
        # else:
        #     # Extract meaningful part of filename
        #     filename = results['filename']
        #     # Try to extract subject ID or use filename
        #     subject_match = re.search(r'sub(\w+)', filename)
        #     if subject_match:
        #         label = f"sub{subject_match.group(1)}"
        #     else:
        #         label = os.path.splitext(filename)[0]
        filename = os.path.basename(results['filename'])
        label = filename[:9]
        
        labels.append(label)
        
        # Use file-level statistics (means of means) if available, otherwise frame-level
        if (results['file_level_rms_mean'] is not None and 
            results['file_level_rms_std'] is not None):
            rms_mean = results['file_level_rms_mean']
            rms_std = results['file_level_rms_std']
        else:
            rms_mean = results['frame_level_rms_mean']
            rms_std = results['frame_level_rms_std']
        
        if (results['file_level_max_mean'] is not None and 
            results['file_level_max_std'] is not None):
            max_mean = results['file_level_max_mean']
            max_std = results['file_level_max_std']
        else:
            max_mean = results['frame_level_max_mean']
            max_std = results['frame_level_max_std']
        
        # Create box plot data
        if rms_mean is not None and rms_std is not None:
            rms_box_data.append(create_box_data(rms_mean, rms_std))
        else:
            print(f"Warning: No RMS data found for {label}")
            rms_box_data.append(create_box_data(0, 0))
        
        if max_mean is not None and max_std is not None:
            max_box_data.append(create_box_data(max_mean, max_std))
        else:
            print(f"Warning: No Max Error data found for {label}")
            max_box_data.append(create_box_data(0, 0))
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figure_size)
    
    # Plot RMS Error box plots
    bp1 = ax1.bxp(rms_box_data, patch_artist=True)
    ax1.set_title('RMS Error Statistics\n(Box: mean ± 0.5σ, Whiskers: mean ± σ)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('RMS Error (meters)', fontsize=11)
    ax1.set_xticklabels(labels)
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3)
    # ax1.tick_params(axis='x', rotation=45)
    
    # Color the RMS boxes
    colors_rms = plt.cm.Blues(np.linspace(0.4, 0.8, len(rms_box_data)))
    for patch, color in zip(bp1['boxes'], colors_rms):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Plot Max Error box plots
    bp2 = ax2.bxp(max_box_data, patch_artist=True)
    ax2.set_title('Max Error Statistics\n(Box: mean ± 0.5σ, Whiskers: mean ± σ)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Max Error (meters)', fontsize=11)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='x', rotation=45)
    
    # Color the Max Error boxes
    colors_max = plt.cm.Reds(np.linspace(0.4, 0.8, len(max_box_data)))
    for patch, color in zip(bp2['boxes'], colors_max):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Add summary text
    fig.suptitle('IK Errors Comparison Across Datasets', fontsize=14, fontweight='bold')
    
    # Add info text
    info_text = f"Datasets: {len(all_results)} | "
    total_files = sum(r['files_processed'] for r in all_results if r['files_processed'])
    total_frames = sum(r['total_frames'] for r in all_results if r['total_frames'])
    info_text += f"Total Files: {total_files} | Total Frames: {total_frames}"
    
    fig.text(0.5, 0.02, info_text, ha='center', fontsize=10, style='italic')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, bottom=0.1)
    
    # Save plot if filename provided
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)
    else:
        output_directory = os.path.dirname(txt_files[0])

    if output_filename:
        # Determine output directory from first input file
        # output_directory = os.path.dirname(txt_files[0])
        output_path = os.path.join(output_directory, output_filename)
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    else:
        # Auto-generate filename
        output_directory = os.path.dirname(txt_files[0])
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        auto_filename = f"ik_errors_comparison_{timestamp}.png"
        output_path = os.path.join(output_directory, auto_filename)
        
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    
    # Display summary statistics
    print(f"\nSummary of {len(all_results)} datasets:")
    print("Dataset\t\t\tRMS Mean±Std\t\tMax Mean±Std")
    print("-" * 70)
    for i, (results, label) in enumerate(zip(all_results, labels)):
        rms_mean = results.get('file_level_rms_mean') or results.get('frame_level_rms_mean', 0)
        rms_std = results.get('file_level_rms_std') or results.get('frame_level_rms_std', 0)
        max_mean = results.get('file_level_max_mean') or results.get('frame_level_max_mean', 0)
        max_std = results.get('file_level_max_std') or results.get('frame_level_max_std', 0)
        
        print(f"{label:<20}\t{rms_mean:.6f}±{rms_std:.6f}\t{max_mean:.6f}±{max_std:.6f}")
    
    if show_plot:
        plt.show()
    
    return fig, (ax1, ax2)



if __name__ == "__main__":

    # Example usage
    data_directory = '/home/lisca/biomec/data/20250306_dataset/results_pool'
    output_dir = '/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/results'
    model_process = osim.ModelProcessor('/home/lisca/biomec/src/athlete/athlete/moco/example_muscle_adjusted_model.osim')
    ik_error_plot_directory = os.path.join(data_directory, 'ik_marker_rmse.png')

    # Analyze all subjects in directory
    results_all = parse_ik_errors_batch(
        directory=data_directory, 
        subject_id="04",
        output_filename="subject04_marker_errors.txt"
    )
    
    # Analyze specific subject only
    results_subject = parse_ik_errors_batch(
        directory=data_directory,
        output_filename="subject00_marker_errors.txt"
    )
    
    # Custom output filename
    results_custom = parse_ik_errors_batch(
        directory=data_directory, 
        subject_id="03",
        output_filename="subject03_marker_errors.txt"
    )
    
    # Access results programmatically
    print(f"Mean RMS Error: {results_all['overall_statistics']['mean_rms_error']:.6f}")
    print(f"Top 3 problematic markers: {[marker for marker, count in results_all['top_3_markers']]}")
    
    fig, axes = plot_ik_errors_boxplot([
        os.path.join(data_directory, "subject03_marker_errors.txt"),
        os.path.join(data_directory, "subject04_marker_errors.txt"),
        os.path.join(data_directory, "subject00_marker_errors.txt")], 
        output_filename=ik_error_plot_directory
        )

    
    model = model_process.process()
    model.initSystem()

    # muscle_mechanics = dict()
    # solution_path = 'example_muscle_driven_46.sto' # os.path.join(data_directory, 'sub03_strifc_cond00000_speed02000_0001_0001_mi_emg.sto') 
    # solution = osim.MocoTrajectory(solution_path)
    # neg_force_assessment = solution_path.replace('.sto', '_neg_forces.csv')

    # muscle_mechanics = calc_muscle_mechanics(model, solution)
    # test_path = os.path.join(output_dir, 'muscle_mechanics.sto')
    # get_negative_muscle_forces(model, solution, output_path=neg_force_assessment)
    
    # osim.STOFileAdapter.write(muscle_mechanics, test_path)
    
    # test_negative_forces_path = os.path.join(output_dir, 'negative_muscle_forces.txt')
    # osim.STOFileAdapter.write(neg_forces, test_negative_forces_path)

    # Analyze a specific coordinate

    for coord in ['knee_angle_r', 'hip_flexion_r', 'ankle_angle_r', 'mtp_angle_r']:
        print(f"Analyzing coordinate: {coord}")
        analyze_single_coordinate(data_directory, coordinate_name=coord, output_dir=output_dir, compute_errors=True)
    
    # analyze_single_coordinate(data_directory, coordinate_name='knee_angle_r', output_dir=output_dir, compute_errors=True)
    
    # Analyze all coordinates
    # analyze_all_coordinates(data_dir, output_dir)
    
    # Get available coordinates
    available_coords = get_available_coordinates(data_directory)
    print(f"Available coordinates: {available_coords}")