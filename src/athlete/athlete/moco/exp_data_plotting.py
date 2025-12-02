import os, sys, argparse
import retrieve_results
import validation as validate
import os
import warnings
from typing import List, Optional, Tuple
from  itertools import product

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
# 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute group averages for biomechanics data")
    parser.add_argument('-n', '--dry-run', action='store_true', help="dry run")
    parser.add_argument('--generate-plots', action='store_true', help="generate plots in addition to data processing")
    args = parser.parse_args()

    # Define directories
    input_dir = "/home/lisca/biomec/data/20250306_dataset/results_pool/"
    plot_dir = os.path.join(input_dir, 'plots/experimental_data_plots/')
    avg_inputs = os.path.join(input_dir, 'averages/')
    models_dir = os.path.join(input_dir, '../')
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(avg_inputs, exist_ok=True)

    # Define subjects and file types to process
    subjects = ["03", "04"]
    # file_types = [
    #     {"content": "emg.sto", "name": "EMG"},
    #     {"content": "ik_move.sto", "name": "Kinematics"}, 
    #     {"content": "grf_move.sto", "name": "Ground Reaction Forces"}
    # ]
    file_type_strings = ["emg.sto", "ik_move.sto", "grf_move.sto"]
    
    # Define analysis parameters - extract these from your existing files
    strike_types = ["fc", "ft", "rc", "rt"]  # Add more if you have different strike patterns: ["fc", "ff", "rf"]
    speeds = ["02000", "03000"]     # Add more if you have different speeds: ["02000", "01500", "02500"]
    
    # Optional: Define specific columns to plot for each file type
    columns_config = {
        "emg.sto": ['right_gastrocnemius_medialis', 'tibialis_anterior', 'right_soleus'],
        "ik_move.sto": ['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
        "grf_move.sto": ['ground_1_force_vx', 'ground_1_force_vy', 'ground_2_force_vx', 'ground_2_force_vy']
    }
    emg_map = {
        'r_sol': 'right_soleus',
        'r_gm': 'right_gastrocnemius_medialis', 
        'r_lm': 'right_gastrocnemius_lateralis'
    }
    
    print("=" * 80)
    print("STARTING BATCH PROCESSING OF GROUP AVERAGES")
    print("=" * 80)
    print(f"Input directory: {input_dir}")
    print(f"Plot directory: {plot_dir}")
    print(f"Output directory: {avg_inputs}")
    print(f"Subjects: {subjects}")
    print(f"Generate plots: {args.generate_plots}")
    print("=" * 80)
    
    # Track processing results
    successful_processing = []
    failed_processing = []
    
    for subject, speed, strike, file_type in product(subjects, speeds, strike_types, file_type_strings):

        print(f"Processing: {subject}, {speed}, {strike}, {file_type} this iteration")
        try:
            fig, ax = retrieve_results.compute_group_averages(
                                input_directory=input_dir,
                                output_directory=plot_dir,
                                sto_output_directory=avg_inputs,
                                file_content_type=file_type,
                                strike_type=strike, # list order goes fc, ft, rc, rt
                                speed=speed,
                                subject_id=subject,
                                generate_plot=True,
                                columns_to_plot=columns_config.get(file_type)
                            )
            plt.close(fig)
        except Exception as e:
            print(f" Error: {e}")

    
    #%% REPORT IK SOLUTION RELIABILITY FIRST

    # ik_error_plot_directory = os.path.join(input_dir, 'validation_plots/')

    # Analyze all subjects in directory
    results_all = validate.parse_ik_errors_batch(
        directory=input_dir, 
        subject_id="04",
        output_filename=os.path.join(plot_dir, "subject04_marker_errors.txt")
    )
    
    # Analyze specific subject only
    results_subject = validate.parse_ik_errors_batch(
        directory=input_dir,
        output_filename=os.path.join(plot_dir, "subject00_marker_errors.txt")
    )
    
    # Custom output filename
    results_custom = validate.parse_ik_errors_batch(
        directory=input_dir, 
        subject_id="03",
        output_filename=os.path.join(plot_dir, "subject03_marker_errors.txt")
    )
    
    # Access results programmatically
    print(f"Mean RMS Error: {results_all['overall_statistics']['mean_rms_error']:.6f}")
    print(f"Top 3 problematic markers: {[marker for marker, count in results_all['top_3_markers']]}")
    
    fig, axes = validate.plot_ik_errors_boxplot([
        os.path.join(plot_dir, "subject03_marker_errors.txt"),
        os.path.join(plot_dir, "subject04_marker_errors.txt"),
        os.path.join(plot_dir, "subject00_marker_errors.txt")], 
        output_filename = os.path.join(plot_dir, 'ik_marker_rmse.png')
        )

    
    
    
    #%% COMPARISON PLOTS - Two Conditions Side by Side

    # Define the conditions to compare
    comparison_configs = [
        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'file_type': 'ik_move.sto',
            'columns': ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'mtp_angle_r'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'ik_strike_comp_{column}.png' 
        },
        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'file_type': 'emg.sto',
            'columns': ['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_tibialis_anterior', 'right_soleus', 'right_peroneus_longus'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'emg_strike_comp_{column}.png' 
        },
        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'file_type': 'grf_move.sto',
            'columns': ['ground_2_force_vy', 'ground_1_force_vy'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'grf_strike_comp_{column}.png' 
        },

        {
            'name': 'FC_vs_FT',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFC 3.0m/s'},
            'condition2': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFT 3.0m/s'},
            'file_type': 'ik_move.sto',
            'columns': ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'mtp_angle_r'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'ffs_ik_lbs_comp_{column}.png' 
        },
        {
            'name': 'FC_vs_FT',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFC 3.0m/s'},
            'condition2': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFT 3.0m/s'},
            'file_type': 'emg.sto',
            'columns': ['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_tibialis_anterior', 'right_soleus', 'right_peroneus_longus'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'ffs_emg_lbs_comp_{column}.png' 
        },
        
        {
            'name': 'RC_vs_RT',
            'condition1': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFC 3.0m/s'},
            'condition2': {'strike_type': 'rt', 'speed': '03000', 'label': 'RFT 3.0m/s'},
            'file_type': 'ik_move.sto',
            'columns': ['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r', 'mtp_angle_r'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'rfs_ik_lbs_comp_{column}.png' 
        },
        {
            'name': 'RC_vs_RT',
            'condition1': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFC 3.0m/s'},
            'condition2': {'strike_type': 'rt', 'speed': '03000', 'label': 'RFT 3.0m/s'},
            'file_type': 'emg.sto',
            'columns': ['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_tibialis_anterior', 'right_soleus', 'right_peroneus_longus'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'rfs_emg_lbs_comp_{column}.png' 
        }
    ]

    # Loop through each comparison configuration
    for config in comparison_configs:
        print(f"\n{'='*80}")
        print(f"CREATING COMPARISON: {config['name']}")
        print(f"{'='*80}")
        
        file_type = config['file_type']
        columns = config['columns']

        output_plot_dir = config.get('plot_dir', plot_dir) #just dumps into generic plots upper folder if not classified
        
        filename_template = config.get('filename_template', '{comparison_name}_{column}.png')

        
        # Loop through each column - ONE PLOT PER COLUMN
        for column in columns:
            print(f"\n  Creating plot for: {column}")
            
            try:
                # Create base plot with CONDITION 1
                cond1 = config['condition1']
                print(f"    Adding condition 1: {cond1['label']}")
                
                fig, ax = retrieve_results.compute_group_averages(
                    input_directory=input_dir,
                    output_directory=config['plot_dir'],
                    sto_output_directory=avg_inputs,
                    file_content_type=file_type,
                    strike_type=cond1['strike_type'],
                    speed=cond1['speed'],
                    subject_id=None,  # All subjects
                    generate_plot=True,
                    columns_to_plot=[column],
                    plot_title=f"{column} - {cond1['label']} vs {config['condition2']['label']}",
                    label=cond1['label']
                )
                
                # Add CONDITION 2 to the same plot
                cond2 = config['condition2']
                print(f"    Adding condition 2: {cond2['label']}")
                
                fig, ax = retrieve_results.compute_group_averages(
                    input_directory=input_dir,
                    output_directory=config['plot_dir'],
                    sto_output_directory=avg_inputs,
                    file_content_type=file_type,
                    strike_type=cond2['strike_type'],
                    speed=cond2['speed'],
                    subject_id=None,  # All subjects
                    generate_plot=True,
                    columns_to_plot=[column],
                    existing_fig_ax=(fig, ax),  # Add to existing plot
                    plot_title=f"{column} - {cond1['label']} vs {cond2['label']}",
                    auto_color=True,
                    label=cond2['label']
                )
                
                # Save the comparison plot
                filename = filename_template.format(
                    comparison_name=config['name'],
                    column=column.replace('/', '_'),
                    cond1_strike=cond1['strike_type'],
                    cond2_strike=cond2['strike_type'],
                    cond1_speed=cond1['speed'],
                    cond2_speed=cond2['speed'],
                    cond1_label=cond1['label'].replace(' ', '_').replace('/', '_'),
                    cond2_label=cond2['label'].replace(' ', '_').replace('/', '_')
                )

                # plot_filename = f"{config['name']}_{column.replace('/', '_')}.png"
                plot_path = os.path.join(plot_dir, filename)
                fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                print(f"    💾 Saved: {filename}")
                
            except Exception as e:
                print(f"    ❌ Error creating plot for {column}: {e}")
                import traceback
                traceback.print_exc()
                continue
 
    # Copy over extrnal loads .xml files for the average GRFs
    source_directories = [os.path.join(models_dir, 'sub03/ses20250306/'),
                        os.path.join(models_dir, 'sub04/ses20250524/')]
    
    for subject_id in subjects:
        print(f"\nProcessing GRF files for subject: sub{subject_id}")
        
        for strike_type in strike_types:
            for speed in speeds:
                print(f"  Condition: strike={strike_type}, speed={speed}")
                
                # Find all matching GRF XML files for this condition
                pattern_prefix = f"sub{subject_id}_stri{strike_type}_cond00000_speed{speed}"
                matching_files = []
                
                try:
                    for dir in source_directories:
                        if not os.path.exists(dir):
                            print(f"    ⚠️  Source directory not found: {dir}")
                            continue

                        try:
                            all_files = os.listdir(dir)
                        except PermissionError:
                            print(f"    ⚠️  Permission denied for: {dir}")
                            continue

                        for filename in all_files:
                            full_path = os.path.join(dir, filename)
                            if (filename.startswith(pattern_prefix) and 
                                filename.endswith('_grf_move.xml') and
                                os.path.isfile(full_path)):
                                matching_files.append((filename, dir, full_path))

                    if len(matching_files) == 0:
                        print(f"    ⚠️  No GRF XML files found in {dir}")
                        failed_processing.append(f"GRF Copy: {pattern_prefix} - No files found")
                        continue

                    # Sort by filename only (not by directory)
                    matching_files.sort(key=lambda x: x[0])

                    
                    if len(matching_files) < 20:
                        print(f"    ⚠️  Only {len(matching_files)} files found (expected at least 20)")
                        # Use the last available file if less than 20
                        file_to_copy, source_dir_used, source_path = matching_files[-1]
                        print(f"    Using last available file: {file_to_copy}")
                    else:
                        # Use the 20th file (index 19)
                        file_to_copy, source_dir_used, source_path = matching_files[19]
                        print(f"    Found {len(matching_files)} files, using 20th: {file_to_copy}")
                    
                    # Define source and destination paths
                    source_path = os.path.join(source_dir_used, file_to_copy)
                    end_filename = f"sub{subject_id}_stri{strike_type}_cond00000_speed{speed}_grf_move_averages.xml"
                    dest_path = os.path.join(avg_inputs, end_filename)
                    
                    # Copy the file
                    if not args.dry_run:
                        import shutil
                        shutil.copy2(source_path, dest_path)
                        # shutil.copy2(source_path, dest_path)
                        print(f"    ✅ Copied to: {end_filename}")
                        successful_processing.append(f"GRF Copy: {end_filename}")
                    
                                        
                        # Update XML references inline
                        try:
                            # Read the XML file
                            with open(dest_path, 'r') as f:
                                xml_content = f.read()
                            
                            # Define the new reference filenames
                            new_datafile = f"sub{subject_id}_stri{strike_type}_cond00000_speed{speed}_grf_move_averages.sto"
                            new_kinematics = f"sub{subject_id}_stri{strike_type}_cond00000_speed{speed}_ik_move_averages.sto"
                            
                            # Update <datafile> reference using regex
                            import re
                            xml_content = re.sub(
                                r'<datafile>(.*?)</datafile>',
                                f'<datafile>{new_datafile}</datafile>',
                                xml_content
                            )
                            
                            # Update <external_loads_model_kinematics_file> reference
                            xml_content = re.sub(
                                r'<external_loads_model_kinematics_file>(.*?)</external_loads_model_kinematics_file>',
                                f'<external_loads_model_kinematics_file>{new_kinematics}</external_loads_model_kinematics_file>',
                                xml_content
                            )
                            
                            # Write the modified content back
                            with open(dest_path, 'w') as f:
                                f.write(xml_content)
                            
                            print(f"    ✅ Updated XML references:")
                            print(f"       datafile: {new_datafile}")
                            print(f"       kinematics: {new_kinematics}")
                            
                        except Exception as e:
                            print(f"    ⚠️  Warning: Could not update XML references: {e}")
                            # File was copied but references weren't updated
                        
                        successful_processing.append(f"GRF Copy: {end_filename}")
                    else:
                        print(f"    🔄 DRY RUN: Would copy {file_to_copy} to {end_filename}")
                        print(f"    🔄 DRY RUN: Would update XML references:")
                        print(f"       datafile: sub{subject_id}_stri{strike_type}_cond00000_speed{speed}_grf_move_averages.sto")
                        print(f"       kinematics: sub{subject_id}_stri{strike_type}_cond00000_speed{speed}_ik_move_averages.sto")
                        successful_processing.append(f"GRF Copy: {end_filename} (dry run)")

                except Exception as e:
                    print(f"    ❌ Error copying file: {e}")
                    failed_processing.append(f"GRF Copy: {pattern_prefix} - {e}")
                    import traceback
                    traceback.print_exc()


