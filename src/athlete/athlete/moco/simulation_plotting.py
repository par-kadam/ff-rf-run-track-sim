import os, sys, argparse
import retrieve_results
import matplotlib as plt
import os
import warnings
from typing import List, Optional, Tuple
import opensim as osim
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
    plot_dir = os.path.join(input_dir, 'plots/')
    avg_inputs = os.path.join(input_dir, 'averages/')
    models_dir = '/home/lisca/biomec/data/20250306_dataset/'
    model_file = '/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim'
    moco_model = '/home/lisca/biomec/src/athlete/athlete/moco/example_muscle_adjusted_model.osim'
    moco_solutions_dir = os.path.join(input_dir, 'averages/end_results/')
    
    # Define subjects and file types to process
    subjects = ["03", "04"]
    subject_colors = {
    "03": '#d62728',  # Red
    "04": '#9467bd',  # Purple
}
    file_types = [
        {"content": "emg.sto", "name": "EMG"},
        # {"content": "ik_move.sto", "name": "Kinematics"}, 
        # {"content": "grf_move.sto", "name": "Ground Reaction Forces"}
    ]
    
    # Define analysis parameters - extract these from your existing files
    strike_types = ["fc", "ft", "rc", "rt"]  # Add more if you have different strike patterns: ["fc", "ff", "rf"]
    speeds = ["02000", "03000"]     # Add more if you have different speeds: ["02000", "01500", "02500"]
    
    # Optional: Define specific columns to plot for each file type
    columns_config = {
        # "emg.sto": ['right_soleus', 'right_gastrocnemius_lateralis', 'right_gastrocnemius_medialis', 'right_tibialis_anterior'],
        "ik_move.sto": ['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
        # "grf_move.sto": ['ground_1_force_vy', 'ground_2_force_vy']  # Plot all available columns
    }

    # Define the conditions to compare
    solution_configs = [
        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'moco_tendon_strike_comp_3ms_con_{column}.png' 
        },
        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation', '/forceset/tibant_r/activation', '/forceset/perlong_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'moco_muscle_strike_comp_3ms_con_{column}.png' 
        },
        

        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rt', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'moco_tendon_strike_comp_3ms_treat_{column}.png' 
        },
        {
            'name': 'FFS_vs_RFS_3ms',
            'condition1': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'condition2': {'strike_type': 'rt', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation', '/forceset/tibant_r/activation', '/forceset/perlong_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'foot_strike_plots/'), 
            'filename_template': 'moco_muscle_strike_comp_3ms_treat_{column}.png' 
        },        

        {
            'name': 'FFC_vs_FFT_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFC 3.0m/s'},
            'condition2': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFT 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'moco_tendon_lbs_comp_3ms_{column}.png' 
        },
        {
            'name': 'FFC_vs_FFT_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFC 3.0m/s'},
            'condition2': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFT 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation', '/forceset/tibant_r/activation', '/forceset/perlong_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'moco_muscle_lbs_comp_3ms_{column}.png' 
        },

        {
            'name': 'RFC_vs_RFT_3ms',
            'condition1': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFC 3.0m/s'},
            'condition2': {'strike_type': 'rt', 'speed': '03000', 'label': 'RFT 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'moco_tendon_lbs_comp_3ms_{column}.png' 
        },
        {
            'name': 'RFC_vs_RFT_3ms',
            'condition1': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFC 3.0m/s'},
            'condition2': {'strike_type': 'rt', 'speed': '03000', 'label': 'RFT 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation', '/forceset/tibant_r/activation', '/forceset/perlong_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'lbs_comp_plots/'), 
            'filename_template': 'moco_muscle_lbs_comp_3ms_{column}.png' 
        },


        {
            'name': 'FFS_2_vs_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '02000', 'label': 'FFS 2.0m/s'},
            'condition2': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'running_speed_plots/'), 
            'filename_template': 'moco_tendon_speed_ffs_control_{column}.png' 
        },
        {
            'name': 'FFS_2_vs_3ms',
            'condition1': {'strike_type': 'fc', 'speed': '02000', 'label': 'FFS 2.0m/s'},
            'condition2': {'strike_type': 'fc', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'running_speed_plots/'), 
            'filename_template': 'moco_muscle_speed_ffs_control_{column}.png' 
        },
        {
            'name': 'FFS_2_vs_3ms',
            'condition1': {'strike_type': 'ft', 'speed': '02000', 'label': 'FFS 2.0m/s'},
            'condition2': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'running_speed_plots/'), 
            'filename_template': 'moco_tendon_speed_ffs_treat_{column}.png' 
        },
        {
            'name': 'FFS_2_vs_3ms',
            'condition1': {'strike_type': 'ft', 'speed': '02000', 'label': 'FFS 2.0m/s'},
            'condition2': {'strike_type': 'ft', 'speed': '03000', 'label': 'FFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'running_speed_plots/'), 
            'filename_template': 'moco_muscle_speed_ffs_treat_{column}.png' 
        },

        {
            'name': 'RFS_2_vs_3ms',
            'condition1': {'strike_type': 'rc', 'speed': '02000', 'label': 'RFS 2.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/normalized_tendon_force', '/forceset/gasmed_r/normalized_tendon_force', '/forceset/soleus_r/normalized_tendon_force'],
            'plot_dir': os.path.join(plot_dir, 'running_speed_plots/'), 
            'filename_template': 'moco_tendon_speed_rfs_control_{column}.png' 
        },
        {
            'name': 'RFS_2_vs_3ms',
            'condition1': {'strike_type': 'rc', 'speed': '02000', 'label': 'RFS 2.0m/s'},
            'condition2': {'strike_type': 'rc', 'speed': '03000', 'label': 'RFS 3.0m/s'},
            'columns': ['/forceset/gaslat_r/activation', '/forceset/gasmed_r/activation', '/forceset/soleus_r/activation'],
            'plot_dir': os.path.join(plot_dir, 'running_speed_plots/'), 
            'filename_template': 'moco_muscle_speed_rfs_control_{column}.png' 
        }
    ]

    # for config in solution_configs:
    #     columns = config['columns']

    #     output_plot_dir = config.get('plot_dir', plot_dir)

    #     filename_template = config.get('filename_template', '{comparison_name}_{column}.png')

    #     for column in columns:
    #         print(f"\n Creating plot for: {column}")

    #         try:
    #             cond1 = config['condition1']
    #             print(f"    Adding condition 1: {cond1['label']}")

    #             fig, ax = retrieve_results.average_moco_solution_columns(
    #                 input_directory=moco_solutions_dir,
    #                 output_directory=config['plot_dir'],
    #                 # file_suffix='_20res_02mesh_muscle_mechanics.csv',
    #                 file_suffix='_20res_02mesh.sto',
    #                 generate_plot=True,
    #                 columns_to_average=[column],
    #                 strike_type = cond1['strike_type'], 
    #                 speed = cond1['speed'],
    #                 plot_title=f"{column} - {cond1['label']} vs {config['condition2']['label']}",
    #                 verbose=False,
    #                 label=cond1['label']
    #             )

    #             cond2 = config['condition2']
    #             print(f"    Adding condition 2: {cond2['label']}")
                
    #             fig, ax = retrieve_results.average_moco_solution_columns(
    #                 input_directory=moco_solutions_dir,
    #                 output_directory=config['plot_dir'],
    #                 # file_suffix='__20res_02mesh_muscle_mechanics.csv',
    #                 file_suffix='_20res_02mesh.sto',
    #                 generate_plot=True,
    #                 columns_to_average=[column],
    #                 strike_type = cond2['strike_type'], 
    #                 speed = cond2['speed'],
    #                 plot_title=f"{column} - {cond1['label']} vs {cond2['label']}",
    #                 existing_fig_ax=(fig, ax),
    #                 verbose=False,
    #                 label=cond2['label']
    #             )

    #             # Save comparison
    #             filename = filename_template.format(
    #                 comparison_name=config['name'],
    #                 column=column.replace('/', '_'),
    #                 cond1_strike=cond1['strike_type'],
    #                 cond2_strike=cond2['strike_type'],
    #                 cond1_speed=cond1['speed'],
    #                 cond2_speed=cond2['speed'],
    #                 cond1_label=cond1['label'].replace(' ', '_').replace('/', '_'),
    #                 cond2_label=cond2['label'].replace(' ', '_').replace('/', '_')
    #             )

    #             # plot_filename = f"{config['name']}_{column.replace('/', '_')}.png"
    #             plot_path = os.path.join(config['plot_dir'], filename)
    #             fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    #             print(f"    💾 Saved: {plot_path}")
    #             plt.close(fig)

    #         except Exception as e:
    #             print(f"    ❌ Error creating plot for {column}: {e}")
    #             import traceback
    #             traceback.print_exc()
    #             continue

    

    #%% IK VALIDATION PLOTS   
    # ======================================================================================
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
    
    for subject_id in subjects:
        print(f"\n{'='*60}")
        print(f"PROCESSING SUBJECT: sub{subject_id}")
        print(f"{'='*60}")
        
        for strike_type in strike_types:
            for speed in speeds:
                print(f"\n--- Condition: strike={strike_type}, speed={speed} ---")
                
                for file_info in file_types:
                    file_type = file_info["content"]
                    file_name = file_info["name"]
                    
                    print(f"\nProcessing {file_name} ({file_type})...")
                    
                    try:
                        # Check if files exist for this combination
                        test_pattern = f"sub{subject_id}_stri{strike_type}_*_speed{speed}_*_{file_type}"
                        matching_files = []
                        
                        if os.path.exists(input_dir):
                            for filename in os.listdir(input_dir):
                                if (f"sub{subject_id}" in filename and 
                                    f"stri{strike_type}" in filename and 
                                    f"speed{speed}" in filename and 
                                    file_type in filename):
                                    matching_files.append(filename)
                        
                        if not matching_files:
                            print(f"  ⚠️  No files found matching pattern: {test_pattern}")
                            failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - No files found")
                            continue
                        
                        print(f"  📁 Found {len(matching_files)} matching files")
                        
                        columns_to_plot = columns_config.get(file_type, None)

                        # Process the group averages
                        if not args.dry_run:
                            if columns_to_plot:
                                print(f" Creating individual plts for {len(columns_to_plot)} columns")
                                for column in columns_to_plot:
                                    print(f"Plotting column: {column}")

                                    # fig, ax = retrieve_results.compute_group_averages(
                                    #     input_directory=input_dir,
                                    #     output_directory=plot_dir,
                                    #     sto_output_directory=avg_inputs,
                                    #     file_content_type='emg.sto',
                                    #     strike_type=strike_type,
                                    #     speed=speed,
                                    #     subject_id=f'sub{subject_id}',
                                    #     generate_plot=False,
                                    #     # columns_to_plot=columns_config['ik_move.sto']
                                    #     # plot_title=f"Subject {subject_id} - {file_name} (Strike: {strike_type}, Speed: {speed})"
                                    # )
                                    fig, ax = retrieve_results.compute_group_averages(
                                        input_directory=input_dir,
                                        output_directory=plot_dir,
                                        sto_output_directory=avg_inputs,
                                        file_content_type=file_type,
                                        strike_type=strike_type,
                                        speed=speed,
                                        subject_id=None,
                                        generate_plot=True,
                                        columns_to_plot=[column],
                                        plot_title=f"{column} - Strike: {strike_type}, Speed: {speed}m/s",
                                        # plot_filename=f"sub{subject_id}_{strike_type}_{speed}_{column.replace('/', '_')}.png"
                                    )

                                    for subject_id in subjects:
                                        print(f"adding sub{subject_id} moco data...")

                                        # Step 2: Find corresponding Moco solution file
                                        moco_pattern = f"sub{subject_id}_stri{strike_type}_*_speed{speed}_*_20res_02mesh.sto"
                                        moco_file = None
                                        
                                        if os.path.exists(moco_solutions_dir):
                                            for filename in os.listdir(moco_solutions_dir):
                                                if (f"sub{subject_id}" in filename and 
                                                    f"stri{strike_type}" in filename and 
                                                    f"speed{speed}" in filename and 
                                                    "_20res_02mesh.sto" in filename):
                                                    moco_file = os.path.join(moco_solutions_dir, filename)
                                                    break
                                        
                                        # Step 3: Add Moco solution data to the same plot
                                        if moco_file and os.path.exists(moco_file):
                                            print(f"      Adding Moco solution: {os.path.basename(moco_file)}")
                                            
                                            # Map column name to Moco coordinate path if needed
                                            # Adjust this mapping based on your coordinate naming
                                            # moco_column_map = {
                                            #     'hip_flexion_r': '/jointset/hip_r/hip_flexion_r/value',
                                            #     'knee_angle_r': '/jointset/walker_knee_r/knee_angle_r/value',
                                            #     'ankle_angle_r': '/jointset/ankle_r/ankle_angle_r/value',
                                            #     'mtp_angle_r': '/jointset/mtp_r/mtp_angle_r/value',
                                            # }
                                            activation_column_map = {
                                                'right_soleus': '/forceset/soleus_r/activation',
                                                'right_gastrocnemius_lateralis': '/forceset/gaslat_r/activation', 
                                                'right_gastrocnemius_medialis': '/forceset/gasmed_r/activation',
                                                'right_tibialis_anterior': '/forceset/tibant_r/activation'
                                            }

                                            activation_column_map = {
                                                'right_soleus': '/forceset/soleus_r/activation',
                                                'right_gastrocnemius_lateralis': '/forceset/gaslat_r/activation', 
                                                'right_gastrocnemius_medialis': '/forceset/gasmed_r/activation',
                                                'right_tibialis_anterior': '/forceset/tibant_r/activation'
                                            }
                                            
                                            # Get the Moco column name (use mapping or original)
                                            moco_column = activation_column_map.get(column, column)
                                            print(f'adding column {moco_column} from solution')
                                            subject_color = subject_colors.get(subject_id, '#000000') # defaults to black
                                            
                                            try:
                                                fig, ax = retrieve_results.plot_moco_solution_columns(
                                                    moco_solution_file=moco_file,
                                                    columns_to_plot=[moco_column],
                                                    model_file=model_file,
                                                    subject_id=f'sub{subject_id}',
                                                    existing_fig_ax=(fig, ax),  # Add to existing plot
                                                    line_style='-',  # solid line for Moco
                                                    line_width=2,
                                                    color=subject_color,
                                                    normalize_time=True,
                                                    time_as_percentage=True,
                                                    show_legend=False,
                                                    legend_location='upper left',
                                                    xlabel = 'Gait Cycle (%)',
                                                    ylabel = f'{column} angle (deg)'
                                                )
                                                print(f"      ✅ Added Moco data to plot")
                                            except Exception as e:
                                                print(f"      ⚠️  Could not add Moco data: {e}")
                                        else:
                                            print(f"      ⚠️  No Moco solution found matching: {moco_pattern}")
                                        
                                    # Step 4: Save the combined plot
                                    # if args.generate_plots:
                                    plot_filename = f"sub{subject_id}_{strike_type}_{speed}_{column.replace('/', '_')}_combined.png"
                                    plot_path = os.path.join(plot_dir, 'validation_plots/', plot_filename)
                                    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
                                    print(f"      💾 Saved combined plot: {plot_filename}")
                                
                                    print(f"    ✅ Created plot for {column} in {plot_path}")
    

                                        # fig_kin, ax_kin = retrieve_results.compute_group_averages(
            #                             input_directory=input_dir,
            #                             output_directory=plot_dir,
            #                             sto_output_directory=plot_dir,
            #                             file_content_type='grf_move.sto',
            #                             strike_type=strike_types[2], # list order goes fc, ft, rc, rt
            #                             speed=speeds[0],
            #                             # subject_id='sub00',
            #                             generate_plot=True,
            #                             columns_to_plot=columns_config['grf_move.sto'],
            #                             plot_title=f"FFS vs. RFS 2m/s Ground Reactions",
            #                             plot_filename='ffs_rfs_0200_con_grf.png',
            #                             existing_fig_ax = (fig, ax)
            #                         )
                            
                            print(f"Created plot for {column}")
                            print(f"  ✅ Successfully processed {file_name} for sub{subject_id}")
                            successful_processing.append(f"sub{subject_id}_{strike_type}_{speed}_ik_move.sto")
                        else:
                            print(f"  🔄 DRY RUN: Would process {file_name} for sub{subject_id}")
                            successful_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} (dry run)")
                        
                    except ValueError as e:
                        if "No matching files found" in str(e):
                            print(f"  ⚠️  {e}")
                            failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - {e}")
                        else:
                            print(f"  ❌ ValueError processing {file_name} for sub{subject_id}: {e}")
                            failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - {e}")
                    except Exception as e:
                        print(f"  ❌ Error processing {file_name} for sub{subject_id}: {e}")
                        failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - {e}")
                        continue

    #%% Plot MUSCLE TENDON COMPARISONS  
    # fig, ax = retrieve_results.plot_muscle_mechanics_bars(
    # csv_directory=moco_solutions_dir,
    # suffix='_20res_02mesh_muscle_mechanics.csv',
    # column_to_plot= ['/forceset/soleus_r|tendon_force', '/forceset/gaslat_r|tendon_force', '/forceset/gasmed_r|tendon_force'],
    # subject_ids=['03', '04'],
    # speed_strike_pair=[['fc', '02000'], 
    #                    ['ft', '02000'],
    #                    ['fc', '03000'],
    #                    ['ft', '03000']],
    # plot_color='blue',
    # ylabel='Peak Tendon Force (N)',
    # title='Peak Tendon Force',
    # output_file=os.path.join(plot_dir, 'lbs_comp_plots/tendon_force_comp_ffs.png')
    # )

    # fig, ax = retrieve_results.plot_muscle_mechanics_bars(
    # csv_directory=moco_solutions_dir,
    # suffix='_20res_02mesh_muscle_mechanics.csv',
    # column_to_plot= ['/forceset/soleus_r|tendon_force', '/forceset/gaslat_r|tendon_force', '/forceset/gasmed_r|tendon_force'],
    # subject_ids=['03', '04'],
    # speed_strike_pair=[['rc', '02000'], 
    #                    ['rt', '02000'],
    #                    ['rc', '03000'],
    #                    ['rt', '03000']],
    # plot_color='blue',
    # ylabel='Peak Tendon Force (N)',
    # title='Peak Tendon Force',
    # output_file=os.path.join(plot_dir, 'lbs_comp_plots/tendon_force_comp_rfs.png')
    # )



    # fig, ax = retrieve_results.plot_muscle_mechanics_bars(
    # csv_directory=moco_solutions_dir,
    # suffix='_20res_02mesh_muscle_mechanics.csv',
    # column_to_plot= ['/forceset/soleus_r|tendon_strain', '/forceset/gaslat_r|tendon_strain', '/forceset/gasmed_r|tendon_strain'],
    # subject_ids=['03', '04'],
    # speed_strike_pair=[['fc', '02000'], 
    #                    ['ft', '02000'],
    #                    ['fc', '03000'],
    #                    ['ft', '03000']],
    # plot_color='blue',
    # ylabel='Peak Tendon Strain (%)',
    # title=None,
    # output_file=os.path.join(plot_dir, 'lbs_comp_plots/tendon_strain_comp_ffs.png')
    # )

    # fig, ax = retrieve_results.plot_muscle_mechanics_bars(
    # csv_directory=moco_solutions_dir,
    # suffix='_20res_02mesh_muscle_mechanics.csv',
    # column_to_plot= ['/forceset/soleus_r|tendon_strain', '/forceset/gaslat_r|tendon_strain', '/forceset/gasmed_r|tendon_strain'],
    # subject_ids=['03', '04'],
    # speed_strike_pair=[['rc', '02000'], 
    #                    ['rt', '02000'],
    #                    ['rc', '03000'],
    #                    ['rt', '03000']],
    # plot_color='blue',
    # ylabel='Peak Tendon Strain (%)',
    # title=None,
    # output_file=os.path.join(plot_dir, 'lbs_comp_plots/tendon_strain_comp_rfs.png')
    # )

    # fig, ax = retrieve_results.plot_muscle_mechanics_bars(
    # csv_directory=moco_solutions_dir,
    # suffix='_20res_02mesh_muscle_mechanics.csv',
    # column_to_plot= ['/forceset/soleus_r|tendon_force', '/forceset/gaslat_r|tendon_force', '/forceset/gasmed_r|tendon_force'],
    # subject_ids=['03', '04'],
    # speed_strike_pair=[['rc', '02000'], 
    #                    ['fc', '02000'],
    #                    ['rc', '03000'],
    #                    ['fc', '03000']],
    # plot_color='blue',
    # ylabel='Peak Tendon Force (N)',
    # title='Peak Tendon Force',
    # output_file=os.path.join(plot_dir, 'foot_strike_plots/tendon_force_comp_con.png')
    # )

    # #%% JOINT MOMENT BREAKDOWNS
    # # # The directory containing all your .sto solution files
    # # # Use 'r' before the string to handle Windows paths correctly
    # input_directory = moco_solutions_dir

    # # The common suffix of the files you want to process
    # # This ensures we only run on "solution" files and not others
    # solution_suffix = "_20res_02mesh.sto" 

    # # The list of coordinate paths for your function
    # coord_paths = ['/jointset/hip_r/hip_flexion_r',
    #             '/jointset/walker_knee_r/knee_angle_r',
    #             '/jointset/ankle_r/ankle_angle_r',
    #             '/jointset/mtp_r/mtp_angle_r'
    # ]

    # # The directory where all generated plots will be saved
    # output_dir = os.path.join(plot_dir, 'validation_plots/')

    # # Directory containing model files
    # models_base_dir = '/home/lisca/biomec/data/20250306_dataset/'

    # def find_model_for_solution(solution_filename, models_dir):
    #     """
    #     Find the model file that matches the first 5 characters of the solution filename.
        
    #     Parameters:
    #     -----------
    #     solution_filename : str
    #         The solution filename (e.g., 'sub03_strifc_cond00000_speed02000_20res_02mesh.sto')
    #     models_dir : str
    #         Base directory containing subject folders with model files
            
    #     Returns:
    #     --------
    #     str or None
    #         Full path to the matching model file, or None if not found
    #     """
    #     # Extract first 5 characters (e.g., 'sub03')
    #     subject_id = solution_filename[:5]
        
    #     # Look for subject directory
    #     subject_dir = None
    #     for item in os.listdir(models_dir):
    #         if item.startswith(subject_id) and os.path.isdir(os.path.join(models_dir, item)):
    #             subject_dir = os.path.join(models_dir, item)
    #             break
        
    #     if not subject_dir:
    #         print(f"  Warning: Could not find directory for {subject_id}")
    #         return None
        
    #     # Look for session directory
    #     session_dir = None
    #     for ses_item in os.listdir(subject_dir):
    #         if 'ses' in ses_item and os.path.isdir(os.path.join(subject_dir, ses_item)):
    #             session_dir = os.path.join(subject_dir, ses_item)
    #             break
        
    #     if not session_dir:
    #         print(f"  Warning: Could not find session directory in {subject_dir}")
    #         return None
        
    #     # Look for model file matching the subject ID
    #     for model_file in os.listdir(session_dir):
    #         if (model_file.startswith(subject_id) and 
    #             model_file.endswith('_muscle_adjusted_model.osim') and
    #             os.path.isfile(os.path.join(session_dir, model_file))):
    #             model_path = os.path.join(session_dir, model_file)
    #             print(f"  Found model: {os.path.basename(model_path)}")
    #             return model_path
        
    #     print(f"  Warning: Could not find model file for {subject_id} in {session_dir}")
    #     return None




    # # ==========================================================================
    # # SCRIPT EXECUTION ---

    # """
    # Finds all solution files and runs the plotting function on them.
    # """
    
    # # Ensure the output directory exists
    # os.makedirs(output_dir, exist_ok=True)
    # print(f"Output directory set to: {output_dir}\n")

    # print("Starting batch plotting...")

    # # Get all files in the input directory
    # # try:
    # all_files = os.listdir(input_directory)
    # # except FileNotFoundError:
    # #     print(f"Error: Input directory not found at '{input_directory}'")
    # #     print("Please check the 'input_directory' variable.")
    # #     return # Stop the script
    # # except Exception as e:
    # #     print(f"An unexpected error occurred while reading the directory: {e}")
    # #     return

    # # Filter for the solution files
    # solution_files = [f for f in all_files 
    #                   if f.endswith(solution_suffix) and 
    #                   os.path.isfile(os.path.join(input_directory, f))]

    # if not solution_files:
    #     print(f"Warning: No files found in '{input_directory}' with suffix '{solution_suffix}'")
    # else:
    #     print(f"Found {len(solution_files)} solution file(s) to process.\n")

    # loaded_models = {}

    # # Loop through each found solution file
    # for filename in solution_files:
    #     # Construct the full, absolute path to the file
    #     full_solution_path = os.path.join(input_directory, filename)
        
    #     # Find the appropriate model for this solution
    #     model_path = find_model_for_solution(filename, models_base_dir)
    
    #     if model_path is None:
    #         print(f"  Skipping {filename} - no matching model found")
    #         continue
        
    #     # Check if we've already loaded this model (to save time)
    #     subject_id = filename[:5]
    #     if subject_id not in loaded_models:
    #         print(f"  Loading model for {subject_id}...")
    #         try:
    #             model = osim.Model(model_path)
    #             modelProcessor = osim.ModelProcessor(model)
    #             model = modelProcessor.process()
    #             model.initSystem()
    #             loaded_models[subject_id] = model
    #             print(f"  Model loaded and initialized successfully")
    #         except Exception as e:
    #             print(f"  Error loading model: {e}")
    #             print(f"  Skipping {filename}")
    #             continue
    #     else:
    #         print(f"  Using cached model for {subject_id}")
    #         model = loaded_models[subject_id]
        
    #     # Load the solution trajectory
    #     try:
    #         solution_trajectory = osim.MocoTrajectory(full_solution_path)
    #     except Exception as e:
    #         print(f"  Error loading solution trajectory: {e}")
    #         print(f"  Skipping {filename}")
    #         continue
        
    #     # Call your plotting function
    #     try:
    #         retrieve_results.plot_joint_moment_constituents(
    #             solution=solution_trajectory,
    #             file_name=filename,
    #             coord_paths=coord_paths,
    #             output_file=output_dir,
    #             model=model
    #         )
    #         retrieve_results.plot_joint_fiber_forces(
    #             model, 
    #             solution_trajectory, 
    #             coord_paths=coord_paths,
    #             file_name=filename, output_dir=output_dir, 
    #             min_moment_arm_threshold=0.00001,
    #             normalize_time=True
    #         )
    #         print(f"Successfully processed {filename}")
    #     except Exception as e:
    #         print(f"Error processing {filename}: {e}")
    #         continue

    # print("\nBatch plotting complete.")