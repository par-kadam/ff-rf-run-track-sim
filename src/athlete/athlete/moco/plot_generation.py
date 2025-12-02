import os, sys, argparse
import retrieve_results

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
    
    # Define subjects and file types to process
    subjects = ["03", "04"]
    file_types = [
        {"content": "emg.sto", "name": "EMG"},
        {"content": "ik_move.sto", "name": "Kinematics"}, 
        {"content": "grf_move.sto", "name": "Ground Reaction Forces"}
    ]
    
    # Define analysis parameters - extract these from your existing files
    strike_types = ["fc", "ft", "rc", "rt"]  # Add more if you have different strike patterns: ["fc", "ff", "rf"]
    speeds = ["02000", "03000"]     # Add more if you have different speeds: ["02000", "01500", "02500"]
    
    # Optional: Define specific columns to plot for each file type
    columns_config = {
        "emg.sto": ['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_soleus'],
        "ik_move.sto": ['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
        "grf_move.sto": None  # Plot all available columns
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
    
    # %%
    # ANALYSIS SECTION: plot illustrating differences between strike types for confirmation at same speeds (x4 EACH FILE TYPE = 12)
    # plot illustrating differences between control and treatment conditions in each condition (x4 EACH FILE TYPE = 12)
    #^^*Will have to shrink/zoom into stance phase only, then apply velocity calculations per coordinate

    # RFS vs FFS
    # fig, ax = retrieve_results.compute_group_averages(
    #                             input_directory=input_dir,
    #                             output_directory=plot_dir,
    #                             sto_output_directory=plot_dir,
    #                             file_content_type='ik_move.sto',
    #                             strike_type=strike_types[0], # list order goes fc, ft, rc, rt
    #                             speed=speeds[0],
    #                             # subject_id='sub00',
    #                             generate_plot=True,
    #                             columns_to_plot=columns_config['ik_move.sto'],
    #                             plot_title=f"FFS vs. RFS 2m/s Kinematics"
    #                         )
    # fig_kin, ax_kin = retrieve_results.compute_group_averages(
    #                             input_directory=input_dir,
    #                             output_directory=plot_dir,
    #                             sto_output_directory=plot_dir,
    #                             file_content_type='ik_move.sto',
    #                             strike_type=strike_types[2], # list order goes fc, ft, rc, rt
    #                             speed=speeds[0],
    #                             # subject_id='sub00',
    #                             generate_plot=True,
    #                             columns_to_plot=columns_config['ik_move.sto'],
    #                             plot_title=f"FFS vs. RFS 2m/s Kinematics",
    #                             plot_filename='ffs_rfs_0200_con_ik.png',
    #                             existing_fig_ax = (fig, ax)
    #                         )
    
    # fig_emg, ax_emg = retrieve_results.compute_group_averages(
    #                             input_directory=input_dir,
    #                             output_directory=plot_dir,
    #                             sto_output_directory=plot_dir,
    #                             file_content_type='emg.sto',
    #                             strike_type=strike_types[0], # list order goes fc, ft, rc, rt
    #                             speed=speeds[0],
    #                             # subject_id='sub00',
    #                             generate_plot=True,
    #                             columns_to_plot=columns_config['emg.sto'],
    #                             plot_title=f"FFS vs. RFS 2m/s EMG"
    #                         )
    # fig_emg_comp, ax_emg_comp = retrieve_results.compute_group_averages(
    #                             input_directory=input_dir,
    #                             output_directory=plot_dir,
    #                             sto_output_directory=plot_dir,
    #                             file_content_type='emg.sto',
    #                             strike_type=strike_types[2], # list order goes fc, ft, rc, rt
    #                             speed=speeds[0],
    #                             # subject_id='sub00',
    #                             generate_plot=True,
    #                             columns_to_plot=columns_config['emg.sto'],
    #                             plot_title=f"FFS vs. RFS 2m/s EMG",
    #                             plot_filename='ffs_rfs_0200_con_emg.png',
    #                             existing_fig_ax = (fig_emg, ax_emg)
    #                         )
    
    # fig, ax = retrieve_results.compute_group_averages(
    #                             input_directory=input_dir,
    #                             output_directory=plot_dir,
    #                             sto_output_directory=plot_dir,
    #                             file_content_type='grf_move.sto',
    #                             strike_type=strike_types[0], # list order goes fc, ft, rc, rt
    #                             speed=speeds[0],
    #                             # subject_id='sub00',
    #                             generate_plot=True,
    #                             columns_to_plot=columns_config['grf_move.sto'],
    #                             plot_title=f"FFS vs. RFS 2m/s Ground Reactions"
    #                         )
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

    #%% SECTION FOR LOOPING THROUGH SUBJECTS, CREATING CONDITION-SPECIFIC PLOTS 
    # ====================================================================
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
                        
                        # Process the group averages
                        if not args.dry_run:
                            fig, ax = retrieve_results.compute_group_averages(
                                input_directory=input_dir,
                                output_directory=plot_dir,
                                sto_output_directory=avg_inputs,
                                file_content_type='emg.sto',
                                strike_type=strike_type,
                                speed=speed,
                                subject_id=f'sub{subject_id}',
                                generate_plot=False,
                                # columns_to_plot=columns_config['ik_move.sto']
                                # plot_title=f"Subject {subject_id} - {file_name} (Strike: {strike_type}, Speed: {speed})"
                            )
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
    
    