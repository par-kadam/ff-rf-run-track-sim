from inverse_kinematics_emg import solveMocoInverseKinematicsEMG
from troubleshoot_post_addspring import muscleDrivenStateTracking
import retrieve_results
import os
import sys
import argparse


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
    
    # for subject_id in subjects:
    #     print(f"\n{'='*60}")
    #     print(f"PROCESSING SUBJECT: sub{subject_id}")
    #     print(f"{'='*60}")
        
    #     for strike_type in strike_types:
    #         for speed in speeds:
    #             print(f"\n--- Condition: strike={strike_type}, speed={speed} ---")
                
    #             for file_info in file_types:
    #                 file_type = file_info["content"]
    #                 file_name = file_info["name"]
                    
    #                 print(f"\nProcessing {file_name} ({file_type})...")
                    
    #                 try:
    #                     # Check if files exist for this combination
    #                     test_pattern = f"sub{subject_id}_stri{strike_type}_*_speed{speed}_*_{file_type}"
    #                     matching_files = []
                        
    #                     if os.path.exists(input_dir):
    #                         for filename in os.listdir(input_dir):
    #                             if (f"sub{subject_id}" in filename and 
    #                                 f"stri{strike_type}" in filename and 
    #                                 f"speed{speed}" in filename and 
    #                                 file_type in filename):
    #                                 matching_files.append(filename)
                        
    #                     if not matching_files:
    #                         print(f"  ⚠️  No files found matching pattern: {test_pattern}")
    #                         failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - No files found")
    #                         continue
                        
    #                     print(f"  📁 Found {len(matching_files)} matching files")
                        
    #                     # Process the group averages
    #                     if not args.dry_run:
    #                         fig, ax = retrieve_results.compute_group_averages(
    #                             input_directory=input_dir,
    #                             output_directory=plot_dir,
    #                             sto_output_directory=avg_inputs,
    #                             file_content_type=file_type,
    #                             strike_type=strike_type,
    #                             speed=speed,
    #                             subject_id=subject_id,
    #                             generate_plot=False,
    #                             # columns_to_plot=columns_config.get(file_type, None),
    #                             # plot_title=f"Subject {subject_id} - {file_name} (Strike: {strike_type}, Speed: {speed})"
    #                         )
                            
    #                         print(f"  ✅ Successfully processed {file_name} for sub{subject_id}")
    #                         successful_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type}")
    #                     else:
    #                         print(f"  🔄 DRY RUN: Would process {file_name} for sub{subject_id}")
    #                         successful_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} (dry run)")
                        
    #                 except ValueError as e:
    #                     if "No matching files found" in str(e):
    #                         print(f"  ⚠️  {e}")
    #                         failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - {e}")
    #                     else:
    #                         print(f"  ❌ ValueError processing {file_name} for sub{subject_id}: {e}")
    #                         failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - {e}")
    #                 except Exception as e:
    #                     print(f"  ❌ Error processing {file_name} for sub{subject_id}: {e}")
    #                     failed_processing.append(f"sub{subject_id}_{strike_type}_{speed}_{file_type} - {e}")
    #                     continue
    
    # STEP 2: Run first iteration of tracking problem using Moco Inverse Kinematics as first guess
    print(f"\n🚀 STEP 2: RUNNING MOCO INVERSE KINEMATICS...")
    
    # Create a list of speeds first to iterate through
    periodic_speed_list = [float(speed) / 1000 for speed in speeds]

    # Track solution guesses (persists across subjects for warm starts)
    solution_guess_history = {}
    
    for subject_id in subjects:
        print(f"\n{'='*60}")
        print(f"RUNNING MOCO FOR SUBJECT: sub{subject_id}")
        print(f"{'='*60}")
        
        # Find subject directory and session
        subject_dir = None
        session_dir = None
        
        for item in os.listdir(models_dir):
            if f"sub{subject_id}" in item and os.path.isdir(os.path.join(models_dir, item)):
                subject_dir = os.path.join(models_dir, item)
                # Find session directory
                for ses_item in os.listdir(subject_dir):
                    if 'ses' in ses_item and os.path.isdir(os.path.join(subject_dir, ses_item)):
                        session_dir = os.path.join(subject_dir, ses_item)
                        break
                break
        
        if not subject_dir or not session_dir:
            print(f"  ❌ Could not find subject directory for sub{subject_id}")
            failed_processing.append(f"Moco: sub{subject_id} - Directory not found")
            continue
        
        # Find model and muscle paths
        path_model = None
        path_musclepaths = None
        
        for f in os.listdir(session_dir):
            if f.endswith('_SCALED_addbio_free_subt_mtp.osim'):
                path_model = os.path.join(session_dir, f)
            if f.endswith('_FunctionBasedPathSet.xml'):
                path_musclepaths = os.path.join(session_dir, f)
        
        if not path_model:
            print(f"  ❌ Model file not found for sub{subject_id}")
            failed_processing.append(f"Moco: sub{subject_id} - Model file not found")
            continue
        
        print(f"  📁 Model: {os.path.basename(path_model)}")
        if path_musclepaths:
            print(f"  📁 Muscle paths: {os.path.basename(path_musclepaths)}")
        else:
            print(f"  ⚠️  Muscle paths file not found (proceeding without)")

        
        
        for i, strike_type in enumerate(strike_types):
            for j, speed_str in enumerate(speeds):
                current_speed = periodic_speed_list[j]
                    
                print(f"\n--- Processing: strike={strike_type}, speed={speed_str} ---")
                
                # Define file paths for averaged data
                base_filename = f"sub{subject_id}_stri{strike_type}_cond00000_speed{speed_str}"
                
                path_kinematics = os.path.join(avg_inputs, f"{base_filename}_ik_move_averages.sto")
                path_emg = os.path.join(avg_inputs, f"{base_filename}_emg_averages.sto")
                path_grf = os.path.join(avg_inputs, f"{base_filename}_grf_move_averages.xml")  # Note: might need .sto extension
                path_states = os.path.join(avg_inputs, f"{base_filename}_rad_averages.sto")
                path_guess = os.path.join(avg_inputs, 'guesses/', f"{base_filename}_mi_emg_solution.sto")
                # path_solution = os.path.join(avg_inputs, 'guesses/', f"{base_filename}_mi_emg_solution.sto")
                
                path_solution = os.path.join(avg_inputs, 'guesses/', f"{base_filename}_50res_1mesh.sto")
                path_log = os.path.join(avg_inputs, 'guesses/', f"{base_filename}_50res_1mesh.log")
                
                # Check if required input files exist
                required_files = {
                    'kinematics': path_kinematics,
                    'emg': path_emg,
                    'grf': path_grf
                }
                
                # missing_files = []
                # for file_type, file_path in required_files.items():
                #     if not os.path.exists(file_path):
                #         # Try alternative extensions for GRF
                #         if file_type == 'grf' and not os.path.exists(file_path):
                #             alt_grf = file_path.replace('.xml', '.sto')
                #             if os.path.exists(alt_grf):
                #                 path_grf = alt_grf
                #             else:
                #                 missing_files.append(file_type)
                #         else:
                #             missing_files.append(file_type)
                
                # if missing_files:
                #     print(f"  ⚠️  Missing required files: {missing_files}")
                #     failed_processing.append(f"Moco: {base_filename} - Missing files: {missing_files}")
                #     continue
                
                # Check if solution already exists
                if os.path.exists(path_solution) and os.path.getsize(path_solution) > 0:
                    print(f"  ✅ Solution already exists: {os.path.basename(path_solution)}")
                    # Save as potential guess for next iteration
                    solution_guess_history[subject_id] = path_solution
                    successful_processing.append(f"Moco: {base_filename} - Already exists")
                    continue
                
                # Generate states file if it doesn't exist
                if not os.path.exists(path_states):
                    print(f"  🔄 Generating states file: {os.path.basename(path_states)}")
                    try:
                        if not args.dry_run:
                            # Import the function (assuming it's available)
                            from osimFunctions import kinematicsToStates
                            kinematicsToStates(
                                kinematicsFileName=path_kinematics,
                                osimModelFileName=path_model,
                                outputFileName=path_states,
                                inDegrees=True,
                                outDegrees=False
                            )
                            print(f"    ✅ States file generated")
                        else:
                            print(f"    🔄 DRY RUN: Would generate states file")
                    except Exception as e:
                        print(f"    ❌ Error generating states file: {e}")
                        failed_processing.append(f"Moco: {base_filename} - States generation failed: {e}")
                        continue
                
                # Get solution guess (previous solution for this subject if available)
                path_solution_guess = solution_guess_history.get(subject_id, None)
                if path_solution_guess:
                    print(f"  🎯 Using solution guess: {os.path.basename(path_solution_guess)}")
                else:
                    print(f"  🎯 No solution guess available (cold start)")
                
                # Run Moco Inverse Kinematics
                print(f"  🚀 Running Moco Inverse Kinematics...")
                print(f"    Model: {os.path.basename(path_model)}")
                print(f"    Kinematics: {os.path.basename(path_kinematics)}")
                print(f"    EMG: {os.path.basename(path_emg)}")
                print(f"    GRF: {os.path.basename(path_grf)}")
                print(f"    States: {os.path.basename(path_states)}")
                print(f"    Solution: {os.path.basename(path_solution)}")
                
                

                try:
                    if not args.dry_run:
                        # Create placeholder file to mark as in-progress
                        with open(path_solution, 'a'):
                            pass
                    
                    
                            
                            # Run the Moco optimization
                            # solveMocoInverseKinematicsEMG(
                            #     path_model=path_model,
                            #     path_musclepaths=path_musclepaths,
                            #     path_grf=path_grf,
                            #     path_kinematics=path_kinematics,
                            #     path_states=path_states,
                            #     time_final=None,
                            #     time_initial=None,
                            #     force_reserve=20,
                            #     mesh_interval=0.02,
                            #     path_solution_guess=path_solution_guess,
                            #     path_solution=path_solution,
                            #     path_emg=path_emg
                            # )
                            if strike_type in ['ft', 'rt']:
                                muscleDrivenStateTracking(
                                path_model          = path_model,
                                path_musclepaths    = path_musclepaths,
                                path_kinematics     = path_kinematics,
                                path_grf            = path_grf,
                                path_log            = path_log,
                                time_initial        = None,
                                time_final          = None,
                                is_guess            = False,
                                control_condition   = False,
                                reserves            = 50.0,
                                mesh_interval       = 0.1,
                                path_solution       = path_solution, 
                                guess_file          = path_guess,   
                                enable_periodicity=False,
                                periodic_values=False,
                                periodic_speeds = current_speed,
                                periodic_actuators=False,
                                periodic_coordinates=None)
                            
                            elif strike_type in ['fc', 'rc']:
                                muscleDrivenStateTracking(
                                path_model          = path_model,
                                path_musclepaths    = path_musclepaths,
                                path_kinematics     = path_kinematics,
                                path_grf            = path_grf,
                                path_log            = path_log,
                                time_initial        = None,
                                time_final          = None,
                                is_guess            = False,
                                control_condition   = True,
                                reserves            = 50.0,
                                mesh_interval       = 0.1,
                                path_solution       = path_solution, 
                                guess_file          = path_guess,   
                                enable_periodicity=False,
                                periodic_values=False,
                                periodic_speeds = current_speed,
                                periodic_actuators=False,
                                periodic_coordinates=None)

                            # Verify solution was created successfully
                            if os.path.exists(path_solution) and os.path.getsize(path_solution) > 0:
                                print(f" Moco optimization completed successfully")
                                # Save as guess for next iteration
                                solution_guess_history[subject_id] = path_solution
                                successful_processing.append(f"Moco: {base_filename} - Completed")
                            else:
                                print(f"* Moco optimization failed - no solution file generated")
                                failed_processing.append(f"Moco: {base_filename} - No solution generated")
                    else:
                        print(f"    🔄 DRY RUN: Would run Moco optimization")
                        successful_processing.append(f"Moco: {base_filename} - Dry run")
                        
                except Exception as e:
                    print(f"    ❌ Error during Moco optimization: {e}")
                    failed_processing.append(f"Moco: {base_filename} - {e}")
                    # Clean up placeholder file
                    if os.path.exists(path_solution) and os.path.getsize(path_solution) == 0:
                        os.remove(path_solution)
                    continue
    

# input_dir = "/home/lisca/biomec/data/20250306_dataset/results_pool/"
# plot_dir = os.path.join(input_dir, 'plots/')
# avg_inputs = os.path.join(input_dir, 'averages/')
# model_path = "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim"
# # model_path = "/home/lisca/biomec/data/20250306_dataset/sub04/ses20250524/sub04_SCALED_addbio_free_subt_mtp.osim"

# # Create the average inputs first 

# # First call - creates new figure
# fig_emg, ax_emg = retrieve_results.compute_group_averages(
#     input_directory=input_dir,
#     output_directory=plot_dir,
#     sto_output_directory=avg_inputs,
#     file_content_type="emg.sto",
#     strike_type="fc",
#     speed="03000",
#     generate_plot=False,
#     # columns_to_plot=['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
#     columns_to_plot=['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_soleus'],
#     subject_id="03"  # Optional
# )

# fig_ik, ax_ik = retrieve_results.compute_group_averages(
#     input_directory=input_dir,
#     output_directory=plot_dir,
#     sto_output_directory=avg_inputs,
#     file_content_type="grf_move.sto",
#     strike_type="fc",
#     speed="03000",
#     generate_plot=False,
#     # columns_to_plot=['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
#     # columns_to_plot=['ground_1_force_vy', 'ground_2_force_vy', 'ground_1_force_vx', 'ground_2_force_vx'],
#     subject_id="03"  # Optional
# )

# fig_grf, ax_grf = retrieve_results.compute_group_averages(
#     input_directory=input_dir,
#     output_directory=plot_dir,
#     sto_output_directory=avg_inputs,
#     file_content_type="ik_move.sto",
#     strike_type="fc",
#     speed="03000",
#     generate_plot=False,
#     # columns_to_plot=['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
#     # columns_to_plot=['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_soleus'],
#     subject_id="03"  # Optional
# )

# path_data = '/home/lisca/biomec/data/20250306_dataset/'

#     subjects = os.listdir(path_data)
#     # Keep only the directories having the sting 'sub' in their names:
#     # - filter out everything else.
#     subjects = [
#         subject for subject in subjects \
#         if 'sub03' in subject and os.path.isdir(os.path.join(path_data, subject))]
#     subjects.sort()

#     # Initialize the MocoInverse parameters.
#     path_emg              = None
#     path_grf              = None
#     path_cycle_kinematics = None
    
#     # Changes only when a new solution is found. Persists between subjects.
#     path_solution_guess   = None
#     path_solution         = None

#     for sub in subjects:
#         sessions = os.listdir(os.path.join(path_data, sub))
#         # Keep only the directories having the sting 'ses' in their names:
#         # - filter out everything else.
#         sessions = [
#             session for session in sessions \
#             if 'ses' in session and os.path.isdir(os.path.join(path_data, sub, session))]
#         sessions.sort()

#         for ses in sessions:

#             # Extract the model and model's muscle paths corresponding
#             # to the current session.
#             path_model       = None
#             path_musclepaths = None  #  "subject_walk_scaled_FunctionBasedPathSet.xml"
#             for f in os.listdir(os.path.join(path_data, sub, ses)):

#                 if f.endswith('_SCALED_addbio_free_subt_mtp.osim'):
#                     path_model = os.path.join(path_data, sub, ses, f)

#                 if f.endswith('_FunctionBasedPathSet.xml'):
#                     path_musclepaths = os.path.join(path_data, sub, ses, f)

#                 # if f.endswith('grf_move.xml'):
#                 #     path_grf = os.path.join(path_data, sub, ses, f)

#             file_cycles_kinematics = os.listdir(os.path.join(path_data, 'results_pool/'))
#             # Keep only the files having the sting '_ik.sto', and '_emg.sto' in their names:
#             # - filter out everything else.
#             file_cycles_kinematics = [
#                 file_cycle_kinematics for file_cycle_kinematics in file_cycles_kinematics \
#                 if '_ik_move.sto' in file_cycle_kinematics and os.path.isfile(os.path.join(path_data, 'results_pool/', file_cycle_kinematics))]
#             file_cycles_kinematics.sort()

#             # file_cycles_grfs = os.listdir(os.path.join(path_data, sub, ses))
#             # # Keep only the files having the sting '_ik.sto', and '_emg.sto' in their names:
#             # # - filter out everything else.
#             # file_cycles_grfs = [
#             #     file_cycles_grfs for file_cycles_grfs in file_cycles_grfs \
#             #     if '_grf_move.xml' in file_cycles_grfs and os.path.isfile(os.path.join(path_data, sub, ses, file_cycles_grfs))]
#             # file_cycles_grfs.sort()

#             file_count = 0

#             for file_cycle_kinematics in file_cycles_kinematics:
#                 if file_count >=10:
#                     break
#                 file_cycle_stem = file_cycle_kinematics[:-12]
#                 # file_grfs_stem = file_cycles_grfs[:-13]

#                 path_emg              = os.path.join(path_data, 'results_pool/', file_cycle_stem + '_emg.sto')
#                 path_states           = os.path.join(path_data, 'results_pool/', file_cycle_stem + '_rad.sto')
#                 path_grf              = os.path.join(path_data, sub, ses, file_cycle_stem + '_grf_move.xml')
#                 path_cycle_kinematics = os.path.join(path_data, 'results_pool/', file_cycle_kinematics)
#                 path_solution         = os.path.join(path_data, sub, ses, 'guesses/initial/', file_cycle_stem + '_mi_emg.sto')
                


#                 if os.path.exists(path_solution):

#                     print(f'\u001b[32m{path_solution} already exits -> skipping this gait cycle ...\u001b[0m')

#                     # Save it as the last solution, which will be used
#                     # as guess for the next optimization.
#                     if  os.path.getsize(path_solution) > 0 :
#                         path_solution_guess = path_solution
#                         print(f'\u001b[32musing solution as guess: {path_solution_guess}\u001b[0m')

#                     continue
#                 else:
#                     print(f'\u001b[33m{path_solution} does not exist -> running OpenSim Moco on this gait cycle ...\u001b[0m')

#                 if not args.dry_run:

#                     # Run the optimizer only if no --dry-run is given as argument.
#                     print(f'running the optimizer on:\n')
#                     print(f'model:          {path_model}')
#                     print(f'muscle_paths:   {path_musclepaths}')
#                     print(f'grf:            {path_grf}')
#                     print(f'kinematics:     {path_cycle_kinematics}')
#                     print(f'solution_guess: {path_solution_guess}')
#                     print(f'solution:       {path_solution}')

#                     # Mark the trial which is processed by creating a
#                     # place holder for its future solution.
#                     with open(path_solution, 'a'):
#                         pass


# path_model='/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim'
# path_musclepaths='/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_FunctionBasedPathSet.xml'
# path_grf='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub03_strifc_cond00000_speed03000_grf_move_averages.xml'
# path_kinematics='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub03_strifc_cond00000_speed03000_ik_move_averages.sto'
# path_states='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub03_strifc_cond00000_speed03000_rad_averages.sto'
# path_solution='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub03_strifc_cond00000_speed03000_mi_emg_averages.sto'
# path_emg='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub03_strifc_cond00000_speed03000_emg_averages.sto'

# solveMocoInverseKinematicsEMG(
#     path_model=path_model,
#     path_musclepaths=path_musclepaths,
#     path_grf=path_grf,
#     path_kinematics=path_kinematics,
#     path_states=path_states,
#     time_final=None,
#     time_initial=None,
#     force_reserve=20,
#     mesh_interval=0.01,
#     path_solution_guess=None,
#     path_solution=path_solution,
#     path_emg=path_emg
# )