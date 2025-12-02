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
    avg_inputs = os.path.join(input_dir, 'averages/')
    os.makedirs(avg_inputs, exist_ok=True)
    
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
    
    # Track processing results
    successful_processing = []
    failed_processing = []

    
    # Run first iteration of tracking problem using Moco Inverse Kinematics as first guess
    print(f"\n RUNNING MOCO INVERSE KINEMATICS...")
    
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
                # path_guess = os.path.join(avg_inputs, 'guesses/', f"{base_filename}_50res_1mesh.sto")

                path_solution = os.path.join(avg_inputs, 'end_results/', f"{base_filename}_20res_02mesh.sto")
                if not os.path.exists(path_solution):
                    os.makedirs(path_solution, exist_ok=True)
                path_log = os.path.join(avg_inputs, 'end_results/', f"{base_filename}_20res_02mesh.log")
                
                # Check if required input files exist
                required_files = {
                    'kinematics': path_kinematics,
                    'emg': path_emg,
                    'grf': path_grf
                }
                
                
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
                            solveMocoInverseKinematicsEMG(
                                path_model=path_model,
                                path_musclepaths=path_musclepaths,
                                path_grf=path_grf,
                                path_kinematics=path_kinematics,
                                path_states=path_states,
                                time_final=None,
                                time_initial=None,
                                force_reserve=20,
                                mesh_interval=0.02,
                                path_solution_guess=path_solution_guess,
                                path_solution=path_guess,
                                path_emg=path_emg
                            )
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
                                reserves            = 20.0,
                                mesh_interval       = 0.02,
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
                                reserves            = 20.0,
                                mesh_interval       = 0.02,
                                path_solution       = path_solution, 
                                guess_file          = path_guess,   
                                enable_periodicity=False,
                                periodic_values=False,
                                periodic_speeds = current_speed,
                                periodic_actuators=False,
                                periodic_coordinates=None)

                            # Verify solution was created successfully
                            if os.path.exists(path_log) and os.path.getsize(path_log) > 0:
                                print(f" Moco optimization completed successfully")
                                # Save as guess for next iteration
                                solution_guess_history[subject_id] = path_log
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
                    if os.path.exists(path_log) and os.path.getsize(path_log) == 0:
                        os.remove(path_log)
                    continue
    
