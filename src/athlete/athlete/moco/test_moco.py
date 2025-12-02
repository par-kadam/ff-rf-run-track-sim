# %% Import packages

import opensim as osim
import argparse
import osimFunctions as helper
import os
import pickle
import numpy as np
import time
import re
import shutil
from scipy.interpolate import interp1d
import sys
import matplotlib.pyplot as plt
from moco_tracker import MocoTracker
# import seaborn as sns
# import pandas as pd
import warnings
warnings.simplefilter(action = 'ignore', category = FutureWarning)
    
#Create a dictionary of the coordinate tasks originally used in Hamner & Delp
rraTasks = {'pelvis_tx': 1.0e1, 'pelvis_ty': 1.0e1, 'pelvis_tz': 1.0e1,
            'pelvis_tilt': 7.5e2, 'pelvis_list': 2.5e2, 'pelvis_rotation': 5.0e1,
            'hip_flexion_r': 7.5e1, 'hip_adduction_r': 5.0e1, 'hip_rotation_r': 1.0e1,
            'knee_angle_r': 1.0e1, 'ankle_angle_r': 1.0e1, 'mtp_angle_r': 1.0e1,
            'hip_flexion_l': 7.5e1, 'hip_adduction_l': 5.0e1, 'hip_rotation_l': 1.0e1,
            'knee_angle_l': 1.0e1, 'ankle_angle_l': 1.0e1, 'mtp_angle_l': 1.0e1,
            'lumbar_extension': 7.5e1, 'lumbar_bending': 5.0e1, 'lumbar_rotation': 2.5e1,
            'arm_flex_r': 1.0e0, 'arm_add_r': 1.0e0, 'arm_rot_r': 1.0e0,
            'elbow_flex_r': 1.0e0, 'pro_sup_r': 1.0e0,
            'arm_flex_l': 1.0e0, 'arm_add_l': 1.0e0, 'arm_rot_l': 1.0e0,
            'elbow_flex_l': 1.0e0, 'pro_sup_l': 1.0e0
            }

#Create a dictionary of the optimal forces for actuators used in Hamner & Delp
rraActuators = {'pelvis_tx': 0.01, 'pelvis_ty': 0.01, 'pelvis_tz': 0.01,
                'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1,
                'hip_flexion_r': 1000, 'hip_adduction_r': 1000, 'hip_rotation_r': 1000,
                'knee_angle_r': 1000, 'ankle_angle_r': 1000, 'mtp_angle_r': 1000,
                'hip_flexion_l': 1000, 'hip_adduction_l': 1000, 'hip_rotation_l': 1000,
                'knee_angle_l': 1000, 'ankle_angle_l': 1000, 'mtp_angle_l': 1000,
                'lumbar_extension': 1000, 'lumbar_bending': 1000, 'lumbar_rotation': 1000,
                'arm_flex_r': 500, 'arm_add_r': 500, 'arm_rot_r': 500,
                'elbow_flex_r': 500, 'pro_sup_r': 500,
                'arm_flex_l': 500, 'arm_add_l': 500, 'arm_rot_l': 500,
                'elbow_flex_l': 500, 'pro_sup_l': 500
                }

#Create a dictionary for actuator limits used in Hamner & Delp
rraLimits = {'pelvis_tx': 0.01, 'pelvis_ty': 0.01, 'pelvis_tz': 0.01,
             'pelvis_tilt': 10000, 'pelvis_list': 10000, 'pelvis_rotation': 10000,
             'hip_flexion_r': 1, 'hip_adduction_r': 1, 'hip_rotation_r': 1,
             'knee_angle_r': 1, 'ankle_angle_r': 1, 'mtp_angle_r': 1,
             'hip_flexion_l': 1, 'hip_adduction_l': 1, 'hip_rotation_l': 1,
             'knee_angle_l': 1, 'ankle_angle_l': 1, 'mtp_angle_l': 1,
             'lumbar_extension': 1, 'lumbar_bending': 1, 'lumbar_rotation': 1,
             'arm_flex_r': 1, 'arm_add_r': 1, 'arm_rot_r': 1,
             'elbow_flex_r': 1, 'pro_sup_r': 1,
             'arm_flex_l': 1, 'arm_add_l': 1, 'arm_rot_l': 1,
             'elbow_flex_l': 1, 'pro_sup_l': 1
             }

#Create a dictionary for kinematic boundary limits (+/- to max and min)
kinematicLimits = {'pelvis_tx': 0.2, 'pelvis_ty': 0.1, 'pelvis_tz': 0.2,
                   'pelvis_tilt': np.deg2rad(10), 'pelvis_list': np.deg2rad(10), 'pelvis_rotation': np.deg2rad(10),
                   'hip_flexion_r': np.deg2rad(10), 'hip_adduction_r': np.deg2rad(5), 'hip_rotation_r': np.deg2rad(5),
                   'knee_angle_r': np.deg2rad(15), 'ankle_angle_r': np.deg2rad(10), 'mtp_angle_r': np.deg2rad(10),
                   'hip_flexion_l': np.deg2rad(10), 'hip_adduction_l': np.deg2rad(5), 'hip_rotation_l': np.deg2rad(5),
                   'knee_angle_l': np.deg2rad(15), 'ankle_angle_l': np.deg2rad(10), 'mtp_angle_l': np.deg2rad(10),
                   'lumbar_extension': np.deg2rad(10), 'lumbar_bending': np.deg2rad(5), 'lumbar_rotation': np.deg2rad(5),
                   'arm_flex_r': np.deg2rad(5), 'arm_add_r': np.deg2rad(5), 'arm_rot_r': np.deg2rad(5),
                   'elbow_flex_r': np.deg2rad(10), 'pro_sup_r': np.deg2rad(5),
                   'arm_flex_l': np.deg2rad(5), 'arm_add_l': np.deg2rad(5), 'arm_rot_l': np.deg2rad(5),
                   'elbow_flex_l': np.deg2rad(10), 'pro_sup_l': np.deg2rad(5)
                   }

#Create a list of markers to set as fixed in the generic model
fixedMarkers = ['LASI', 'LPSI', 'LTHI', 'LKNE', 'LTIB', 'LANK', 'LHEE', 'LTOE', 'LMT5', 'LMAnkle', 'LSH1', 'LSH2', 'LSH3', 'LTH1', 'LTH2', 'LTH3', 
                'LClavicle', 'RASI', 'RPSI', 'RTHI', 'RKNE', 'RTIB', 'RANK', 'RHEE', 'RTOE', 'RMT5', 'RMAnkle', 'RSH1', 'RSH2', 'RSH3', 'RTH1', 'RTH2', 'RTH3',
                'LAPX1', 'RAPX1', 'LShoulder', 'LBiceps', 'LElbow', 'LMElbow', 'LForearm', 'LWrist','RShoulder', 'RBiceps', 'RElbow', 'RMElbow', 'RForearm', 'RWrist']


# from moco_tracker import run_moco_tracking

# def runMoco(
#         path_model          = None,
#         path_grf_xml        = None,
#         path_kinematics     = None,
#         path_states         = None,
#         time_initial        = None,
#         time_final          = None,
#         force_reserve       = 1.0,
#         mesh_interval       = 0.10,
#         path_solution_guess = None,
#         path_solution       = None,
#         path_solution_kinematics = None,
#         path_emg            = None,
#         path_log            = None,
#         path_setup_file     = None):
#     """
#     Modified runMoco function using moco_tracker.py
#     """
    
#     # Convert your existing parameters to moco_tracker format
#     output_dir = os.path.dirname(path_solution)
    
#     # Extract cycle name from solution path for configuration naming
#     cycle_name = os.path.basename(path_solution).replace('_moco.sto', '')
    
#     # Define custom weights based on your needs
#     # You can adjust these or make them parameters
#     weights = {
#         'state_tracking_weight': 50,
#         'control_weight': 25,
#         'grf_tracking_weight': 7500,
#         'torso_orientation_weight': 10,
#         'feet_orientation_weight': 10,
#         'horizontal_translation_weight': 100,
#         'aux_deriv_weight': 1000,
#         'acceleration_weight': 1
#     }
    
#     try:
#         # Call the moco_tracker function
#         result = run_moco_tracking(
#             model_path=path_model,
#             coordinates_path=path_states,  # Using states instead of kinematics
#             extloads_path=path_grf_xml,
#             initial_time=time_initial,
#             final_time=time_final,
#             walking_speed=1.25,  # You might want to make this a parameter
#             mesh_interval=mesh_interval,
#             output_dir=output_dir,
#             config_name=cycle_name,
#             weights=weights,
#             guess_path=path_solution_guess,
#             # Additional parameters matching your original solver settings
#             convergence_tolerance=1e-2,
#             constraint_tolerance=1e-2,
#             max_iterations=1000
#         )
        
#         # Convert moco_tracker results to match your original return format
#         return {
#             'solved': result['solved'],
#             'run_time': result['solve_time'],
#             'iterations': result['iterations'],
#             'solution_path': result['solution_path'],
#             'states_path': result['states_path'],
#             'kinematics_path': path_solution_kinematics  # Will be created below
#         }
        
#     except Exception as e:
#         print(f"Error in moco_tracker: {e}")
#         return {
#             'solved': False,
#             'run_time': 0,
#             'iterations': 0,
#             'solution_path': path_solution,
#             'states_path': path_states,
#             'kinematics_path': path_solution_kinematics,
#             'error': str(e)
#         }
    
def runMoco(
        path_model          = None,
        path_grf_xml        = None,
        path_kinematics     = None,
        path_states         = None,
        time_initial        = None,
        time_final          = None,
        force_reserve       = 1.0,
        mesh_interval       = 0.0010,
        path_solution_guess = None,
        path_solution       = None,
        path_solution_kinematics = None,
        path_emg            = None,
        path_log            = None,
        path_setup_file     = None):

    """
    Run Moco tracking optimization for a single motion cycle.
    Parameters:
    -----------
    path_model: str
        Path to the OpenSim model (.osim)
    path_grf_xml: str
        Path to the external loads XML file for ground reaction forces
    path_kinematics: str
        Path to the kinematic data (.sto/.mot)
    time_initial: float, optional
        Start time for analysis (if None, use the start time from kinematics)
    time_final: float, optional
        End time for analysis (if None, use the end time from kinematics)
    force_reserve: float
        Reserve actuator optimal force scaling
    mesh_interval: float
        Time interval for the optimization mesh
    path_solution_guess: str, optional
        Path to a previous solution to use as an initial guess
    path_solution: str
        Path to write the solution file
    path_emg: str, optional
        Path to EMG data for tracking
    path_states: str, optional
        Path to write the states output file
    path_setup_file: str, optional
        Path to write the Moco setup file
    path_log: str, optional
        Path to write the log file
    """

    # %% Run the standard Moco driven residual reduction approach

    #Add in opensim logger
    osim.Logger.removeFileSink()
    osim.Logger.addFileSink(path_log)

    #Create a generic tracking tool to manipulate for the 3 cycles
    mocoTrack = osim.MocoTrack()
    mocoTrack.setName('mocoResidualReduction')
    
    # Construct a ModelProcessor and set it on the tool.
    modelProcessor = osim.ModelProcessor(path_model)
    modelProcessor.append(osim.ModOpAddExternalLoads(path_grf_xml))
    # modelProcessor.append(osim.ModOpRemoveMuscles())
    
    #Process model to edit
    mocoModel = modelProcessor.process()
    
    #Add in torque actuators that replicate the RRA actuators
    mocoModel = helper.addTorqueActuators(osimModel = mocoModel,
                                            optForces = rraActuators,
                                            controlLimits = rraLimits)
    
    #Set model in tracking tool
    mocoTrack.setModel(osim.ModelProcessor(mocoModel))
    
    #Construct a table processor to append to the tracking tool for kinematics
    #The kinematics can't be filtered here with the operator as it messes with
    #time stamps in a funky way. This however has already been done in the 
    #conversion to state coordinates
    tableProcessor = osim.TableProcessor(path_states)
    mocoTrack.setStatesReference(tableProcessor)
    
    #Create a dictionary to set kinematic bounds
    #Create this based on maximum and minimum values in the kinematic data
    #plus/minus some generic values
    
    #Load the kinematics file as a table
    ikTable = osim.TimeSeriesTable(path_states)

        # Get time bounds from the kinematics if not provided
    if time_initial is None:
        time_initial = ikTable.getIndependentColumn()[0]
    if time_final is None:
        time_final = ikTable.getIndependentColumn()[-1]
    
    # Set time bounds in the tracking tool
    mocoTrack.set_initial_time(time_initial)
    mocoTrack.set_final_time(time_final)

    #Create the bounds dictionary
    kinematicBounds = {}
    #Loop through the coordinates
    for coord in kinematicLimits.keys():
        try:
            #Get the coordinate path
            coordPath = mocoModel.updCoordinateSet().get(coord).getAbsolutePathString()+'/value'
            #Set bounds in dictionary
            kinematicBounds[coord] = [
                ikTable.getDependentColumn(coordPath).to_numpy().min() - kinematicLimits[coord],
                ikTable.getDependentColumn(coordPath).to_numpy().max() + kinematicLimits[coord]
                ]
        except:
            continue

    #Set the global states tracking weight in the tracking problem
    mocoTrack.set_states_global_tracking_weight(1)
    
    #Set tracking tool to apply states to guess
    mocoTrack.set_apply_tracked_states_to_guess(False)
    
    #Provide the setting to ignore unused columns in kinematic data
    mocoTrack.set_allow_unused_references(True)
    
    #Set Moco to track the speed derivatives from kinematic data
    #This is probably a fair comparison to RRA given that it tracks accelerations
    #The case study test of this looked like it would help with smoothness as well
    mocoTrack.set_track_reference_position_derivatives(True)
    
    #Set tracking mesh interval time
    mocoTrack.set_mesh_interval(mesh_interval) #### note that this is likely a different time step to RRA
    
    #Set the coordinate reference task weights to match RRA

    #Create weight set for state tracking
    stateWeights = osim.MocoWeightSet()
    
    #Set a scaling factor for tracking the speeds
    speedsTrackingScale = 0.01
    
    #Loop through coordinates to apply weights
    for coordInd in range(mocoModel.updCoordinateSet().getSize()):
        
        #Get name and absolute path to coordinate
        coordName = mocoModel.updCoordinateSet().get(coordInd).getName()
        coordPath = mocoModel.updCoordinateSet().get(coordInd).getAbsolutePathString()
    
        #If a task weight is provided, add it in
        if coordName in list(rraTasks.keys()):
            #Append state into weight set
            #Track the coordinate value
            stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/value',
                                                        rraTasks[coordName]))
            #Don't track the Coordinate speed
            stateWeights.cloneAndAppend(osim.MocoWeight(f'{coordPath}/speed',
                                                        rraTasks[coordName] * speedsTrackingScale))
            
    #Add state weights to the tracking tool
    mocoTrack.set_states_weight_set(stateWeights)

    #Initialise the Moco study
    study = mocoTrack.initialize()
    problem = study.updProblem()
    
    #Set the parameters for the regularization term on MocoTrack problem
    #(minimize squared excitations)
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.001)
            
    #Lock time bounds to the IK data
    problem.setTimeBounds(time_initial, time_final)
    
    # Set kinematic bounds
    for coordInd in range(mocoModel.updCoordinateSet().getSize()):
        coordName = mocoModel.updCoordinateSet().get(coordInd).getName()
        if coordName in kinematicBounds:
            coordPath = mocoModel.updCoordinateSet().get(coordInd).getAbsolutePathString() + '/value'
            problem.setStateInfo(
                coordPath,
                [kinematicBounds[coordName][0], kinematicBounds[coordName][1]]
            )

    # Configure solver
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.set_optim_max_iterations(1000)
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_convergence_tolerance(1e-2)
    
    # Set initial guess if provided
    if path_solution_guess is not None and os.path.exists(path_solution_guess):
        solver.setGuessFile(path_solution_guess)
    
    # Reset problem (required for solver settings to take effect)
    solver.resetProblem(problem)

    # Save setup file
    study.printToXML(path_setup_file)
    
    # Solve the optimization problem
    print(f"Running Moco optimization from {time_initial} to {time_final}...")
    start_time = time.time()
    solution = study.solve()
    run_time = round(time.time() - start_time, 2)
    print(f"Optimization completed in {run_time} seconds")
    
    # Check solution status
    if solution.isSealed():
        solution.unseal()
        print("Warning: Solution is sealed, optimization may not have converged")
        solved = False
    else:
        print(f"Optimization successful in {solution.getNumIterations()} iterations")
        solved = True
    
    # Write solution to file
    solution.write(path_solution)
    
    # Export states
    osim.STOFileAdapter().write(solution.exportToStatesTable(), path_states)
    
    # Close logger
    osim.Logger.removeFileSink()
    
    # Return key results
    return {
        'solved': solved,
        'run_time': run_time,
        'iterations': solution.getNumIterations(),
        'solution_path': path_solution,
        'states_path': path_solution,
        'kinematics_path': path_solution_kinematics
    }



if __name__ == "__main__":
    
    # # Configure output directory
    # output_dir = os.path.dirname(path_solution)
    # os.makedirs(output_dir, exist_ok=True)
    
    # # # Create a unique output identifier from the solution path
    # # file_stem = os.path.basename(path_solution).split('.')[0]
    
    # # Convert kinematics to states for Moco
    # if path_states is None:
    #     path_states = os.path.join(output_dir, f"{file_stem}_states.sto")
        
    # This problem solves in about 5 minutes.

    parser = argparse.ArgumentParser(description="Run the OpenSim Moco optimizer")
    parser.add_argument('-n', '--dry-run', action='store_true', help="dry run")
    args = parser.parse_args()

    path_data = '/home/lisca/biomec/data/20240904_dataset/'

    subjects = os.listdir(path_data)
    # Keep only the directories having the sting 'sub' in their names:
    # - filter out everything else.
    subjects = [
        subject for subject in subjects
        if 'sub' in subject and os.path.isdir(os.path.join(path_data, subject))]
    subjects.sort()

    moco_results = {}

    path_solution_guess = None

    # Initialize the MocoInverse parameters.
    path_emg              = None
    path_grf_xml          = None
    path_cycle_kinematics = None
    path_states           = None

    # Changes only when a new solution is found. Persists between subjects.
    path_solution_guess   = None
    path_solution         = None

    path_solution_kinematics= None
    path_log              = None
    

    for sub in subjects:
        moco_results[sub] = {}

        sessions = os.listdir(os.path.join(path_data, sub))
        # Keep only the directories having the sting 'ses' in their names:
        # - filter out everything else.
        sessions = [
            session for session in sessions \
            if 'ses' in session and os.path.isdir(os.path.join(path_data, sub, session))]
        sessions.sort()

        for ses in sessions:
            moco_results[sub][ses] = {}
            session_dir = os.path.join(path_data, sub, ses)

            # Extract the model and model's muscle paths corresponding
            # to the current session.
            path_model       = None
            path_musclepaths = None  #  "subject_walk_scaled_FunctionBasedPathSet.xml"
            for f in os.listdir(os.path.join(path_data, sub, ses)):

                if f.endswith('_SCALED.osim'):
                    path_model = os.path.join(path_data, sub, ses, f)

                if f.endswith('_FunctionBasedPathSet.xml'):
                    path_musclepaths = os.path.join(path_data, sub, ses, f)

            if path_model is None:
                print(f"No model file found for {sub}/{ses}, skipping...")
                continue

            # Keep only the files having the sting '_ik.sto', '_grf.xml', and '_emg.sto' in their names:
            # - filter out everything else.
            file_cycles_kinematics = [
                file_cycle_kinematics for file_cycle_kinematics in os.listdir(session_dir) \
                if '_ik.sto' in file_cycle_kinematics and os.path.isfile(os.path.join(session_dir, file_cycle_kinematics))]
            file_cycles_kinematics.sort()

            # Process each kinematic file (gait cycle)
            for file_cycle_kinematics in file_cycles_kinematics:
                # Extract the file stem to build related file paths
                file_cycle_stem = file_cycle_kinematics[:-7]  # Remove '_ik.sto'
                
                # Define all file paths needed for this cycle
                path_cycle_kinematics = os.path.join(session_dir, file_cycle_kinematics)
                path_grf_xml = os.path.join(session_dir, f"{file_cycle_stem}_grf.xml")
                path_emg = os.path.join(session_dir, f"{file_cycle_stem}_emg.sto")
                path_solution = os.path.join(session_dir, f"{file_cycle_stem}_moco.sto")

                path_states = os.path.join(session_dir, f"{file_cycle_stem}_rad.sto")
                path_log = os.path.join(session_dir, f"{file_cycle_stem}_moco.log")
                path_setup = os.path.join(session_dir, f"{file_cycle_stem}_moco_setup.omoco")

                skip_preprocessing = (
                    os.path.exists(path_grf_xml) and 
                    os.path.exists(path_states) and 
                    os.path.getsize(path_grf_xml) > 0 and 
                    os.path.getsize(path_states) > 0
                )
                
                # Check if solution already exists
                if os.path.exists(path_solution) and os.path.getsize(path_solution) > 0:
                    print(f"{path_solution} already exists, skipping...")
                    # Update solution guess for next run
                    path_solution_guess = path_solution
                    continue
     
                required_inputs_exist = (
                    os.path.exists(path_grf_xml) and os.path.getsize(path_grf_xml) > 0 and
                    os.path.exists(path_states) and os.path.getsize(path_states) > 0
                )
                                
                print(f"Processing: {sub}/{ses}/{file_cycle_kinematics}")
                print(f"  Model: {path_model}")
                print(f"  Kinematics: {path_cycle_kinematics}")
                print(f"  GRF: {path_grf_xml}")
                print(f"  EMG: {path_emg}")
                print(f"  Solution: {path_solution}")
                
                # Create placeholder file to mark as in-progress
                if not args.dry_run:
                    with open(path_solution, 'a'):
                        pass
                    
                    # Run Moco optimization
                    try:
                        if not os.path.exists(path_states) or os.path.getsize(path_states) == 0:
                            print("Creating states file from kinematics...")
                            helper.kinematicsToStates(
                                kinematicsFileName=path_cycle_kinematics,
                                outputFileName=path_states, 
                                osimModelFileName=path_model,
                                inDegrees=True, outDegrees=False, filtFreq=None
                            )
                        else:
                            print("States file already exists, skipping creation...")

                        results = runMoco(
                            path_model=path_model,
                            path_grf_xml=path_grf_xml,
                            path_kinematics=path_cycle_kinematics,
                            path_states=path_states,
                            path_emg=path_emg,
                            time_initial=None,
                            time_final=None,
                            path_solution_guess=path_solution_guess,
                            path_solution=path_solution,
                            path_solution_kinematics=path_solution_kinematics,
                            path_setup_file=path_setup,
                            path_log=path_log,
                            force_reserve=1.0,
                            mesh_interval=0.010
                        )
                        
                        
                        # Convert states back to kinematics
                        path_solution_kinematics = os.path.join(session_dir, f"{file_cycle_stem}_mocokin.sto")
                        helper.process_sto_file(
                            path_solution, path_model,
                            path_solution_kinematics,
                            inDegrees=False, outDegrees=True)
                        # Store results
                        moco_results[sub][ses][file_cycle_stem] = results

                        # Update solution guess if successful
                        if results['solved']:
                            path_solution_guess = path_solution
                            
                    except Exception as e:
                        print(f"Error processing {file_cycle_stem}: {e}")
                        # Optionally remove the placeholder file on error
                        if os.path.exists(path_solution) and os.path.getsize(path_solution) == 0:
                            os.remove(path_solution)
