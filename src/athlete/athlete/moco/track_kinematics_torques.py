# -------------------------------------------------------------------------- #
# OpenSim Moco: exampleMocoTrack.py                                          #
# -------------------------------------------------------------------------- #
# Copyright (c) 2019 Stanford University and the Authors                     #
#                                                                            #
# Author(s): Nicholas Bianco                                                 #
#                                                                            #
# Licensed under the Apache License, Version 2.0 (the "License") you may     #
# not use this file except in compliance with the License. You may obtain a  #
# copy of the License at http://www.apache.org/licenses/LICENSE-2.0          #
#                                                                            #
# Unless required by applicable law or agreed to in writing, software        #
# distributed under the License is distributed on an "AS IS" BASIS,          #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   #
# See the License for the specific language governing permissions and        #
# limitations under the License.                                             #
# -------------------------------------------------------------------------- #

# This example features two different tracking problems solved using the
# MocoTrack tool. 
#  - The first problem demonstrates the basic usage of the tool interface
#    to solve a torque-driven marker tracking problem. 
#  - The second problem shows how to customize a muscle-driven state tracking 
#    problem using more advanced features of the tool interface.
# 
# See the README.txt next to this file for more information.

import os
import opensim as osim
import numpy as np
from osimFunctions import kinematicsToStates

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Geometry")

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Models/RajagopalModel/Geometry")

path_opensim_core = '/home/lisca/biomec/src/opensim-core/'
path_data = '/home/lisca/biomec/data/20250306_dataset/'

name_problem = "torque_driven_marker_tracking"

# path_model = os.path.join(
#     path_data,
#     'sub01_SCALED.osim')
# path_grf = os.path.join(
#     path_data,
#     'sub01_strifs_cond00000_speedp0300_0001_0001_grf.xml')
# path_markers = os.path.join(
#     path_data,
#     'sub01_strifs_cond00000_speedp0300_0001_0001_trc.trc')
# path_states = os.path.join(
#     path_data,
#     'sub01_strifs_cond00000_speedp0300_0001_0001_ik.sto')

# path_guess    = None
# path_solution = os.path.join(path_data, 'sub01_strifs_cond00000_speedp0300_tdmts.sto')

def torqueDrivenMarkerTracking(        
        path_model          = None,
        path_grf            = None,
        path_markers        = None,
        path_kinematics     = None,
        path_states         = None, 
        time_initial        = None,
        time_final          = None,
        force_reserve       = 1.0,
        mesh_interval       = 0.0010,
        path_solution       = None):

    # Create and name an instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName(name_problem)

    # Construct a ModelProcessor and add it to the tool. ModelProcessors
    # accept a base model and allow you to easily modify the model by appending
    # ModelOperators. Operations are performed in the order that they are
    # appended to the model.
    # Create the base Model by passing in the model file.
    modelProcessor = osim.ModelProcessor(path_model)

    # Add ground reaction external loads in lieu of a ground-contact model.
    modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))

    # Remove all the muscles in the model's ForceSet.
    modelProcessor.append(osim.ModOpRemoveMuscles())
    # Add CoordinateActuators to the model degrees-of-freedom. This ignores the 
    # pelvis coordinates which already have residual CoordinateActuators.
    modelProcessor.append(osim.ModOpAddReserves(force_reserve))
    track.setModel(modelProcessor)

    # Use this convenience function to set the MocoTrack markers reference
    # directly from a TRC file. By default, the markers data is filtered at
    # 6 Hz and if in millimeters, converted to meters.
    # track.setMarkersReferenceFromTRC(path_markers)

    kinematicsToStates(kinematicsFileName=path_kinematics,
                    osimModelFileName=path_model,
                    outputFileName=path_states,
                    inDegrees=True,
                    outDegrees=False,
                    excludeColumns=['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz'])

    kin_table = osim.TimeSeriesTable(path_kinematics)
    time = kin_table.getIndependentColumn()
    states_table_processor = osim.TableProcessor(path_states)
    states_table_processor.append(osim.TabOpUseAbsoluteStateNames())

    track.setStatesReference(states_table_processor)

    # There is marker data in the 'marker_trajectories.trc' associated with
    # model markers that no longer exists (i.e. markers on the arms). Set this
    # flag to avoid an exception from being thrown.
    track.set_allow_unused_references(True)

    # Increase the global marker tracking weight, which is the weight
    # associated with the internal MocoMarkerTrackingCost term.
    track.set_markers_global_tracking_weight(10)

    # Increase the tracking weights for individual markers in the data set 
    # placed on bony landmarks compared to markers located on soft tissue.
    track.set_track_reference_position_derivatives(True)
    # track.set_apply_tracked_states_to_guess(True)

    state_weights = osim.MocoWeightSet()
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_r/value", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_r/speed", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_r/value", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_r/speed", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_r/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_r/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_r/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_r/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_r/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_r/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_r/value", 5))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_r/speed", 5))

    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_l/value", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_l/speed", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_l/value", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_l/speed", 20))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_l/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_l/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_l/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_l/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_l/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_l/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_l/value", 5))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_l/speed", 5))

    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tilt/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tilt/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_list/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_list/speed", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_rotation/value", 10))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_rotation/speed", 10))

    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tx/value", 30))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tx/speed", 30))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_ty/value", 30))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_ty/speed", 30))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tz/value", 30))
    state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tz/speed", 30))    
    track.set_states_weight_set(state_weights)

    # markerWeights = osim.MocoWeightSet()
    # markerWeights.cloneAndAppend(osim.MocoWeight("RASI", 20))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LASI", 20))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RPSI", 20))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LPSI", 20))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RTHI", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RTH1", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RTH2", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RTH3", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTHI", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTH1", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTH2", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTH3", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RTIB", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RSH1", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RSH2", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTIB", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTH1", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTH2", 5))
    
    # markerWeights.cloneAndAppend(osim.MocoWeight("RKNE", 10))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RANK", 10))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RHEE", 10))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RMT5", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RTOE", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LKNE", 10))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LANK", 10))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LHEE", 10))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LMT5", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LTOE", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RD1M", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LD1M", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LSTL", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LLCA", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RSTL", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RLCA", 5))

    # markerWeights.cloneAndAppend(osim.MocoWeight("C7", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RShoulder", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LShoulder", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LClavicle", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RBiceps", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LBiceps", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RMElbow", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RElbow", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LMElbow", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LElbow", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RForearm", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LForearm", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("RWrist", 5))
    # markerWeights.cloneAndAppend(osim.MocoWeight("LWrist", 5))
    # track.set_markers_weight_set(markerWeights)

    # Initial time, final time, and mesh interval. The number of mesh points
    # used to discretize the problem is computed internally using these values.
    if time_initial is not None:
        track.set_initial_time(time_initial)
    else:
        track.set_initial_time(time[0])

    if time_final is not None:
        track.set_final_time(time_final)
    else:
        track.set_final_time(time[-1])

    track.set_mesh_interval(mesh_interval)


    # Instead of calling solve(), call initialize() to receive a pre-configured
    # MocoStudy object based on the settings above. Use this to customize the
    # problem beyond the MocoTrack interface.
    study = track.initialize()
    problem = study.updProblem()

    speed_bounds = [-20, 20]
    problem.setStateInfoPattern('/jointset/.*/speed', speed_bounds, [], [])

    torque_bounds = [-500, 500]
    problem.setStateInfoPattern('/forceset/.*', torque_bounds, [], [])

    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    # solver.set_optim_constraint_tolerance(1e-4)
    # solver.set_optim_convergence_tolerance(1e-2)
    # solver.set_num_mesh_intervals(
    #     int(np.round((time_final - time_initial) / mesh_interval)))
    # solver.set_multibody_dynamics_mode('explicit')
    # solver.set_minimize_implicit_auxiliary_derivatives(True)
    # solver.set_implicit_auxiliary_derivatives_weight(aux_deriv_weight / 6.0)
    # solver.set_implicit_auxiliary_derivative_bounds(
    #     osim.MocoBounds(-100.0, 100.0))
    solver.set_parallel(28)
    # solver.set_parameters_require_initsystem(False)
    # solver.set_optim_max_iterations()
    # solver.set_scale_variables_using_bounds(True)
    # solver.set_optim_finite_difference_scheme('forward')

    # Solve! Use track.solve() to skip visualizing.
    solution = track.solve()
    solution.write(path_solution)

    # study.visualize(solution)

import argparse 

if __name__ == "__main__":
    # This problem solves in about 5 minutes.

    parser = argparse.ArgumentParser(description="Run the OpenSim Moco optimizer")
    parser.add_argument('-n', '--dry-run', action='store_true', help="dry run")
    args = parser.parse_args()

    path_data = '/home/lisca/biomec/data/20250306_dataset/'

    subjects = os.listdir(path_data)
    # Keep only the directories having the sting 'sub' in their names:
    # - filter out everything else.
    subjects = [
        subject for subject in subjects \
        if 'sub03' in subject and os.path.isdir(os.path.join(path_data, subject))]
    subjects.sort()

    path_grf              = None
    path_cycle_kinematics = None
    
    # Changes only when a new solution is found. Persists between subjects.
    path_solution_guess   = None
    path_solution         = None

    for sub in subjects:
        sessions = os.listdir(os.path.join(path_data, sub))
        # Keep only the directories having the sting 'ses' in their names:
        # - filter out everything else.
        sessions = [
            session for session in sessions \
            if 'ses' in session and os.path.isdir(os.path.join(path_data, sub, session))]
        sessions.sort()

        for ses in sessions:

            # Extract the model and model's muscle paths corresponding
            # to the current session.
            path_model       = None
            path_musclepaths = None  #  "subject_walk_scaled_FunctionBasedPathSet.xml"
            for f in os.listdir(os.path.join(path_data, sub, ses)):

                if f.endswith('_SCALED.osim'):
                    path_model = os.path.join(path_data, sub, ses, f)

                if f.endswith('_FunctionBasedPathSet.xml'):
                    path_musclepaths = os.path.join(path_data, sub, ses, f)

            file_cycles_kinematics = os.listdir(os.path.join(path_data, sub, ses))
            # Keep only the files having the sting '_ik.sto', '_grf.xml', and '_emg.sto' in their names:
            # - filter out everything else.
            file_cycles_kinematics = [
                file_cycle_kinematics for file_cycle_kinematics in file_cycles_kinematics \
                if '_trc.trc' in file_cycle_kinematics and os.path.isfile(os.path.join(path_data, sub, ses, file_cycle_kinematics))]
            file_cycles_kinematics.sort()

            for file_cycle_kinematics in file_cycles_kinematics:
                file_cycle_stem = file_cycle_kinematics[:-8]

                path_grf              = os.path.join(path_data, sub, ses, file_cycle_stem + '_grf.xml')
                path_markers          = os.path.join(path_data, sub, ses, file_cycle_kinematics)
                path_states           = os.path.join(path_data, sub, ses, file_cycle_stem + '_rad.sto')
                path_kinematics       = os.path.join(path_data, sub, ses, file_cycle_stem + '_ik.sto')
                path_solution         = os.path.join(path_data, sub, ses, file_cycle_stem + '_td_kin.sto')
                


                if os.path.exists(path_solution):

                    print(f'\u001b[32m{path_solution} already exits -> skipping this gait cycle ...\u001b[0m')

                    # Save it as the last solution, which will be used
                    # as guess for the next optimization.
                    if  os.path.getsize(path_solution) > 0 :
                        path_solution_guess = path_solution
                        print(f'\u001b[32musing solution as guess: {path_solution_guess}\u001b[0m')

                    continue
                else:
                    print(f'\u001b[33m{path_solution} does not exist -> running OpenSim Moco on this gait cycle ...\u001b[0m')

                if not args.dry_run:

                    # Run the optimizer only if no --dry-run is given as argument.
                    print(f'running the optimizer on:\n')
                    print(f'model:          {path_model}')
                    print(f'muscle_paths:   {path_musclepaths}')
                    print(f'grf:            {path_grf}')
                    print(f'kinematics:     {path_cycle_kinematics}')
                    print(f'solution_guess: {path_solution_guess}')
                    print(f'solution:       {path_solution}')

                    # Mark the trial which is processed by creating a
                    # place holder for its future solution.
                    with open(path_solution, 'a'):
                        pass

                    torqueDrivenMarkerTracking(
                        path_model          = path_model,
                        path_grf            = path_grf,
                        path_markers        = path_markers,
                        path_kinematics     = path_kinematics,
                        path_states         = path_states,
                        time_initial        = None,
                        time_final          = None,
                        force_reserve       = 100.0,
                        mesh_interval       = 0.010,
                        path_solution       = path_solution
                    )