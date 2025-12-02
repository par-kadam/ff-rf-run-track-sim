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
import pandas as pd
import argparse 
from osimFunctions import kinematicsToStates, addTorqueActuators, compute_sto_stdevs
from retrieve_results import analyze_moco_ik_tracking_errors
import re

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Geometry")

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Models/RajagopalModel/Geometry")

path_opensim_core = '/home/lisca/biomec/src/opensim-core/'
path_data = '/home/lisca/biomec/data/20250306_dataset/'

name_problem = "residual_reduction_marker_tracking"

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

#Create a dictionary of the optimal forces for actuators used in Hamner & Delp
rraActuators = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
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
rraLimits = {'pelvis_tx': 10000, 'pelvis_ty': 10000, 'pelvis_tz': 10000,
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
                   'knee_angle_r': np.deg2rad(15), 'ankle_angle_r': np.deg2rad(10), 'mtp_angle_r': np.deg2rad(5),
                   'hip_flexion_l': np.deg2rad(10), 'hip_adduction_l': np.deg2rad(5), 'hip_rotation_l': np.deg2rad(5),
                   'knee_angle_l': np.deg2rad(15), 'ankle_angle_l': np.deg2rad(10), 'mtp_angle_l': np.deg2rad(5),
                   'lumbar_extension': np.deg2rad(10), 'lumbar_bending': np.deg2rad(5), 'lumbar_rotation': np.deg2rad(5),
                   'arm_flex_r': np.deg2rad(5), 'arm_add_r': np.deg2rad(5), 'arm_rot_r': np.deg2rad(5),
                   'elbow_flex_r': np.deg2rad(10), 'pro_sup_r': np.deg2rad(5),
                   'arm_flex_l': np.deg2rad(5), 'arm_add_l': np.deg2rad(5), 'arm_rot_l': np.deg2rad(5),
                   'elbow_flex_l': np.deg2rad(10), 'pro_sup_l': np.deg2rad(5)
                   }


def modify_trc_markers_with_accumulated_distance(input_file=None, output_file=None):
   
    """
    Read TRC file using TRCFileAdapter and modify x-coordinates of all markers
    by adding accumulated distance based on trial speed and time.
    
    This function creates a TRC file that is fully compatible with OpenSim's
    setMarkersReferenceFromTRC() function.
    
    Parameters:
    -----------
    trc_file_path : str
        Path to the input .trc file
    output_file_path : str, optional
        Path to save the modified TRC file. If None, creates a new file with '_modified' suffix
    
    Returns:
    --------
    output_file_path : str
        Path to the modified TRC file (compatible with setMarkersReferenceFromTRC)
    trial_speed : float
        Extracted trial speed from filename
    modified_table : osim.TimeSeriesTableVec3
        Modified marker data table
    """

    # Read TRC file
    marker_table = osim.TimeSeriesTableVec3(input_file)
    
    # Get table information
    num_frames = marker_table.getNumRows()
    num_markers = marker_table.getNumColumns()
    marker_names = marker_table.getColumnLabels()
    times = marker_table.getIndependentColumn()

    # get speed from name string alone
    file_name = os.path.basename(input_file)
    speed_parse = re.search('speed(\d+)', file_name)
    if speed_parse:
        speed_str = speed_parse.group(1)
        trial_speed = float(speed_str) / 1000
        print(f"Extracted trial speed: {trial_speed} m/s from pattern '{speed_parse.group(0)}' in filename: {file_name}")
    else:
        raise ValueError("No speed pattern found")
    # if len(file_name) < 31:
    #     raise ValueError(f"Filename too short to extract speed from 31st character: {file_name}")
    # try:
    #     # assuming character at pos 31 denotes speed x.0 m/s
    #     speed_char = file_name[32]
    #     trial_speed = float(speed_char)
    #     print(f"Extracted trial speed: {trial_speed} from filename: {file_name}")
    # except (ValueError, IndexError) as e:
    #     raise ValueError(f"Could not extract trial speed from 31st character of filename: {e}")
    
    # original_metadata = marker_table.getTableMetaData()
    # print(f"Processing {num_frames} frames with {num_markers} markers")

    # # Verify essential metadata exists
    # required_keys = ['DataRate', 'Units', 'NumMarkers', 'NumFrames']
    # for key in required_keys:
    #     if not original_metadata.hasKey(key):
    #         print(f"Warning: Missing required metadata key: {key}")
    
    # Iterate through all data points
    for frame_idx in range(num_frames):
        current_time = times[frame_idx]
        
        dx = current_time * trial_speed

        for marker_idx in range(num_markers):
            marker_name = marker_names[marker_idx]

            marker_column = marker_table.updDependentColumnAtIndex(marker_idx)

            row_data = marker_table.getRowAtIndex(frame_idx)
            position = marker_column.getElt(frame_idx, 0)
            
            # Extract coordinates
            x, y, z = position.get(0), position.get(1), position.get(2)


            # Apply modification to X coordinate
            modified_x = x + dx
            new_position = osim.Vec3(modified_x, y, z)

            # updRowAtIndex to get a mutable reference and modify
            # mutable_row = marker_table.updRowAtIndex(frame_idx)
            # mutable_row.setElt(0, marker_idx, new_position)
            # Update position within column
            marker_column.updElt(frame_idx, 0).set(0, modified_x)
            marker_column.updElt(frame_idx, 0).set(1, y)
            marker_column.updElt(frame_idx, 0).set(2, z)
            
            # Print details for first few frames and first marker for debugging
            if frame_idx < 3 and marker_idx == 0:
                print(f"Frame {frame_idx}, Time {current_time:.3f}s, {marker_name}:")
                print(f"  Offset = {trial_speed} * {current_time:.3f} = {dx:.3f}")
                print(f"  X: {x:.2f} -> {modified_x:.2f}")

    # The table should already have the correct metadata, but let's verify
    # if not marker_table.getTableMetaData().hasKey('DataRate'):
    #     print("Adding missing DataRate metadata")
    #     marker_table.addTableMetaDataString('DataRate', '1000')  # Default value
    
    # if not marker_table.getTableMetaData().hasKey('Units'):
    #     print("Adding missing Units metadata")
    #     marker_table.addTableMetaDataString('Units', 'm')  # Default value

    # Write modified data
    adapter = osim.TRCFileAdapter()
    adapter.write(marker_table, output_file)
    print(f"Modified data written to {output_file}")
    

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
    # modelProcessor.append(osim.ModOpAddReserves(force_reserve))
    model = modelProcessor.process()
    model = addTorqueActuators(osimModel=model,
                               optForces=rraActuators,
                               controlLimits=rraLimits)
    
    track.setModel(osim.ModelProcessor(model))

    markers_out = path_solution.replace('.trc', '_trc_move.trc')
    #Forward translate all markers?
    # path_markers = modify_trc_markers_with_accumulated_distance(path_markers, output_file=markers_out)
    # Use this convenience function to set the MocoTrack markers reference
    # directly from a TRC file. By default, the markers data is filtered at
    # 6 Hz and if in millimeters, converted to meters.
    track.setMarkersReferenceFromTRC(path_markers)

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
    track.set_states_global_tracking_weight(1)

    # Increase the tracking weights for individual markers in the data set 
    # placed on bony landmarks compared to markers located on soft tissue.
    track.set_track_reference_position_derivatives(True)
    # track.set_apply_tracked_states_to_guess(True)

    state_weights = osim.MocoWeightSet()
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_r/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_r/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_r/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_r/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_r/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_r/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/walker_knee_r/knee_angle_r/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/walker_knee_r/knee_angle_r/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/patellofemoral_r/knee_angle_r_beta/value", 5))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/patellofemoral_r/knee_angle_r_beta/speed", 5))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_r/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_r/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_r/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_r/speed", 10))

    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_l/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_l/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_l/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_l/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_l/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_l/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/walker_knee_l/knee_angle_l/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/walker_knee_l/knee_angle_l/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/patellofemoral_l/knee_angle_l_beta/value", 5))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/patellofemoral_l/knee_angle_l_beta/speed", 5))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_l/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_l/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_l/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_l/speed", 10))

    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tilt/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tilt/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_list/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_list/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_rotation/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_rotation/speed", 10))

    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tx/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tx/speed", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_ty/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_ty/speed", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tz/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tz/speed", 30))    
    
    temp_stdev_path = path_solution.replace('.sto', 'stdev.csv')
    compute_sto_stdevs(stoFileName=path_states, outputCsvFileName=temp_stdev_path)

    ref_table_processor = osim.TableProcessor(osim.TimeSeriesTable(path_states))
    ref_table_processor.append(osim.TabOpUseAbsoluteStateNames())

    # Create state tracking weights using dict above
    #-------------------------------------------
    # states_table = osim.TableProcessor(converted_states)
    # model = modelProcessor.process()
    model.initSystem()
    model_coordinate_paths = ref_table_processor.process(model).getColumnLabels()

    ##  STILL HAVE TO GENERATE CSV STDEVS FILE IN COORDINATE CONVERSION
    coordinates_std = pd.read_csv(temp_stdev_path, index_col=0)
    state_weights = osim.MocoWeightSet()
    for valuePath in model_coordinate_paths:
        speed_path = valuePath.replace('/value', '/speed')
        value_weight = 1.0
        speed_weight = 1.0

        for name in coordinates_std.index:
            stdev = coordinates_std.loc[name][0]
            denom = 1 * stdev

            if name in valuePath:
                if 'lumbar' in name:
                    # if 'torso_orientation_weight' in weights:
                    #     value_weight = 1.0/denom
                    #     speed_weight = 0.01/denom
                    # else:
                    value_weight = 1/denom
                    speed_weight = 0.01/denom

                elif('beta' in name or
                        'subtalar' in name or 
                        'wrist' in name):
                        value_weight = 0.0
                        speed_weight = 0.0

                elif ('ankle' in name or
                        'mtp' in name):
                    value_weight = 2.0/denom
                    speed_weight = 0.02/denom

                elif 'pelvis' in name:
                    if 'pelvis_ty' in name:
                        value_weight = 0.0/denom
                        speed_weight = 0.02/denom
                    else:
                        value_weight = 1.0/denom
                        speed_weight = 0.01/denom

                else:
                    value_weight = 1.0/denom
                    speed_weight = 0.01/denom

        state_weights.cloneAndAppend(
            osim.MocoWeight(valuePath, value_weight))
        state_weights.cloneAndAppend(
            osim.MocoWeight(speed_path, speed_weight))
    
    track.set_states_weight_set(state_weights)

    

    markerWeights = osim.MocoWeightSet()
    markerWeights.cloneAndAppend(osim.MocoWeight("RASI", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LASI", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("RPSI", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LPSI", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("RTHI", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RTH1", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RTH2", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RTH3", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTHI", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTH1", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTH2", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTH3", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RTIB", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RSH1", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RSH2", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTIB", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTH1", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTH2", 50))
    
    markerWeights.cloneAndAppend(osim.MocoWeight("RKNE", 100))
    markerWeights.cloneAndAppend(osim.MocoWeight("RANK", 100))
    markerWeights.cloneAndAppend(osim.MocoWeight("RHEE", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("RMT5", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("RTOE", 500))
    markerWeights.cloneAndAppend(osim.MocoWeight("LKNE", 100))
    markerWeights.cloneAndAppend(osim.MocoWeight("LANK", 100))
    markerWeights.cloneAndAppend(osim.MocoWeight("LHEE", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LMT5", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LTOE", 500))
    markerWeights.cloneAndAppend(osim.MocoWeight("RD1M", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LD1M", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LSTL", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("LLCA", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("RSTL", 200))
    markerWeights.cloneAndAppend(osim.MocoWeight("RLCA", 200))

    markerWeights.cloneAndAppend(osim.MocoWeight("C7", 500))
    markerWeights.cloneAndAppend(osim.MocoWeight("RShoulder", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LShoulder", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LClavicle", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RBiceps", 5))
    markerWeights.cloneAndAppend(osim.MocoWeight("LBiceps", 5))
    markerWeights.cloneAndAppend(osim.MocoWeight("RMElbow", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RElbow", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LMElbow", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LElbow", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("RForearm", 5))
    markerWeights.cloneAndAppend(osim.MocoWeight("LForearm", 5))
    markerWeights.cloneAndAppend(osim.MocoWeight("RWrist", 50))
    markerWeights.cloneAndAppend(osim.MocoWeight("LWrist", 50))
    track.set_markers_weight_set(markerWeights)

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

    # Define contact sphere force names for right foot
    # Adjust these based on your specific model's contact sphere names
    forceNamesRightFoot = [
        'forceset/contactHeel_r',
        'forceset/contactLateralRearfoot_r', 
        'forceset/contactLateralMidfoot_r',
        'forceset/contactLateralToe_r',
        'forceset/contactMedialToe_r',
        'forceset/contactMedialMidfoot_r'
    ]
    
    # Define contact sphere force names for left foot
    forceNamesLeftFoot = [
        'forceset/contactHeel_l',
        'forceset/contactLateralRearfoot_l',
        'forceset/contactLateralMidfoot_l', 
        'forceset/contactLateralToe_l',
        'forceset/contactMedialToe_l',
        'forceset/contactMedialMidfoot_l'
    ]

    # Instead of calling solve(), call initialize() to receive a pre-configured
    # MocoStudy object based on the settings above. Use this to customize the
    # problem beyond the MocoTrack interface.
    study = track.initialize()
    problem = study.updProblem()

    # speed_bounds = [-40, 40]
    # problem.setStateInfoPattern('/jointset/.*/speed', speed_bounds, [], [])

    # torque_bounds = [-1500, 1500]
    # problem.setStateInfoPattern('/forceset/.*', torque_bounds, [], [])

    # Insert control effort goal
    effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(0.001)

    # Add contact tracking
    # contact_tracking = osim.MocoContactTrackingGoal('contact', 7500)
    # contact_tracking.setExternalLoadsFile(path_grf)
    # right_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'ExternalForce_2', ['/bodyset/calcn_r', '/bodyset/toes_r'])
    # left_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'ExternalForce_1', ['/bodyset/calcn_l', '/bodyset/toes_l'])
    # contact_tracking.addContactGroup(right_contact_group)
    # contact_tracking.addContactGroup(left_contact_group)
    # contact_tracking.setNormalizeTrackingError(True)
    # contact_tracking.setEnabled(True)

    # problem.addGoal(contact_tracking)

    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_convergence_tolerance(1e-2)
    # solver.set_num_mesh_intervals(
    #     int(np.round((time_final - time_initial) / mesh_interval)))
    # solver.set_multibody_dynamics_mode('explicit')
    # solver.set_minimize_implicit_auxiliary_derivatives(True)
    # solver.set_implicit_auxiliary_derivatives_weight(aux_deriv_weight / 6.0)
    # solver.set_implicit_auxiliary_derivative_bounds(
    #     osim.MocoBounds(-100.0, 100.0))
    solver.set_optim_max_iterations(10000)
    solver.set_parallel(28)
    # solver.set_parameters_require_initsystem(False)
    
    # solver.set_scale_variables_using_bounds(True)
    # solver.set_optim_finite_difference_scheme('forward')
    solver.resetProblem(problem)

    # Solve! Use track.solve() to skip visualizing.
    solution = study.solve()
    solution.write(path_solution)

    # study.visualize(solution)
    path_problem = path_solution.replace('.sto', '_setup.omoco')
    study.printToXML(path_problem)


    model = modelProcessor.process()
    external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)
    external_loads_file_path = path_solution.replace('.sto', '_external_loads.sto')
    osim.STOFileAdapter.write(external_loads, external_loads_file_path)

    path_error_report = path_solution.replace('.sto', '_errs.sto')
    analyze_moco_ik_tracking_errors(moco_solution_file=path_solution, 
                                    ik_kinematics_file=path_states, 
                                    output_file=path_error_report, 
                                    model_file=path_model,
                                    omoco_setup_file=path_problem)

    force_paths_ankle_muscles = [
    '/forceset/gasmed_r',
    '/forceset/gaslat_r', 
    '/forceset/soleus_r',
    '/forceset/tibant_r',
    '/forceset/gasmed_l',
    '/forceset/gaslat_l',
    '/forceset/soleus_l',
    '/forceset/tibant_l']

    # Get all forceset/reserve values from the solution
    solution_trajectory = osim.MocoTrajectory(path_solution)
    net_moments = study.calcGeneralizedForces(solution_trajectory, force_paths_ankle_muscles)
    net_moments_file_path = path_solution.replace('.sto', '_net_moments.sto')
    osim.STOFileAdapter.write(net_moments, net_moments_file_path)


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

                if f.endswith('_SCALED_addbiomech.osim'):
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
                path_solution         = os.path.join(path_data, sub, ses, file_cycle_stem + '_td_marks.sto')
                


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
                        mesh_interval       = 0.2,
                        path_solution       = path_solution
                    )