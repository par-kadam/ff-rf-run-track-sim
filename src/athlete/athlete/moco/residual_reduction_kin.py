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

"""
This variation is built to run on treadmill stationary data; marker data has been included for tracking
"""
import os
import opensim as osim
import numpy as np
from osimFunctions import kinematicsToStates, addTorqueActuators
import re

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Geometry")

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Models/RajagopalModel/Geometry")

path_opensim_core = '/home/lisca/biomec/src/opensim-core/'
path_data = '/home/lisca/biomec/data/20250306_dataset/'

name_problem = "torque_driven_state_tracking"

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

def set_state_bounds(coord_bounds=None, reference_states_path=None, custom_bounds=None):

    # Define angular radian & velocity values for bounds (in radians for angles, meters for translations)
    pi = np.pi
    coord_bounds = {
        # Pelvis translations (meters)
        '/jointset/ground_pelvis/pelvis_tx/value': [-5, 10], 
        '/jointset/ground_pelvis/pelvis_ty/value': [0, 1.5],
        '/jointset/ground_pelvis/pelvis_tz/value': [-2, 2],
        
        # Pelvis rotations radians)
        '/jointset/ground_pelvis/pelvis_tilt/value': [-90*pi/180, 90*pi/180],
        '/jointset/ground_pelvis/pelvis_list/value': [-90*pi/180, 90*pi/180],
        '/jointset/ground_pelvis/pelvis_rotation/value': [-90*pi/180, 90*pi/180],
        
        # Hip joints radians)
        '/jointset/hip_l/hip_flexion_l/value': [-30*pi/180, 120*pi/180],
        '/jointset/hip_r/hip_flexion_r/value': [-30*pi/180, 120*pi/180],
        '/jointset/hip_l/hip_adduction_l/value': [-30*pi/180, 30*pi/180],
        '/jointset/hip_r/hip_adduction_r/value': [-30*pi/180, 30*pi/180],
        '/jointset/hip_l/hip_rotation_l/value': [-45*pi/180, 45*pi/180],
        '/jointset/hip_r/hip_rotation_r/value': [-45*pi/180, 45*pi/180],
        
        # Knee joints (radians)
        '/jointset/walker_knee_l/knee_angle_l/value':  [0, 120*pi/180],
        '/jointset/walker_knee_r/knee_angle_r/value':  [0, 120*pi/180],
        
        # Ankle joints (radians)
        '/jointset/ankle_l/ankle_angle_l/value':  [-40*pi/180, 30*pi/180],
        '/jointset/ankle_r/ankle_angle_r/value':  [-40*pi/180, 30*pi/180],
        
        # MTP joints  radians) - if present
        '/jointset/mtp_l/mtp_angle_l/value':  [-30*pi/180, 60*pi/180],
        '/jointset/mtp_r/mtp_angle_r/value':  [-30*pi/180, 60*pi/180],
        
        # Lumbar joint  radians) - if present
        '/jointset/back/lumbar_extension/value':  [-30*pi/180, 30*pi/180],
        '/jointset/back/lumbar_bending/value':  [-30*pi/180, 30*pi/180],
        '/jointset/back/lumbar_rotation/value':  [-30*pi/180, 30*pi/180],
        
        # Pelvis translation speeds  m/s)
        '/jointset/ground_pelvis/pelvis_tx/speed':  [-5, 5],  # Forward walking speed
        '/jointset/ground_pelvis/pelvis_ty/speed':  [-3, 3],      # Vertical speed
        '/jointset/ground_pelvis/pelvis_tz/speed':  [-2, 2],     # Lateral speed
    }
    
    reference_data = {}
    if reference_states_path:
        try:
            reference_table = osim.TimeSeriesTable(reference_states_path)
            for coord_path in coord_bounds.keys():
                if reference_table.hasColumn(coord_path):
                    coord_data = reference_table.getDependentColumn(coord_path).to_numpy()
                    reference_data[coord_path] = {
                        'initial': float(coord_data[0]),
                        'final': float(coord_data[-1])
                    }
            print(f"Loaded reference data from: {reference_states_path}")
            print(f"reference has {len(reference_data)} coordinates")
        except Exception as e:
            print(f"Error loading reference file: {e}")
            reference_data = {}

    final_state_bounds = {}

    all_coords = set(coord_bounds.keys())
    if custom_bounds:
            all_coords.update(custom_bounds.keys())

    for coord_path in all_coords:
        if custom_bounds and coord_path in custom_bounds:
            custom_coord_bounds = custom_bounds[coord_path]

            if isinstance(custom_coord_bounds, dict):
                #Dictionary format required: (trajectory_bounds, initial_bounds, final_bounds)
                trajectory_bounds = custom_coord_bounds.get('trajectory', coord_bounds[coord_path], [-10, 10])
                initial_bounds = custom_coord_bounds.get('initial', [])
                final_bounds = custom_coord_bounds.get('final', [])
            elif isinstance(custom_coord_bounds, (tuple, list)) and len(custom_coord_bounds) >= 1:
                #Tuple/list format: (trajectory_bounds, initial_bounds, final_bounds)
                trajectory_bounds = custom_coord_bounds[0]
                initial_bounds = custom_coord_bounds[1] if len(custom_coord_bounds) > 1 else []
                final_bounds = custom_coord_bounds[2] if len(custom_coord_bounds) > 2 else []
            else:
                # If custom bounds are not in the expected format, use default
                trajectory_bounds = coord_bounds.get(coord_path, [-10, 10])
                initial_bounds = []
                final_bounds = []

            if initial_bounds is None:
                initial_bounds = []
            if final_bounds is None:
                final_bounds = []

            final_state_bounds[coord_path] = (trajectory_bounds, initial_bounds, final_bounds)

        else:
            # Otherwise use automatic extraction from reference data
            trajectory_bounds = coord_bounds.get(coord_path, [-10, 10])

            if reference_data and coord_path in reference_data:
                # Reference data for initial/final bounds
                init_val = reference_data[coord_path]['initial']
                final_val = reference_data[coord_path]['final']
                final_state_bounds[coord_path] = (trajectory_bounds, [init_val], [final_val])
            else:
                # If no reference data provided
                final_state_bounds[coord_path] = (trajectory_bounds, [], [])
    
    return final_state_bounds


#Create a dictionary of the optimal forces for actuators used in Hamner & Delp
rraActuators = {'pelvis_tx': 1, 'pelvis_ty': 1, 'pelvis_tz': 1,
                'pelvis_tilt': 1, 'pelvis_list': 1, 'pelvis_rotation': 1,
                'hip_flexion_r': 1000, 'hip_adduction_r': 1000, 'hip_rotation_r': 1000,
                'knee_angle_r': 1000, 'ankle_angle_r': 1000,
                'hip_flexion_l': 1000, 'hip_adduction_l': 1000, 'hip_rotation_l': 1000,
                'knee_angle_l': 1000, 'ankle_angle_l': 1000,
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
             'knee_angle_r': 1, 'ankle_angle_r': 1,
             'hip_flexion_l': 1, 'hip_adduction_l': 1, 'hip_rotation_l': 1,
             'knee_angle_l': 1, 'ankle_angle_l': 1,
             'lumbar_extension': 1, 'lumbar_bending': 1, 'lumbar_rotation': 1,
             'arm_flex_r': 1, 'arm_add_r': 1, 'arm_rot_r': 1,
             'elbow_flex_r': 1, 'pro_sup_r': 1,
             'arm_flex_l': 1, 'arm_add_l': 1, 'arm_rot_l': 1,
             'elbow_flex_l': 1, 'pro_sup_l': 1
             }

#Create a dictionary of the coordinate tasks originally used in Hamner & Delp
rraTasks = {'pelvis_tx': 2.5e2, 'pelvis_ty': 1.0e2, 'pelvis_tz': 2.5e1,
            'pelvis_tilt': 7.5e2, 'pelvis_list': 2.5e2, 'pelvis_rotation': 5.0e1,
            'hip_flexion_r': 7.5e1, 'hip_adduction_r': 5.0e1, 'hip_rotation_r': 1.0e1,
            'knee_angle_r': 1.0e1, 'ankle_angle_r': 1.0e2,
            'hip_flexion_l': 7.5e1, 'hip_adduction_l': 5.0e1, 'hip_rotation_l': 1.0e1,
            'knee_angle_l': 1.0e1, 'ankle_angle_l': 1.0e2,
            'lumbar_extension': 7.5e1, 'lumbar_bending': 5.0e1, 'lumbar_rotation': 2.5e1,
            'arm_flex_r': 1.0e0, 'arm_add_r': 1.0e0, 'arm_rot_r': 1.0e0,
            'elbow_flex_r': 1.0e0, 'pro_sup_r': 1.0e0,
            'arm_flex_l': 1.0e0, 'arm_add_l': 1.0e0, 'arm_rot_l': 1.0e0,
            'elbow_flex_l': 1.0e0, 'pro_sup_l': 1.0e0, 'mtp_angle_r': 1.0e2, 'mtp_angle_l': 1.0e2
            }

#Create a dictionary for kinematic boundary limits (+/- to max and min)
kinematicLimits = {'pelvis_tx': 0.2, 'pelvis_ty': 0.1, 'pelvis_tz': 0.2,
                   'pelvis_tilt': np.deg2rad(10), 'pelvis_list': np.deg2rad(10), 'pelvis_rotation': np.deg2rad(10),
                   'hip_flexion_r': np.deg2rad(10), 'hip_adduction_r': np.deg2rad(5), 'hip_rotation_r': np.deg2rad(5),
                   'knee_angle_r': np.deg2rad(15), 'ankle_angle_r': np.deg2rad(10),
                   'hip_flexion_l': np.deg2rad(10), 'hip_adduction_l': np.deg2rad(5), 'hip_rotation_l': np.deg2rad(5),
                   'knee_angle_l': np.deg2rad(15), 'ankle_angle_l': np.deg2rad(10),
                   'lumbar_extension': np.deg2rad(10), 'lumbar_bending': np.deg2rad(5), 'lumbar_rotation': np.deg2rad(5),
                   'arm_flex_r': np.deg2rad(5), 'arm_add_r': np.deg2rad(5), 'arm_rot_r': np.deg2rad(5),
                   'elbow_flex_r': np.deg2rad(10), 'pro_sup_r': np.deg2rad(5),
                   'arm_flex_l': np.deg2rad(5), 'arm_add_l': np.deg2rad(5), 'arm_rot_l': np.deg2rad(5),
                   'elbow_flex_l': np.deg2rad(10), 'pro_sup_l': np.deg2rad(5)
                   }

def torqueDrivenStateTracking(        
        path_model          = None,
        path_grf            = None,
        path_markers        = None,
        path_kinematics     = None,
        path_states         = None, 
        time_initial        = None,
        time_final          = None,
        force_reserve       = 1.0,
        mesh_interval       = 0.0010,
        path_solution       = None,
        path_guess          = None):

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

    model = modelProcessor.process()
    # model = addTorqueActuators(osimModel=model,
    #                            optForces=rraActuators,
    #                            controlLimits= rraLimits)
    

    track.setModel(osim.ModelProcessor(model))

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

    kin_table = osim.TimeSeriesTable(path_states)
    time = kin_table.getIndependentColumn()
    states_table_processor = osim.TableProcessor(path_states)
    # states_table_processor.append(osim.TabOpUseAbsoluteStateNames())

    track.setStatesReference(states_table_processor)

    coord_bounds = {}
    for coord, limit_value in kinematicLimits.items():
        for coord_ind in range(model.updCoordinateSet().getSize()):
            coord_name = model.updCoordinateSet().get(coord_ind).getName()
            coord_path = model.updCoordinateSet().get(coord_ind).getAbsolutePathString()+'/value'
            
            if coord in coord_path:
                    
                print(f"current coord:", coord, "and its limit differences:", limit_value)
                # coord_path = model.updCoordinateSet().get('pelvis_tz').getAbsolutePathString()+'value'
                print(f"The full coordinate path string:", coord_path)
                coord_data = kin_table.getDependentColumn(coord_path).to_numpy()

                coord_bounds[coord] = [coord_data.min() - limit_value,
                                    coord_data.max() + limit_value]

    # There is marker data in the 'marker_trajectories.trc' associated with
    # model markers that no longer exists (i.e. markers on the arms). Set this
    # flag to avoid an exception from being thrown.
    track.set_allow_unused_references(True)

    # Increase the global marker tracking weight, which is the weight
    # associated with the internal MocoMarkerTrackingCost term.
    track.set_states_global_tracking_weight(10)

    # Increase the tracking weights for individual markers in the data set 
    # placed on bony landmarks compared to markers located on soft tissue.
    track.set_track_reference_position_derivatives(True)
    track.set_apply_tracked_states_to_guess(True)

    state_weights = osim.MocoWeightSet()

    speed_tracking_scale = 0.1

    for coord in range(model.updCoordinateSet().getSize()):

        coord_name = model.updCoordinateSet().get(coord).getName()
        coord_path = model.updCoordinateSet().get(coord).getAbsolutePathString()

        if coord_name in list(rraTasks.keys()):
            state_weights.cloneAndAppend(osim.MocoWeight(f'{coord_path}/value', 
                                                         rraTasks[coord_name]))
            state_weights.cloneAndAppend(osim.MocoWeight(f'{coord_path}/speed', 
                                                         rraTasks[coord_name] * speed_tracking_scale))
            
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_r/value", 20))
    
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_r/value", 20))

    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_r/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_r/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_r/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_r/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_l/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_l/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_l/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_l/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_l/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_l/value", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tilt/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_list/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_rotation/value", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tx/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_ty/value", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tz/value", 30))

    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_r/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_r/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_r/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_r/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_r/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_r/speed", 5))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_flexion_l/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_adduction_l/speed", 20))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/hip_rotation_l/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/knee_angle_l/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/ankle_angle_l/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/hip_r/mtp_angle_l/speed", 5))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tilt/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_list/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_rotation/speed", 10))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tx/speed", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_ty/speed", 30))
    # state_weights.cloneAndAppend(osim.MocoWeight("/jointset/ground_pelvis/pelvis_tz/speed", 30))    

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
        final_index = len(time) - 1
        track.set_final_time(time[final_index])

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

    # speed_bounds = [-20, 20]
    # problem.setStateInfoPattern('/jointset/.*/speed', speed_bounds, [], [])

    # torque_bounds = [-1500, 1500]
    # problem.setStateInfoPattern('/forceset/.*', torque_bounds, [], [])
 

    effort = osim.MocoControlGoal().safeDownCast(problem.updGoal('control_effort'))
    effort.setWeight(10)

    # state_bounds = set_state_bounds(coord_bounds=None, reference_states_path=path_states)
    
    # for path in state_bounds:
    #     traj_bounds = state_bounds[path][0]
    #     init_bounds = state_bounds[path][1] if len(state_bounds[path]) > 1 else []
    #     final_bounds = state_bounds[path][2] if len(state_bounds[path]) > 2 else []
    #     problem.setStateInfo(path, traj_bounds, init_bounds, final_bounds)

    # kinematic bounds set in problem
    for coord in range(model.updCoordinateSet().getSize()):
        if model.updCoordinateSet().get(coord).getName() in list(coord_bounds.keys()):
            coord_name = model.updCoordinateSet().get(coord).getName()
            coord_path = model.updCoordinateSet().get(coord).getAbsolutePathString()+'/value'

            problem.setStateInfo(coord_path, 
                                 [coord_bounds[coord_name][0], coord_bounds[coord_name][1]])
   

    model = modelProcessor.process()
    model.initSystem()
    force_set = model.getForceSet()
    # force_set = model.updForceSet()
    # remove_actuators = ['reserve_pelvis_tx', 'reserve_pelvis_ty']
    remove_model_indices = []
    # for name in remove_actuators:
    #     try:
    #         index = force_set.getIndex(name, 0)
    #         remove_model_indices.append(index)
    #     except:
    #         print(f"Actuator not found: {name}")
    # # Remove in reverse order to void index shifting
    # for index in sorted(remove_model_indices, reverse=True):
    #     success = force_set.remove(index)
    #     if not success:
    #         print(f"Failed to remove actuator at index {index}")
    # print(f"Remaining actuators: {force_set.getSize()}")
    
    # for i in range(force_set.getSize()):
    #     force_path = force_set.get(i).getAbsolutePathString()
    #     if 'pelvis_t' in str(force_path):
    #         # effort.setWeightForControl(force_path, 500)
    #         remove_model_indices.append(i)

    # # Remove indices in reverse order since .remove() function expects indices for forcesets
    # for index in sorted(remove_model_indices, reverse=True):
    #     try:
    #         success = force_set.remove(index)
    #         if success:
    #             print(f"Removed actuator at index {index}")
    #     except Exception as e:
    #         print(f"Failed to remove actuator at index {index}: {e}")

    # Add contact tracking
    contact_tracking = osim.MocoContactTrackingGoal('contact', 7500)
    contact_tracking.setExternalLoadsFile(path_grf)
    right_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'ExternalForce_2', ['/bodyset/calcn_r', '/bodyset/toes_r'])
    left_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'ExternalForce_1', ['/bodyset/calcn_l', '/bodyset/toes_l'])
    contact_tracking.addContactGroup(right_contact_group)
    contact_tracking.addContactGroup(left_contact_group)
    contact_tracking.setNormalizeTrackingError(True)
    contact_tracking.setEnabled(True)

    # problem.addGoal(contact_tracking)


    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.set_optim_max_iterations(10000)
    solver.set_parallel(28)
    
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_convergence_tolerance(1e-2)
    solver.resetProblem(problem)
    # solver.set_num_mesh_intervals(
    #     int(np.round((time_final - time_initial) / mesh_interval)))
    # solver.set_multibody_dynamics_mode('explicit')
    # solver.set_minimize_implicit_auxiliary_derivatives(True)
    # solver.set_implicit_auxiliary_derivatives_weight(aux_deriv_weight / 6.0)
    # solver.set_implicit_auxiliary_derivative_bounds(
    #     osim.MocoBounds(-100.0, 100.0))
    
    # solver.set_parameters_require_initsystem(False)
    
    # solver.set_scale_variables_using_bounds(True)
    # solver.set_optim_finite_difference_scheme('forward')

    guess_file = path_guess
    input_guess = osim.MocoTrajectory(guess_file)
    guess = solver.createGuess()
    # guess.randomizeAdd()
    states_table = input_guess.exportToStatesTable()
    control_table = input_guess.exportToControlsTable()
    state_labels = states_table.getColumnLabels()
    control_labels = control_table.getColumnLabels()

    columns_to_remove = []
    pattern = re.compile(r'/forceset/(.+?)(?!/activation)$')

    for label in state_labels:
        if 'forceset' in label and 'reserve' not in label:
            states_table.removeColumn(label)

    for label in control_labels:
        if 'reserve' in label and 'reserve' not in label:
            control_table.removeColumn(label)

    for label in state_labels:
        if 'normalized_tendon_force' in label:
            states_table.removeColumn(label)

    for label in state_labels:
        if 'subtalar' in label:
            state_labels.removeColumn(label)
        else:
            continue
    for label in control_labels:
        if 'subtalar' in label:
            control_labels.removeColumn(label)
        else:
            continue

    for label in state_labels:
        if pattern.search(label):
            columns_to_remove.append(label)

    for label in control_labels:
        if pattern.search(label):
            columns_to_remove.append(label)

    for label in columns_to_remove:
        if states_table.hasColumn(label):
            states_table.removeColumn(label)
        if control_table.hasColumn(label):
            control_table.removeColumn(label)
                
    guess.insertControlsTrajectory(control_table, True)    
    guess.insertStatesTrajectory(states_table, True)

    solver.setGuess(guess)

    # Solve! Use track.solve() to skip visualizing.
    solution = study.solve()
    solution.write(path_solution)

    # study.visualize(solution)

    model = modelProcessor.process()
    external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)
    external_loads_file_path = path_solution.replace('.sto', '_external_loads.sto')
    osim.STOFileAdapter.write(external_loads, external_loads_file_path)

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
    path_guess            = None
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
                path_solution         = os.path.join(path_data, sub, ses, file_cycle_stem + '_td_states.sto')
                
                path_guess            = os.path.join(path_data, sub, ses, file_cycle_stem + '_mi.sto')
                


                if os.path.exists(path_solution):

                    print(f'\u001b[32m{path_solution} already exits -> skipping this gait cycle ...\u001b[0m')

                    # Save it as the last solution, which will be used
                    # as guess for the next optimization.
                    # if  os.path.getsize(path_solution) > 0 :
                    #     path_solution_guess = path_solution
                    #     print(f'\u001b[32musing solution as guess: {path_solution_guess}\u001b[0m')

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
                    print(f'guess:          {path_guess}')
                    print(f'solution:       {path_solution}')

                    # Mark the trial which is processed by creating a
                    # place holder for its future solution.
                    with open(path_solution, 'a'):
                        pass

                    torqueDrivenStateTracking(
                        path_model          = path_model,
                        path_grf            = path_grf,
                        path_markers        = path_markers,
                        path_kinematics     = path_kinematics,
                        path_states         = path_states,
                        time_initial        = None,
                        time_final          = None,
                        force_reserve       = 500.0,
                        mesh_interval       = 0.20,
                        path_solution       = path_solution,
                        path_guess          = path_guess
                    )