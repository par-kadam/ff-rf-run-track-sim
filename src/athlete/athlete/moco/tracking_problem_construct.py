import os
import opensim as osim
import numpy as np
from osimFunctions import kinematicsToStates, compute_sto_stdevs, generate_guess_file, statesToKinematics
import retrieve_results as get_results
import pandas as pd
import datetime
import sys
from datetime import datetime
import validation as validate

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Geometry")


#Create a dictionary for kinematic boundary limits (+/- to max and min)
kinematicLimits = {'pelvis_tx': 0.2, 'pelvis_ty': 0.1, 'pelvis_tz': 0.2,
                   'pelvis_tilt': np.deg2rad(10), 'pelvis_list': np.deg2rad(10), 'pelvis_rotation': np.deg2rad(10),
                   'hip_flexion_r': np.deg2rad(10), 'hip_adduction_r': np.deg2rad(5), 'hip_rotation_r': np.deg2rad(5),
                   'knee_angle_r': np.deg2rad(15), 'ankle_angle_r': np.deg2rad(10), 'subtalar_angle_r': np.deg2rad(10), 'mtp_angle_r': np.deg2rad(10),
                   'hip_flexion_l': np.deg2rad(10), 'hip_adduction_l': np.deg2rad(5), 'hip_rotation_l': np.deg2rad(5),
                   'knee_angle_l': np.deg2rad(15), 'ankle_angle_l': np.deg2rad(10), 'subtalar_angle_l': np.deg2rad(10), 'mtp_angle_l': np.deg2rad(10),
                   'lumbar_extension': np.deg2rad(10), 'lumbar_bending': np.deg2rad(5), 'lumbar_rotation': np.deg2rad(5),
                   'arm_flex_r': np.deg2rad(5), 'arm_add_r': np.deg2rad(5), 'arm_rot_r': np.deg2rad(5),
                   'elbow_flex_r': np.deg2rad(10), 'pro_sup_r': np.deg2rad(5),
                   'arm_flex_l': np.deg2rad(5), 'arm_add_l': np.deg2rad(5), 'arm_rot_l': np.deg2rad(5),
                   'elbow_flex_l': np.deg2rad(10), 'pro_sup_l': np.deg2rad(5)
                   }

def set_state_bounds(coord_bounds=None, reference_states_path=None, custom_bounds=None, kinematicLimits=None):

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

                    coord_name = coord_path.split('/')[-2]
                    if kinematicLimits and coord_name in kinematicLimits:
                        # Get new upper and lower bounds using inputted deg2rad vals from kinematicLimits dict
                        min_val = coord_data.min() - 2*kinematicLimits[coord_name]
                        max_val = coord_data.max() + 2*kinematicLimits[coord_name]
                        # Update coord_bounds with calculated values for computing upper & lower bounds
                        coord_bounds[coord_path] = [min_val, max_val]

                        print(f"Updated bounds for {coord_name}: [{min_val:.4f}, {max_val:.4f}]")
                    else:
                        coord_bounds[coord_path] = coord_bounds[coord_path]

                    # Storing reference data for initial/final bounds as their unchanged vlaues
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



# def muscleDrivenStateTracking(
#                         path_model          = path_model,
#                         path_musclepaths    = path_musclepaths,
#                         path_kinematics     = path_kinematics,
#                         path_grf            = path_grf,
#                         path_log            = None,  
#                         time_initial        = None,
#                         time_final          = None,
#                         reserves            = 100.0,
#                         mesh_interval       = 0.010,
#                         guess_file          = None,   
#                         path_solution       = path_solution, 
#                         is_guess            = True,
#                         control_condition   = True,
#                         enable_periodicity=False,
#                         periodic_values=False,
#                         periodic_speeds=False,
#                         periodic_actuators=False,
#                         periodic_coordinates=None):
def muscleDrivenStateTracking(
                        path_model          = None,
                        path_musclepaths    = None,
                        path_kinematics     = None,
                        path_grf            = None,
                        path_log            = None,  
                        time_initial        = None,
                        time_final          = None,
                        reserves            = 100.0,
                        mesh_interval       = 0.010,
                        guess_file          = None,   
                        path_solution       = None, 
                        is_guess            = True,
                        control_condition   = True,
                        enable_periodicity=False,
                        periodic_values=False,
                        periodic_speeds=False,
                        periodic_actuators=False,
                        periodic_coordinates=None):


    osim.Logger.removeFileSink()
    osim.Logger.addFileSink(path_log)

     

    # Construct a ModelProcessor and set it on the tool. The default
    # muscles in the model are replaced with optimization-friendly
    # DeGrooteFregly2016Muscles, and adjustments are made to the default muscle
    # parameters.
    model = osim.Model(path_model)
    state = model.initSystem()
    mass = model.getTotalMass(state)

    coord_set = model.updCoordinateSet()
    for coord_name in ['subtalar_angle', 'mtp_angle']:
        for side in ['_l', '_r']:
            coord = coord_set.get(f'{coord_name}{side}')
            coord.set_locked(False)

    # =============================================================================
    # TORQUE ACTUATORS
    # =============================================================================
    

    coordNames = ['arm_flex', 'arm_add', 'arm_rot', 'elbow_flex', 'pro_sup']
    strengths = [0.5, 0.5, 0.5, 0.5, 0.5] # [1.0, 1.0, 1.0, 1.0, 1.0]
    for coordName, strength in zip(coordNames, strengths):
        for side in ['_l', '_r']:
            actu = osim.ActivationCoordinateActuator()
            actu.set_coordinate(f'{coordName}{side}')
            actu.setName(f'torque_{coordName}{side}')
            actu.setOptimalForce(strength*mass)
            actu.setMinControl(-1.0)
            actu.setMaxControl(1.0)
            model.addForce(actu)

    # Then add in lumbar actuators
    coordNames = ['lumbar_extension', 'lumbar_bending', 'lumbar_rotation']
    for coordName in coordNames:
        actu = osim.ActivationCoordinateActuator()
        actu.set_coordinate(f'{coordName}')
        actu.setName(f'torque_{coordName}')
        actu.setOptimalForce(mass)
        actu.setMinControl(-1.0)
        actu.setMaxControl(1.0)
        model.addForce(actu)

    # =============================================================================
    # JOINT STIFFNESS MODIFICATIONS
    # =============================================================================

    lumbar_stiffness_scale = 1.0
        
    # Stiffness values for lumbar joint
    stiffnesses = [1.0, 1.5, 0.5] # in N-m/rad*kg
    for coordName, stiffness in zip(coordNames, stiffnesses):
        sgf = osim.SpringGeneralizedForce(coordName)
        sgf.setName(f'passive_stiffness_{coordName}')
        sgf.setStiffness(lumbar_stiffness_scale * stiffness * mass)
        sgf.setViscosity(2.0)
        model.addForce(sgf)

    
    # Custom ModOp to scale lumbar stiffness would go here
    # This is a placeholder - actual implementation would require custom ModOp
    # modelProcessor.append(osim.ModOpScaleLumbarStiffness(lumbar_stiffness_scale))

    # =============================================================================
    # RESERVE ACTUATORS
    # =============================================================================
    
    # reserve_strength = 50
        
    coordNames = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz', 
                'pelvis_list', 'pelvis_tilt', 'pelvis_rotation',
                'hip_adduction_r', 'hip_rotation_r', 'hip_flexion_r',
                'hip_adduction_l', 'hip_rotation_l', 'hip_flexion_l',
                'knee_angle_r', 'ankle_angle_r', 
                'knee_angle_l', 'ankle_angle_l', 
                'mtp_angle_r', 'mtp_angle_l',
                'subtalar_angle_r', 'subtalar_angle_l']
    
    for coordName in coordNames:
        actu = osim.CoordinateActuator()
        actu.set_coordinate(coordName)
        actu.setName(f'reserve_{coordName}')
        scale = 1.0
        if ((coordName == 'pelvis_tx') or
                (coordName == 'pelvis_ty') or
                (coordName == 'pelvis_tz')):
            scale = 0.1 
        actu.setOptimalForce(scale * reserves)
        actu.setMinControl(-1.0)
        actu.setMaxControl(1.0)
        model.addForce(actu)



    model.finalizeConnections()


    # =============================================================================
    # ADD IN SPRING STIFFNESS COMPONENTS
    # =============================================================================
    for side in ['_r', '_l']:
        spring = osim.SpringGeneralizedForce()
        spring.setName(f'mtp_spring{side}')
        spring.set_coordinate(f'mtp_angle{side}')
        if control_condition:
            spring.setStiffness(2.5)  # N-m/rad
            spring.setViscosity(2.0)  # Damping coefficient
            spring.setRestLength(0.0)  # Rest length
        else:
            spring.setStiffness(10.5)  # N-m/rad
            spring.setViscosity(2.0)  # Damping coefficient
            spring.setRestLength(0.0)  # Rest length

        model.addForce(spring)
        # model.addComponent(spring)
    model.finalizeConnections()

    initial_forceset = model.getForceSet()
    weld_joints = list()
    remove_actuators = list()
    for side in ['_l', '_r']:
        # weld_joints.append(f"subtalar{side}")
        weld_joints.append(f"radius_hand{side}")
        # weld_joints.append(f"wrist_flex{side}")
        # weld_joints.append(f"wrist_dev{side}")
        remove_actuators.append(f"wrist_flex{side}")
        remove_actuators.append(f"wrist_dev{side}")

    for i in range(initial_forceset.getSize() - 1, -1, -1):
        force = initial_forceset.get(i).getName()
        if force in remove_actuators:
            initial_forceset.remove(i)
            print(f"Removed: {force}")

    modelProcessor = osim.ModelProcessor(model)
    modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weld_joints))
    modelProcessor.append(osim.ModOpAddReserves(reserves))
    
    # modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))
    # Use a function-based representation for the muscle paths. This is
    # recommended to speed up convergence, but if you would like to use
    # the original GeometryPath muscle wrapping instead, simply comment out
    # this line. To learn how to create a set of function-based paths for
    # your model, see the example 'examplePolynomialPathFitter.py'.
    modelProcessor.append(osim.ModOpReplacePathsWithFunctionBasedPaths(
            path_musclepaths))
    
    

    
    # # =============================================================================
    # # MUSCLE MODELING MODIFICATIONS
    # # =============================================================================

    # Enable tendon compliance for ankle plantarflexors
    muscle_name_files = "muscle_names_and_paths.txt"
    model = modelProcessor.process()
    model.initSystem()
    muscles = model.updMuscles()
    with open(muscle_name_files, 'w') as f:
        f.write("Muscle Names and full paths\n")
        f.write("="*50 + "\n\n")
        f.write(f"Total number of muscles: {muscles.getSize()}\n\n")

        for imusc in np.arange(muscles.getSize()):
            muscle = osim.DeGrooteFregly2016Muscle.safeDownCast(muscles.get(int(imusc)))
            muscName = muscle.getName()
            muscPath = muscle.getAbsolutePathString()

            f.write(f"{imusc+1:3d}: {muscName}\n")
            f.write(f"  Path: {muscPath}\n\n")

            #**REFER TO ARNOLD ET AL. FOR SETTING CORRECT FIBER PARAMETERS 
            if ('gas' in muscName) or ('soleus' in muscName):
                muscle.set_ignore_tendon_compliance(False)
                muscle.set_tendon_strain_at_one_norm_force(0.10)
                muscle.set_passive_fiber_strain_at_one_norm_force(2.0)

    model.finalizeConnections()
    modelProcessor = osim.ModelProcessor(model)
    
    modelProcessor.append(
        osim.ModOpUseImplicitTendonComplianceDynamicsDGF())
    
    model = modelProcessor.process()
    model.initSystem()


    # Construct a TableProcessor of the coordinate data and pass it to the 
    # tracking tool. TableProcessors can be used in the same way as
    # ModelProcessors by appending TableOperators to modify the base table.
    # A TableProcessor with no operators, as we have here, simply returns the
    # base table.
    converted_states = 'converted_states.sto'
    state_stdevs = 'state_stds.csv'

    kinematicsToStates(kinematicsFileName=path_kinematics,
                       osimModelFileName=path_model,
                       outputFileName=converted_states,
                       inDegrees=True,
                       outDegrees=False,
                       excludeColumns=['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz'],
                       include_speeds=True,
                       include_accelerations=False
                       )
    
    track_coords_table = osim.TimeSeriesTable(converted_states)

    # labels_prescribed = track_coords_table.getColumnLabels()
    # # write list of coordinates that have been locked in model adjustments above that need to be removed from reference kinematics
    exclude_coords = ['wrist_dev_l', 'wrist_flex_l', 'wrist_dev_r', 'wrist_flex_r']
    # # Find columns to keep
    # keep_cols = []
    # for label in labels_prescribed:
    #     #check for excluded coordinates in current label
    #     should_exclude = False
    #     for coord in exclude_coords:
    #         if coord in label:
    #             should_exclude = True
    #             break
    #     if not should_exclude:
    #         keep_cols.append(label)
        
    # # Generate new table with kept cols only
    # filtered_table = osim.TimeSeriesTable()
    # filtered_table.setColumnLabels(keep_cols)
    # # Write data into said table
    # for i in range(track_coords_table.getNumRows()):
    #     time = track_coords_table.getIndependentColumn()[i]
    #     row = osim.RowVector(len(keep_cols), 0.0)
    #     for j, label in enumerate(keep_cols):
    #         original_index = list(labels_prescribed).index(label)
    #         row[j] = track_coords_table.getRowAtIndex(i)[original_index]
      
    #     filtered_table.appendRow(time, row)

    # # return filtered_table


    time =  track_coords_table.getIndependentColumn()

    nrows =  track_coords_table.getNumRows()
    pelvis_tx =  track_coords_table.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value')
    pelvis_tx_new = osim.Vector(nrows, 0)
    for i in np.arange(nrows):
        dx = (time[int(i)] - time[0]) * periodic_speeds
        pelvis_tx_new[int(i)] = pelvis_tx[int(i)] + dx

        pelvis_tx_upd =  track_coords_table.updDependentColumn('/jointset/ground_pelvis/pelvis_tx/value')
        for i in np.arange(nrows):
            pelvis_tx_upd[int(i)] = pelvis_tx_new[int(i)]

    for side in ['_l', '_r']:
        track_coords_table.removeColumn(f'/jointset/radius_hand{side}/wrist_dev{side}/value')
        print(f"effectively removed wrist_dev{side} from track_coords_table")
        track_coords_table.removeColumn(f'/jointset/radius_hand{side}/wrist_flex{side}/value')
        print(f"effectively removed wrist_flex{side} from track_coords_table")
        track_coords_table.removeColumn(f'/jointset/radius_hand{side}/wrist_dev{side}/speed')
        track_coords_table.removeColumn(f'/jointset/radius_hand{side}/wrist_flex{side}/speed')
    
    adjusted_coord_file = converted_states.replace('.sto', '_filt.sto')
    osim.STOFileAdapter.write(track_coords_table, adjusted_coord_file)

    ref_table_processor = osim.TableProcessor(track_coords_table)
    ref_table_processor.append(osim.TabOpLowPassFilter(6))
    ref_table_processor.append(osim.TabOpUseAbsoluteStateNames())


    scale = 1
    # weights = {
    # 'state_tracking_weight':    1 * scale,
    # 'control_weight':           1 * scale,
    # 'grf_tracking_weight':      1 * scale,
    # 'torso_orientation_weight': 1 * scale,
    # 'feet_orientation_weight':  1 * scale,
    # 'control_tracking_weight':  1 * scale, #*
    # 'aux_deriv_weight':         1 * scale,
    # 'acceleration_weight':      1 * scale,
    # }

    weights = {
    'state_tracking_weight':    50 * scale* 2,
    'control_weight':           25 * scale,
    'grf_tracking_weight':      7500 * scale,
    'torso_orientation_weight': 10 * scale,
    'feet_orientation_weight':  10 * scale,
    'control_tracking_weight':  0 * scale, 
    'aux_deriv_weight':         1000 * scale,
    'acceleration_weight':      1 * scale,
    }

    compute_sto_stdevs(stoFileName=adjusted_coord_file, outputCsvFileName=state_stdevs)

    # Create state tracking weights using dict above
    #-------------------------------------------
    # states_table = osim.TableProcessor(converted_states)
    model = modelProcessor.process()
    model.initSystem()
    model_coordinate_paths = ref_table_processor.process(model).getColumnLabels()

    ##  STILL HAVE TO GENERATE CSV STDEVS FILE IN COORDINATE CONVERSION
    coordinates_std = pd.read_csv(state_stdevs, index_col=0)
    stds_df = coordinates_std[~coordinates_std.index.isin(exclude_coords)]
    # stds_df = coordinates_std.drop(exclude_coords, errors='ignore')

    print(f"updating stdev coordinates: {stds_df}")
    state_weights = osim.MocoWeightSet()
    for valuePath in model_coordinate_paths:
        speed_path = valuePath.replace('/value', '/speed')
        value_weight = 1.0
        speed_weight = 1.0
        
        for name in stds_df.index:
            if name not in exclude_coords:
                stdev = stds_df.loc[name][0]
                denom = 1 * stdev

                if name in valuePath:
                    if 'lumbar' in name:
                        if 'torso_orientation_weight' in weights:
                            value_weight = 1.0/denom
                            speed_weight = 0.0/denom
                        else:
                            value_weight = 1/denom
                            speed_weight = 0.01/denom

                    elif('beta' in name or
                            # 'subtalar' in name or 
                            'wrist' in name):
                            value_weight = 0.0
                            speed_weight = 0.0

                    elif ('ankle' in name or 
                        #   'hip' in name or
                        'knee' in name or
                        'mtp' in name):
                        value_weight = 4.0/denom
                        speed_weight = 0.03/denom

                    elif 'pelvis' in name:
                        if 'pelvis_ty' in name:
                            value_weight = 0.0/denom
                            speed_weight = 0.01/denom
                        else:
                            value_weight = 1.0/denom
                            speed_weight = 0.01/denom

                    # elif 'subtalar' in name:
                    #     value_weight = 1.0
                    #     speed_weight = 0.0

                    else:
                        value_weight = 1.0/denom
                        speed_weight = 0.01/denom

        state_weights.cloneAndAppend(
            osim.MocoWeight(valuePath, value_weight))
        state_weights.cloneAndAppend(
            osim.MocoWeight(speed_path, speed_weight))


    numForces = 0
    for actu in model.getComponentsList():
        try:
            if (actu.getConcreteClassName().endswith('Muscle') or 
                actu.getConcreteClassName().endswith('Actuator')):
                numForces += 1
        except AttributeError: 
            continue

    # This setting allows extra data columns contained in the states
    # reference that don't correspond to model coordinates.
        # Create and name an instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName("muscle_driven_state_tracking")
    track.setModel(modelProcessor)
    track.set_allow_unused_references(True)
    track.setStatesReference(ref_table_processor)

    track.set_apply_tracked_states_to_guess(True)
    track.set_states_weight_set(state_weights)
    

    

    # Since there is only coordinate position data in the states references,
    # this setting is enabled to fill in the missing coordinate speed data using
    # the derivative of splined position data.
    track.set_track_reference_position_derivatives(True)
    if is_guess:
        track.set_states_global_tracking_weight(weights['state_tracking_weight'] * 2 / model.getNumCoordinates())
    else:
        track.set_states_global_tracking_weight(weights['state_tracking_weight'] / model.getNumCoordinates())
    track.set_control_effort_weight(weights['control_weight'] / numForces)
    
    # Initial time, final time, and mesh interval.
    coordinates_table = osim.TimeSeriesTable(adjusted_coord_file)
    time = coordinates_table.getIndependentColumn()
    
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

    # Get a reference to the MocoControlCost that is added to every MocoTrack
    # problem by default.
    problem = study.updProblem()

    # MTP toe spring stiffness for optimizaiton
    # =========================================
    # mtp_stiffness_param = osim.MocoParameter()
    # mtp_stiffness_param.setName('mtp_spring_stiffness')
    # mtp_stiffness_param.appendComponentPath("/forceset/mtp_spring_r")
    # mtp_stiffness_param.appendComponentPath("/forceset/mtp_spring_l")
    # mtp_stiffness_param.setPropertyName('stiffness')
    # mtp_stiffness_param.setBounds(osim.MocoBounds(1, 100))  # in N.m/rad
    # problem.addParameter(mtp_stiffness_param)

    # spring_stiffness_bounds = osim.MocoBounds(0, 100) # in N.m/rad
    # problem.addParameter('mtp_spring_stiffness', 
    #                      "/forceset/mtp_spring_r", 
    #                      'stiffness', 
    #                      spring_stiffness_bounds)
    
    # spring_rest_bounds = osim.MocoBounds(-0.2, 0.2) # in N.m/rad
    # problem.addParameter('mtp_spring_rest_lengths', 
    #                      "/forceset/mtp_spring_r", 
    #                      'rest_length', 
    #                      spring_rest_bounds)
    
    # spring_viscosity_bounds = osim.MocoBounds(0, 10) # in N.m.s/rad
    # problem.addParameter('mtp_spring_damping', 
    #                      "/forceset/mtp_spring_r", 
    #                      'viscosity', 
    #                      spring_viscosity_bounds)

    # spring_stiffness_bounds = osim.MocoBounds(0, 100) # in N.m/rad
    # problem.addParameter('mtp_spring_stiffness', 
    #                      "/forceset/mtp_spring_l", 
    #                      'stiffness', 
    #                      spring_stiffness_bounds)
    
    # spring_rest_bounds = osim.MocoBounds(-0.2, 0.2) # in N.m/rad
    # problem.addParameter('mtp_spring_rest_lengths', 
    #                      "/forceset/mtp_spring_l", 
    #                      'rest_length', 
    #                      spring_rest_bounds)
    
    # spring_viscosity_bounds = osim.MocoBounds(0, 10) # in N.m.s/rad
    # problem.addParameter('mtp_spring_damping', 
    #                      "/forceset/mtp_spring_l", 
    #                      'viscosity', 
    #                      spring_viscosity_bounds)    

    
    coord_bounds = {}
    reference_data = {}
    state_bounds = set_state_bounds(coord_bounds=None, reference_states_path=adjusted_coord_file, kinematicLimits=None)

    # Add in State bounds now
    # =======================
    speed_bounds = [-50.0, 50.0]
    musc_bounds = [0.01, 1.0]
    torque_bounds = [-1.0, 1.0]
    tendon_bounds = [0.01, 1.8]
    problem.setStateInfoPattern('/jointset/.*/speed', speed_bounds, [], [])
    problem.setControlInfoPattern('/forceset/.*', musc_bounds, [], [])
    problem.setControlInfoPattern('/forceset/torque.*', torque_bounds, [], [])
    problem.setControlInfoPattern('/forceset/reserve.*', torque_bounds, [], [])
    problem.setStateInfoPattern('.*/activation', musc_bounds, [], [])
    problem.setStateInfoPattern('.*torque.*/activation', torque_bounds, [], [])
    problem.setStateInfoPattern('.*reserve.*/activation', torque_bounds, [], [])
    problem.setStateInfoPattern('.*tendon_force', tendon_bounds, [], [])
    for path in state_bounds:
        traj_bounds = state_bounds[path][0]
        init_bounds = state_bounds[path][1] if len(state_bounds[path]) > 1 else []
        final_bounds = state_bounds[path][2] if len(state_bounds[path]) > 2 else []
        problem.setStateInfo(path, traj_bounds, init_bounds, final_bounds)


    # Modify control effort goal
    effort = osim.MocoControlGoal.safeDownCast(problem.updGoal("control_effort"))
    effort.setDivideByDisplacement(True)
    # effort.setWeight(weights['control_weight'] / numForces)
    # effort.setWeight(0.1)
    effort.setEnabled(True)

    average_speed = osim.MocoAverageSpeedGoal()
    average_speed.setName('avg_speed')
    average_speed.set_desired_average_speed(periodic_speeds)
    problem.addGoal(average_speed)

    # Put a large weight on the pelvis CoordinateActuators, which act as the
    # residual, forces which we would like to keep as small
    # as possible.
    model = modelProcessor.process()
    model.initSystem()
    forceSet = model.getForceSet()
    for i in range(forceSet.getSize()):
        forcePath = forceSet.get(i).getAbsolutePathString()
        if 'pelvis' in str(forcePath):
            effort.setWeightForControl(forcePath, 100)
        if not is_guess:
            if 'reserve' in str(forcePath):
                effort.setWeightForControl(forcePath, 50)


    # Add torso Orientation Goal
    exp_states = osim.TableProcessor(adjusted_coord_file).process()
    exp_state_labels = exp_states.getColumnLabels()
    torso_table = osim.TimeSeriesTable(exp_states)
    state_vars = model.getStateVariableNames()
    activation_vec = osim.Vector(torso_table.getNumRows(), 0.02)
    for isv in np.arange(state_vars.getSize()):
        state_var = state_vars.get(int(isv))
        if not state_var in exp_state_labels:
            torso_table.appendColumn(state_var, activation_vec)

    # for state_label in exp_state_labels:
    #     # exp_state_label = exp_state_labels.get(int(isv))
    #     if not state_label in state_vars:
    #         torso_table.removeColumn(state_label)

    # torso_table.removeColumn('/forceset/addbrev_r')
    torso_orientation_goal = osim.MocoOrientationTrackingGoal('torso_orientation_goal', 
                                                                weights['torso_orientation_weight'])
    torso_orientation_goal.setStatesReference(osim.TableProcessor(torso_table))
    paths = osim.StdVectorString()
    paths.push_back('/bodyset/torso')
    torso_orientation_goal.setFramePaths(paths)
    torso_orientation_goal.setEnabled(True)
    
    problem.addGoal(torso_orientation_goal)

    torso_speed_goal = osim.MocoAngularVelocityTrackingGoal('torso_angular_velocity_goal', 
                                                            0.1 * weights['torso_orientation_weight'])
    torso_speed_goal.setStatesReference(osim.TableProcessor(torso_table))
    paths = osim.StdVectorString()
    paths.push_back('/bodyset/torso')
    torso_speed_goal.setFramePaths(paths)
    torso_speed_goal.setEnabled(True)

    problem.addGoal(torso_speed_goal)

    # Foot Orientation goal
    # exp_state_table = expiremental_states.process()
    exp_state_labels = exp_states.getColumnLabels()
    foot_table = osim.TimeSeriesTable(exp_states)
    state_vars = model.getStateVariableNames()
    activation_vec = osim.Vector(foot_table.getNumRows(), 0.02)
    for isv in np.arange(state_vars.getSize()):
        state_var = state_vars.get(int(isv))
        if not state_var in exp_state_labels:
            foot_table.appendColumn(state_var, activation_vec)


    foot_orientation_goal = osim.MocoOrientationTrackingGoal('feet_orientation_goal', 
                                                                weights['feet_orientation_weight'])
    foot_orientation_goal.setStatesReference(osim.TableProcessor(foot_table))
    paths = osim.StdVectorString()
    paths.push_back('/bodyset/calcn_r')
    paths.push_back('/bodyset/calcn_l')
    foot_orientation_goal.setFramePaths(paths)
    foot_orientation_goal.setEnabled(True)

    problem.addGoal(foot_orientation_goal)

    foot_speed_goal = osim.MocoAngularVelocityTrackingGoal('feet_speed_goal', 
                                                                0.01 * weights['feet_orientation_weight'])
    foot_speed_goal.setStatesReference(osim.TableProcessor(foot_table))
    paths = osim.StdVectorString()
    paths.push_back('/bodyset/calcn_r')
    paths.push_back('/bodyset/calcn_l')
    foot_speed_goal.setFramePaths(paths)
    foot_speed_goal.setEnabled(True)

    problem.addGoal(foot_speed_goal)



    # This is supposed tohelp with preventing the model from starting off in a completely co-contracted state    
    problem.addGoal(osim.MocoInitialActivationGoal('initial_activation'))

    # Distance constraint for paths to prevent intersecting bodies
    # ==============================================

    distance_constraint = osim.MocoFrameDistanceConstraint()
    distance_constraint.setName('distance_constraint')
    foot_distance = 0.05 # min separation in METERS
    distance_constraint.addFramePair(
        osim.MocoFrameDistanceConstraintPair(
            'bodyset/calcn_l', 'bodyset/calcn_r', foot_distance, np.inf))
    distance_constraint.addFramePair(
        osim.MocoFrameDistanceConstraintPair(
            'bodyset/toes_l', 'bodyset/toes_r', foot_distance, np.inf))
    distance_constraint.addFramePair(
        osim.MocoFrameDistanceConstraintPair(
            'bodyset/calcn_l', 'bodyset/toes_r', foot_distance, np.inf))
    distance_constraint.addFramePair(
        osim.MocoFrameDistanceConstraintPair(
            'bodyset/toes_l', 'bodyset/calcn_r', foot_distance, np.inf))
    
    distance_constraint.setProjection('vector')
    distance_constraint.setProjectionVector(osim.Vec3(0, 0, 1))
    problem.addPathConstraint(distance_constraint)

    # Periodicity Goal
    # ================
    if enable_periodicity:
        print("**periodicity goal for joint coords only so far...")
            
        periodic = osim.MocoPeriodicityGoal('periodicity')
        coord_count = 0
        
        for coord in model.getComponentsList():
            if not type(coord) is osim.Coordinate: continue
            coordName = coord.getName()
            coordValue = coord.getStateVariableNames().get(0)
            coordSpeed = coord.getStateVariableNames().get(1)

            if 'beta' in coordName: continue
            # if periodic_coordinates is not None:
            #     includeCoord = False
            #     for includedCoord in periodic_coordinates:
            #         if includedCoord in coordName:
            #             includeCoord = True
            #     if not includeCoord: continue

            if 'pelvis_tx' in coordName:
                if periodic_speeds:
                    periodic.addStatePair(osim.MocoPeriodicityGoalPair(coordSpeed))
                    coord_count += 1
            else:
                if periodic_values:
                    periodic.addStatePair(osim.MocoPeriodicityGoalPair(coordValue))
                    coord_count += 1
                if periodic_speeds:
                    periodic.addStatePair(osim.MocoPeriodicityGoalPair(coordSpeed))
                    coord_count += 1

        if periodic_actuators:
            muscle_count = 0
            tendon_count = 0
            fiber_length_count = 0
            actuator_count = 0
            for actu in model.getComponentsList():
                # if not type(actu) is osim.Muscle: continue
            
                # state_var_names = actu.getStateVariableNames().get(0)
                # state_list = [state_var_names.get(i) for i in range(state_var_names.size())]

                is_actuator = actu.getConcreteClassName().endswith('Actuator')
                is_muscle = actu.getConcreteClassName().endswith('Muscle')
                # is_tendon_force = actu.getConcreteClassName().endswith('normalized_tendon_force')
                # has_activation_state = actu.getConcreteClassName().endswith('activation')
                
                if not (is_actuator or is_muscle):
                    continue
                if 'reserve' in actu.getAbsolutePathString(): continue
                if 'torque' in actu.getAbsolutePathString(): continue

                # if(not (actu.getConcreteClassName().endswith('Actuator') or 
                #         actu.getConcreteClassName().endswith('Muscle'+'/activation'))): continue

                if is_muscle:
                    try:
                        activ_state = actu.getStateVariableNames().get(0)    
                        # act_path = actu.getAbsolutePathString() + '/activation'
                        periodic.addStatePair(osim.MocoPeriodicityGoalPair(activ_state))
                        muscle_count += 1

                        # Check for fiber length state variable
                        for i in range(actu.getStateVariableNames().getSize()):
                            state_name = actu.getStateVariableNames().get(i)
                            if 'fiber_length' in state_name.lower():
                                periodic.addStatePair(osim.MocoPeriodicityGoalPair(state_name))
                                fiber_length_count += 1
                                print(f"Added fiber length for {actu.getName()}: {state_name}")
                                break

                        if actu.getStateVariableNames().getSize() > 1:
                            tendon_state = actu.getStateVariableNames().get(1)
                            # tendon_path = actu.getAbsolutePathString() + '/normalized_tendon_force'
                            periodic.addStatePair(osim.MocoPeriodicityGoalPair(tendon_state))
                            tendon_count += 1

                    except Exception as e:
                        print(f'error adding muscle {actu.getName()}: {e}')
                        continue
                    
                elif is_actuator:
                    try:
                        if actu.getStateVariableNames().getSize() > 0:
                            activ_state = actu.getStateVariableNames().get(0)
                            periodic.statePair(osim.MocoPeriodicityGoalPair(activ_state))
                            actuator_count += 1

                        control_path = actu.getAbsolutePathString()
                        periodic.addControlPair(osim.MocoPeriodicityGoalPair(control_path))
                        actuator_count += 1
                    except Exception as e:
                        print(f'Error adding actuator control {actu.getName()}: {e}')
                        continue

                # fiber_len_path = actu.getAbsolutePathString() + '/fiber_length'
                # periodic.addStatePair(osim.MocoPeriodicityGoalPair(fiber_len_path))
                # periodic.addStatePair(osim.MocoPeriodicityGoalPair(actu.getStateVariableNames().get(1)))


                # if is_actuator:
                #     periodic.addStatePair(osim.MocoPeriodicityGoalPair(actu.getStateVariableNames().get(0)))
                #     periodic.addStatePair(osim.MocoPeriodicityGoalPair(actu.getAbsolutePathString()))

                # if (actu.getConcreteClassName().endswith('Muscle') and 
                #     actu.getStateVariableNames().getSize() > 1):
                #     periodic.addStatePair(osim.MocoPeriodicityGoalPair(actu.getStateVariableNames().get(1)))

        print(f"coordinates added to periodicity constraint: {coord_count}")
        print(f"muscles added to periodicity constraint: {muscle_count}")
        print(f"tendons added to periodicity constraint: {tendon_count}")
        print(f"fiber lengths added to periodicity constraint: {fiber_length_count}")
        print(f"other actuators added to periodicity constraint: {actuator_count}")
        problem.addGoal(periodic)

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
    
    # Add in Contact Tracking
    # =======================
    contact_tracking = osim.MocoContactTrackingGoal('contact', weights['grf_tracking_weight'])
    contact_tracking.setExternalLoadsFile(path_grf)
    right_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'ExternalForce_2', ['/bodyset/calcn_r', '/bodyset/toes_r'])
    left_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'ExternalForce_1', ['/bodyset/calcn_l', '/bodyset/toes_l'])
    contact_tracking.addContactGroup(right_contact_group)
    contact_tracking.addContactGroup(left_contact_group)
    contact_tracking.setNormalizeTrackingError(True)
    contact_tracking.setEnabled(True)
        
    # Add the goal to the problem
    problem.addGoal(contact_tracking)

    
    print("Contact tracking goal added successfully")
    print(f"  - Right foot contact spheres: {len(forceNamesRightFoot)}")
    print(f"  - Left foot contact spheres: {len(forceNamesLeftFoot)}")

    # Update the solver problem and tolerances.
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    solver.set_optim_convergence_tolerance(1e-2)
    solver.set_optim_constraint_tolerance(1e-2)
    solver.set_optim_max_iterations(5000)
    # solver.set_num_mesh_intervals(int(np.round((time[-1] - time[0]) / mesh_interval)))
    solver.set_multibody_dynamics_mode('explicit')
    solver.set_parallel(28)
    # solver.set_variables_using_bounds(True)
    solver.set_scale_variables_using_bounds(True)
    solver.set_optim_finite_difference_scheme('forward')

    solver.set_minimize_implicit_auxiliary_derivatives(True)
    solver.set_implicit_auxiliary_derivatives_weight(weights['aux_deriv_weight'] / 6.0)
    solver.set_implicit_auxiliary_derivative_bounds(osim.MocoBounds(-100, 100))

    # solver.set_parameters_require_initsystem(False)
    # solver.set_verbosity(1)
    solver.set_optim_solver("ipopt")
    # solver.set_optim_ipopt_print_level(3)
    # solver.optim_print_frequency(1)
    # solver.set_write_solution("yes")
    
    
    
    input_guess = osim.MocoTrajectory(guess_file)

    states_table = input_guess.exportToStatesTable()
    control_table = input_guess.exportToControlsTable()
    state_labels = states_table.getColumnLabels()
    control_labels = control_table.getColumnLabels()

    for label in state_labels:
        if 'reserve' in label:
            states_table.removeColumn(label)

    for label in control_labels:
        if 'reserve' in label:
            control_table.removeColumn(label)


    guess = solver.createGuess()
    guess.insertControlsTrajectory(control_table, True)    
    guess.insertStatesTrajectory(states_table, True)
    # guess.randomizeAdd()
    solver.setGuess(guess)

    # Print the problem setup
    solution_dir = os.path.dirname(path_solution)
    solution_basename = os.path.basename(path_solution)

    problem_name = solution_basename.replace('.sto', '_spring.omoco')
    
    print(f"Problem directory: {problem_name}")

    path_problem = os.path.join(solution_dir, problem_name)

    # path_problem = path_solution.replace('.sto', '_spring.omoco')


    study.printToXML(path_problem)


    # Solve!
    solution = study.solve()
    solution.write(path_solution)

    osim.Logger.setLevel(osim.Logger.Level_Info)
    log_file = '_example_muscle_driven_output.log'
    osim.Logger.addFileSink(log_file)



    # Also grab the model output that was used too
    end_model_path = path_model.replace('.osim', '_moco_adjusted.osim')
    model.printToXML(end_model_path)

    osim.Logger.removeFileSink()
    osim.Logger.setLevel(osim.Logger.Level_Off)

    # Retrieve simulated GRFs generated by each specific contact sphere
    solution = osim.MocoTrajectory(path_solution)

    force_paths_ankle_muscles = [
    '/forceset/gasmed_r',
    '/forceset/gaslat_r', 
    '/forceset/soleus_r',
    '/forceset/tibant_r',
    '/forceset/gasmed_l',
    '/forceset/gaslat_l',
    '/forceset/soleus_l',
    '/forceset/tibant_l'
    ]
    
    # Get net moments from solution
    solution_trajectory = osim.MocoTrajectory(path_solution)
    force_paths = osim.StdVectorString()
    force_paths.push_back(path_grf)
    # net_moments = calc_net_joint_moments(model, solution_trajectory)
    net_moments = study.calcGeneralizedForces(solution_trajectory, force_paths_ankle_muscles)
    net_moments_file = path_solution.replace('.sto', '_net_moments.sto')
    osim.STOFileAdapter.write(net_moments, net_moments_file)

    path_comp = path_solution.replace('.sto', '_joint_traj_errors.sto')
    get_results.analyze_moco_ik_tracking_errors(moco_solution_file=path_solution, ik_kinematics_file=path_kinematics, output_file=path_comp, 
                                    model_file=path_model, omoco_setup_file=path_problem)

    external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)
    osim.STOFileAdapter.write(external_loads, path_solution.replace('.sto', '_simulated_loads.sto'))


    path_spring_solution = path_solution.replace('.sto', '_spring_forces_mtp.sto')
    get_results.write_mtp_spring_forces(solution=solution, model=model, output_file=path_spring_solution)

    neg_force_check = path_solution.replace('.sto', '_neg_muscles.csv')
    validate.get_negative_muscle_forces(solution=solution, model=model, output_path=neg_force_check)

    # Store the muscle data from the solution
    muscle_data = validate.calc_muscle_mechanics(model, solution)
    muscle_data_file = path_solution.replace('.sto', '_muscle_mechanics.csv')
    osim.STOFileAdapter.write(muscle_data, muscle_data_file)


    main_coords = ['/jointset/hip_r/hip_flexion_r',
                '/jointset/walker_knee_r/knee_angle_r',
                '/jointset/ankle_r/ankle_angle_r',
                '/jointset/mtp_r/mtp_angle_r']
    moments_file = path_solution.replace('.sto', '_leg_net_moments.png')
    get_results.plot_joint_moment_constituents(model, solution, coord_paths=main_coords, 
                                               muscle_paths=None, output_file=moments_file)

    contacts_file = path_solution.replace('.sto', '_contact_forces.sto')
    cop_file = path_solution.replace('.sto', '_cop_output.sto')

    contact_forces, copTable = get_results.get_solution_contacts(model, solution)
    osim.STOFileAdapter.write(contact_forces, contacts_file)
    osim.STOFileAdapter.write(copTable, cop_file)



    # Visualize the solution.
    # study.visualize(solution)


if __name__ == "__main__":

    # path_opensim_core = "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/"
    path_opensim_core = "/home/lisca/biomec/data/20250306_dataset/"

    name_problem = 'track_muscle_driven_running'

    path_model = os.path.join(
        path_opensim_core,
        "sub04/ses20250524/sub04_SCALED_addbio_free_subt_mtp.osim")
        # "sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim")
        # '2015_match_markers_and_physics.osim')
    path_grf = os.path.join(
        path_opensim_core,
        'results_pool/averages/sub04_strirt_cond00000_speed03000_grf_move_averages.xml')
    path_markers = os.path.join(
        path_opensim_core,
        'sub03/ses20250306/sub04_strirt_cond00000_speed03000_0001_0003_trc.trc')
    path_kinematics = os.path.join(
        path_opensim_core,
        'results_pool/averages/sub04_strirt_cond00000_speed03000_ik_move_averages.sto')
    prescr_states = os.path.join(
        path_opensim_core,
        'results_pool/averages/sub04_strirt_cond00000_speed03000_rad_averages.sto')
    path_musclepaths = os.path.join(
        path_opensim_core,
        'sub04/ses20250524/sub04_FunctionBasedPathSet.xml')
        # 'sub03/ses20250306/sub03_FunctionBasedPathSet.xml')


    guess_file          = os.path.join(path_opensim_core, 'results_pool/averages/guesses/sub04_strirt_cond00000_speed03000_mi_emg_solution.sto') 
    # guess_file          = os.path.join(path_opensim_core, 'results_pool/averages/sub03_strifc_cond00000_speed02000_2res_02mesh.sto') 
    # guess_file          = os.path.join(path_opensim_core, 'results_pool/averages/guesses/sub04_strift_cond00000_speed03000_50res_1mesh.sto') 
    # guess_file          = 'example_muscle_driven_test_50b_new_model_skip.sto'
    # path_solution       = os.path.join(path_opensim_core, 'results_pool/averages/guesses/sub04_strirt_cond00000_speed03000_20res_02mesh.sto')
    # path_log            = os.path.join(path_opensim_core, 'results_pool/averages/guesses/sub04_strirt_cond00000_speed03000_20res_02mesh.log')

    path_solution       = 'sub04_strirt_cond00000_speed03000_20res_02mesh.sto'
    path_log            = 'sub04_strirt_cond00000_speed03000_20res_02mesh.log'


    muscleDrivenStateTracking(
                            path_model          = path_model,
                            path_musclepaths    = path_musclepaths,
                            path_kinematics     = path_kinematics,
                            path_grf            = path_grf,
                            path_log            = path_log,
                            time_initial        = None,
                            time_final          = None,
                            reserves            = 200.0,
                            mesh_interval       = 0.2,
                            is_guess            = False,
                            path_solution       = path_solution, 
                            guess_file          = guess_file,
                            control_condition   = False,
                            enable_periodicity=False,
                            periodic_values=False,
                            periodic_speeds=3.0,
                            periodic_actuators=False,
                            periodic_coordinates=None)