
import numpy as np
import opensim as osim
import pandas as pd
import os
import datetime
from scipy.interpolate import interp1d
from PIL import Image
# import utilities
# from utilities import simtk2numpy, numpy2simtk
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from osimFunctions import kinematicsToStates, kinematicsToStatesForMoco
import re
import warnings
from typing import Optional, List, Tuple, Dict, Callable, Any
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
from functools import wraps



def calc_negative_muscle_force_base(model, solution):
    """
    ** model has to be a modelProcessor().process() object, solution must be an osim.MocoTrajectory() object
    """
    model. initSystem()
    outputs = osim.analyze(model, solution.exportToStatesTable(),
                           solution.exportToControlsTable(), ['.*\|tendon_force'])
    def simtkmin(simtkvec):
        lowest = np.inf
        for i in range(simtkvec.size()):
            if simtkvec[i] < lowest:
                lowest = simtkvec[i]
            return lowest
        
        negforces = list()
        muscle_names = list()
        for imusc in range(model.getMuscles().getSize()):
            musc = model.updMuscles().get(imusc)
            max_iso = musc.get_max_isometric_force()
            force = musc.getDependentColumn(
                musc.getAbsolutePathString() + "|tendon_force")
            neg = simtkmin(force) / max_iso

            if neg < 0:
                negforces.append(neg)
                muscle_names.append(musc.getName())
                print(f'    {musc.getName()}: {neg} F_iso')
        if len(negforces) == 0:
            print('No negative forces')
        else:
            imin = np.argmin(negforces)
            print(f'Largest negative force: {muscle_names[imin]}  '
                  f'with {negforces[imin]} F_iso')
        return min([0] + negforces)
    
def calc_muscle_mechanics(model, solution):
    output_list = list()
    for output in ['normalized_fiber_length', 'normalized_fiber_velocity', 
                       'active_fiber_force', 'passive_fiber_force', 
                       'tendon_force', 'activation', 'cos_pennation_angle', 
                       'active_force_length_multiplier', 
                       'force_velocity_multiplier', 'passive_force_multiplier',
                       'tendon_length', 'tendon_strain']:
        for imusc in range(model.getMuscles().getSize()):
            musc = model.updMuscles().get(imusc)
            output_list.append(f'.*{musc.getName()}.*\|{output}')

    for output in ['length', 'lengthening_speed']:
        output_list.append(f'.*\|{output}')

    outputs = osim.analyze(model, solution.exportToStatesTable(), 
        solution.exportToControlsTable(), output_list)
    
    return outputs

def calc_negative_muscle_forces(model, solution):
    model.initSystem()
    print(f'Negative force report for {model.getName()}')
    
    return calc_negative_muscle_force_base(model, solution)
    
def create_contact_sphere_force_table(model, solution):
    """
    Create two tables of contact sphere forces and CoP location from the solution
    """
    model.initSystem()

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

    # external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)
    # osim.STOFileAdapter.write(external_loads, os.path.join('example_muscle_driven' + '_simulated_external_loads.sto'))





    force_names = ['forceset/contactHeel',
                      'forceset/contactLateralRearfoot',
                      'forceset/contactLateralMidfoot',
                      'forceset/contactLateralToe',
                      'forceset/contactMedialToe',
                      'forceset/contactMedialMidfoot'
                      ]
    
    force_labels = ['heel',
                       'lat_rear',
                       'lat_mid', 
                       'lat_toe', 
                       'med_toe', 
                       'med_mid',
                       ]
    
    sphere_names = ['contactgeometryset/heel',
                       'contactgeometryset/lateralRearfoot',
                       'contactgeometryset/lateralMidfoot',
                       'contactgeometryset/lateralToe',
                       'contactgeometryset/medialToe',
                       'contactgeometryset/medialMidfoot',
                       ]
    
    labels = osim.StdVectorString()
    for iside, side in enumerate(['_r', '_l']):
        # zipped = zip(force_labels, sphere_names)
        for force_label in force_labels:
            for suffix in ['_force_v', '_force_p', '_torque_']:
                labels.append(f'{force_label}{side}{suffix}')

    external_force_table = osim.TimeSeriesTableVec3()
    cop_table = osim.TimeSeriesTableVec3()        
    states_trajectory = solution.exportToStatesTrajectory(model)
    # time_vector = states_trajectory.getIndependentColumn()
    # num_states = len(time_vector)
    num_states = states_trajectory.getSize()
    

    # previous_time = None
    # processed_count = 0
    time_points = []
    for istate in range(num_states):
        # Access each state object
        state = states_trajectory.get(istate)
        # Get the time for the current state and append to the list
        time_points.append(state.getTime())
    print(time_points)

    for istate in range(num_states):
        state = states_trajectory.get(istate)
        # model.realizeVelocity(state)
        model.realizeAcceleration(state)
        # current_time = float(time_vector[istate])
        current_time = state.getTime()


        # if previous_time is not None and abs(current_time - previous_time) < 1e-10:
        #     # Skip this state if the time is the same as the previous state
        #     continue

        # create state from the tabel data
        # state = model.initSystem()

        # state_names = states_trajectory.getColumnLabels()
        # for j in range(len(state_names)):
        #     state_name = state_names[j]
        #     state_value = states_table.getDependentColumn(state_name)[istate]
        #     try:
        #         model.setStateVariableValue(state, state_name, state_value)
        #     except:
        #         # skip if state variable is not found in current model
        #         pass

        # model.realizeDynamics(state)


        
        row = osim.RowVectorVec3(3*2*len(force_names)) # 6 contact points * 3 values (force, position, torque)
        cop_row = osim.RowVectorVec3(2)
        
        # zipped = zip(force_names, force_labels, sphere_names)
        for iside, side in enumerate(['_r', '_l']):

            cop = osim.Vec3(0)
            cop_torque_sum_x = 0
            cop_torque_sum_z = 0
            cop_force_sum = 0

            offset = 3 * iside * len(sphere_names)
            # reset zipped iterator per side
            # zipped = zip(force_names, force_labels, sphere_names)
            
            # print(f"Zipped list confirmed to contain {len(list(zipped))} tuples.")
            
            for i, (force_name, force_label, sphere_name) in enumerate(
                zip(force_names, force_labels, sphere_names)):

                force = osim.Vec3(0)
                torque = osim.Vec3(0)

                try:
                    force_obj = osim.Force.safeDownCast(model.getComponent(f'{force_name}{side}'))
                    force_vals = force_obj.getRecordValues(state)

                    force[0] = force_vals.get(0)
                    force[1] = force_vals.get(1)
                    force[2] = force_vals.get(2)
                    torque[0] = force_vals.get(3)
                    torque[1] = force_vals.get(4)
                    torque[2] = force_vals.get(5)
                
                    # Print the retrieved values for debugging
                    print(f"Force: [{force[0]:.4f}, {force[1]:.4f}, {force[2]:.4f}]")
                    print(f"Torque: [{torque[0]:.4f}, {torque[1]:.4f}, {torque[2]:.4f}]")

                except Exception as e:
                    print(f"Error getting force for {force_name}{side}: {e}")
                    continue

                # Get sphere position
                try:

                    sphere = osim.ContactSphere.safeDownCast(model.getComponent(f'{sphere_name}{side}'))
                    frame = sphere.getFrame()
                    position = frame.getPositionInGround(state)
                    location = sphere.get_location()
                    location_in_ground = frame.expressVectorInGround(state, location)

                    position[0] = position[0] + location_in_ground[0]
                    position[1] = position[1] + location_in_ground[1]
                    position[2] = position[2] + location_in_ground[2]
                except Exception as e:
                    print(f"Error getting position for {sphere_name}{side}: {e}")
                    position = osim.Vec3(0, 0, 0)
                
                # Store in row
                row[3 * i + offset] = force
                row[3 * i + 1 + offset] = position
                row[3 * i + 2 + offset] = torque

                force_y = abs(force[1])
                if force_y > 0.1:
                    cop_force_sum += force_y
                    cop_torque_sum_x += force_y * position[0]
                    cop_torque_sum_z += force_y * position[2]

                    print(f"Time: {state.getTime():.4f}, Side: {side}, Contact Point: {force_label}{side}")
                    print(f"Force at point: {force}")
                    print(f"Position at point: {position}")
                    print(f"Current summed force: {cop_force_sum:.4f}")
                    print(f"Current summed torque (X): {cop_torque_sum_x:.4f}")
                    print(f"Current summed torque (Z): {cop_torque_sum_z:.4f}")
                
                # for suffix in ['_force_v', '_force_p', '_torque_']:
                #     labels.append(f'{force_label}{side}{suffix}')

            if np.abs(cop_force_sum) > 0.1:
                cop[0] = cop_torque_sum_x / cop_force_sum
                cop[2] = cop_torque_sum_z / cop_force_sum
                print(f"  CoP{side}: [{cop[0]:.3f}, {cop[1]:.3f}, {cop[2]:.3f}], Total Force: {cop_force_sum:.2f}N")
            else:
                print(f"  CoP{side}: No contact (force sum = {cop_force_sum:.2f}N)")

            cop_row[iside] = cop

        external_force_table.appendRow(current_time, row)
        cop_table.appendRow(current_time, cop_row)

    external_force_table.setColumnLabels(labels)

    labels = osim.StdVectorString()
    labels.append('2_ground_force_p')
    labels.append('1_ground_force_p')
    cop_table.setColumnLabels(labels)

    suffixes = osim.StdVectorString()
    suffixes.append('x')
    suffixes.append('y')
    suffixes.append('z')

    return external_force_table.flatten(suffixes), cop_table.flatten(suffixes)

def get_solution_contacts(model, solution):
    print(f"contact sphere forces for {solution}")
    model.initSystem()
    return create_contact_sphere_force_table(model, solution)

# # Retrieve simulated GRFs generated by each specific contact sphere
# model = osim.Model('/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED.osim')
# solution = osim.MocoTrajectory('example_muscle_driven_21.sto')
# contact_forces, copTable = create_contact_sphere_force_table(model, solution)
# osim.STOFileAdapter.write(contact_forces, os.path.join('example_muscle_driven_21.sto' + 'contact_sphere_forces.sto'))
# osim.STOFileAdapter.write(contact_forces, os.path.join('example_muscle_driven_21.sto' + 'cop_table.sto'))

def calc_net_joint_moments(model, solution):
    """
    Calculate net joint moments from a Moco solution using calcGeneralizedForces().
    
    Parameters:
    -----------
    model : osim.Model
        OpenSim model used in the simulation
    solution : osim.MocoTrajectory  
        Moco solution trajectory

    coord_paths : list of str denoting paths of each coordinate to include
 
    Returns:
    --------
    osim.TimeSeriesTable
        Table containing net joint moments for all coordinates
    """
    
    # Initialize the model
    model.initSystem()
    

    import tempfile
    with tempfile.TemporaryDirectory() as temp_dirname:
        id_tool = osim.InverseDynamicsTool()
        modelID = osim.Model(model)
        id_tool.setModel(modelID)
        
        table = solution.exportToStatesTable()
        labels = list(table.getColumnLabels())
        import re
        for label in range(len(labels)):
            labels[label] = labels[label].replace('/value', '')
            labels[label] = re.sub('/jointset/(.*?)/', '', labels[label])
        table.setColumnLabels(labels)

        storage = osim.Storage()
        sto_labels = osim.ArrayStr()
        sto_labels.append('time')
        for label in labels:
            sto_labels.append(label)
        storage.setColumnLabels(sto_labels)

        times = table.getIndependentColumn()
        for itime in np.arange(table.getNumRows()):
            row_vec = osim.RowVector(table.getRowAtIndex(int(itime)))
            storage.append(times[itime], row_vec.transpose())

        id_tool.setCoordinateValues(storage)
        excluded_forces = osim.ArrayStr()
        excluded_forces.append('ACTUATORS')
        id_tool.setExcludedForces(excluded_forces)
        id_output = 'joint_moment_residuals.sto'
        id_tool.setResultsDir(temp_dirname)
        id_tool.setOutputGenForceFileName(id_output)
        id_tool.run()

        net_moments_time_series = osim.TimeSeriesTable(os.path.join(temp_dirname, id_output))
    return net_moments_time_series
    
def plot_joint_moment_constituents(model, solution, coord_paths=None,
                                 file_name=None, muscle_paths=None, coordact_paths=[], 
                                 output_file=None):
    
    """
    solution is a moco trajectory processed solution iwth osim.MocoTrajectory()
    coord_paths : list of str denoting paths of each coordinate to include
    muscle_paths : optional list of str for getting, plotting tendon forces for full paths;
    coordact_paths : Additional list of optional str that would contribute to net moment that isn't already in muscle set
    
    Outputs: 
    A matplotlib figure object
    """

    
    if not muscle_paths:
        muscle_paths = list()
        for muscle in model.getMuscleList():
            muscle_paths.append(muscle.getAbsolutePathString())
    num_coord_acts = len(coordact_paths)
    num_coords = len(coord_paths)
    num_muscles = len(muscle_paths)

    net_moments = calc_net_joint_moments(model, solution)
    time = solution.getTimeMat()
    states_trajectory = solution.exportToStatesTrajectory(model)
    
    fig = plt.figure(figsize=(8.5, 11))
    tendon_forces = np.empty((len(time), num_muscles))
    for imusc, muscle_path in enumerate(muscle_paths):
        muscle = model.getComponent(muscle_path)
        for itime in range(len(time)):
            state = states_trajectory.get(itime)
            model.realizeDynamics(state)
            tendon_forces[itime, imusc] = muscle.getTendonForce(state)

    coordact_moments = np.empty((len(time), num_coord_acts))
    for iact, coord_act_path in enumerate(coordact_paths):
        coord_act = model.getComponent(coord_act_path)
        for itime in range(len(time)):
            state = states_trajectory.get(itime)
            model.realizeDynamics(state)
            coordact_moments[itime, iact] = coord_act.getActuation(state)

    for icoord, coord_path in enumerate(coord_paths):
        coord = model.getComponent(coord_path)
        label = os.path.split(coord_path)[-1] + '_moment'

        #* Get moments as defined by prior dictionary ealrier in function
        net_moment = simtk2numpy(net_moments.getDependentColumn(label))

        muscle_moment_arms = np.empty((len(time), num_muscles))
        for imusc, muscle_path in enumerate(muscle_paths):
            muscle = model.getComponent(muscle_path)
            for itime in range(len(time)):
                state = states_trajectory.get(itime)
                muscle_moment_arms[itime, imusc] = muscle.computeMomentArm(state, coord)

        axis = fig.add_subplot(num_coords, 1, icoord + 1)
        net_integrated_moment = np.trapz(np.abs(net_moment), x=time)
        sum_displayed_actuators = np.zeros_like(time)
        for imusc, muscle_path in enumerate(muscle_paths):
            if np.any(muscle_moment_arms[:, imusc]) > 0.00001:
                current_moment = tendon_forces[:, imusc] * muscle_moment_arms[:, imusc]
                moment_integral = np.trapz(np.abs(current_moment), x=time)
                if moment_integral > 0.01 * net_integrated_moment:
                    axis.plot(time, current_moment, label=muscle_path)
                    sum_displayed_actuators += current_moment

        for ext_act, coordact_path in enumerate(coordact_paths):
                current_moment = coordact_moments[:, ext_act]
                axis.plot(time, current_moment, label=coordact_path)
                sum_displayed_actuators += current_moment

        axis.plot(time, sum_displayed_actuators, label='sum actuators shown', color='gray', linewidth=2)
        
        axis.plot(time, net_moment, label='sum actuators shown', color='black', linewidth=2)

        axis.set_title(coord_path)
        axis.set_ylabel('Moment (Nm)')
        axis.legend(frameon=False, bbox_to_anchor=(1, 1), loc='upper left', ncol=2)
        axis.tick_params(axis='both')
    axis.set_xlabel('Time (s)')

    fig.tight_layout()
    if output_file:
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if file_name is not None:
            filename = f'joint_moment_constituents_{file_name}.png'
        else:
            # os.makedirs(output_file, exist_ok=True)
            filename = f'joint_moment_constituents_{solution}.png'
        filepath = os.path.join(output_file, filename)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f'Saved figure to {filepath}')
    return fig
        




def analyze_moco_ik_tracking_errors(moco_solution_file=None, ik_kinematics_file=None, output_file=None, 
                                   model_file=None, omoco_setup_file=None, include_weighted_quotients=True):
    """
    Compute tracking errors between Moco solution and IK kinematics, calculating RMS errors
    and identifying maximum error coordinates at each time point.
    
    Both input files are converted to degrees using kinematicsToStates() before comparison.
    If omoco_setup_file is provided, weighted error quotients are also computed and included.
    
    Parameters:
    -----------
    moco_solution_file : str
        Path to the Moco solution .sto file (typically in radians)
    ik_kinematics_file : str
        Path to the IK kinematics .sto file (typically in degrees)
    output_file : str
        Path for the output .sto file containing error analysis
    model_file : str, optional
        Path to the OpenSim model file (.osim) - required for unit conversion
        If None, will attempt to infer from common naming patterns
    omoco_setup_file : str, optional
        Path to the .omoco setup file for extracting coordinate weights
        If provided, weighted error quotients (error/weight) will be included
        
    Returns:
    --------
    dict
        Dictionary containing error statistics and file paths
    """

    def load_sto_file(filename):
        """Load STO file and return pandas DataFrame"""
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the line with 'endheader'
        endheader_idx = None
        for i, line in enumerate(lines):
            if 'endheader' in line.lower():
                endheader_idx = i
                break
        
        if endheader_idx is None:
            raise ValueError(f"No 'endheader' found in file: {filename}")
        
        # Headers are in the line after 'endheader'
        header_line = lines[endheader_idx + 1]
        headers = [h.strip() for h in header_line.split('\t') if h.strip()]
        
        # Data starts two lines after 'endheader'
        data_start_idx = endheader_idx + 2
        
        # Parse the data
        data_dict = {header: [] for header in headers}
        
        for i in range(data_start_idx, len(lines)):
            line = lines[i].strip()
            if not line:  # Skip empty lines
                continue
            
            values = line.split('\t')
            if len(values) >= len(headers):
                for j, header in enumerate(headers):
                    try:
                        data_dict[header].append(float(values[j]))
                    except (ValueError, IndexError):
                        # Handle non-numeric data or missing values
                        data_dict[header].append(np.nan)
        
        return pd.DataFrame(data_dict)
    
    def extract_coordinate_name(column_name):
        """
        Extract coordinate name from different column naming conventions.
        Handles both direct coordinate names and path-style names.
        """
        # Remove common prefixes and suffixes
        name = column_name.replace('/jointset/', '').replace('/value', '').replace('/speed', '')
        
        # Handle path-style naming (e.g., '/jointset/ground_pelvis/pelvis_tx/value')
        if '/' in name:
            parts = name.split('/')
            # Usually the coordinate name is the last part or second to last
            if len(parts) >= 2:
                return parts[-1]  # Take the last part as coordinate name
        
        return name
    
    def match_coordinates(moco_cols, ik_cols):
        """
        Bi-directional coordinate matching between Moco solution and IK data.
        Tries both directions:
        1. Extract coordinate names from IK paths and search in Moco paths
        2. Extract coordinate names from Moco paths and search in IK paths
        
        Excludes columns containing 'speed', 'beta', or 'acceleration' in their name.
        Returns dictionary mapping Moco columns to IK columns.
        """
        coord_mapping = {}
            
        excluded_patterns = ['speed', 'beta', 'acceleration']
        
        # Filter out non-coordinate columns for both datasets
        moco_coords = []
        ik_coords = []
        
        for col in moco_cols:
            if col.lower() != 'time' and not any(pattern in col.lower() for pattern in excluded_patterns):
                moco_coords.append(col)
        
        for col in ik_cols:
            if col.lower() != 'time' and not any(pattern in col.lower() for pattern in excluded_patterns):
                ik_coords.append(col)
        
        print(f"Debug: Found {len(moco_coords)} Moco coordinates and {len(ik_coords)} IK coordinates")
        print(f"Debug: Sample Moco coordinates: {moco_coords[:3] if moco_coords else 'None'}")
        print(f"Debug: Sample IK coordinates: {ik_coords[:3] if ik_coords else 'None'}")
        
        # Common coordinate name mappings for known variations
        name_mappings = {
            'pro_sup_l': ['pronation_supination_l', 'radioulnar_l', 'pro_sup_l'],
            'pro_sup_r': ['pronation_supination_r', 'radioulnar_r', 'pro_sup_r'],
            'lumbar_extension': ['back_extension', 'torso_extension', 'spine_extension'],
            'lumbar_bending': ['back_bending', 'torso_bending', 'spine_bending'],
            'lumbar_rotation': ['back_rotation', 'torso_rotation', 'spine_rotation'],
            'subtalar_angle_l': ['subtalar_l', 'subtalar_angle_l'],
            'subtalar_angle_r': ['subtalar_r', 'subtalar_angle_r'],
            'mtp_angle_l': ['mtp_l', 'mtp_angle_l'],
            'mtp_angle_r': ['mtp_r', 'mtp_angle_r'],
            # Add more mappings as needed
        }
        
        def try_match_coordinate(coord_name, target_paths, source_path, direction):
            """Helper function to try matching a coordinate name against target paths"""
            best_match = None
            
            # Method 1: Direct coordinate name matching
            for target_path in target_paths:
                target_coord_name = extract_coordinate_name(target_path)
                
                if coord_name.lower() == target_coord_name.lower():
                    best_match = target_path
                    print(f"  {direction} Direct match: {coord_name} -> {target_path}")
                    return best_match
            
            # Method 2: Search for coordinate name within target paths
            for target_path in target_paths:
                if coord_name.lower() in target_path.lower():
                    # Additional check: make sure it's a coordinate name, not just partial match
                    path_parts = target_path.lower().split('/')
                    if coord_name.lower() in path_parts:
                        best_match = target_path
                        print(f"  {direction} Path component match: {coord_name} -> {target_path}")
                        return best_match
            
            # Method 3: Fuzzy matching
            for target_path in target_paths:
                path_parts = target_path.lower().split('/')
                
                # Try matching with common variations
                coord_variations = [
                    coord_name.lower(),
                    coord_name.lower().replace('_', ''),
                    coord_name.lower().replace('-', '_'),
                ]
                
                for variation in coord_variations:
                    for part in path_parts:
                        # Check for exact match of variation in path part
                        if variation == part:
                            best_match = target_path
                            print(f"  {direction} Fuzzy exact match: {coord_name} ({variation}) -> {target_path}")
                            return best_match
                        # Check for partial match within path part
                        elif (variation in part and 
                            len(variation) > 3 and  # Avoid very short matches
                            abs(len(variation) - len(part)) <= 2):  # Similar length
                            best_match = target_path
                            print(f"  {direction} Fuzzy partial match: {coord_name} ({variation} in {part}) -> {target_path}")
                            return best_match
            
            # Method 4: Common coordinate name mappings
            potential_names = name_mappings.get(coord_name.lower(), [coord_name.lower()])
            
            for potential_name in potential_names:
                for target_path in target_paths:
                    if potential_name in target_path.lower():
                        path_parts = target_path.lower().split('/')
                        if any(potential_name in part for part in path_parts):
                            best_match = target_path
                            print(f"  {direction} Mapped match: {coord_name} -> {potential_name} -> {target_path}")
                            return best_match
            
            return None
        
        # DIRECTION 1: Extract coordinate names from IK paths and search in Moco paths
        print(f"\n=== Direction 1: Matching IK coordinates to Moco paths ===")
        for ik_path in ik_coords:
            # Skip if already matched
            if ik_path in coord_mapping.values():
                continue
                
            ik_coord_name = extract_coordinate_name(ik_path)
            print(f"Processing IK path '{ik_path}' -> coordinate name '{ik_coord_name}'")
            
            # Find available Moco paths (not already matched)
            available_moco = [m for m in moco_coords if m not in coord_mapping.keys()]
            
            best_match = try_match_coordinate(ik_coord_name, available_moco, ik_path, "IK->Moco")
            
            if best_match:
                coord_mapping[best_match] = ik_path
            else:
                print(f"  No match found for IK coordinate: {ik_coord_name} (from path: {ik_path})")
        
        # DIRECTION 2: Extract coordinate names from Moco paths and search in IK paths
        print(f"\n=== Direction 2: Matching Moco coordinates to IK paths ===")
        for moco_path in moco_coords:
            # Skip if already matched
            if moco_path in coord_mapping.keys():
                continue
                
            moco_coord_name = extract_coordinate_name(moco_path)
            print(f"Processing Moco path '{moco_path}' -> coordinate name '{moco_coord_name}'")
            
            # Find available IK paths (not already matched)
            available_ik = [i for i in ik_coords if i not in coord_mapping.values()]
            
            best_match = try_match_coordinate(moco_coord_name, available_ik, moco_path, "Moco->IK")
            
            if best_match:
                coord_mapping[moco_path] = best_match
            else:
                print(f"  No match found for Moco coordinate: {moco_coord_name} (from path: {moco_path})")
        
        print(f"\n=== Final Matching Results ===")
        print(f"Total matches found: {len(coord_mapping)}")
        
        return coord_mapping
    
    print(f"Loading Moco solution: {moco_solution_file}")
    print(f"Loading IK kinematics: {ik_kinematics_file}")
    
    # Extract coordinate weights if setup file is provided
    coordinate_weights = {}
    if omoco_setup_file and include_weighted_quotients:
        try:
            print(f"Extracting coordinate weights from: {omoco_setup_file}")
            coordinate_weights = extract_coordinate_weights_from_omoco(omoco_setup_file, verbose=False)
            print(f"Extracted {len(coordinate_weights)} coordinate weights")
        except Exception as e:
            print(f"Warning: Could not extract weights from setup file: {e}")
            print("Continuing without weighted error calculations...")
    

    # Try to find model file if not provided
    if model_file is None:
        # Try common patterns to find the model file
        base_dir = os.path.dirname(moco_solution_file)
        potential_model_files = [
            os.path.join(base_dir, f) for f in os.listdir(base_dir) 
            if f.endswith('_SCALED.osim') or f.endswith('.osim')
        ]
        
        if potential_model_files:
            model_file = potential_model_files[0]
            print(f"Auto-detected model file: {model_file}")
        else:
            # Try parent directory
            parent_dir = os.path.dirname(base_dir)
            if os.path.exists(parent_dir):
                potential_model_files = [
                    os.path.join(parent_dir, f) for f in os.listdir(parent_dir) 
                    if f.endswith('_SCALED.osim') or f.endswith('.osim')
                ]
                if potential_model_files:
                    model_file = potential_model_files[0]
                    print(f"Auto-detected model file in parent dir: {model_file}")
    
    if model_file is None or not os.path.exists(model_file):
        raise ValueError(f"Model file not found. Please provide a valid model_file parameter.")
    
    print(f"Using model file: {model_file}")
    
    # Create temporary files for unit-converted data
    temp_dir = os.path.dirname(output_file)
    # os.makedirs(temp_dir, exist_ok=True)
    
    moco_temp_file = os.path.join(temp_dir, "temp_moco_degrees.sto")
    ik_temp_file = os.path.join(temp_dir, "temp_ik_degrees.sto")
    
    try:
        # Convert Moco solution from radians to degrees
        print("Converting Moco solution to degrees...")
        kinematicsToStates(
            kinematicsFileName=moco_solution_file,
            osimModelFileName=model_file,
            outputFileName=moco_temp_file,
            inDegrees=False,    # Moco solutions are typically in radians
            outDegrees=True,     # Convert to degrees
            include_speeds=False,
            include_accelerations=False
        )
        
        # Convert IK file to ensure it's in degrees (may already be, but this ensures consistency)
        print("Converting IK kinematics to degrees...")
        kinematicsToStates(
            kinematicsFileName=ik_kinematics_file,
            osimModelFileName=model_file,
            outputFileName=ik_temp_file,
            inDegrees=True,     # input ik files are expected to be in radians
            outDegrees=True,     # Keep in degrees
            include_speeds=False,
            include_accelerations=False
        )
        ik_temp_file = ik_kinematics_file
        
        # Now load the converted files
        moco_data = load_sto_file(moco_temp_file)
        ik_data = load_sto_file(ik_temp_file)
        
    except Exception as e:
        print(f"Error during unit conversion: {e}")
        print("Attempting to load files without conversion...")
        # Fallback: load original files without conversion
        moco_data = load_sto_file(moco_solution_file)
        ik_data = load_sto_file(ik_kinematics_file)
        print("Warning: Files loaded without unit conversion. Results may be inaccurate if units differ.")
    
    # finally:
    #     # Clean up temporary files
    #     for temp_file in [moco_temp_file, ik_temp_file]:
    #         if os.path.exists(temp_file):
    #             try:
    #                 os.remove(temp_file)
    #             except:
    #                 pass  # Ignore cleanup errors
    
    # Ensure time columns exist
    moco_time_col = 'time' if 'time' in moco_data.columns else moco_data.columns[0]
    ik_time_col = 'time' if 'time' in ik_data.columns else ik_data.columns[0]
    
    print(f"Moco solution time range: {moco_data[moco_time_col].min():.3f} to {moco_data[moco_time_col].max():.3f}")
    print(f"IK data time range: {ik_data[ik_time_col].min():.3f} to {ik_data[ik_time_col].max():.3f}")
    
    # Match coordinate columns between datasets
    coord_mapping = match_coordinates(moco_data.columns, ik_data.columns)
    
    if not coord_mapping:
        raise ValueError("No matching coordinates found between Moco solution and IK data")
    
    print(f"Found {len(coord_mapping)} matching coordinates:")
    for moco_col, ik_col in coord_mapping.items():
        print(f"  {moco_col} <-> {ik_col}")
    
    # Interpolate IK data to match Moco solution time points
    moco_time = moco_data[moco_time_col].values
    ik_time = ik_data[ik_time_col].values
    
    # Find overlapping time range
    time_start = max(moco_time.min(), ik_time.min())
    time_end = min(moco_time.max(), ik_time.max())
    
    # Filter Moco data to overlapping time range
    moco_mask = (moco_time >= time_start) & (moco_time <= time_end)
    moco_time_filtered = moco_time[moco_mask]
    
    print(f"Analysis time range: {time_start:.3f} to {time_end:.3f}")
    print(f"Number of Moco time points: {len(moco_time_filtered)}")
    
    # Create interpolation functions for IK data
    ik_interp_funcs = {}
    for moco_col, ik_col in coord_mapping.items():
        if not np.all(np.isnan(ik_data[ik_col])):
            # Create interpolation function (linear interpolation)
            ik_interp_funcs[moco_col] = interp1d(
                ik_time, ik_data[ik_col], 
                kind='linear', bounds_error=False, fill_value='extrapolate'
            )
    
    # Calculate errors for each coordinate at each time point
    error_data = {'time': moco_time_filtered}
    coordinate_errors = {}
    coordinate_weights_used = {}
    
    for moco_col in coord_mapping.keys():
        if moco_col in ik_interp_funcs:
            # Get Moco values at filtered time points
            moco_values = moco_data[moco_col].values[moco_mask]
            
            # Get interpolated IK values
            ik_values_interp = ik_interp_funcs[moco_col](moco_time_filtered)
            
            # Calculate error (Moco - IK)
            error_values = moco_values - ik_values_interp
            
            coordinate_errors[moco_col] = error_values
            error_data[f'{moco_col}_error'] = error_values
            error_data[f'{moco_col}_moco'] = moco_values
            error_data[f'{moco_col}_ik'] = ik_values_interp
            
            # Store weight info for later summary calculation
            if coordinate_weights:
                coord_name = extract_coordinate_name(moco_col)
                
                # Find matching weight (exact match first, then partial)
                weight = None
                if coord_name in coordinate_weights:
                    weight = coordinate_weights[coord_name]
                else:
                    # Try partial matching
                    for weight_key, weight_value in coordinate_weights.items():
                        if weight_key in coord_name or coord_name in weight_key:
                            weight = weight_value
                            break
                
                if weight is not None and weight > 0:
                    coordinate_weights_used[moco_col] = weight
                    print(f"Found weight for {coord_name}: {weight}")
                else:
                    print(f"No weight found for {coord_name}, excluding from weighted summary")
    
    # gather all coordinates from error file that shouldn't be considered in RMS values
    def should_exclude_coordinate(coord):
        coord_name = extract_coordinate_name(coord)
        excluded_coords = ['subtalar_angle_r','subtalar_angle_l', 'arm_flex_r', 'arm_flex_l', 'pro_sup_r', 'pro_sup_l']
        return coord_name in excluded_coords

    # Calculate RMS error at each time point
    rms_errors = []
    max_error_coords = []
    max_error_values = []
    
    for i in range(len(moco_time_filtered)):
        # Get all coordinate errors at this time point
        time_point_errors = []
        time_point_coord_names = []
        time_point_abs_errors = []
        
        for coord, errors in coordinate_errors.items():
            if not np.isnan(errors[i]):
                time_point_errors.append(errors[i] ** 2)
                time_point_coord_names.append(coord)
                time_point_abs_errors.append(abs(errors[i]))
        
        # Calculate regular RMS error
        if time_point_errors:
            rms_error = np.sqrt(np.mean(time_point_errors))
            rms_errors.append(rms_error)
            
            # Find coordinate with maximum absolute error
            max_error_idx = np.argmax(time_point_abs_errors)
            max_error_coords.append(time_point_coord_names[max_error_idx])
            max_error_values.append(time_point_abs_errors[max_error_idx])
        else:
            rms_errors.append(np.nan)
            max_error_coords.append('')
            max_error_values.append(np.nan)
    
    # Add RMS and max error data
    error_data['rms_error_all_coords'] = rms_errors
    error_data['max_error_coordinate'] = max_error_coords
    error_data['max_error_value'] = max_error_values
    
    # Create the summary row with weighted error quotients
    if coordinate_weights_used and include_weighted_quotients:
        print(f"\nCreating weighted error quotient summary row...")
        
        # Calculate cumulative errors and weighted quotients for summary row
        summary_row_data = {}
        
        # Time column for summary row (use a distinctive value like -999 or the word "SUMMARY")
        summary_row_data['time'] = [-999.0]  # Distinctive time value for summary row
        
        # Initialize other non-coordinate columns for summary row
        summary_row_data['rms_error_all_coords'] = [np.nan]
        summary_row_data['max_error_coordinate'] = ['WEIGHTED_SUMMARY']
        summary_row_data['max_error_value'] = [np.nan]
        
        total_weighted_quotient = 0.0
        num_weighted_coords = 0
        coordinate_quotients = {}  # Store for ranking
        
        for coord, errors in coordinate_errors.items():
            # Calculate cumulative absolute error for this coordinate
            cumulative_error = np.nansum(np.abs(errors))
            
            # Add cumulative error to summary row
            error_col = f'{coord}_error'
            moco_col = f'{coord}_moco'
            ik_col = f'{coord}_ik'
            
            summary_row_data[error_col] = [cumulative_error]
            summary_row_data[moco_col] = [np.nan]  # Not meaningful for summary
            summary_row_data[ik_col] = [np.nan]    # Not meaningful for summary
            
            # Calculate weighted quotient if weight is available
            if coord in coordinate_weights_used and not should_exclude_coordinate(coord):
                weight = coordinate_weights_used[coord]
                weighted_quotient = cumulative_error / weight
                
                total_weighted_quotient += weighted_quotient
                num_weighted_coords += 1
                
                coord_name = extract_coordinate_name(coord)
                coordinate_quotients[coord_name] = weighted_quotient
                
                print(f"  {coord_name}: Cumulative_error={cumulative_error:.4f}, Weight={weight:.1f}, Quotient={weighted_quotient:.6f}")
        
        # Calculate overall weighted quotient for summary
        if num_weighted_coords > 0:
            overall_weighted_quotient = total_weighted_quotient / num_weighted_coords
            summary_row_data['max_error_value'] = [overall_weighted_quotient]
            print(f"  Overall weighted error quotient: {overall_weighted_quotient:.6f}")
        
        # Append summary row to the main error data
        for key in error_data.keys():
            if key in summary_row_data:
                # Convert to list if it's a numpy array, then extend
                if isinstance(error_data[key], np.ndarray):
                    error_data[key] = error_data[key].tolist()
                error_data[key].extend(summary_row_data[key])
            else:
                # For any columns we missed, add NaN to summary row
                if isinstance(error_data[key], np.ndarray):
                    error_data[key] = error_data[key].tolist()
                error_data[key].append(np.nan)
        
        print(f"Added weighted error quotient summary row at time = -999")
        
        # Create top 5 summary section
        if coordinate_quotients:
            print(f"\nCreating top 5 weighted error summary section...")
            
            # Sort coordinates by weighted quotient (highest first)
            sorted_quotients = sorted(coordinate_quotients.items(), key=lambda x: x[1], reverse=True)
            top_5_coords = sorted_quotients[:5]
            
            coord_rms_values = {}
            for coord, errors in coordinate_errors.items():
                coord_name = extract_coordinate_name(coord)
                coord_rms = np.sqrt(np.nanmean(errors**2))
                coord_rms_values[coord_name] = coord_rms
            else: 
                # mark excluded coords
                coord_rms_values[coord_name] = np.nan
                print(f" Excluding {coord_name} from rms calculations")

            # Add separator rows and top 5 summary
            separator_time_start = -998.0
            
            # Add separator row
            separator_row_data = {}
            separator_row_data['time'] = [separator_time_start]
            separator_row_data['rms_error_all_coords'] = [np.nan]
            separator_row_data['max_error_coordinate'] = ['TOP_5_SUMMARY_START']
            separator_row_data['max_error_value'] = [np.nan]
            
            # Initialize all other columns for separator
            for key in error_data.keys():
                if key not in separator_row_data:
                    separator_row_data[key] = [np.nan]
            
            # Add separator to main data
            for key in error_data.keys():
                if isinstance(error_data[key], list):
                    error_data[key].extend(separator_row_data[key])
                else:
                    error_data[key] = error_data[key].tolist()
                    error_data[key].extend(separator_row_data[key])
            
            # Add top 5 rows with their rms
            for rank, (coord_name, quotient) in enumerate(top_5_coords, 1):
                top5_row_data = {}
                top5_row_data['time'] = [separator_time_start + rank]  # -997, -996, -995, -994, -993
                
                coord_rms = coord_rms_values.get(coord_name, np.nan)
                top5_row_data['rms_error_all_coords'] = [coord_rms]  # Rank number
                
                top5_row_data['max_error_coordinate'] = [coord_name]  # Coordinate name
                top5_row_data['max_error_value'] = [quotient]  # Weighted quotient value
                
                # Initialize all other columns as NaN for top 5 rows
                for key in error_data.keys():
                    if key not in top5_row_data:
                        top5_row_data[key] = [np.nan]
                
                # Add to main data
                for key in error_data.keys():
                    error_data[key].extend(top5_row_data[key])
                
                print(f"  Rank {rank}: {coord_name} = {quotient:.6f}")
    
    # Create output DataFrame
    error_df = pd.DataFrame(error_data)
    
    # Write to .sto file
    write_sto_file(error_df, output_file)
    
    # Exclude any unwanted coordinates here
    filtered_errors = []
    for coord, errors in coordinate_errors.items():
        if not should_exclude_coordinate(coord):
            filtered_errors.extend(errors)
    # Calculate summary statistics
    overall_rms = np.sqrt(np.nanmean([err**2 for errors in coordinate_errors.values() for err in errors]))
    
    coord_rms_stats = {}
    for coord, errors in coordinate_errors.items():
        if not should_exclude_coordinate(coord):
            coord_rms = np.sqrt(np.nanmean(errors**2))
            coord_max_abs = np.nanmax(np.abs(errors))
            coord_mean_abs = np.nanmean(np.abs(errors))
            cumulative_abs_error = np.nansum(np.abs(errors))
            
            stats = {
                'rms_error': coord_rms,
                'max_abs_error': coord_max_abs,
                'mean_abs_error': coord_mean_abs,
                'cumulative_abs_error': cumulative_abs_error, 
                'excluded': False
            }
        else:
            # still compute stats for excluded coords, but just mark them
            coord_rms = np.sqrt(np.nanmean(errors**2))
            coord_max_abs = np.nanmax(np.abs(errors))
            coord_mean_abs = np.nanmean(np.abs(errors))
            cumulative_abs_error = np.nansum(np.abs(errors))

            stats = {
                'rms_error': coord_rms,
                'max_abs_error': coord_max_abs,
                'mean_abs_error': coord_mean_abs,
                'cumulative_abs_error': cumulative_abs_error, 
                'excluded': True
            }
            
        # Add weighted statistics if available
        if coord in coordinate_weights_used and include_weighted_quotients:
            weight = coordinate_weights_used[coord]
            weighted_cumulative_quotient = coord_rms / weight
            
            stats.update({
                'weight': weight,
                'weighted_cumulative_quotient': weighted_cumulative_quotient
            })
        
        coord_rms_stats[coord] = stats
        
    # Print summary
    print(f"\nError Analysis Summary:")
    print(f"Overall RMS error across all coordinates: {overall_rms:.4f}")
    print(f"Mean RMS error across time: {np.nanmean(rms_errors):.4f}")
    print(f"Max RMS error across time: {np.nanmax(rms_errors):.4f}")
    
    if coordinate_weights_used and include_weighted_quotients:
        weighted_quotients = [stats['weighted_cumulative_quotient'] for stats in coord_rms_stats.values() 
                            if 'weighted_cumulative_quotient' in stats]
        if weighted_quotients:
            print(f"Mean weighted cumulative quotient: {np.mean(weighted_quotients):.6f}")
    
    print(f"\nPer-coordinate statistics:")
    for coord, stats in coord_rms_stats.items():
        coord_name = extract_coordinate_name(coord)
        base_stats = f"RMS={stats['rms_error']:.4f}, Cumulative_abs={stats['cumulative_abs_error']:.4f}"
        
        if stats.get('excluded', False):
            exclusion_note = '(Excluded from rankings)'
        else:
            exclusion_note = ''

        if 'weight' in stats:
            weighted_stats = f", Weight={stats['weight']:.1f}, Cumulative_quotient={stats['weighted_cumulative_quotient']:.6f}"
            print(f"  {coord_name}: {base_stats}{weighted_stats}")
        else:
            print(f"  {coord_name}: {base_stats}")
    
    result_dict = {
        'overall_rms_error': overall_rms,
        'time_series_rms_errors': rms_errors,
        'coordinate_statistics': coord_rms_stats,
        'output_file': output_file,
        'num_matched_coordinates': len(coord_mapping),
        'time_range': (time_start, time_end),
        'num_time_points': len(moco_time_filtered)
    }
    
    # Add weighted results if available
    if coordinate_weights_used and include_weighted_quotients:
        result_dict.update({
            'coordinate_weights_used': coordinate_weights_used,
            'num_weighted_coordinates': len(coordinate_weights_used)
        })
    
    return result_dict

def compute_weighted_error_quotient(moco_solution_file, ik_kinematics_file, coordinate_weights, 
                                   model_file=None, temp_dir=None, verbose=True):
    """
    Compute a single weighted error quotient by dividing total coordinate errors by their Moco weights.
    
    This function uses analyze_moco_ik_tracking_errors() to get detailed error analysis, then computes
    a weighted error quotient that normalizes tracking performance by coordinate importance.
    
    Parameters:
    -----------
    moco_solution_file : str
        Path to the Moco solution .sto file
    ik_kinematics_file : str
        Path to the IK kinematics .sto file
    coordinate_weights : dict
        Dictionary mapping coordinate names to their Moco tracking weights
        Example: {'pelvis_tilt': 10.0, 'hip_flexion_r': 1.0, 'knee_angle_r': 1.0}
    model_file : str, optional
        Path to the OpenSim model file (.osim)
    temp_dir : str, optional
        Directory for temporary files (defaults to same as solution file)
    verbose : bool
        Whether to print detailed output
        
    Returns:
    --------
    dict
        Dictionary containing:
        - 'weighted_error_quotient': The final quotient value
        - 'total_weighted_error': Sum of (coordinate_error / coordinate_weight)
        - 'total_weight': Sum of all coordinate weights used
        - 'coordinate_quotients': Individual quotients per coordinate
        - 'detailed_results': Full results from analyze_moco_ik_tracking_errors
    """
    
    if temp_dir is None:
        temp_dir = os.path.dirname(moco_solution_file)
    
    # Create temporary output file for detailed analysis
    temp_error_file = os.path.join(temp_dir, "temp_detailed_errors.sto")
    
    try:
        # Run the detailed error analysis
        if verbose:
            print("Running detailed error analysis...")
        
        detailed_results = analyze_moco_ik_tracking_errors(
            moco_solution_file=moco_solution_file,
            ik_kinematics_file=ik_kinematics_file,
            output_file=temp_error_file,
            model_file=model_file
        )
        
        if verbose:
            print(f"Found {detailed_results['num_matched_coordinates']} matching coordinates")
        
        # Extract coordinate statistics
        coord_stats = detailed_results['coordinate_statistics']
        
        # Calculate weighted error quotients for each coordinate
        coordinate_quotients = {}
        total_weighted_error = 0.0
        total_weight = 0.0
        used_coordinates = []
        missing_weights = []
        
        for coord_column, stats in coord_stats.items():
            # Extract coordinate name from column (remove path prefixes)
            coord_name = extract_coordinate_name(coord_column)
            
            # Find matching weight (try exact match first, then partial match)
            weight = None
            
            # Try exact match
            if coord_name in coordinate_weights:
                weight = coordinate_weights[coord_name]
            else:
                # Try partial matching (useful for path-style names)
                for weight_key, weight_value in coordinate_weights.items():
                    if weight_key in coord_name or coord_name in weight_key:
                        weight = weight_value
                        if verbose:
                            print(f"Matched '{coord_name}' with weight key '{weight_key}'")
                        break
            
            if weight is not None and weight > 0:
                # Calculate error quotient for this coordinate
                coord_rms_error = stats['rms_error']
                coord_quotient = coord_rms_error / weight
                
                coordinate_quotients[coord_name] = {
                    'rms_error': coord_rms_error,
                    'weight': weight,
                    'quotient': coord_quotient,
                    'column_name': coord_column
                }
                
                total_weighted_error += coord_quotient
                total_weight += weight
                used_coordinates.append(coord_name)
                
                if verbose:
                    print(f"  {coord_name}: RMS={coord_rms_error:.4f}, Weight={weight:.1f}, Quotient={coord_quotient:.6f}")
            
            else:
                missing_weights.append(coord_name)
                if verbose:
                    print(f"  {coord_name}: No weight found - excluding from quotient calculation")
        
        # Calculate final weighted error quotient
        if total_weight > 0:
            weighted_error_quotient = total_weighted_error / len(used_coordinates) if used_coordinates else 0
        else:
            weighted_error_quotient = float('inf')
            print("Warning: No valid weights found, cannot compute quotient")
        
        # Summary
        if verbose:
            print(f"\nWeighted Error Quotient Summary:")
            print(f"  Used coordinates: {len(used_coordinates)}")
            print(f"  Missing weights: {len(missing_weights)}")
            print(f"  Total weighted error: {total_weighted_error:.6f}")
            print(f"  Average weighted error quotient: {weighted_error_quotient:.6f}")
            
            if missing_weights:
                print(f"  Coordinates without weights: {missing_weights}")
        
        return {
            'weighted_error_quotient': weighted_error_quotient,
            'total_weighted_error': total_weighted_error,
            'total_weight': total_weight,
            'coordinate_quotients': coordinate_quotients,
            'used_coordinates': used_coordinates,
            'missing_weights': missing_weights,
            'detailed_results': detailed_results
        }
        
    finally:
        # Clean up temporary error file
        if os.path.exists(temp_error_file):
            try:
                os.remove(temp_error_file)
            except:
                pass

def extract_coordinate_name(column_name):
    """
    Extract coordinate name from different column naming conventions.
    Enhanced version that handles more naming patterns.
    """
    # Remove common prefixes and suffixes
    name = column_name.replace('/jointset/', '').replace('/value', '').replace('/speed', '')
    
    # Handle path-style naming (e.g., '/jointset/ground_pelvis/pelvis_tx/value')
    if '/' in name:
        parts = name.split('/')
        # Usually the coordinate name is the last part
        if len(parts) >= 1:
            return parts[-1]
    
    return name

def extract_coordinate_weights_from_omoco(omoco_setup_file, verbose=True):
    """
    Extract coordinate tracking weights from a Moco setup (.omoco) XML file.
    
    Parameters:
    -----------
    omoco_setup_file : str
        Path to the .omoco setup file generated by study.printToXML()
    verbose : bool
        Whether to print detailed information about found weights
        
    Returns:
    --------
    dict
        Dictionary mapping coordinate names to their tracking weights
    """
    
    if not os.path.exists(omoco_setup_file):
        raise FileNotFoundError(f"Moco setup file not found: {omoco_setup_file}")
    
    try:
        import xml.etree.ElementTree as ET
        
        if verbose:
            print(f"Parsing Moco setup file: {omoco_setup_file}")
        
        tree = ET.parse(omoco_setup_file)
        root = tree.getroot()
        
        coordinate_weights = {}
        
        # Look for MocoStateTrackingGoal (this is where coordinate weights are stored)
        state_tracking_goals = root.findall('.//MocoStateTrackingGoal')
        
        if not state_tracking_goals:
            print("Warning: No MocoStateTrackingGoal found in setup file")
            return coordinate_weights
        
        for goal_idx, goal in enumerate(state_tracking_goals):
            goal_name = goal.get('name', f'StateTracking_{goal_idx}')
            
            if verbose:
                print(f"\nFound MocoStateTrackingGoal: {goal_name}")
            
            # Look for state_weights/MocoWeightSet structure
            state_weights_elements = goal.findall('.//MocoWeightSet[@name="state_weights"]')
            
            for weight_set in state_weights_elements:
                # Find the objects container and then MocoWeight elements
                objects_container = weight_set.find('objects')
                if objects_container is not None:
                    weights = objects_container.findall('MocoWeight')
                    
                    for weight_elem in weights:
                        # Get the state path and weight value
                        state_path = weight_elem.get('name', '')
                        weight_value_elem = weight_elem.find('weight')
                        
                        if weight_value_elem is not None:
                            weight_value = weight_value_elem.text
                        else:
                            weight_value = '1.0'
                        
                        try:
                            weight_value = float(weight_value)
                        except (ValueError, TypeError):
                            if verbose:
                                print(f"  Warning: Invalid weight value '{weight_value}' for {state_path}")
                            continue
                        
                        if state_path and '/value' in state_path:
                            # Extract coordinate name from state path
                            # Example: "/jointset/ground_pelvis/pelvis_tilt/value" -> "pelvis_tilt"
                            coord_name = extract_coordinate_name_from_path(state_path)
                            
                            if coord_name:
                                coordinate_weights[coord_name] = weight_value
                                
                                if verbose:
                                    print(f"  {coord_name}: {weight_value}")
        
        # Also check for alternative weight structures (fallback)
        if not coordinate_weights:
            if verbose:
                print("No weights found in standard structure, trying alternative formats...")
            
            # Look for any MocoWeight elements anywhere in the file
            all_weights = root.findall('.//MocoWeight')
            for weight_elem in all_weights:
                state_path = weight_elem.get('name', '')
                weight_value_elem = weight_elem.find('weight')
                
                if weight_value_elem is not None and '/value' in state_path:
                    try:
                        weight_value = float(weight_value_elem.text)
                        coord_name = extract_coordinate_name_from_path(state_path)
                        if coord_name:
                            coordinate_weights[coord_name] = weight_value
                            if verbose:
                                print(f"  Alternative format - {coord_name}: {weight_value}")
                    except (ValueError, TypeError):
                        continue
        
        if verbose:
            print(f"\nExtracted {len(coordinate_weights)} coordinate weights:")
            for coord, weight in sorted(coordinate_weights.items()):
                print(f"  {coord}: {weight}")
        
        if not coordinate_weights:
            print("Warning: No coordinate weights found in the setup file")
            print("This might indicate the file format is different or weights are stored elsewhere")
        
        return coordinate_weights
        
    except ET.ParseError as e:
        raise ValueError(f"Error parsing XML file {omoco_setup_file}: {e}")
    except Exception as e:
        raise ValueError(f"Error extracting weights from {omoco_setup_file}: {e}")

def extract_coordinate_name_from_path(state_path):
    """
    Extract coordinate name from OpenSim state path.
    
    Examples:
    "/jointset/ground_pelvis/pelvis_tilt/value" -> "pelvis_tilt"
    "/jointset/hip_r/hip_flexion_r/value" -> "hip_flexion_r"
    """
    
    # Remove /value suffix
    path = state_path.replace('/value', '').replace('/speed', '')
    
    # Split by / and get parts
    parts = path.split('/')
    
    # The coordinate name is typically the last non-empty part
    for part in reversed(parts):
        if part and part != 'jointset':
            return part
    
    return None

def auto_compute_weighted_error_quotient(moco_solution_file, ik_kinematics_file, 
                                       omoco_setup_file, model_file=None, 
                                       temp_dir=None, verbose=True):
    """
    Compute weighted error quotient using coordinate weights automatically extracted 
    from the Moco setup file.
    
    Parameters:
    -----------
    moco_solution_file : str
        Path to the Moco solution .sto file
    ik_kinematics_file : str
        Path to the IK kinematics .sto file
    omoco_setup_file : str
        Path to the .omoco setup file generated by study.printToXML()
    model_file : str, optional
        Path to the OpenSim model file (.osim)
    temp_dir : str, optional
        Directory for temporary files
    verbose : bool
        Whether to print detailed output
        
    Returns:
    --------
    dict
        Dictionary containing weighted error quotient results and extracted weights
    """
    
    if verbose:
        print("=== Automatic Weighted Error Quotient Analysis ===")
        print("Step 1: Extracting coordinate weights from Moco setup file...")
    
    # Extract coordinate weights from the setup file
    try:
        coordinate_weights = extract_coordinate_weights_from_omoco(omoco_setup_file, verbose=verbose)
        
        if not coordinate_weights:
            raise ValueError("No coordinate weights found in the setup file")
            
    except Exception as e:
        print(f"Error extracting weights: {e}")
        raise
    
    if verbose:
        print(f"\nStep 2: Computing weighted error quotient...")
    
    # Use the extracted weights to compute the error quotient
    results = compute_weighted_error_quotient(
        moco_solution_file=moco_solution_file,
        ik_kinematics_file=ik_kinematics_file,
        coordinate_weights=coordinate_weights,
        model_file=model_file,
        temp_dir=temp_dir,
        verbose=verbose
    )
    
    # Add the extracted weights to the results
    results['extracted_weights'] = coordinate_weights
    results['omoco_setup_file'] = omoco_setup_file
    
    return results


def plot_moco_solution_columns(moco_solution_file, 
                               columns_to_plot,
                               model_file=None,
                               subject_id=None,
                               existing_fig_ax=None,
                               normalize_time=True,
                               time_as_percentage=False,
                               plot_title=None,
                               output_file=None,
                               line_style='-',
                               line_width=2,
                               color=None,
                               label_prefix=None,
                               show_legend=True,
                               legend_location='upper left',
                               grid=True,
                               xlabel=None,
                               ylabel=None,
                               figsize=(12, 8)):
    """
    Plot specified columns from a Moco solution .sto file, with optional superimposition
    on an existing figure.
    
    Parameters:
    -----------
    moco_solution_file : str
        Path to the Moco solution .sto file
    columns_to_plot : list of str
        List of column names to plot (e.g., ['hip_flexion_r', 'knee_angle_r'])
    existing_fig_ax : tuple, optional
        Existing (figure, axes) tuple to add plots to. If None, creates new figure
    normalize_time : bool, optional (default=True)
        If True, normalizes time to start at 0
    time_as_percentage : bool, optional (default=False)
        If True, converts time to percentage of cycle (0-100%)
    plot_title : str, optional
        Custom title for the plot
    output_file : str, optional
        Path to save the plot image
    line_style : str, optional (default='-')
        Matplotlib line style (e.g., '-', '--', '-.', ':')
    line_width : float, optional (default=2)
        Width of plot lines
    color : str or list, optional
        Color(s) for the lines. If None, uses default colormap
    label_prefix : str, optional
        Prefix to add to legend labels (useful when comparing multiple files)
    show_legend : bool, optional (default=True)
        Whether to display the legend
    legend_location : str, optional (default='upper left')
        Location for the legend
    grid : bool, optional (default=True)
        Whether to show grid lines
    xlabel : str, optional
        Custom x-axis label
    ylabel : str, optional
        Custom y-axis label
    figsize : tuple, optional (default=(12, 8))
        Figure size (width, height) in inches (only used if creating new figure)
        
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects for the plot
        
    Examples:
    ---------
    # Create a new plot
    fig, ax = plot_moco_solution_columns(
        'solution1.sto',
        columns_to_plot=['hip_flexion_r', 'knee_angle_r']
    )
    
    # Add to existing plot for comparison
    fig, ax = plot_moco_solution_columns(
        'solution2.sto',
        columns_to_plot=['hip_flexion_r', 'knee_angle_r'],
        existing_fig_ax=(fig, ax),
        line_style='--',
        label_prefix='Trial 2'
    )
    
    plt.show()
    """
    
    # Load the .sto file
    print(f"Loading Moco solution from: {moco_solution_file}")
    temp_converted_file = 'temp_degrees_solution.sto'
    # kinematicsToStates(kinematicsFileName=moco_solution_file,
    #                    osimModelFileName=model_file,
    #                    outputFileName=temp_converted_file,
    #                    inDegrees=False,
    #                    outDegrees=True,
    #                    include_speeds=False,
    #                    include_accelerations=False)
    
    # file_path = moco_solution_file[1] if isinstance(moco_solution_file, tuple) else moco_solution_file

    kinematicsToStatesForMoco(mocoSolutionFileName=moco_solution_file, outputFileName=temp_converted_file,
                              inDegrees=False, outDegrees=True, verbose=True)

    # Read the .sto file
    with open(temp_converted_file, 'r') as f:
        lines = f.readlines()
    
    # Find the line with 'endheader'
    endheader_idx = None
    for i, line in enumerate(lines):
        if 'endheader' in line.lower():
            endheader_idx = i
            break
    
    if endheader_idx is None:
        raise ValueError(f"No 'endheader' found in file: {temp_converted_file}")
    
    # Headers are in the line after 'endheader'
    header_line = lines[endheader_idx + 1]
    headers = [h.strip() for h in header_line.split('\t') if h.strip()]
    
    # Data starts two lines after 'endheader'
    data_start_idx = endheader_idx + 2
    
    # Parse the data
    data_dict = {header: [] for header in headers}
    
    for i in range(data_start_idx, len(lines)):
        line = lines[i].strip()
        if not line:  # Skip empty lines
            continue
        
        values = line.split('\t')
        if len(values) >= len(headers):
            for j, header in enumerate(headers):
                try:
                    data_dict[header].append(float(values[j]))
                except (ValueError, IndexError):
                    data_dict[header].append(np.nan)
    
    df = pd.DataFrame(data_dict)
    
    # Get time column
    time_col = 'time' if 'time' in df.columns else df.columns[0]
    time_data = df[time_col].values
    
    # Normalize time if requested
    if normalize_time:
        time_data = time_data - time_data[0]
    
    # Convert to percentage if requested
    if time_as_percentage:
        time_duration = time_data[-1] - time_data[0]
        time_data = (time_data / time_duration) * 100
    
    # Validate columns exist
    available_columns = [col for col in columns_to_plot if col in df.columns]
    missing_columns = [col for col in columns_to_plot if col not in df.columns]
    
    if missing_columns:
        print(f"Warning: The following columns were not found in the file: {missing_columns}")
    
    if not available_columns:
        raise ValueError("None of the specified columns were found in the file")
    
    print(f"Plotting columns: {available_columns}")
    
    # Create or use existing figure
    if existing_fig_ax is not None:
        fig, ax = existing_fig_ax
        print("Adding to existing figure")
    else:
        fig, ax = plt.subplots(figsize=figsize)
        print("Creating new figure")
    
    # Generate colors if not provided
    if color is None:
        colors = plt.cm.tab10(np.linspace(0, 1, len(available_columns)))
    elif isinstance(color, str):
        colors = [color] * len(available_columns)
    else:
        colors = color
    
    # Plot each column
    for i, col in enumerate(available_columns):
        y_data = df[col].values
        
        # Create label with subject ID
        if i == 0 and subject_id:
            label = subject_id
        else:
            label = None
        
        # label = ' - '.join(label_parts)
        print(f"DEBUG: Created label: '{label}'")
        
        # label = subject_id

        # Plot the data
        ax.plot(time_data, y_data, 
                linestyle=line_style, 
                linewidth=line_width,
                color=colors[i] if hasattr(colors[i], '__iter__') and not isinstance(colors[i], str) else colors[i],
                label=label)
    
    # Customize plot
    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        if time_as_percentage:
            ax.set_xlabel('Gait Cycle (%)')
        else:
            ax.set_xlabel('Time (s)')
    
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel('Value')
    
    if plot_title:
        ax.set_title(plot_title)
    
    if show_legend:
        ax.legend(loc=legend_location, bbox_to_anchor=(1.05, 1) if legend_location == 'upper left' else None)
    
    if grid:
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot if output file specified
    if output_file:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")
    
    return fig, ax


def write_mtp_spring_forces(solution = None, model = None, output_file=None):

    """
    *Takes in solution that has to already passed through osim.MocoTrajectory(), processed model through osim.ModelProcessor().process()
    ^ Then outputs file into a given path
    """
    # model_proc = osim.ModelProcessor(model_file)
    # model = model_proc.process()
    # solution = osim.MocoTrajectory(moco_solution)
    forceset = model.getForceSet()
    # Initialize force dictionary
    spring_forces = {}
    
    for i in range(forceset.getSize()):
        force = forceset.get(i)
        if "mtp_spring" in force.getName():
            # Store index w/ force name as the key
            spring_forces[force.getName()] = i
            print(f"Found spring force: {force.getName()} at index {i}")

    spring_r = forceset.get(spring_forces["mtp_spring_r"])
    spring_l = forceset.get(spring_forces["mtp_spring_l"])
    spring_array_r_vector = osim.Vector(solution.getNumTimes(), 0.0)
    spring_array_l_vector = osim.Vector(solution.getNumTimes(), 0.0)
    states_problem_trajectory = solution.exportToStatesTrajectory(model)
    state = model.initSystem()

    # Use timeseries table to retrieve time vector
    # solution_time_series = osim.TimeSeriesTable(solution)
    solution_table = solution.exportToStatesTable()
    time = solution_table.getIndependentColumn()
    
    for i in range(solution.getNumTimes()):
        # getting state object at current time step
        current_state = states_problem_trajectory.get(i)
        # Realizing the system to velocity stage is needed for forces.
        model.realizeVelocity(current_state)

        # Using the spring_forces dictionary, store force values in pre-initialized vector objects
        # Get single force values from each spring at current state
        spring_array_r = spring_r.getRecordValues(current_state).get(0)
        spring_array_l = spring_l.getRecordValues(current_state).get(0)
        spring_array_r_vector.set(i, spring_array_r)
        spring_array_l_vector.set(i, spring_array_l)

    # given that you already have the time vector given from timeseries operation above, establish a matrix for holding two columns of data
    spring_matrix = osim.Matrix(solution.getNumTimes(), 2)
    # Populate matrix with spring force data
    for i in range(solution.getNumTimes()):
        spring_matrix.set(i, 0, spring_array_l_vector.get(i))
        spring_matrix.set(i, 1, spring_array_r_vector.get(i))
        
    # create vector of column labels
    column_labels = osim.StdVectorString()
    column_labels.append("mtp_spring_l")
    column_labels.append("mtp_spring_r")

    # create NEW time series table and write to file
    spring_table = osim.TimeSeriesTable(time, spring_matrix, column_labels)
    osim.STOFileAdapter.write(spring_table, output_file)
    print(f"Spring force data written to: {output_file}")
    # Just for debugging in output window
    print(spring_array_l_vector)
    print(spring_array_r_vector)


def write_sto_file(dataframe, filename):
    """
    Write a pandas DataFrame to OpenSim .sto file format
    """
    with open(filename, 'w') as f:
        # Write header
        f.write(f"{filename}\n")
        f.write("version=1\n")
        f.write(f"nRows={len(dataframe)}\n")
        f.write(f"nColumns={len(dataframe.columns)}\n")
        f.write("inDegrees=yes\n")
        f.write("endheader\n")
        
        # Write column headers
        f.write('\t'.join(dataframe.columns) + '\n')
        
        # Write data
        for _, row in dataframe.iterrows():
            formatted_row = []
            for col in dataframe.columns:
                value = row[col]
                if isinstance(value, str):
                    formatted_row.append(value)
                elif pd.isna(value):
                    formatted_row.append('0.0')
                else:
                    formatted_row.append(f"{value:.6f}")
            f.write('\t'.join(formatted_row) + '\n')
    
    print(f"Error analysis results written to: {filename}")

def simtk2numpy(simtk_vector):
    """
    Convert an OpenSim SimTK Vector to a NumPy array
    """
    # Instantiate an array the same size as the input vec
    # array = np.empty(simtk_vector.size())
    # for i in range(simtk_vector.size()):
    #     array[i] = simtk_vector.get(i)
    # return array
    return np.array([simtk_vector[i] for i in range(simtk_vector.size())])
    
def mass_normalize_grf(df: pd.DataFrame,
                       model_path: str,
                       grf_columns: Optional[List[str]] = None) -> pd.DataFrame:
    
    """
    Normalize GRF data by body mass from OpenSim model.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing GRF data
    model_path : str
        Path to OpenSim model file (.osim)
    grf_columns : list of str, optional
        Specific GRF column names to normalize. If None, auto-detects GRF columns
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with normalized GRF columns (values divided by body mass)
    """

    try:
        # Load model and get total mass
        model = osim.Model(model_path)
        state = model.initSystem()
        total_mass = model.getTotalMass(state)
        
        print(f"Model total mass: {total_mass:.2f} kg")
        
    except Exception as e:
        raise RuntimeError(f"Error loading OpenSim model: {e}")
    
    # Auto-detect GRF columns if not provided
    if grf_columns is None:
        grf_columns = []
        grf_patterns = [
            'ground_force', 'grf', 'force_vx', 'force_vy', 'force_vz',
            'force_px', 'force_py', 'force_pz', 'moment_x', 'moment_y', 'moment_z',
            'torque_x', 'torque_y', 'torque_z'
        ]
        
        for col in df.columns:
            col_lower = col.lower()
            if any(pattern in col_lower for pattern in grf_patterns):
                grf_columns.append(col)
        
        if grf_columns:
            print(f"Auto-detected GRF columns: {grf_columns}")
        else:
            print("Warning: No GRF columns detected automatically")
            return df.copy()  # Return unchanged if no GRF columns found
    
    # Validate that specified columns exist
    missing_columns = [col for col in grf_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Specified GRF columns not found in data: {missing_columns}")
    
    # Create copy of dataframe for normalization
    df_normalized = df.copy()
    
    # Normalize GRF columns by body mass
    for col in grf_columns:
        df_normalized[col] = df[col] / total_mass
        print(f"Normalized column: {col}")
    
    print(f"Normalized {len(grf_columns)} GRF columns by body mass ({total_mass:.2f} kg)")
    
    return df_normalized

def mass_normalize_grf_by_mass(df: pd.DataFrame,
                              subject_mass: float,
                              grf_columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Normalize GRF data by a specific subject mass value.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing GRF data
    subject_mass : float
        Subject mass in kg
    grf_columns : list of str, optional
        Specific GRF column names to normalize. If None, auto-detects GRF columns
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with normalized GRF columns (values divided by body mass)
    """
    
    if subject_mass <= 0:
        raise ValueError(f"Invalid subject mass: {subject_mass}")
    
    # Auto-detect GRF columns if not provided
    if grf_columns is None:
        grf_columns = []
        grf_patterns = [
            'ground_force', 'grf', 'force_vx', 'force_vy', 'force_vz',
            'force_px', 'force_py', 'force_pz', 'moment_x', 'moment_y', 'moment_z',
            'torque_x', 'torque_y', 'torque_z'
        ]
        
        for col in df.columns:
            col_lower = col.lower()
            if any(pattern in col_lower for pattern in grf_patterns):
                grf_columns.append(col)
        
        if grf_columns:
            print(f"Auto-detected GRF columns: {grf_columns}")
        else:
            print("Warning: No GRF columns detected automatically")
            return df.copy()  # Return unchanged if no GRF columns found
    
    # Validate that specified columns exist
    missing_columns = [col for col in grf_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Specified GRF columns not found in data: {missing_columns}")
    
    # Create copy of dataframe for normalization
    df_normalized = df.copy()
    
    # Normalize GRF columns by body mass
    for col in grf_columns:
        df_normalized[col] = df[col] / subject_mass
        print(f"Normalized column {col} by mass {subject_mass} kg")
    
    print(f"Normalized {len(grf_columns)} GRF columns by subject mass ({subject_mass} kg)")
    
    return df_normalized

def normalize_grf_by_manual_body_weight(df: pd.DataFrame, 
                                       body_weight: float,
                                       grf_columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Normalize GRF data by manually provided body weight.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame containing GRF data
    body_weight : float
        Body weight in kg
    grf_columns : list of str, optional
        Specific GRF column names to normalize. If None, auto-detects GRF columns
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with normalized GRF columns (values divided by body weight)
    """
    
    if body_weight <= 0:
        raise ValueError(f"Body weight must be positive, got: {body_weight}")
    
    # Auto-detect GRF columns if not provided
    if grf_columns is None:
        grf_columns = []
        grf_patterns = [
            'ground_force', 'grf', 'force_vx', 'force_vy', 'force_vz',
            'force_px', 'force_py', 'force_pz', 'moment_x', 'moment_y', 'moment_z',
            'torque_x', 'torque_y', 'torque_z'
        ]
        
        for col in df.columns:
            col_lower = col.lower()
            if any(pattern in col_lower for pattern in grf_patterns):
                grf_columns.append(col)
        
        if grf_columns:
            print(f"Auto-detected GRF columns: {grf_columns}")
        else:
            print("Warning: No GRF columns detected automatically")
            return df.copy()  # Return unchanged if no GRF columns found
    
    # Validate that specified columns exist
    missing_columns = [col for col in grf_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Specified GRF columns not found in data: {missing_columns}")
    
    # Create copy of dataframe for normalization
    df_normalized = df.copy()
    
    # Normalize GRF columns by body weight
    for col in grf_columns:
        df_normalized[col] = df[col] / body_weight
        print(f"Normalized column: {col} by body weight ({body_weight} kg)")
    
    return df_normalized

def compute_group_averages(input_directory: str, 
                          file_content_type: str, 
                          strike_type: str, 
                          speed: str, 
                          apply_mass_norm: bool = False,  # Add missing parameter
                          sub_masses: Optional[Dict[str, float]] = None,  # Add missing parameter
                          model_path: Optional[str] = None,
                          columns_to_plot: Optional[List[str]] = None,
                          grf_columns: Optional[List[str]] = None,  # Fix type annotation
                          subject_id: Optional[str] = None,
                          plot_title: Optional[str] = None,
                          plot_filename: Optional[str] = None,
                          sto_output_directory: Optional[str] = None,
                          output_directory: Optional[str] = None,
                          generate_plot: bool = True,
                          verbose: bool = True,
                          existing_fig_ax: Optional[Tuple] = None,
                          auto_color: bool = True, 
                          label: Optional[str] = None) -> Tuple[Optional[plt.Figure], Optional[plt.Axes]]:
    """
    Computes average values across groups of files that contain the same type of data.
    
    Parameters:
    -----------
    input_directory : str
        Directory containing the input files
    file_content_type : str
        Type of file content to process (e.g., "grf_move.sto")
    strike_type : str
        Strike type to filter files by
    speed : str
        Speed to filter files by
    apply_mass_norm : bool, optional
        If True, normalize GRF values by subject mass before averaging
    sub_masses : dict, optional
        Dictionary mapping subject IDs to their masses (kg)
        Required when apply_mass_norm=True
    model_path : str, optional
        Path to OpenSim model file (.osim) - used for fallback mass if sub_masses not provided
    grf_columns : list of str, optional
        Specific GRF column names to normalize. If None, auto-detects GRF columns
    subject_id : str, optional
        If specified, only files with this subject ID will be processed
    plot_title : str, optional
        Custom title for the plot
    sto_output_directory : str, optional
        Directory to save averaged .sto files
    output_directory : str, optional
        Directory to save plots
    generate_plot : bool, optional
        Whether to generate plots
    existing_fig_ax : tuple, optional
        Existing (figure, axes) tuple to add plots to
        
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects for the plot
    """
    
    def load_sto_file(filepath: str) -> pd.DataFrame:
        """Load STO file and return pandas DataFrame"""
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Find the line with 'endheader'
        endheader_idx = None
        for i, line in enumerate(lines):
            if 'endheader' in line.lower():
                endheader_idx = i
                break
        
        if endheader_idx is None:
            raise ValueError(f"No 'endheader' found in file: {filepath}")
        
        # Headers are in the line after 'endheader'
        header_line = lines[endheader_idx + 1]
        headers = [h.strip() for h in header_line.split('\t') if h.strip()]
        
        # Data starts two lines after 'endheader'
        data_start_idx = endheader_idx + 2
        
        # Parse the data
        data_dict = {header: [] for header in headers}
        
        for i in range(data_start_idx, len(lines)):
            line = lines[i].strip()
            if not line:  # Skip empty lines
                continue
            
            values = line.split('\t')
            if len(values) >= len(headers):
                for j, header in enumerate(headers):
                    try:
                        data_dict[header].append(float(values[j]))
                    except (ValueError, IndexError):
                        data_dict[header].append(np.nan)
        
        return pd.DataFrame(data_dict)
    
    def write_sto_file(df: pd.DataFrame, filename: str) -> None:
        """Write pandas DataFrame to STO file format"""
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
        
        with open(filename, 'w') as f:
            # Write header information
            f.write("OpenSimSTO version=1\n")
            f.write("name=averaged_data\n")
            f.write(f"datacolumns={len(df.columns)}\n")
            f.write(f"datarows={len(df)}\n")
            # make sure you're using percentage values of cycle as first col vals
            f.write(f"range={df.iloc[:, 0].min()} {df.iloc[:, 0].max()}\n")
            f.write("inDegrees=yes\n")
            f.write("endheader\n")
            
            # Write column headers
            f.write("\t".join(df.columns) + "\n")
            
            # Write data rows
            df.to_csv(f, sep='\t', index=False, header=False, float_format='%.8f')
    
    def normalize_to_1000_points(df: pd.DataFrame) -> pd.DataFrame:
        """Normalize data to 1000 points using interpolation"""
        if len(df) == 1000:
            return df
        
        # Create new time vector with 1000 points
        time_original = df['time'].values
        time_new = np.linspace(time_original.min(), time_original.max(), 1000)
        
        # Initialize new dataframe
        df_normalized = pd.DataFrame({'time': time_new})
        
        # Interpolate each column (except time)
        for col in df.columns:
            if col != 'time':
                # Remove NaN values for interpolation
                mask = ~np.isnan(df[col])
                if mask.sum() > 1:  # Need at least 2 points for interpolation
                    interp_func = interp1d(time_original[mask], df[col][mask], 
                                         kind='linear', bounds_error=False, 
                                         fill_value='extrapolate')
                    df_normalized[col] = interp_func(time_new)
                else:
                    # If too few valid points, fill with NaN
                    df_normalized[col] = np.full(1000, np.nan)
        
        return df_normalized
    
    def extract_file_attributes(filename: str) -> Dict[str, str]:
        """Extract attributes from filename based on the specified patterns"""
        attributes = {}
        
        # Extract subject ID (2 characters after 'sub')
        sub_match = re.search(r'sub(\w{2})', filename)
        if sub_match:
            attributes['subject'] = sub_match.group(1)
        
        # Extract strike type (2 characters after 'stri')
        stri_match = re.search(r'stri(\w{2})', filename)
        if stri_match:
            attributes['strike'] = stri_match.group(1)
        
        # Extract speed (5 characters after 'speed')
        speed_match = re.search(r'speed(\w{5})', filename)
        if speed_match:
            attributes['speed'] = speed_match.group(1)
        
        return attributes
    
        # Normalize subject_id input - handle both "03" and "sub03" formats
    if subject_id is not None:
        if subject_id.startswith('sub') and len(subject_id) > 3:
            # Extract just the ID part from "sub03" -> "03"
            subject_id = subject_id[3:5]  # Take characters 3-4 (0-indexed)
            if verbose:
                print(f"Normalized subject_id to: '{subject_id}'")
        elif len(subject_id) != 2:
            raise ValueError(f"subject_id must be 2 characters (e.g., '03') or format 'sub03', got: '{subject_id}'")
    
    # Validation for mass normalization
    if apply_mass_norm:
        if sub_masses is None and model_path is None:
            raise ValueError("apply_mass_norm=True requires either sub_masses dictionary or model_path")
        if 'grf' not in file_content_type.lower():
            print("Warning: apply_mass_norm=True but file_content_type doesn't contain 'grf'")
    
    # Find matching files with detailed debugging
    if not os.path.exists(input_directory):
        raise FileNotFoundError(f"Input directory does not exist: {input_directory}")
    
    all_files = [f for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f))]
    
    if verbose:
        print(f"\nFile matching criteria:")
        print(f"  Directory: {input_directory}")
        print(f"  Content type: '{file_content_type}'")
        print(f"  Strike type: '{strike_type}'")
        print(f"  Speed: '{speed}'")
        print(f"  Subject ID: {subject_id if subject_id else 'ALL'}")
        print(f"  Total files in directory: {len(all_files)}")
    
    # Step-by-step filtering with debugging
    matching_files = []
    content_matches = []
    attribute_debug = []
    
    for filename in all_files:
        # Step 1: Content type filter
        if file_content_type not in filename:
            continue
        content_matches.append(filename)
        
        # Step 2: Extract attributes
        attrs = extract_file_attributes(filename)
        
        # Step 3: Apply filters
        strike_match = attrs.get('strike') == strike_type
        speed_match = attrs.get('speed') == speed
        
        # Subject filter: if subject_id is None, accept all; otherwise must match
        if subject_id is None:
            subject_match = True
        else:
            subject_match = attrs.get('subject') == subject_id
        
        all_criteria_match = strike_match and speed_match and subject_match
        
        attribute_debug.append({
            'filename': filename,
            'extracted_attrs': attrs,
            'strike_match': strike_match,
            'speed_match': speed_match, 
            'subject_match': subject_match,
            'overall_match': all_criteria_match
        })
        
        if all_criteria_match:
            filepath = os.path.join(input_directory, filename)
            matching_files.append(filepath)
    
    if verbose:
        print(f"\nStep 1 - Files containing '{file_content_type}': {len(content_matches)}")
        
        if len(content_matches) > 0 and len(content_matches) <= 20:
            print("Files with matching content type:")
            for debug_info in attribute_debug:
                filename = debug_info['filename']
                attrs = debug_info['extracted_attrs']
                matches = debug_info
                print(f"  {filename}")
                print(f"    Extracted: subject='{attrs.get('subject')}', strike='{attrs.get('strike')}', speed='{attrs.get('speed')}'")
                print(f"    Matches: strike={matches['strike_match']}, speed={matches['speed_match']}, subject={matches['subject_match']}")
                print(f"    Overall: {matches['overall_match']}")
        elif len(content_matches) > 20:
            print(f"  Too many files to show details ({len(content_matches)} files)")
            
            # Show summary of available attributes
            subjects_found = set()
            strikes_found = set()
            speeds_found = set()
            
            for debug_info in attribute_debug:
                attrs = debug_info['extracted_attrs']
                if attrs.get('subject'):
                    subjects_found.add(attrs['subject'])
                if attrs.get('strike'):
                    strikes_found.add(attrs['strike'])
                if attrs.get('speed'):
                    speeds_found.add(attrs['speed'])
            
            print(f"  Available subjects: {sorted(subjects_found)}")
            print(f"  Available strikes: {sorted(strikes_found)}")
            print(f"  Available speeds: {sorted(speeds_found)}")
    
    if not matching_files:
        error_msg = f"No matching files found for strike={strike_type}, speed={speed}, content={file_content_type}"
        if subject_id:
            error_msg += f", subject={subject_id}"
        
        if verbose and len(content_matches) > 0:
            print(f"\nTroubleshooting: Found {len(content_matches)} files with correct content type but wrong attributes")
            
            # Show what values are actually available
            subjects_found = set()
            strikes_found = set()
            speeds_found = set()
            
            for debug_info in attribute_debug:
                attrs = debug_info['extracted_attrs']
                if attrs.get('subject'):
                    subjects_found.add(attrs['subject'])
                if attrs.get('strike'):
                    strikes_found.add(attrs['strike'])
                if attrs.get('speed'):
                    speeds_found.add(attrs['speed'])
            
            print(f"Available combinations in {file_content_type} files:")
            print(f"  Subjects: {sorted(subjects_found)}")
            print(f"  Strike types: {sorted(strikes_found)}")
            print(f"  Speeds: {sorted(speeds_found)}")
        
        raise ValueError(error_msg)
    
    if verbose:
        print(f"\nFinal result: {len(matching_files)} files match all criteria:")
        for f in matching_files:
            print(f"  {os.path.basename(f)}")

    
    # Load and process files
    normalized_dfs = []
    file_durations = []
    excluded_files = []
    min_duration = 0.65
    max_duration = 0.9

    for filepath in matching_files:
        try:
            df = load_sto_file(filepath)

            # Now normalize to 1000 time points
            # df_norm = normalize_to_1000_points(df)
            # normalized_dfs.append(df_norm)
            
            # Get this file's duration
            file_duration = df['time'].iloc[-1] - df['time'].iloc[0]
            # file_durations.append(file_duration)
            
            if file_duration < min_duration or file_duration > max_duration:
                excluded_files.append({
                    'filename': os.path.basename(filepath),
                    'duration': file_duration,
                    'reason': f"Duration: {file_duration} outside range"
                })
                continue

            if verbose:
                print(f"Valid duration: {os.path.basename(filepath)} = {file_duration:.3f}s")
            
            # Apply mass normalization BEFORE time normalization if requested
            if apply_mass_norm and file_content_type == 'grf_move.sto':
            
                # Extract subject ID from filename for mass lookup
                filename = os.path.basename(filepath)
                attrs = extract_file_attributes(filename)
                file_subject_id = attrs.get('subject')
                
                if file_subject_id and sub_masses and file_subject_id in sub_masses:
                    # Use subject-specific mass from dictionary
                    subject_mass = sub_masses[file_subject_id]
                    if verbose:
                        print(f"Using subject {file_subject_id} mass: {subject_mass} kg for {filename}")
                    
                    # Apply mass normalization to this file
                    df = mass_normalize_grf_by_mass(df, subject_mass, grf_columns)
                    
                elif model_path and os.path.exists(model_path):
                    # Fallback to model-based mass
                    if verbose:
                        print(f"Subject mass not found in sub_masses, using model mass for {filename}")
                    df = mass_normalize_grf(df, model_path, grf_columns)
                    
                else:
                    print(f"Warning: Cannot apply mass normalization for {filename} - no mass available")

            df_norm = normalize_to_1000_points(df)
            normalized_dfs.append(df_norm)
            file_durations.append(file_duration)

            if verbose:
                print(f"Processed: {os.path.basename(filepath)} (duration: {file_duration:.3f}s)")


            
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            excluded_files.append({
                'filename': os.path.basename(filepath),
                'duration': file_duration,
                'reason': f"Duration: {file_duration} outside range"
            })
            continue
    
    if not normalized_dfs:
        raise ValueError("No files were successfully processed")
    
    # Store average duration from all files being processed in function call
    average_duration = np.mean(file_durations)
    print(f"\nFile range duration: {average_duration:.3f}s (range: {min(file_durations):.3f}s - {max(file_durations):.3f}s)")
    print(f"Number of analyzed files: {len(file_durations)}")
    # Guarantee the x-axis/time values have 1000 points
    time_points = np.linspace(0, average_duration, 1000)
    
    # Compute averages and standard deviations
    all_columns = set()
    for df in normalized_dfs:
        all_columns.update(df.columns)
    all_columns = sorted(list(all_columns))
    
    # Initialize result arrays (for cycle vals, 1000 points from value 0 until 100)
    gait_cycle_percentage = np.linspace(0, 100, 1000)
    
    averages = {}
    std_devs = {}
    
    for col in all_columns:
        if col == 'time':
            # Replace original time with calculated avg duration, with cycle percentage included
            averages['time'] = time_points
            std_devs['time'] = np.zeros_like(time_points)
            averages['gait_cycle_percent'] = gait_cycle_percentage
            std_devs['gait_cycle_percent'] = np.zeros_like(gait_cycle_percentage)
        else:
            # Collect data for this column from all files
            col_data = []
            for df in normalized_dfs:
                if col in df.columns:
                    col_data.append(df[col].values)
                else:
                    col_data.append(np.full(1000, np.nan))
            
            # Convert to array (n_files x 1000)
            col_data = np.array(col_data)
            
            # Compute mean and std along axis 0 (across files)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                averages[col] = np.nanmean(col_data, axis=0)
                std_devs[col] = np.nanstd(col_data, axis=0, ddof=1)
    
    column_order = ['time', 'gait_cycle_percent']
    remain_columns = [col for col in sorted(all_columns) if col != 'time']
    column_order.extend(remain_columns)

    # Create average DataFrame
    avg_df = pd.DataFrame()
    # avg_df['time'] = time_points
    # avg_df['gait_cycle_percent'] = gait_cycle_percentage
    for col in column_order:
        if col in averages: # Skips over original time col
            avg_df[col] = averages[col]
    
    if subject_id is not None:
        subject_prefix = f"sub{subject_id}"
    else:
        subject_prefix = "sub00"

    if output_directory is None:
        plot_output_dir = input_directory
    else:
        plot_output_dir = output_directory
        os.makedirs(plot_output_dir, exist_ok=True)

    if sto_output_directory is None:
        sto_output_dir = plot_output_dir
    else:
        sto_output_dir = sto_output_directory
        os.makedirs(sto_output_dir,exist_ok=True)

    # # Save average .sto file
    # output_filename = f"{subject_prefix}_stri{strike_type}_cond00000_speed{speed}_{file_content_type.replace('.sto', '')}_averages.sto"
    # output_path = os.path.join(sto_output_dir, output_filename)
    # write_sto_file(avg_df, output_path)
    # print(f"Average data saved to: {output_path}")

    if sto_output_dir != input_directory: # copy to direct to additional path
        # The files written here are for moco inputs and therefore don't need percent values in them
        avg_df_input = avg_df.drop(columns=['gait_cycle_percent'], errors='ignore')

        output_filename_no_percent = f"{subject_prefix}_stri{strike_type}_cond00000_speed{speed}_{file_content_type.replace('.sto', '')}_averages.sto"
        output_path_no_percent = os.path.join(sto_output_dir, output_filename_no_percent)
        write_sto_file(avg_df_input, output_path_no_percent)
        print(f"Average data (without gait_cycle_percent) saved to: {output_path_no_percent}")
    else:
        print("sto_output_directory is same as input_directory - single file saved with complete data")
    

    # if apply_mass_norm:
    #     for filepath in matching_files:
    #         df = load_sto_file(filepath)
    #         for mass in sub_masses.values():
    #             try:
    #                 # body_mass = sub_masses[mass]
    #                 avg_df_normalized = normalize_grf_by_manual_body_weight(df, mass, grf_columns)

    #                 # Save the data frame to a file
    #                 normalized_filename = f"{subject_prefix}_stri{strike_type}_cond00000_speed{speed}_averages_normalized.sto"
    #                 normalized_path = os.path.join(sto_output_dir, normalized_filename)
    #                 write_sto_file(avg_df_normalized, normalized_path)
    #                 print(f"Mass-normalized average data saved to: {normalized_path}")

    #                 # Save normalized averages file to sto_output_directory (without gait_cycle_percent)
    #                 if sto_output_dir != input_directory:  # Only create second copy if directories are different
    #                     avg_df_normalized_no_gait = avg_df_normalized.drop(columns=['gait_cycle_percent'], errors='ignore')
                        
    #                     normalized_filename_no_gait = f"{subject_prefix}_stri{strike_type}_cond00000_speed{speed}_averages_normalized.sto"
    #                     normalized_path_no_gait = os.path.join(sto_output_dir, normalized_filename_no_gait)
    #                     write_sto_file(avg_df_normalized_no_gait, normalized_path_no_gait)
    #                     print(f"Mass-normalized average data (without gait_cycle_percent) saved to: {normalized_path_no_gait}")
    #                 else:
    #                     print("sto_output_directory is same as input_directory - single normalized file saved with complete data")
                    
    #             except Exception as e:
    #                 print(f"Error during mass normalization: {e}")
    #                 print("Continuing with non-normalized data only...")

    if generate_plot:
        avg_df_input = avg_df
        
        output_filename_plot = f"{subject_prefix}_stri{strike_type}_cond00000_speed{speed}_{file_content_type.replace('.sto', '')}_averages.sto"
        output_filename_plot_path = os.path.join(sto_output_dir, output_filename_plot)
        write_sto_file(avg_df_input, output_filename_plot_path)
        print(f"Average data saved to without cycle percents: {output_filename_plot_path}")

    else:
        avg_df_input = avg_df.drop(columns=['gait_cycle_percent'], errors='ignore')
        print("Plot generation disabled")
        return None,None

    # Create plot
    if existing_fig_ax is not None:
        fig, ax = existing_fig_ax
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
    
    # Specify columns to plot
    if columns_to_plot is not None:
        plot_columns = []
        for col in columns_to_plot:
            if col in averages and col != 'time' and col != 'gait_cycle_percent':
                if not np.all(np.isnan(averages[col])):
                    plot_columns.append(col)
                else:
                    print(f"Warning: column '{col}' contains only NaN vals, skipping")
            elif col not in averages:
                print(f"Warning: Column '{col}' not found in data, skipping")
    else:
        # Plots all cols except for time and gait_cycle_percent
        plot_columns = [col for col in all_columns
                        if col not in ['time', 'gait_cycle_percent'] and 
                        col in averages and not np.all(np.isnan(averages[col]))]
        
    if not plot_columns:
        print("Warning: No valid columns to plot")
        return None, None if not generate_plot else (fig, ax)
    print(f'Plotting columns: {plot_columns}')

    # Generate x-axis vals as percentage of cycle
    x_axis_values = gait_cycle_percentage / 1.0

    if auto_color and existing_fig_ax is not None:
        # Get existing line colors from the plot
        existing_colors = []
        for line in ax.get_lines():
            existing_colors.append(line.get_color())
        
        # Convert existing colors to RGB tuples for comparison
        # import matplotlib.colors as mcolors
        existing_rgb = []
        for color in existing_colors:
            try:
                rgb = mcolors.to_rgb(color)
                existing_rgb.append(rgb)
            except:
                pass
        
        # Define color palette to choose from
        color_palette = [
            '#d62728',  # Red
            '#ff7f0e',  # Orange
            '#2ca02c',  # Green
            '#1f77b4',  # Blue
            '#9467bd',  # Purple
            '#8c564b',  # Brown
        ]
        
        # Find colors that haven't been used yet
        available_colors = []
        for color in color_palette:
            color_rgb = mcolors.to_rgb(color)
            # Check if this color is sufficiently different from existing colors
            is_different = True
            for existing in existing_rgb:
                # Calculate color distance
                distance = sum((c1 - c2)**2 for c1, c2 in zip(color_rgb, existing))
                if distance < 0.1:  # Threshold for "too similar"
                    is_different = False
                    break
            if is_different:
                available_colors.append(color)
        
        # If we've used all colors, cycle through them again
        if not available_colors:
            available_colors = color_palette
        
        # Assign colors from available pool
        colors = [available_colors[i % len(available_colors)] for i in range(len(plot_columns))]
        
        print(f"Auto-assigned colors (existing plot detected): {colors}")
    
    else:
        # Defaults to original color map
        colors = plt.cm.Set3(np.linspace(0, 1, len(plot_columns)))
    

    # Plot the selected data
    for i, col in enumerate(plot_columns):
        # if col in averages and not np.all(np.isnan(averages[col])):
        mean_vals = averages[col]
        std_vals = std_devs[col]
        
        if label is not None:
            # Uses custom label if provided
            # If plotting multiple columns at once, append column name
            if len(plot_columns) > 1:
                mean_label = f"{label} - {col}"
                std_label = f"{label} - {col} (+/-2σ)"
            else:
                # Single column instance
                mean_label = label
                std_label = f"{label} (±2σ)"
        else:
            # Default labelling
            mean_label = f"{col} (mean)"
            std_label = f"{label} (±2σ)"

        # Plot mean line
        ax.plot(x_axis_values, mean_vals, linestyle='--', color=colors[i], 
                linewidth=2, label=mean_label)
        
        # Plot standard deviation bands (mean ± 1 std)
        upper_bound = mean_vals + (std_vals)
        lower_bound = mean_vals - (std_vals)
        
        ax.fill_between(x_axis_values, lower_bound, upper_bound, 
                        color=colors[i], alpha=0.3, 
                        # label=f"{col} (+/-2σ)"
                        )
    
    # Customize plot
    ax.set_xlabel('cycle %')
    ax.set_ylabel('Value')
    
    if plot_title:
        ax.set_title(plot_title)
    else:
        title_parts = [f"Strike: {strike_type}", f"Speed: {speed}"]
        # if subject_id:
        #     title_parts.insert(0, f"Subject: {subject_id}")
        # title_parts.append(f"Content: {file_content_type}")
        # ax.set_title("Group Averages - " + ", ".join(title_parts))
    
    # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    if plot_filename is not None:
        plot_filename = plot_filename
    else:
        plot_filename = f"{subject_prefix}_stri{strike_type}_cond00000_speed{speed}_{file_content_type.replace('.sto', '')}_plot.png"
    
    plot_path = os.path.join(plot_output_dir, plot_filename)
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_path}")
    
    return fig, ax

def average_moco_solution_columns(input_directory: str,
                                   file_suffix: str,
                                   columns_to_average: List[str],
                                   strike_type: Optional[str] = None, 
                                   speed: Optional[str] = None,
                                   output_directory: Optional[str] = None,
                                   output_filename: Optional[str] = None,
                                   generate_plot: bool = True,
                                   plot_title: Optional[str] = None,
                                   existing_fig_ax: Optional[Tuple] = None,
                                   color_list: Optional[List] = None,
                                   verbose: bool = True, 
                                   label: Optional[str] = None) -> Tuple[Optional[plt.Figure], Optional[plt.Axes]]:
    """
    Average specified columns from multiple MoCo solution .sto files after time-normalizing to 1000 points.
    
    This function:
    1. Finds all .sto files in the input directory matching the specified suffix
    2. Time-normalizes each file to 1000 data points
    3. Averages the specified columns across all files
    4. Optionally plots the results (with option to superimpose on existing figure)
    5. Saves the averaged data to a .sto file
    
    Parameters:
    -----------
    input_directory : str
        Directory containing the MoCo solution .sto files
    file_suffix : str
        Common suffix to identify files to process (e.g., "solution.sto", "states.sto")
    columns_to_average : list of str
        List of column names to average and optionally plot
    output_directory : str, optional
        Directory to save output files. If None, uses input_directory
    output_filename : str, optional
        Custom filename for output .sto file. If None, generates from suffix
    generate_plot : bool, optional (default=True)
        Whether to generate plots of the averaged data
    plot_title : str, optional
        Custom title for the plot. If None, generates automatic title
    existing_fig_ax : tuple, optional
        Existing (figure, axes) tuple to add plots to for superimposing
    verbose : bool, optional (default=True)
        Whether to print detailed processing information
        
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects for the plot, or (None, None) if generate_plot=False
        
    Example:
    --------
    # First call - creates new figure
    fig, ax = average_moco_solution_columns(
        input_directory="/path/to/solutions/",
        file_suffix="solution.sto",
        columns_to_average=['hip_flexion_r', 'knee_angle_r', 'ankle_angle_r'],
        generate_plot=True
    )
    
    # Second call - superimpose on existing figure
    fig, ax = average_moco_solution_columns(
        input_directory="/path/to/solutions/",
        file_suffix="solution.sto",
        columns_to_average=['hip_flexion_l', 'knee_angle_l', 'ankle_angle_l'],
        existing_fig_ax=(fig, ax)
    )
    """
    
    # Helper function to load .sto file
    def load_sto_file(filepath: str) -> pd.DataFrame:
        """Load STO file and return pandas DataFrame"""
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Find the line with 'endheader'
        endheader_idx = None
        for i, line in enumerate(lines):
            if 'endheader' in line.lower():
                endheader_idx = i
                break
        
        if endheader_idx is None:
            raise ValueError(f"No 'endheader' found in file: {filepath}")
        
        # Headers are in the line after 'endheader'
        header_line = lines[endheader_idx + 1]
        headers = [h.strip() for h in header_line.split('\t') if h.strip()]

        headers = [h.replace('|', '/') for h in headers]
        
        # Data starts two lines after 'endheader'
        data_start_idx = endheader_idx + 2
        
        # Parse the data
        data_dict = {header: [] for header in headers}
        
        for i in range(data_start_idx, len(lines)):
            line = lines[i].strip()
            if not line:  # Skip empty lines
                continue
            
            values = line.split('\t')
            if len(values) >= len(headers):
                for j, header in enumerate(headers):
                    try:
                        data_dict[header].append(float(values[j]))
                    except (ValueError, IndexError):
                        data_dict[header].append(np.nan)
        
        return pd.DataFrame(data_dict)
    
    def write_sto_file(df: pd.DataFrame, filename: str) -> None:
        """Write pandas DataFrame to STO file format"""
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
        
        with open(filename, 'w') as f:
            # Write header information
            f.write("OpenSimSTO version=1\n")
            f.write("name=averaged_moco_solution\n")
            f.write(f"datacolumns={len(df.columns)}\n")
            f.write(f"datarows={len(df)}\n")
            f.write(f"range={df.iloc[:, 0].min()} {df.iloc[:, 0].max()}\n")
            f.write("inDegrees=yes\n")
            f.write("endheader\n")
            
            # Write column headers
            f.write("\t".join(df.columns) + "\n")
            
            # Write data rows
            df.to_csv(f, sep='\t', index=False, header=False, float_format='%.8f')
    
    def normalize_to_1000_points(df: pd.DataFrame) -> pd.DataFrame:
        """Normalize data to 1000 points using interpolation"""
        if len(df) == 1000:
            return df
        
        # Create new time vector with 1000 points
        time_original = df['time'].values
        time_new = np.linspace(time_original.min(), time_original.max(), 1000)
        
        # Initialize new dataframe
        df_normalized = pd.DataFrame({'time': time_new})
        
        # Interpolate each column (except time)
        for col in df.columns:
            if col != 'time':
                # Remove NaN values for interpolation
                mask = ~np.isnan(df[col])
                if mask.sum() > 1:  # Need at least 2 points for interpolation
                    interp_func = interp1d(time_original[mask], df[col][mask], 
                                         kind='linear', bounds_error=False, 
                                         fill_value='extrapolate')
                    df_normalized[col] = interp_func(time_new)
                else:
                    # If too few valid points, fill with NaN
                    df_normalized[col] = np.full(1000, np.nan)
        
        return df_normalized
    
    all_files = [f for f in os.listdir(input_directory) if os.path.isfile(os.path.join(input_directory, f))]
    
    # Filter by suffix
    matching_files = [f for f in all_files if f.endswith(file_suffix)]
    
    # Filter by strike_type if provided
    if strike_type:
        matching_files = [f for f in matching_files if strike_type in f]

    # Filter by speed if provided
    if speed:
        matching_files = [f for f in matching_files if speed in f]
        
    # Create full paths for the final list of files
    matching_files_paths = [os.path.join(input_directory, f) for f in matching_files]

    if not matching_files_paths:
        filter_msg = [f"suffix='{file_suffix}'"]
        if strike_type: filter_msg.append(f"strike_type='{strike_type}'")
        if speed: filter_msg.append(f"speed='{speed}'")
        raise ValueError(f"No files found with filters ({', '.join(filter_msg)}) in {input_directory}")

    if verbose:
        print(f"\nFound {len(matching_files_paths)} files matching the criteria:")
        for f in matching_files_paths:
            print(f"  {os.path.basename(f)}")


    # Load and process files
    normalized_dfs = []
    file_durations = []
    
    for filepath in matching_files_paths:
        try:
            df = load_sto_file(filepath)
            
            # Get file duration
            file_duration = df['time'].iloc[-1] - df['time'].iloc[0]
            file_durations.append(file_duration)
            
            # Normalize to 1000 points
            df_norm = normalize_to_1000_points(df)
            normalized_dfs.append(df_norm)
            
            if verbose:
                print(f"Processed: {os.path.basename(filepath)} (duration: {file_duration:.3f}s)")
                
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            continue
    
    if not normalized_dfs:
        raise ValueError("No files were successfully processed")
    
    # Calculate average duration
    average_duration = np.mean(file_durations)
    if verbose:
        print(f"\nAverage duration: {average_duration:.3f}s (range: {min(file_durations):.3f}s - {max(file_durations):.3f}s)")
        print(f"Number of processed files: {len(normalized_dfs)}")
    
    # Create time vector for averaged data
    time_points = np.linspace(0, average_duration, 1000)
    
    # Collect all available columns from all files
    all_columns = set()
    for df in normalized_dfs:
        all_columns.update(df.columns)
    all_columns = sorted(list(all_columns))
    
    # Validate requested columns exist
    missing_columns = [col for col in columns_to_average if col not in all_columns]
    if missing_columns:
        print(f"Warning: Requested columns not found in data: {missing_columns}")
        columns_to_average = [col for col in columns_to_average if col in all_columns]
        if not columns_to_average:
            raise ValueError("None of the requested columns exist in the data")
    
    # Compute averages and standard deviations
    averages = {}
    std_devs = {}
    
    for col in all_columns:
        if col == 'time':
            averages['time'] = time_points
            std_devs['time'] = np.zeros_like(time_points)
        else:
            # Collect data for this column from all files
            col_data = []
            for df in normalized_dfs:
                if col in df.columns:
                    col_data.append(df[col].values)
                else:
                    col_data.append(np.full(1000, np.nan))
            
            # Convert to array (n_files x 1000)
            col_data = np.array(col_data)
            
            # Compute mean and std along axis 0 (across files)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                averages[col] = np.nanmean(col_data, axis=0)
                std_devs[col] = np.nanstd(col_data, axis=0, ddof=1)
    
    # Create average DataFrame
    column_order = ['time'] + [col for col in sorted(all_columns) if col != 'time']
    avg_df = pd.DataFrame()
    
    for col in column_order:
        if col in averages:
            avg_df[col] = averages[col]
    
    # Set output directory
    if output_directory is not None:
        os.makedirs(output_directory, exist_ok=True)
    
    # Generate output filename
    if output_filename is None:
        name_parts = ["averaged"]
        if strike_type: name_parts.append(strike_type)
        if speed: name_parts.append(speed)
        name_parts.append(file_suffix.lstrip('_'))
        output_filename = "_".join(name_parts)
    
    output_path = os.path.join(output_directory, output_filename)
    write_sto_file(avg_df, output_path)
    
    if verbose:
        print(f"\nAveraged data saved to: {output_path}")
        
    if not generate_plot:
        if verbose: print("Plot generation disabled")
        return None, None
    
    if existing_fig_ax is not None:
        fig, ax = existing_fig_ax
        if verbose: print("Superimposing plot on existing figure")
    else:
        fig, ax = plt.subplots(figsize=(12, 8))
        
    plot_columns = [col for col in columns_to_average if col in averages and col != 'time' and not np.all(np.isnan(averages[col]))]
    
    # Filter columns to plot (only those that exist and have valid data)
    plot_columns = []
    for col in columns_to_average:
        if col in averages and col != 'time':
            if not np.all(np.isnan(averages[col])):
                plot_columns.append(col)
            else:
                print(f"Warning: column '{col}' contains only NaN values, skipping")
    
    if not plot_columns:
        print("Warning: No valid columns to plot")
        return fig, ax
    
    if verbose:
        print(f'Plotting columns: {plot_columns}')

    # Extract and format subject IDs for the legend
    subject_ids = sorted(list(set([os.path.basename(f)[:5] for f in matching_files])))
    
    if not subject_ids:
        id_string = "Subjects"  # Fallback if no IDs are found
    elif len(subject_ids) <= 3:
        id_string = ", ".join(subject_ids)
    else:
        # For many subjects, create a range like "S0001 to S0005"
        id_string = f"{subject_ids[0]} to {subject_ids[-1]}"
    
    # Helper function to convert colors to RGB for comparison
    def color_to_rgb(color):
        """Convert matplotlib color to RGB tuple."""
        if isinstance(color, str):
            return mcolors.to_rgb(color)
        elif isinstance(color, (list, tuple)) and len(color) >= 3:
            return tuple(color[:3])
        return color
    
    # Get colors already used in the plot (if superimposing)
    used_colors = []
    if existing_fig_ax is not None:
        for line in ax.get_lines():
            used_colors.append(color_to_rgb(line.get_color()))
        for collection in ax.collections:
            if hasattr(collection, 'get_facecolor'):
                fc = collection.get_facecolor()
                if len(fc) > 0:
                    used_colors.append(tuple(fc[0][:3]))
    
    # Define default color list if not provided
    if color_list is None:
        color_list = [
            'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
            'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan',
            'darkblue', 'darkgreen', 'darkred', 'darkorange', 'darkviolet',
            'steelblue', 'forestgreen', 'firebrick', 'goldenrod', 'mediumorchid'
        ]
    
    # Select colors that are not already used
    def colors_are_similar(c1, c2, threshold=0.3):
        """Check if two RGB colors are similar."""
        c1_rgb = color_to_rgb(c1)
        c2_rgb = color_to_rgb(c2)
        return sum((a - b) ** 2 for a, b in zip(c1_rgb, c2_rgb)) ** 0.5 < threshold
    
    available_colors = []
    for candidate_color in color_list:
        candidate_rgb = color_to_rgb(candidate_color)
        is_too_similar = any(colors_are_similar(candidate_rgb, used_rgb) for used_rgb in used_colors)
        if not is_too_similar:
            available_colors.append(candidate_color)
    
    if len(available_colors) < len(plot_columns):
        import colorsys
        for i in range(len(plot_columns) - len(available_colors)):
            hue = (i * 0.618033988749895) % 1.0
            rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
            available_colors.append(rgb)
    
    # Assign colors to columns
    colors = [available_colors[i % len(available_colors)] for i in range(len(plot_columns))]
    
    # Plot each column
    for i, col in enumerate(plot_columns):
        mean_vals = averages[col]
        std_vals = std_devs[col]
        
        # Construct the new label for the legend
        # Start with the ID string and append the optional input label
        base_label = f"{id_string} - {label}" if label else id_string
        
        # If plotting multiple columns, add the column name to differentiate them
        if len(plot_columns) > 1:
            mean_label = f"{base_label} - {col}"
        else:
            mean_label = base_label

        # Plot mean line with the new label
        ax.plot(time_points, mean_vals, color=colors[i], 
                linewidth=2, label=mean_label)
        
        # Plot standard deviation bands (mean ± 1 std)
        upper_bound = mean_vals + std_vals
        lower_bound = mean_vals - std_vals
        
        # The shaded area is not labeled to keep the legend clean
        # ax.fill_between(time_points, lower_bound, upper_bound, 
        #                 color=colors[i], alpha=0.3)
    
    # Customize plot
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Value')
    
    if plot_title:
        ax.set_title(plot_title)
    else:
        ax.set_title(f"Averaged MoCo Solutions - {file_suffix}")
    
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_filename = output_filename.replace('.sto', '_plot.png')
    plot_path = os.path.join(output_directory, plot_filename)
    fig.savefig(plot_path, dpi=300, bbox_inches='tight')
    
    if verbose:
        print(f"Plot saved to: {plot_path}")
    
    return fig, ax

def plot_joint_fiber_forces(model, solution, coord_paths=None,
                           file_name=None, output_dir=None, 
                           min_moment_arm_threshold=0.00001,
                           normalize_time=False):
    """
    Plot summed active and passive fiber forces for muscles spanning specified joints.
    
    Creates subplots showing the sum of active fiber forces and sum of passive fiber forces
    for all muscles that cross each specified coordinate (joint).
    
    Parameters
    ----------
    model : opensim.Model
        OpenSim model (must be initialized)
    solution : osim.MocoTrajectory, osim.StatesTrajectory, or str
        Moco solution trajectory object OR path to .sto file
        If string path provided, will load as TimeSeriesTable and convert to StatesTrajectory
    coord_paths : list of str
        List of coordinate paths to analyze (e.g., ['/jointset/ankle_r/ankle_angle_r'])
    file_name : str, optional
        Base name for the output file
    output_dir : str, optional
        Directory to save the figure
    min_moment_arm_threshold : float, optional
        Minimum moment arm magnitude to consider a muscle as spanning the joint (default: 0.00001 m)
    normalize_time : bool, optional
        If True, normalize time to 0-100% (default: False)
    
    Returns
    -------
    matplotlib.figure.Figure
        The created figure object
        
    Examples
    --------
    >>> model_proc = osim.ModelProcessor('model.osim')
    >>> model = model_proc.process()
    >>> model.initSystem()
    >>> 
    >>> # Option 1: Pass file path (will be loaded automatically)
    >>> fig = plot_joint_fiber_forces(model, 'solution.sto', 
    ...                                coord_paths=['/jointset/ankle_r/ankle_angle_r'])
    >>> 
    >>> # Option 2: Pass trajectory object
    >>> table = osim.TimeSeriesTable('solution.sto')
    >>> solution = osim.StatesTrajectory.createFromStatesTable(model, table)
    >>> fig = plot_joint_fiber_forces(model, solution,
    ...                                coord_paths=['/jointset/ankle_r/ankle_angle_r'])
    """
    
    if coord_paths is None:
        raise ValueError("coord_paths must be provided")
    
    # Handle solution input - convert string path to trajectory if needed
    if isinstance(solution, str):
        # Load from file path
        solution_table = osim.TimeSeriesTable(solution)
        solution, _ = osim.StatesTrajectory.createFromStatesTable(model, solution_table)
    
    # Get muscle paths from model
    muscle_paths = []
    for muscle in model.getMuscleList():
        muscle_paths.append(muscle.getAbsolutePathString())
    
    num_coords = len(coord_paths)
    num_muscles = len(muscle_paths)
    
    # Get time vector
    time = np.array(solution.getTimeMat()).flatten()
    states_trajectory = solution.exportToStatesTrajectory(model)
    
    # Calculate active and passive fiber forces for all muscles
    active_fiber_forces = np.empty((len(time), num_muscles))
    passive_fiber_forces = np.empty((len(time), num_muscles))
    
    for imusc, muscle_path in enumerate(muscle_paths):
        muscle = model.getComponent(muscle_path)
        for itime in range(len(time)):
            state = states_trajectory.get(itime)
            model.realizeDynamics(state)
            active_fiber_forces[itime, imusc] = muscle.getActiveFiberForce(state)
            passive_fiber_forces[itime, imusc] = muscle.getPassiveFiberForce(state)
    
    # Create figure with subplots
    fig, axes = plt.subplots(num_coords, 1, figsize=(10, 3 * num_coords), squeeze=False)
    axes = axes.flatten()
    
    # Process each coordinate
    for icoord, coord_path in enumerate(coord_paths):
        coord = model.getComponent(coord_path)
        coord_name = os.path.split(coord_path)[-1]
        
        # Calculate moment arms for all muscles at this coordinate
        muscle_moment_arms = np.empty((len(time), num_muscles))
        for imusc, muscle_path in enumerate(muscle_paths):
            muscle = model.getComponent(muscle_path)
            for itime in range(len(time)):
                state = states_trajectory.get(itime)
                muscle_moment_arms[itime, imusc] = muscle.computeMomentArm(state, coord)
        
        # Initialize summed forces
        sum_active_forces = np.zeros_like(time)
        sum_passive_forces = np.zeros_like(time)
        
        muscles_included = []
        
        # Sum forces for muscles that span this joint
        for imusc, muscle_path in enumerate(muscle_paths):
            # Check if muscle has non-negligible moment arm at any point
            if np.any(np.abs(muscle_moment_arms[:, imusc]) > min_moment_arm_threshold):
                sum_active_forces += active_fiber_forces[:, imusc]
                sum_passive_forces += passive_fiber_forces[:, imusc]
                muscle_name = os.path.split(muscle_path)[-1]
                muscles_included.append(muscle_name)
        
        # Prepare time axis
        if normalize_time:
            time_axis = np.linspace(0, 100, len(time))
            xlabel = 'Time (% gait cycle)'
        else:
            time_axis = time
            xlabel = 'Time (s)'
        
        # Plot on subplot
        ax = axes[icoord]
        ax.plot(time_axis, sum_active_forces, label='Sum Active Fiber Forces', 
                color='red', linewidth=2)
        ax.plot(time_axis, sum_passive_forces, label='Sum Passive Fiber Forces', 
                color='blue', linewidth=2)
        ax.plot(time_axis, sum_active_forces + sum_passive_forces, 
                label='Total Fiber Forces', color='black', linewidth=2, linestyle='--')
        
        ax.set_title(f'{coord_name}\n({len(muscles_included)} muscles spanning joint)', 
                    fontsize=12, fontweight='bold')
        ax.set_ylabel('Force (N)', fontsize=10)
        # ax.legend(frameon=False, loc='best', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis='both', labelsize=9)
        
        # Only set xlabel on bottom subplot
        if icoord == num_coords - 1:
            ax.set_xlabel(xlabel, fontsize=10)
        
        # Print muscle list for this coordinate
        if muscles_included:
            print(f"\n{coord_name} - Muscles included ({len(muscles_included)}):")
            for muscle_name in sorted(muscles_included):
                print(f"  - {muscle_name}")
    
    fig.tight_layout()
    
    # Save figure if output directory provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        if file_name is not None:
            filename = f'joint_fiber_forces_{file_name}.png'
        else:
            filename = 'joint_fiber_forces.png'
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f'\nSaved figure to {filepath}')
    
    return fig

def plot_muscle_mechanics_bars(csv_directory,
                               suffix='_muscle_mechanics.csv',
                               speed_strike_pair=None,
                               column_to_plot=None,
                               subject_ids=None,
                               plot_color='blue',
                               existing_fig_ax=None,
                               figsize=(12, 6),
                               ylabel='Peak Value',
                               title=None,
                               output_file=None):
    """
    Create bar plots of peak values from muscle mechanics CSV files.
    
    Each bar represents the average peak value across subjects for a specific condition,
    with error bars showing standard error of the mean (SEM).
    
    Parameters:
    -----------
    csv_directory : str
        Directory containing the CSV files
    suffix : str, optional (default='_muscle_mechanics.csv')
        Common suffix for the CSV files to process
    column_to_plot : str, optional
        Column name to extract and plot (e.g., '/forceset/soleus_r/activation')
    subject_ids : list of str, optional
        List of subject IDs to include (e.g., ['03', '04']). If None, includes all found
    plot_color : str, optional (default='blue')
        Base color for the bars. Will use shades of this color for different subjects
    existing_fig_ax : tuple, optional
        Existing (figure, axes) to add bars to. If None, creates new figure
    figsize : tuple, optional (default=(12, 6))
        Figure size if creating new figure
    ylabel : str, optional (default='Peak Value')
        Label for y-axis
    title : str, optional
        Plot title
    output_file : str, optional
        Path to save the plot
        
    Returns:
    --------
    tuple
        (figure, axes) matplotlib objects
        
    Examples:
    ---------
    # Create new plot
    fig, ax = plot_muscle_mechanics_bars(
        csv_directory='path/to/csvs/',
        column_to_plot='/forceset/soleus_r/activation',
        plot_color='blue',
        title='Soleus Activation'
    )
    
    # Add another muscle to same plot
    fig, ax = plot_muscle_mechanics_bars(
        csv_directory='path/to/csvs/',
        column_to_plot='/forceset/gaslat_r/activation',
        plot_color='red',
        existing_fig_ax=(fig, ax)
    )
    """
    
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from scipy import stats
    import os
    
    
    # Find all CSV files with the specified suffix
    csv_files = []
    if os.path.exists(csv_directory):
        for filename in os.listdir(csv_directory):
            if filename.endswith(suffix):
                csv_files.append(os.path.join(csv_directory, filename))
    
    if not csv_files:
        raise ValueError(f"No CSV files found with suffix '{suffix}' in {csv_directory}")
    
    print(f"Found {len(csv_files)} CSV files")
    
    # convert single column into list
    if isinstance(column_to_plot, str):
        column_to_plot = [column_to_plot]
    elif column_to_plot is None:
        raise ValueError("must specify column_to_plot")
    # Parse file names and organize by condition
    conditions = {}  # Key: (strike_type, speed), Value: list of (subject_id, filepath, data)
    
    for filepath in csv_files:
        filename = os.path.basename(filepath)
        
        # Extract subject ID (assumes format: sub03_stri... or just 03_stri...)
        try:
            if 'sub' in filename:
                subject_id = filename.split('sub')[1].split('_')[0]
            else:
                subject_id = filename.split('_')[0]
        except:
            print(f"Warning: Could not extract subject ID from {filename}, skipping")
            continue
        
        # Filter by subject_ids if specified
        if subject_ids is not None and subject_id not in subject_ids:
            continue
        
        # Extract strike type (2 chars after 'stri')
        try:
            strike_idx = filename.index('stri') + 4
            strike_type = filename[strike_idx:strike_idx+2]
        except:
            print(f"Warning: Could not extract strike type from {filename}, skipping")
            continue
        
        # Extract speed (4 chars after 'speed')
        try:
            speed_idx = filename.index('speed') + 5
            speed = filename[speed_idx:speed_idx+5]
        except:
            print(f"Warning: Could not extract speed from {filename}, skipping")
            continue

        if speed_strike_pair is not None:
            current_condition = [strike_type, speed]
            if current_condition not in speed_strike_pair:
                continue
        
        # Read CSV file
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Find 'endheader' line
            header_idx = None
            for i, line in enumerate(lines):
                if 'endheader' in line.lower():
                    header_idx = i
                    break
            
            if header_idx is None:
                print(f"Warning: No 'endheader' found in {filename}, skipping")
                continue
            
            # Column headers are on the line after 'endheader'
            header_line = lines[header_idx + 1]
            headers = [h.strip() for h in header_line.split('\t') if h.strip()]
            
            # Data starts two lines after 'endheader'
            data_start_idx = header_idx + 2
            
            # Parse data
            data_dict = {header: [] for header in headers}
            for i in range(data_start_idx, len(lines)):
                line = lines[i].strip()
                if not line:
                    continue
                values = line.split('\t')
                if len(values) >= len(headers):
                    for j, header in enumerate(headers):
                        try:
                            data_dict[header].append(float(values[j]))
                        except (ValueError, IndexError):
                            data_dict[header].append(np.nan)
                            
            df = pd.DataFrame(data_dict)

            for column in column_to_plot:
                if column not in df.columns:
                    print(f"Warning: Column '{column} not found in {filename}")
                    print(f"Available columns: {list(df.columns)[:9]}...")  # Show first 5
                    continue

                # Get peak value from this file
                peak_value = df[column].max()
                
                # Store data
                condition_key = (strike_type, speed, column)
                if condition_key not in conditions:
                    conditions[condition_key] = []
                conditions[condition_key].append({
                    'subject_id': subject_id,
                    'filepath': filepath,
                    'peak_value': peak_value
                })
                
                print(f"  {filename}: subject={subject_id}, strike={strike_type}, speed={speed}, column={column}, peak={peak_value:.4f}")
                
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue
    
    if not conditions:
        raise ValueError("No valid data found to plot")
    
    # Calculate statistics for each condition
    condition_stats = {}
    for condition_key, data_list in conditions.items():
        peak_values = [d['peak_value'] for d in data_list]
        n_subjects = len(peak_values)
        
        mean_peak = np.mean(peak_values)
        sem_peak = stats.sem(peak_values)  # Standard error of mean
        
        condition_stats[condition_key] = {
            'mean': mean_peak,
            'sem': sem_peak,
            'n': n_subjects,
            'subjects': [d['subject_id'] for d in data_list]
        }
        
        strike_type, speed, column = condition_key
        print(f"\nCondition {strike_type}_{speed}: mean={mean_peak:.4f}, SEM={sem_peak:.4f}, n={n_subjects}")
        
      
    # Group by (strike, speed) for x-axis positions
    unique_conditions = sorted(set((s, sp) for s, sp, _ in condition_stats.keys()))
    
    # Plot bars grouped by column
    if existing_fig_ax:
        fig, ax = existing_fig_ax
        bar_offset = len(ax.patches) // len(unique_conditions) if ax.patches else 0
    else:
        fig, ax = plt.subplots(figsize=figsize)
        bar_offset = 0
    
    bar_width = 0.15
    
    for col_idx, column in enumerate(column_to_plot):
        x_positions = np.arange(len(unique_conditions)) + (bar_offset + col_idx) * bar_width
        
        means = []
        sems = []
        for strike, speed in unique_conditions:
            key = (strike, speed, column)
            if key in condition_stats:
                means.append(condition_stats[key]['mean'])
                sems.append(condition_stats[key]['sem'])
            else:
                means.append(0)
                sems.append(0)
        

       # Create color palette for conditions
    color_palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    n_actuators = len(column_to_plot)
    n_conditions = len(unique_conditions)
    
    # Create bars: each actuator gets a group, conditions within group have different colors
    for col_idx, column in enumerate(column_to_plot):
        # Extract actuator name
        actuator_name = column.split('/')[-1].split('|')[0] if '|' in column else column.split('/')[-1]
        
        # Base position for this actuator's group
        actuator_base_pos = col_idx * (n_conditions * bar_width + 0.2)  # 0.2 = gap between groups
        
        # Plot each condition within this actuator's group
        for cond_idx, (strike, speed) in enumerate(unique_conditions):
            key = (strike, speed, column)
            
            if key in condition_stats:
                mean_val = condition_stats[key]['mean']
                sem_val = condition_stats[key]['sem']
            else:
                mean_val = 0
                sem_val = 0
            
            # Position for this condition's bar within the actuator group
            x_pos = actuator_base_pos + cond_idx * bar_width
            
            # Condition label for legend (only add label on first actuator to avoid duplicates)
            condition_label = f"{strike}_{speed[:2]}"
            
            ax.bar(x_pos, mean_val, bar_width,
                   yerr=sem_val,
                   color=color_palette[cond_idx % len(color_palette)],
                   alpha=0.7,
                   capsize=5,
                   label=condition_label,
                   error_kw={'elinewidth': 2, 'capthick': 2})
    
    # Set x-ticks: one per actuator, centered under each group
    actuator_names = [col.split('/')[-1].split('|')[0] if '|' in col else col.split('/')[-1] 
                      for col in column_to_plot]
    
    tick_positions = [i * (n_conditions * bar_width + 0.2) + (n_conditions - 1) * bar_width / 2 
                      for i in range(n_actuators)]
    
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(actuator_names, rotation=45, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # Creates dict that removes duplicates
    ax.legend(by_label.values(), by_label.keys(), loc='best', title='Conditions')

    if output_file:
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_file}")

    return fig, ax

def write_results_to_xml(results: Dict[str, Any], 
                         output_file: str,
                         root_tag: str = "MocoObjectiveAnalysis",
                         pretty_print: bool = True) -> str:
    """
    Write analysis results dictionary to an XML file.
    
    Parameters
    ----------
    results : dict
        Dictionary containing analysis results
    output_file : str
        Path to the output XML file
    root_tag : str, optional
        Name of the root XML element (default: "MocoObjectiveAnalysis")
    pretty_print : bool, optional
        If True, format XML with indentation (default: True)
    
    Returns
    -------
    str
        Path to the created XML file
    
    Examples
    --------
    >>> results = {'objective_contact': 25.5, 'num_solutions': 10}
    >>> write_results_to_xml(results, 'analysis_results.xml')
    """
    
    def dict_to_xml(parent, data):
        """Recursively convert dictionary to XML elements."""
        if isinstance(data, dict):
            for key, value in data.items():
                # Create sanitized tag name (XML tags can't start with numbers or contain spaces)
                tag_name = str(key).replace(' ', '_').replace('-', '_')
                if tag_name[0].isdigit():
                    tag_name = 'item_' + tag_name
                
                child = ET.SubElement(parent, tag_name)
                dict_to_xml(child, value)
        elif isinstance(data, (list, tuple)):
            for i, item in enumerate(data):
                item_elem = ET.SubElement(parent, f'item_{i}')
                dict_to_xml(item_elem, item)
        else:
            # Convert value to string and set as text
            parent.text = str(data)
    
    # Create root element
    root = ET.Element(root_tag)
    
    # Add timestamp
    timestamp = ET.SubElement(root, 'timestamp')
    timestamp.text = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Add results
    results_elem = ET.SubElement(root, 'results')
    dict_to_xml(results_elem, results)
    
    # Create tree
    tree = ET.ElementTree(root)
    
    # Write to file
    if pretty_print:
        # Use minidom for pretty printing
        xml_str = ET.tostring(root, encoding='unicode')
        dom = minidom.parseString(xml_str)
        pretty_xml = dom.toprettyxml(indent='  ')
        
        with open(output_file, 'w') as f:
            f.write(pretty_xml)
    else:
        tree.write(output_file, encoding='unicode', xml_declaration=True)
    
    return output_file


def export_to_xml(output_file: Optional[str] = None,
                  root_tag: str = "MocoObjectiveAnalysis",
                  pretty_print: bool = True,
                  auto_filename: bool = True):
    """
    Decorator to automatically export function results to XML file.
    
    Parameters
    ----------
    output_file : str, optional
        Path to output XML file. If None and auto_filename=True, 
        generates filename based on function name and timestamp
    root_tag : str, optional
        Name of the root XML element (default: "MocoObjectiveAnalysis")
    pretty_print : bool, optional
        If True, format XML with indentation (default: True)
    auto_filename : bool, optional
        If True and output_file is None, auto-generate filename (default: True)
    
    Returns
    -------
    function
        Decorated function that exports results to XML
    
    Examples
    --------
    >>> @export_to_xml(output_file='my_results.xml')
    ... def my_analysis_function():
    ...     return {'result': 42}
    
    >>> @export_to_xml()  # Auto-generates filename
    ... def another_analysis():
    ...     return {'data': [1, 2, 3]}
    """
    
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Execute the original function
            results = func(*args, **kwargs)
            
            # Determine output filename
            xml_file = output_file
            if xml_file is None and auto_filename:
                timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                xml_file = f"{func.__name__}_{timestamp}.xml"
            
            # Export to XML if filename is provided or auto-generated
            if xml_file:
                try:
                    write_results_to_xml(
                        results=results,
                        output_file=xml_file,
                        root_tag=root_tag,
                        pretty_print=pretty_print
                    )
                    print(f"\nResults exported to XML: {xml_file}")
                except Exception as e:
                    print(f"\nWarning: Failed to export results to XML: {e}")
            
            return results
        
        return wrapper
    return decorator


def export_to_xml(output_file: Optional[str] = None,
                  root_tag: str = "MocoObjectiveAnalysis",
                  pretty_print: bool = True,
                  auto_filename: bool = True):
    """
    Decorator to automatically export function results to XML file.
    
    Parameters
    ----------
    output_file : str, optional
        Path to output XML file. If None and auto_filename=True, 
        generates filename based on function name and timestamp
    root_tag : str, optional
        Name of the root XML element (default: "MocoObjectiveAnalysis")
    pretty_print : bool, optional
        If True, format XML with indentation (default: True)
    auto_filename : bool, optional
        If True and output_file is None, auto-generate filename (default: True)
    
    Returns
    -------
    function
        Decorated function that exports results to XML
    
    Examples
    --------
    >>> @export_to_xml(output_file='my_results.xml')
    ... def my_analysis_function():
    ...     return {'result': 42}
    
    >>> @export_to_xml()  # Auto-generates filename
    ... def another_analysis():
    ...     return {'data': [1, 2, 3]}
    """
    
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Execute the original function
            results = func(*args, **kwargs)
            
            # Determine output filename
            xml_file = output_file
            if xml_file is None and auto_filename:
                timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                xml_file = f"{func.__name__}_{timestamp}.xml"
            
            # Export to XML if filename is provided or auto-generated
            if xml_file:
                try:
                    write_results_to_xml(
                        results=results,
                        output_file=xml_file,
                        root_tag=root_tag,
                        pretty_print=pretty_print
                    )
                    print(f"\nResults exported to XML: {xml_file}")
                except Exception as e:
                    print(f"\nWarning: Failed to export results to XML: {e}")
            
            return results
        
        return wrapper
    return decorator


def parse_objectives_from_sto_file(file_path: str, 
                                    term_labels: Dict[str, str],
                                    verbose: bool = False) -> Dict[str, float]:
    """
    Parse objective values from the header of a MoCo solution STO file.
    
    MoCo solution files store objective information in the file header as metadata.
    This function reads the file and extracts objective term values.
    
    Parameters
    ----------
    file_path : str
        Path to the MoCo solution STO file
    term_labels : dict
        Dictionary mapping keys to objective label names
    verbose : bool
        If True, print parsing information
        
    Returns
    -------
    dict
        Dictionary with objective values for each term
    """
    objective_values = {key: 0.0 for key in term_labels.keys()}
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() == 'endheader':
                break
            
            if not line.strip().startswith('#'):
                continue
            
            # Remove '# ' and parse label=value
            content = line.strip().lstrip('#').strip()
            if '=' not in content:
                continue
            
            parts = content.split('=', 1)
            if len(parts) != 2:
                continue
            
            label_in_file = parts[0].strip()
            value_str = parts[1].strip()
            
            # Exact label match only
            for key, label in term_labels.items():
                if label_in_file == label:
                    try:
                        objective_values[key] = float(value_str)
                    except ValueError:
                        pass
                    break
    
    return objective_values

# def parse_objectives_from_sto_file(file_path: str, 
#                                 term_labels: Dict[str, str],
#                                 verbose: bool = False) -> Dict[str, float]:
#     """
#     Parse objective values from MoCo solution STO file header metadata.
    
#     Parameters
#     ----------
#     file_path : str
#         Path to the MoCo solution STO file
#     term_labels : dict
#         Dictionary mapping keys to objective label names
#     verbose : bool
#         If True, print parsing information
        
#     Returns
#     -------
#     dict
#         Dictionary with objective values for each term
#     """
#     objective_values = {key: 0.0 for key in term_labels.keys()}
    
#     try:
#         # Read the entire table to access metadata
#         table = osim.TimeSeriesTable(file_path)
        
#         # Objective information is stored in table metadata
#         # Try to access each objective term from metadata
#         for key, label in term_labels.items():
#             try:
#                 # Check if this key exists in the table metadata
#                 if table.hasTableMetaDataKey(label):
#                     value_str = table.getTableMetaDataAsString(label)
#                     objective_values[key] = float(value_str)
#                     if verbose:
#                         print(f"    Found {label}: {objective_values[key]}")
#             except Exception as e:
#                 if verbose:
#                     print(f"    Could not extract {label}: {e}")
#                 objective_values[key] = 0.0
    
#     except Exception as e:
#         if verbose:
#             print(f"  Error accessing table metadata from {file_path}: {e}")
    
#     return objective_values


def _parse_objectives_from_sto_header(filepath: str, term_labels: Dict[str, str]) -> Dict[str, float]:
    """Parse objective values from STO file header metadata."""
    objective_values = {key: 0.0 for key in term_labels.keys()}
    
    with open(filepath, 'r') as f:
        for line in f:
            # Stop when we reach the data section
            if line.strip() == 'endheader':
                break

            content = line.strip().lstrip('#').strip()
            if '=' not in content:
                continue
            parts = content.split('=', 1)
            if len(parts) !=2:
                continue
            label_in_file = parts[0].strip()
            value_str = parts[1].strip()
            
            if label_in_file == 'objective':
                try:
                    objective_values['objective'] = float(value_str)
                except ValueError:
                    pass 
                continue

            # Look for objective term lines in header
            for key, label in term_labels.items():
                if key != 'objective' and label_in_file == label:
                    try:
                        objective_values[key] = float(value_str)
                    except ValueError:
                        pass
                    break
    
    return objective_values


def compute_objective_term_percentages(directory: str, 
                                       solution_suffix: str = '.sto',
                                       output_xml: Optional[str] = None,
                                       verbose: bool = True) -> Dict[str, float]:
    """
    Compute average percentage values of MoCo objective terms across multiple solution files.
    
    Parameters
    ----------
    directory : str
        Path to the directory containing MoCo solution files
    solution_suffix : str, optional
        Common suffix shared by all solution files (default: '.sto')
    output_xml : str, optional
        If provided, saves results to an XML file at this path (default: None)
    verbose : bool, optional
        If True, print progress information (default: True)
    
    Returns
    -------
    dict
        Dictionary with average percentage contributions for each objective term
    """
    import xml.etree.ElementTree as ET
    from xml.dom import minidom
    
    # Find solution files
    solution_files = [f for f in os.listdir(directory) if f.endswith(solution_suffix)]
    
    if not solution_files:
        raise FileNotFoundError(
            f"No solution files with suffix '{solution_suffix}' found in {directory}"
        )
    
    if verbose:
        print(f"Found {len(solution_files)} solution files in {directory}")
    
    # Objective term labels
    term_labels = {
        'objective_contact': 'objective_contact',
        'objective_control_effort': 'objective_control_effort',
        'objective_feet_orientation_goal': 'objective_feet_orientation_goal',
        'objective_feet_speed_goal': 'objective_feet_speed_goal',
        'objective_state_tracking': 'objective_state_tracking',
        'objective_torso_angular_velocity_goal': 'objective_torso_angular_velocity_goal',
        'objective_torso_orientation_goal': 'objective_torso_orientation_goal',
        'objective': 'objective'
    }
    
    # Store percentages
    all_percentages = {
        'objective_contact': [],
        'objective_control_effort': [],
        'objective_feet_goals': [],
        'objective_state_tracking': [],
        'objective_torso_goals': []
    }
    
    successful_files = []
    failed_files = []
    
    # Process each solution file
    for sol_file in solution_files:
        sol_path = os.path.join(directory, sol_file)
        
        try:
            # Try getObjectiveTerm first (newer OpenSim versions)
            objective_values = {}
            use_file_parsing = False
            
            try:
                trajectory = osim.MocoTrajectory(sol_path)
                for key, label in term_labels.items():
                    value = trajectory.getObjectiveTerm(label)
                    objective_values[key] = value
            except AttributeError:
                # Fallback to parsing from file header
                use_file_parsing = True
                objective_values = _parse_objectives_from_sto_header(sol_path, term_labels)
            
            if verbose and use_file_parsing:
                print(f"  {sol_file}: Using file header parsing")
            
            total_objective = objective_values['objective']
            
            if total_objective == 0:
                if verbose:
                    print(f"  Warning: Total objective is zero for {sol_file}, skipping...")
                failed_files.append(sol_file)
                continue

            with open('retrieve_results.py', 'r') as f:
                lines = f.readlines()

            # Find line "continue" after "if total_objective == 0"
            insert_pos = None
            for i in range(4000, 4010):
                if 'continue' in lines[i] and 'total_objective == 0' in lines[i-4]:
                    insert_pos = i + 1
                    break

            if insert_pos:
                verification_code = '''            
                        # Verify terms sum to total
                        terms_sum = (objective_values['objective_contact'] +
                                    objective_values['objective_control_effort'] +
                                    objective_values['objective_state_tracking'] +
                                    objective_values['objective_feet_orientation_goal'] +
                                    objective_values['objective_feet_speed_goal'] +
                                    objective_values['objective_torso_angular_velocity_goal'] +
                                    objective_values['objective_torso_orientation_goal'])
                        
                        if verbose and abs(terms_sum - total_objective) > 1e-6:
                            print(f"  WARNING {sol_file}: Terms sum ({terms_sum:.6f}) != Total ({total_objective:.6f})")
                            print(f"    Individual values:")
                            for key, val in objective_values.items():
                                print(f"      {key}: {val:.6e}")
            '''
                
                lines.insert(insert_pos, verification_code)
                
                with open('retrieve_results.py', 'w') as f:
                    f.writelines(lines)
                
                print(f"Inserted verification code at line {insert_pos}")
            else:
                print("Could not find insertion point")
            
            # Calculate percentages
            contact_pct = (objective_values['objective_contact'] / total_objective) * 100
            control_pct = (objective_values['objective_control_effort'] / total_objective) * 100 
            state_tracking_pct = (objective_values['objective_state_tracking'] / total_objective) * 100 
            
            feet_goals_pct = (
                (objective_values['objective_feet_orientation_goal'] + 
                 objective_values['objective_feet_speed_goal']) / total_objective
            ) * 100 
            
            torso_goals_pct = (
                (objective_values['objective_torso_angular_velocity_goal'] + 
                 objective_values['objective_torso_orientation_goal']) / total_objective
            ) * 100 
            
            # Store percentages
            all_percentages['objective_contact'].append(contact_pct)
            all_percentages['objective_control_effort'].append(control_pct)
            all_percentages['objective_feet_goals'].append(feet_goals_pct)
            all_percentages['objective_state_tracking'].append(state_tracking_pct)
            all_percentages['objective_torso_goals'].append(torso_goals_pct)
            
            successful_files.append(sol_file)
            
            if verbose:
                print(f"\n  {sol_file}:")
                print(f"    Contact: {contact_pct:.2f}%")
                print(f"    Control Effort: {control_pct:.2f}%")
                print(f"    Feet Goals: {feet_goals_pct:.2f}%")
                print(f"    State Tracking: {state_tracking_pct:.2f}%")
                print(f"    Torso Goals: {torso_goals_pct:.2f}%")
                print(f"    Total: {total_objective:.6f}")
                
        except Exception as e:
            if verbose:
                print(f"  Error processing {sol_file}: {e}")
            failed_files.append(sol_file)
            continue
    
    # Calculate averages
    if not successful_files:
        raise ValueError(
            f"No solutions could be successfully processed. "
            f"Failed files: {failed_files}"
        )
    
    average_percentages = {}
    std_percentages = {}
    for key, values in all_percentages.items():
        if values:
            average_percentages[key] = np.mean(values)
            std_percentages[key] = np.std(values, ddof=1) if len(values) > 1 else 0.0
        else:
            average_percentages[key] = 0.0
            std_percentages[key] = 0.0
    
    # Add metadata
    average_percentages['num_solutions'] = len(successful_files)
    average_percentages['solution_files'] = successful_files
    average_percentages['std_percentages'] = std_percentages
    
    if verbose:
        print("\n" + "="*60)
        print(f"AVERAGE PERCENTAGES (n={len(successful_files)} solutions):")
        print("="*60)
        print(f"  Contact: {average_percentages['objective_contact']:.2f}% ± {std_percentages['objective_contact']:.2f}%")
        print(f"  Control Effort: {average_percentages['objective_control_effort']:.2f}% ± {std_percentages['objective_control_effort']:.2f}%")
        print(f"  Feet Goals: {average_percentages['objective_feet_goals']:.2f}% ± {std_percentages['objective_feet_goals']:.2f}%")
        print(f"  State Tracking: {average_percentages['objective_state_tracking']:.2f}% ± {std_percentages['objective_state_tracking']:.2f}%")
        print(f"  Torso Goals: {average_percentages['objective_torso_goals']:.2f}% ± {std_percentages['objective_torso_goals']:.2f}%")
        print("="*60)
        
        # Verify total is approximately 100%
        total_pct = sum([
            average_percentages['objective_contact'],
            average_percentages['objective_control_effort'],
            average_percentages['objective_feet_goals'],
            average_percentages['objective_state_tracking'],
            average_percentages['objective_torso_goals']
        ])
        print(f"  Sum of percentages: {total_pct:.2f}%")
        
        if failed_files:
            print(f"\n  Warning: {len(failed_files)} files could not be processed:")
            for f in failed_files:
                print(f"    - {f}")
    
    # Write to XML if requested
    if output_xml:
        root = ET.Element('ObjectiveTermAnalysis')
        
        # Metadata
        metadata = ET.SubElement(root, 'Metadata')
        ET.SubElement(metadata, 'AnalysisDate').text = datetime.datetime.now().isoformat()
        ET.SubElement(metadata, 'Directory').text = directory
        ET.SubElement(metadata, 'SolutionSuffix').text = solution_suffix
        ET.SubElement(metadata, 'NumSolutions').text = str(len(successful_files))
        
        # Average percentages with standard deviations
        averages = ET.SubElement(root, 'AveragePercentages')
        
        contact_elem = ET.SubElement(averages, 'Contact', unit='percent')
        contact_elem.text = f"{average_percentages['objective_contact']:.4f}"
        contact_elem.set('std', f"{std_percentages['objective_contact']:.4f}")
        
        control_elem = ET.SubElement(averages, 'ControlEffort', unit='percent')
        control_elem.text = f"{average_percentages['objective_control_effort']:.4f}"
        control_elem.set('std', f"{std_percentages['objective_control_effort']:.4f}")
        
        feet_elem = ET.SubElement(averages, 'FeetGoals', unit='percent')
        feet_elem.text = f"{average_percentages['objective_feet_goals']:.4f}"
        feet_elem.set('std', f"{std_percentages['objective_feet_goals']:.4f}")
        
        tracking_elem = ET.SubElement(averages, 'StateTracking', unit='percent')
        tracking_elem.text = f"{average_percentages['objective_state_tracking']:.4f}"
        tracking_elem.set('std', f"{std_percentages['objective_state_tracking']:.4f}")
        
        torso_elem = ET.SubElement(averages, 'TorsoGoals', unit='percent')
        torso_elem.text = f"{average_percentages['objective_torso_goals']:.4f}"
        torso_elem.set('std', f"{std_percentages['objective_torso_goals']:.4f}")
        
        total_elem = ET.SubElement(averages, 'Total', unit='percent')
        total_elem.text = f"{total_pct:.4f}"
        
        # Processed solutions
        solutions = ET.SubElement(root, 'ProcessedSolutions')
        for sol_file in successful_files:
            ET.SubElement(solutions, 'Solution').text = sol_file
        
        # Failed files
        if failed_files:
            failed = ET.SubElement(root, 'FailedSolutions')
            for fail_file in failed_files:
                ET.SubElement(failed, 'Solution').text = fail_file
        
        # Pretty print XML
        xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
        
        with open(output_xml, 'w') as f:
            f.write(xml_str)
        
        if verbose:
            print(f"\nResults saved to: {output_xml}")
    
    return average_percentages


# # Example usage
# if __name__ == "__main__":
#     # Example file paths - adjust these to your actual file paths
#     moco_solution = "example_muscle_driven_47.sto"
#     ik_kinematics = 'converted_states.sto' # "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_strifc_cond00000_speed02000_0001_0003_rad.sto"
#     model_file = "example_muscle_adjusted_model.osim"
#     omoco_setup_file = "example_muscle_driven_47_spring.omoco"  # Generated by study.printToXML()
#     output_error_file = "moco_ik_tracking_errors.sto"
    
#     solution = osim.MocoTrajectory(moco_solution)
#     model_proc = osim.ModelProcessor(model_file)
#     model = model_proc.process()
#     model.initSystem()

#     spring_sto = moco_solution.replace('.sto', '_springforces.sto')
#     write_mtp_spring_forces(solution=solution, model=model, output_file=spring_sto)
#     main_coords = ['/jointset/hip_r/hip_flexion_r',
#                 '/jointset/walker_knee_r/knee_angle_r',
#                 '/jointset/ankle_r/ankle_angle_r',
#                 '/jointset/mtp_r/mtp_angle_r']
#     dropoff_folder = './sub03_testing/'
#     moments_file = moco_solution.replace('.sto', '_leg_net_moments.sto')
#     moments_path = os.path.join(dropoff_folder, moments_file)
#     plot_joint_moment_constituents(model, solution, coord_paths=main_coords, 
#                                 muscle_paths=None, output_dir=moments_file)

#     # solution = osim.TimeSeriesTable(moco_solution)
#     # print(f'Solution has {solution.getNumRows()} rows and {solution.getNumColumns()} columns.')

#     # Method 1: Complete analysis with weighted errors in .sto file (RECOMMENDED)
#     try:
#         print("=== Complete Error Analysis with Weighted Errors ===")
#         detailed_results = analyze_moco_ik_tracking_errors(
#             moco_solution_file=moco_solution,
#             ik_kinematics_file=ik_kinematics,
#             output_file=output_error_file,
#             model_file=model_file,
#             omoco_setup_file=omoco_setup_file  # This adds weighted error columns
#         )
        
#         print(f"\nDetailed analysis completed!")
#         print(f"Results saved to: {detailed_results['output_file']}")
        
#         if 'coordinate_weights_used' in detailed_results:
#             print(f"Weighted error columns added for {detailed_results['num_weighted_coordinates']} coordinates")
#             print("The .sto file now contains:")
#             print("  - Raw errors (coordinate_error)")
#             print("  - Weighted errors (coordinate_weighted_error)")
#             print("  - Weighted RMS errors (weighted_rms_error_all_coords)")
#             print("  - Max weighted error tracking")
        
#         # Also compute the overall weighted quotient
#         print("\n=== Computing Overall Weighted Error Quotient ===")
#         auto_results = auto_compute_weighted_error_quotient(
#             moco_solution_file=moco_solution,
#             ik_kinematics_file=ik_kinematics,
#             omoco_setup_file=omoco_setup_file,
#             model_file=model_file,
#             verbose=False  # Less verbose since we already have detailed output
#         )
        
#         print(f"Overall Weighted Error Quotient: {auto_results['weighted_error_quotient']:.6f}")
        
#         # Save comprehensive results
#         import json
#         results_file = "complete_tracking_analysis.json"
        
#         # Combine results (exclude detailed_results to avoid duplication)
#         combined_results = {
#             'weighted_error_quotient': auto_results['weighted_error_quotient'],
#             'coordinate_weights': auto_results['extracted_weights'],
#             'analysis_summary': {
#                 'overall_rms_error': detailed_results['overall_rms_error'],
#                 'num_coordinates': detailed_results['num_matched_coordinates'],
#                 'num_weighted_coordinates': detailed_results.get('num_weighted_coordinates', 0),
#                 'time_range': detailed_results['time_range'],
#                 'num_time_points': detailed_results['num_time_points']
#             },
#             'files': {
#                 'moco_solution': moco_solution,
#                 'ik_kinematics': ik_kinematics,
#                 'setup_file': omoco_setup_file,
#                 'output_errors': output_error_file,
#                 'model_file': model_file
#             }
#         }
        
#         with open(results_file, 'w') as f:
#             json.dump(combined_results, f, indent=2, default=str)
#         print(f"Complete results saved to: {results_file}")
        
#     except Exception as e:
#         print(f"Error in complete analysis: {e}")
#         import traceback
#         traceback.print_exc()
        
#         # Fallback: analysis without weights
#         print("\nTrying analysis without weighted errors...")
#         try:
#             fallback_results = analyze_moco_ik_tracking_errors(
#                 moco_solution_file=moco_solution,
#                 ik_kinematics_file=ik_kinematics,
#                 output_file=output_error_file,
#                 model_file=model_file
#                 # No omoco_setup_file = no weighted errors
#             )
#             print(f"Fallback analysis completed: {fallback_results['output_file']}")
            
#         except Exception as e2:
#             print(f"Fallback analysis also failed: {e2}")

# Example usage
if __name__ == "__main__":
    
    # Example file paths - adjust these to your actual file paths
#     moco_solution = "example_muscle_driven_21_nospring.sto"
#     ik_kinematics = "converted_states.sto"
#     output_error_file = "moco_ik_tracking_errors.sto"
#     model_file = "example_muscle_adjusted_model.osim"
#     omoco_setup_file = "example_muscle_driven_test_11_spring.omoco"  # Generated by study.printToXML()

#     solution_table = osim.TimeSeriesTable(moco_solution)
#     time_vector = solution_table.getIndependentColumn()
#     print('Time vector from solution:')
#     print(time_vector)
    
#     solution = osim.MocoTrajectory(moco_solution)
#     model_proc = osim.ModelProcessor(model_file)
#     model = model_proc.process()
#     # model.initSystem()
    
    
#         # contact_table = create_contact_sphere_force_table(model, solution)

#         # contacts_file = moco_solution.replace('.sto', 'contact_forces.sto')
#         # cop_file = moco_solution.replace('.sto', 'cop_output.sto')

#         # contact_forces, copTable = get_solution_contacts(model, solution)
#         # osim.STOFileAdapter.write(contact_forces, contacts_file)
#         # osim.STOFileAdapter.write(copTable, cop_file)

        
#     except Exception as e:
#         print(f"Error during analysis: {e}")
#         import traceback
#         traceback.print_exc()

# output_path = "moco_ik_tracking_errors.sto"

    # Example function calls showing how to use the function multiple times
    # to create plots with multiple averages
    
    input_dir = "/home/lisca/biomec/data/20250306_dataset/results_pool/"
    plot_dir = os.path.join(input_dir, 'plots/')
    avg_inputs = os.path.join(input_dir, 'averages/')
    # model_path = "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim"
    model_path = "sub03_muscle_adjusted_model.osim"
    
    track_solution_file="/home/lisca/biomec/data/20250306_dataset/results_pool/averages/end_results/sub03_strift_cond00000_speed03000_20res_02mesh.sto", 
    moco_solution_file= '/home/lisca/biomec/data/20250306_dataset/results_pool/averages/guesses/sub03_strifc_cond00000_speed02000_mi_emg_solution.sto'
    ik_kinematics_file="/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub04_strifc_cond00000_speed02000_ik_move_averages.sto",     
    omoco_setup_file='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/guesses/sub04_strifc_cond00000_speed02000_20res_02mesh_spring.omoco', 



    solution = osim.MocoTrajectory(moco_solution_file)
    model_proc = osim.ModelProcessor(model_path)
    model = model_proc.process()
    model.initSystem()

    # solution = osim.MocoTrajectory(moco_solution_file)
    # table = osim.TimeSeriesTable(moco_solution_file)
    # solution = osim.MocoTrajectory.createFromStatesTable(table)

    contacts_file = moco_solution_file.replace('.sto', '_contact_forces.sto')
    cop_file = moco_solution_file.replace('.sto', '_cop_output.sto')

    contact_forces, copTable = get_solution_contacts(model, solution)
    osim.STOFileAdapter.write(contact_forces, contacts_file)
    osim.STOFileAdapter.write(copTable, cop_file)

    main_coords = ['/jointset/hip_r/hip_flexion_r',
                '/jointset/walker_knee_r/knee_angle_r',
                '/jointset/ankle_r/ankle_angle_r',
                '/jointset/mtp_r/mtp_angle_r']
    fibers_file = ('fiber_profile')

    plot_joint_fiber_forces(model, track_solution_file, coord_paths=main_coords,
                           file_name=fibers_file, output_dir='./', 
                           min_moment_arm_threshold=0.00001,
                           normalize_time=True)
    

    results = compute_objective_term_percentages(
    directory='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/end_results/',
    solution_suffix='_20res_02mesh.sto',
    output_xml= os.path.join(input_dir, 'averages/end_results/cost_term_contribution.xml'),
    verbose=True)

    print("\nResults:")
    print(f"Contact: {results['objective_contact']:.2f}%")
    print(f"Control Effort: {results['objective_control_effort']:.2f}%")
    print(f"State Tracking: {results['objective_state_tracking']:.2f}%")
    print(f"Feet Goals: {results['objective_feet_goals']:.2f}%")
    print(f"Torso Goals: {results['objective_torso_goals']:.2f}%")
    print(f"\nProcessed {results['num_solutions']} solutions")

    # # Plot soleus activation
    # fig, ax = plot_muscle_mechanics_bars(
    # csv_directory=os.path.join(avg_inputs, 'end_results/'),
    # suffix='_20res_02mesh_muscle_mechanics.csv',
    # column_to_plot= ['/forceset/soleus_r|tendon_force', '/forceset/gaslat_r|tendon_force'],
    # subject_ids=['03', '04'],
    # speed_strike_pair=[['fc', '02000']],
    # plot_color='blue',
    # ylabel='Peak Activation',
    # title='Soleus Peak Activation by Condition',
    # output_file='z_soleus_activation.png'
    # )

    # First call - creates new figure
    # fig, ax = compute_group_averages(
    #     input_directory=input_dir,
    #     output_directory=plot_dir,
    #     sto_output_directory=avg_inputs,
    #     file_content_type="emg.sto",
    #     strike_type="fc",
    #     speed="02000",
    #     generate_plot=True,
    #     # columns_to_plot=['hip_flexion_r', 'ankle_angle_r', 'knee_angle_r', 'mtp_angle_r'],
    #     columns_to_plot=['right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_soleus'],
    #     subject_id="03"  # Optional
    # )
 
    
    # # Second call - adds to existing figure
    # fig, ax = compute_group_averages(
    #     input_directory=input_dir,
    #     output_directory=plot_dir,
    #     sto_output_directory=input_dir,
    #     file_content_type="ik_move.sto", 
    #     strike_type="ft",
    #     speed="02000",
    #     columns_to_plot=['ankle_angle_r'],
    #     # subject_id="03",
    #     # existing_fig_ax=(fig, ax)  # Add to existing plot
    # )
    
    # plot_moco_solution_columns(moco_solution_file="/home/lisca/biomec/data/20250306_dataset/results_pool/averages/guesses/sub04_strifc_cond00000_speed02000_20res_02mesh.sto", 
    #                            columns_to_plot=['/jointset/ankle_r/ankle_angle_r/value'],
    #                            model_file=model_path,
    #                            subject_id='sub04',
    #                            existing_fig_ax=(fig, ax),
    #                            normalize_time=True,
    #                            time_as_percentage=True,
    #                            plot_title=None,
    #                            output_file='z_comparison.png',
    #                            line_style='-',
    #                            line_width=2,
    #                            color=None,
    #                            label_prefix=None,
    #                            show_legend=True,
    #                            legend_location='upper left',
    #                            grid=True,
    #                            xlabel=None,
    #                            ylabel='ankle plantarflexion angle (degrees)',
    #                            figsize=(12, 8))  
    
    # Third call - adds to existing figure
    # fig, ax = compute_group_averages(
    #     input_directory=input_dir,
    #     output_directory=plot_dir,
    #     sto_output_directory=avg_inputs,
    #     file_content_type="ik_move.sto",
    #     strike_type="rc",
    #     speed="03000",
    #     apply_mass_norm= False,
    #     model_path=model_path,
    #     # grf_columns=['ground_1_force_vy', 'ground_2_force_vy'],
    #     # columns_to_plot=['ground_1_force_vy', 'ground_2_force_vy'],
    #     columns_to_plot=['ankle_angle_r', 'hip_flexion_r', 'knee_angle_r'],
    #     generate_plot=True,
    #     # subject_id="03",
    #     sub_masses= {
    #         "03": 57.8,
    #         # "04": 88.7
    #     }
    # )
    
    # plt.show()
    # analyze_moco_ik_tracking_errors(moco_solution_file="/home/lisca/biomec/data/20250306_dataset/results_pool/averages/guesses/sub04_strifc_cond00000_speed02000_20res_02mesh.sto", 
    #                                 ik_kinematics_file="/home/lisca/biomec/data/20250306_dataset/results_pool/averages/sub04_strifc_cond00000_speed02000_ik_move_averages.sto", 
    #                                 output_file='z_errors.sto',
    #                                 model_file='/home/lisca/biomec/data/20250306_dataset/sub04/ses20250524/sub04_SCALED_addbio_free_subt_mtp.osim',
    #                                 omoco_setup_file='/home/lisca/biomec/data/20250306_dataset/results_pool/averages/guesses/sub04_strifc_cond00000_speed02000_20res_02mesh_spring.omoco', 
    #                                 include_weighted_quotients=False)

