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
from osimFunctions import kinematicsToStates, compute_sto_stdevs, generate_guess_file
import pandas as pd
import datetime
import sys
from datetime import datetime

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Geometry")

path_opensim_core = "/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/"

name_problem = 'track_muscle_driven_running'

path_model = os.path.join(
    path_opensim_core,
    "sub03_SCALED.osim")
path_grf = os.path.join(
    path_opensim_core,
    'sub03_strifc_cond00000_speed02000_0001_0003_grf.xml')
path_markers = os.path.join(
    path_opensim_core,
    'sub03_strifc_cond00000_speed02000_0001_0003_trc.trc')
path_states = os.path.join(
    path_opensim_core,
    'sub03_strifc_cond00000_speed02000_0001_0003_ik.sto')

path_solution_guess = 'muscle_driven_state_tracking_sol_5.sto'
path_solution       = os.path.join('./sub03_testing/', 'muscle_driven_state_tracking_sol_9.sto')


# CONSIDER PUTTING THE FOLLOWING FUNCTIONS IN ANOTHER SCRIPT AND IMPORTING FROM THERE FOR ORGANIZATION:
# ====================================================================================================
def create_contact_sphere_force_table(model, solution):
    
    model.initSystem()
    external_forces_table = osim.TimeSeriesTableVec3()
    cop_table = osim.TimeSeriesTableVec3()
    states_trajectory = solution.exportToStatesTrajectory(model)
    num_states = states_trajectory.getSize()

    forceNames = ['forceset/contactHeel',
                  'forceset/contactLateralRearfoot',
                  'forceset/contactLateralMidfoot',
                  'forceset/contactLateralToe',
                  'forceset/contactMedialToe',
                  'forceset/contactMedialMidfoot'
                ]
    forceLabels = ['heel',
                   'lat_rear',
                   'lat_mid', 
                   'lat_toe', 
                   'med_toe', 
                   'med_mid',
                    ]

    sphereNames = ['contactgeometryset/heel',
                   'contactgeometryset/lateralRearfoot',
                   'contactgeometryset/lateralMidfoot',
                   'contactgeometryset/lateralToe',
                   'contactgeometryset/medialToe',
                   'contactgeometryset/medialMidfoot',
                    ]
    
    for istate in range(num_states):
        state = states_trajectory.get(istate)
        model.realizeVelocity(state)
        row = osim.RowVectorVec3(3*2*len(forceNames))
        cop_row = osim.RowVectorVec3(2)
        labels = osim.StdVectorString()

        for iside, side in enumerate(['_r', '_l']):
            cop = osim.Vec3(0)
            cop_torque_sum_x = 0
            cop_force_sum = 0
            cop_torque_sum_z = 0

            offset = 3 * iside * len(sphereNames)
            zipped = zip(forceNames, forceLabels, sphereNames)

            for i, (forceName, forceLabel, sphereName) in enumerate(zipped):
                force = osim.Vec3(0)
                torque = osim.Vec3(0)

                force_obj = osim.Force.safeDownCast(
                    model.getComponent(f'{forceName}{side}'))
                force_values = force_obj.getRecordValues(state)

                force[0] = force_values.get(0)
                force[1] = force_values.get(1)
                force[2] = force_values.get(2)
                torque[0] = force_values.get(3)
                torque[1] = force_values.get(4)
                torque[2] = force_values.get(5)

                sphere = osim.ContactSphere.safeDownCast(model.getComponent(f'{sphereName}{side}'))
                frame = sphere.getFrame()
                position = frame.getPositionInGround(state)
                location = sphere.get_location()
                locationInGround = frame.expressVectorInGround(state, location)

                position[0] = position[0] + locationInGround[0]
                position[1] = position[1] + locationInGround[1]
                position[2] = position[2] + locationInGround[2]

                row[3*i + offset] = force
                row[3*i + 1 + offset] = position
                row[3*i + 2 + offset] = torque

                #** CHECK THE CALCULATIONS ARE CORRECT FOR THE BELOW:
                cop_force_sum += force[1]
                cop_torque_sum_x += force[1] * position[0]
                cop_torque_sum_z += force[1] * position[2]

                for suffix in ['_force_v', '_force_p', '_torque_']:
                    labels.append(f'{forceLabel}{side}{suffix}')

            if np.abs(cop_force_sum) > 1e-3:
                cop[0] = cop_torque_sum_x / cop_force_sum
                cop[2] = cop_torque_sum_z / cop_force_sum

            cop_row[iside] = cop

        external_forces_table.appendRow(state.getTime(), row)
        cop_table.appendRow(state.getTime(), cop_row)

    external_forces_table.setColumnLabels(labels)

    labels = osim.StdVectorString()
    labels.append('2_ground_force_p')
    labels.append('1_ground_force_p')
    cop_table.setColumnLabels(labels)

    suffixes = osim.StdVectorString()
    suffixes.append('x')
    suffixes.append('y')
    suffixes.append('z')
    return external_forces_table.flatten(suffixes), cop_table



def create_valid_path(path):
    path = f'{path}'
    path = path.replace(' ', '')
    path = path.replace(':', '_')
    path = path.replace('.', '_')
    return path

def get_sol_archive_path(result_path, name):
    now = datetime.datetime.now()
    now.strftime('%Y-%m-%dT%H:%M:%S')
    now = create_valid_path(now)
    return os.path.join(f'{result_path}', 'archive', 
                        f'{name}_{now}.sto')



def setup_logging(log_file_path=None, also_print_to_console=True):
    """
    Set up logging to redirect print statements to a file while optionally keeping console output.
    
    Parameters:
    -----------
    log_file_path : str, optional
        Path to the log file. If None, creates a timestamped file.
    also_print_to_console : bool
        If True, output goes to both file and console. If False, only to file.
    
    Returns:
    --------
    file_handle : file object
        The file handle for the log file (to close later)
    """
    
    # Create log file name if not provided
    if log_file_path is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file_path = f"moco_tracking_log_{timestamp}.txt"
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(log_file_path)), exist_ok=True)
    
    # Open the log file
    log_file = open(log_file_path, 'w', encoding='utf-8')
    
    if also_print_to_console:
        # Create a custom class that writes to both file and console
        class TeeOutput:
            def __init__(self, file1, file2):
                self.file1 = file1
                self.file2 = file2
            
            def write(self, data):
                self.file1.write(data)
                self.file2.write(data)
                self.file1.flush()  # Ensure immediate writing
                self.file2.flush()
            
            def flush(self):
                self.file1.flush()
                self.file2.flush()
        
        # Redirect stdout to both console and file
        sys.stdout = TeeOutput(sys.stdout, log_file)
        sys.stderr = TeeOutput(sys.stderr, log_file)
    else:
        # Redirect stdout only to file
        sys.stdout = log_file
        sys.stderr = log_file
    
    print(f"=== MOCO TRACKING LOG STARTED ===")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Log file: {os.path.abspath(log_file_path)}")
    print(f"Console output: {'Enabled' if also_print_to_console else 'Disabled'}")
    print("=" * 50)
    
    return log_file

def debug_model_processor(model_processor, original_model_path):
    """
    Debug function to analyze the ModelProcessor and its resulting model states.
    Enhanced with detailed force categorization using getConcreteClassName().
    """
    print(f"\n" + "="*80)
    print(f"🔧 DEBUGGING MODEL PROCESSOR")
    print(f"="*80)
    
    try:
        # Process the ModelProcessor to get the final model
        print(f"📋 Original model path: {original_model_path}")
        print(f"🔄 Processing ModelProcessor...")
        
        processed_model = model_processor.process()
        processed_model.initSystem()
        
        print(f"✅ ModelProcessor successfully processed")
        print(f"📊 Final model name: '{processed_model.getName()}'")
        
        # Get all state variable names from the processed model
        allStates = processed_model.getStateVariableNames()
        print(f"📈 Total state variables in processed model: {allStates.getSize()}")
        
        # Categorize the states
        coordinateValues = []
        coordinateSpeeds = []
        muscleActivations = []
        muscleFiberLengths = []
        muscleTendonForces = []
        torqueActivations = []
        reserveActivations = []
        otherStates = []
        
        print(f"\n🔍 Categorizing all {allStates.getSize()} states...")
        
        for i in range(allStates.getSize()):
            stateName = allStates.get(i)
            
            if '/value' in stateName:
                coordinateValues.append(stateName)
            elif '/speed' in stateName:
                coordinateSpeeds.append(stateName)
            elif '/activation' in stateName:
                if 'torque_' in stateName.lower():
                    torqueActivations.append(stateName)
                elif 'reserve_' in stateName.lower():
                    reserveActivations.append(stateName)
                else:
                    muscleActivations.append(stateName)
            elif '/fiber_length' in stateName:
                muscleFiberLengths.append(stateName)
            elif '/normalized_tendon_force' in stateName:
                muscleTendonForces.append(stateName)
            else:
                otherStates.append(stateName)
        
        print(f"\n📊 STATE BREAKDOWN:")
        print(f"   🎯 Coordinate values: {len(coordinateValues)}")
        print(f"   🏃 Coordinate speeds: {len(coordinateSpeeds)}")
        print(f"   💪 Muscle activations: {len(muscleActivations)}")
        print(f"   📏 Muscle fiber lengths: {len(muscleFiberLengths)}")
        print(f"   🔗 Muscle tendon forces: {len(muscleTendonForces)}")
        print(f"   ⚙️  Torque activations: {len(torqueActivations)}")
        print(f"   🔄 Reserve activations: {len(reserveActivations)}")
        print(f"   ❓ Other states: {len(otherStates)}")
        
        # === ENHANCED FORCE CATEGORIZATION SECTION ===
        print(f"\n" + "="*60)
        print(f"🔧 DETAILED FORCE ANALYSIS (getConcreteClassName)")
        print(f"="*60)
        
        # Initialize force categorization dictionaries
        force_categories = {
            'DeGrooteFregly2016Muscle': [],
            'Thelen2003Muscle': [],
            'Millard2012EquilibriumMuscle': [],
            'Other_Muscle': [],
            'CoordinateActuator': [],
            'PointActuator': [],
            'TorqueActuator': [],
            'SpringGeneralizedForce': [],
            'SmoothSphereHalfSpaceForce': [],
            'HuntCrossleyForce': [],
            'ExternalForce': [],
            'Other_Force': []
        }
        
        num_forces = 0
        muscles = []
        torqueActuators = []
        reserveActuators = []
        contactForces = []
        externalForces = []
        otherForces = []
        
        print(f"🔍 Analyzing all force components in model...")
        
        # Iterate through all components and categorize forces
        for actu in processed_model.getComponentsList():
            className = actu.getConcreteClassName()
            componentName = actu.getName()
            componentPath = actu.getAbsolutePathString()
            
            # Check if this is a force-generating component
            if (className.endswith('Muscle') or 
                className.endswith('Actuator') or 
                className.endswith('Force') or
                'Force' in className):
                
                num_forces += 1
                
                # Detailed categorization based on concrete class name
                if 'DeGrooteFregly2016Muscle' in className:
                    force_categories['DeGrooteFregly2016Muscle'].append((componentName, componentPath))
                    muscles.append(componentName)
                elif 'Thelen2003Muscle' in className:
                    force_categories['Thelen2003Muscle'].append((componentName, componentPath))
                    muscles.append(componentName)
                elif 'Millard2012EquilibriumMuscle' in className:
                    force_categories['Millard2012EquilibriumMuscle'].append((componentName, componentPath))
                    muscles.append(componentName)
                elif className.endswith('Muscle'):
                    force_categories['Other_Muscle'].append((componentName, componentPath))
                    muscles.append(componentName)
                elif 'CoordinateActuator' in className:
                    force_categories['CoordinateActuator'].append((componentName, componentPath))
                    if 'torque' in componentName.lower():
                        torqueActuators.append(componentName)
                    elif 'reserve' in componentName.lower():
                        reserveActuators.append(componentName)
                elif 'PointActuator' in className:
                    force_categories['PointActuator'].append((componentName, componentPath))
                    torqueActuators.append(componentName)
                elif 'TorqueActuator' in className:
                    force_categories['TorqueActuator'].append((componentName, componentPath))
                    torqueActuators.append(componentName)
                elif 'SpringGeneralizedForce' in className:
                    force_categories['SpringGeneralizedForce'].append((componentName, componentPath))
                    otherForces.append(componentName)
                elif 'SmoothSphereHalfSpaceForce' in className:
                    force_categories['SmoothSphereHalfSpaceForce'].append((componentName, componentPath))
                    contactForces.append(componentName)
                elif 'HuntCrossleyForce' in className:
                    force_categories['HuntCrossleyForce'].append((componentName, componentPath))
                    contactForces.append(componentName)
                elif 'ExternalForce' in className:
                    force_categories['ExternalForce'].append((componentName, componentPath))
                    externalForces.append(componentName)
                else:
                    force_categories['Other_Force'].append((componentName, componentPath))
                    otherForces.append(componentName)
        
        # Print detailed force breakdown
        print(f"\n📊 FORCE COMPONENT SUMMARY:")
        print(f"   🔄 Total force components: {num_forces}")
        
        print(f"\n💪 MUSCLE TYPES:")
        for muscle_type, components in force_categories.items():
            if 'Muscle' in muscle_type and components:
                print(f"   {muscle_type}: {len(components)}")
                if len(components) <= 5:  # Show details for small lists
                    for name, path in components:
                        print(f"      - {name}")
                else:  # Show first few for large lists
                    for name, path in components[:3]:
                        print(f"      - {name}")
                    print(f"      ... and {len(components) - 3} more")
        
        print(f"\n⚙️  ACTUATOR TYPES:")
        for actuator_type in ['CoordinateActuator', 'PointActuator', 'TorqueActuator']:
            components = force_categories[actuator_type]
            if components:
                print(f"   {actuator_type}: {len(components)}")
                for name, path in components:
                    print(f"      - {name}")
        
        print(f"\n🔗 CONTACT FORCE TYPES:")
        for contact_type in ['SmoothSphereHalfSpaceForce', 'HuntCrossleyForce']:
            components = force_categories[contact_type]
            if components:
                print(f"   {contact_type}: {len(components)}")
                if len(components) <= 10:  # Show details for reasonable lists
                    for name, path in components:
                        print(f"      - {name}")
                else:  # Show first few for large lists
                    for name, path in components[:5]:
                        print(f"      - {name}")
                    print(f"      ... and {len(components) - 5} more")
        
        print(f"\n🌐 OTHER FORCE TYPES:")
        for other_type in ['ExternalForce', 'SpringGeneralizedForce', 'Other_Force']:
            components = force_categories[other_type]
            if components:
                print(f"   {other_type}: {len(components)}")
                for name, path in components:
                    print(f"      - {name}")
        
        # Print categorized lists for backward compatibility
        print(f"\n📋 FORCE CATEGORIES:")
        print(f"   💪 Muscles: {len(muscles)}")
        print(f"   ⚙️  Torque actuators: {len(torqueActuators)}")
        print(f"   🔄 Reserve actuators: {len(reserveActuators)}")
        print(f"   📱 Contact forces: {len(contactForces)}")
        print(f"   🌐 External forces: {len(externalForces)}")
        print(f"   ❓ Other forces: {len(otherForces)}")
        
        if muscles:
            print(f"\n💪 MUSCLES (first 10 of {len(muscles)}):")
            for i, muscle in enumerate(muscles[:10]):
                print(f"   {i+1:2d}: {muscle}")
            if len(muscles) > 10:
                print(f"       ... and {len(muscles) - 10} more")
        
        if torqueActuators:
            print(f"\n⚙️  TORQUE ACTUATORS ({len(torqueActuators)}):")
            for i, actuator in enumerate(torqueActuators):
                print(f"   {i+1:2d}: {actuator}")
        
        if reserveActuators:
            print(f"\n🔄 RESERVE ACTUATORS ({len(reserveActuators)}):")
            for i, actuator in enumerate(reserveActuators):
                print(f"   {i+1:2d}: {actuator}")
        
        if contactForces:
            print(f"\n📱 CONTACT FORCES ({len(contactForces)}):")
            for i, force in enumerate(contactForces):
                print(f"   {i+1:2d}: {force}")
        
        # Summary
        print(f"\n📋 SUMMARY:")
        print(f"   ✅ ModelProcessor created successfully")
        print(f"   📊 Total states: {allStates.getSize()}")
        print(f"   🎯 Coordinates: {len(coordinateValues)} values + {len(coordinateSpeeds)} speeds")
        print(f"   💪 Muscle states: {len(muscleActivations)} activations + {len(muscleTendonForces)} tendon forces")
        print(f"   ⚙️  Actuator states: {len(torqueActivations)} torque + {len(reserveActivations)} reserve")
        print(f"   🔄 Total force components: {num_forces}")
        
        return {
            'total_states': allStates.getSize(),
            'coordinate_values': coordinateValues,
            'coordinate_speeds': coordinateSpeeds,
            'muscle_activations': muscleActivations,
            'tendon_forces': muscleTendonForces,
            'torque_activations': torqueActivations,
            'reserve_activations': reserveActivations,
            'processed_model': processed_model,
            'force_categories': force_categories,
            'num_forces': num_forces
        }
        
    except Exception as e:
        print(f"\n❌ ERROR analyzing ModelProcessor:")
        print(f"   Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return None
    
    finally:
        print(f"="*80)

# =====================================================================================================
def convert_coordinates_to_radians(input_file_path, exclude_columns=None):
    """
    Convert coordinate columns from degrees to radians, excluding specified columns
    
    Parameters:
    -----------
    input_file_path : str
        Path to the input .sto file
    exclude_columns : list
        List of column name substrings to exclude from conversion
        
    Returns:
    --------
    osim.TimeSeriesTable
        Modified table with converted coordinates
    """
    
    if exclude_columns is None:
        exclude_columns = ['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz']
    
    # Load the coordinates table
    # coordinatesTable = osim.TimeSeriesTable(input_file_path)
    
    # Get column labels
    # labels = coordinatesTable.getColumnLabels()

    storage = osim.Storage(input_file_path)
    labels = storage.getColumnLabels()

    # Determine which columns to convert
    for i in range(labels.getSize()):
        label = labels.get(i)
        
        should_exclude = False
        for exclude_col in exclude_columns:
            if exclude_col in label:
                should_exclude = True
                break
        
        if not should_exclude and label != 'time':
            # Get the data column
            state_array = osim.Array()
            storage.getDataColumn(label, state_array)
            
            # Convert each value in the array
            for j in range(state_array.getSize()):
                deg_value = state_array.get(j)
                rad_value = deg_value * (np.pi / 180.0)
                state_array.set(j, rad_value)
            
            # Set the converted data back
            storage.setDataColumn(label, state_array)
    
    # Write to temporary file and reload as TimeSeriesTable
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sto', delete=False) as tmp_file:
        temp_path = tmp_file.name
    
    # Write the modified storage to temp file
    storage.print(temp_path)
    
    # Load as TimeSeriesTable
    coordinatesTable = osim.TimeSeriesTable(temp_path)
    
    # Clean up temp file
    os.unlink(temp_path)
    
    return coordinatesTable

def setup_contact_tracking(problem, grf_file_path, 
                           forceNamesRightFoot=None, forceNamesLeftFoot=None, contact_weight=1000, 
                          normalize_tracking_err=True, projection_type='vector', tracking_enabled=True):
    """
    Setup contact tracking goals for ground reaction forces
    
    Parameters:
    -----------
    problem : osim.MocoProblem
        The Moco problem to add contact tracking to
    grf_file_path : str
        Path to the external loads file (.xml or .mot)
    contact_weight : float
        Weight for contact tracking goal (default: 1000)
    normalize_tracking_err : bool
        Whether to normalize tracking error (default: True)
    projection_type : str
        Type of projection ('vector' or 'plane', default: 'vector')
    
    Returns:
    --------
    osim.MocoProblem
        Problem with contact tracking goal added
    """
    
    print(f"Setting up contact tracking with weight: {contact_weight}")
    
    
    # Contact Tracking
    contact_tracking = osim.MocoContactTrackingGoal('contact', contact_weight)
    contact_tracking.setExternalLoadsFile(grf_file_path)
    right_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesRightFoot, 'ExternalForce_2', ['/bodyset/toes_r'])
    left_contact_group = osim.MocoContactTrackingGoalGroup(forceNamesLeftFoot, 'ExternalForce_1', ['/bodyset/toes_l'])
    contact_tracking.addContactGroup(right_contact_group)
    contact_tracking.addContactGroup(left_contact_group)
    contact_tracking.setNormalizeTrackingError(normalize_tracking_err)
    contact_tracking.setEnabled(tracking_enabled)
        
    # Add the goal to the problem
    problem.addGoal(contact_tracking)
    
    print("Contact tracking goal added successfully")
    print(f"  - Right foot contact spheres: {len(forceNamesRightFoot)}")
    print(f"  - Left foot contact spheres: {len(forceNamesLeftFoot)}")
    
    return problem


def get_model_contact_forces(model_path):
    """
    Inspect the model to identify available contact force names
    
    Parameters:
    -----------
    model_path : str
        Path to the OpenSim model file
        
    Returns:
    --------
    dict
        Dictionary with 'right_foot' and 'left_foot' lists of contact force names
    """
    import opensim as osim
    
    model = osim.Model(model_path)
    model.initSystem()
    
    contact_forces = {'right_foot': [], 'left_foot': [], 'other': []}
    
    forceSet = model.getForceSet()
    print(f"Inspecting {forceSet.getSize()} forces in the model:")
    
    for i in range(forceSet.getSize()):
        force = forceSet.get(i)
        force_name = force.getName()
        force_path = force.getAbsolutePathString()
        
        print(f"  Force {i}: {force_name} ({force.getConcreteClassName()})")
        
        # Categorize contact forces based on naming patterns
        if 'contact' in force_name.lower():
            if '_r' in force_name or 'right' in force_name.lower():
                contact_forces['right_foot'].append(force_path)
            elif '_l' in force_name or 'left' in force_name.lower():
                contact_forces['left_foot'].append(force_path)
            else:
                contact_forces['other'].append(force_path)
    
    print(f"\nFound contact forces:")
    print(f"  Right foot: {len(contact_forces['right_foot'])}")
    print(f"  Left foot: {len(contact_forces['left_foot'])}")
    print(f"  Other: {len(contact_forces['other'])}")
    
    return contact_forces


def setup_contact_tracking_auto(problem, model_path, grf_file_path, 
                               contact_weight=1000, normalize_tracking_err=True):
    """
    Automatically setup contact tracking by inspecting the model for contact forces
    
    Parameters:
    -----------
    problem : osim.MocoProblem
        The Moco problem to add contact tracking to
    model_path : str
        Path to the OpenSim model file
    grf_file_path : str
        Path to the external loads file (.xml or .mot)
    contact_weight : float
        Weight for contact tracking goal (default: 1000)
    normalize_tracking_err : bool
        Whether to normalize tracking error (default: True)
        
    Returns:
    --------
    osim.MocoProblem
        Problem with contact tracking goal added
    """

    
    print("Auto-detecting contact forces from model...")
    
    # Get contact forces from the model
    contact_forces = get_model_contact_forces(model_path)
    
    if not contact_forces['right_foot'] and not contact_forces['left_foot']:
        print("WARNING: No contact forces found in model. Skipping contact tracking.")
        return problem
    
    # Create the contact tracking goal
    contactTracking = osim.MocoContactTrackingGoal('contact_tracking', contact_weight)
    contactTracking.setExternalLoadsFile(grf_file_path)
    contactTracking.setNormalizeTrackingError(normalize_tracking_err)
    
    # Add right foot contact group if forces exist
    if contact_forces['right_foot']:
        rightContactGroup = osim.MocoContactTrackingGoalGroup(
            contact_forces['right_foot'], 'ExternalForce_2', ['/bodyset/toes_r'])
        rightContactGroup.setNormalizeTrackedForceMagnitude(True)
        contactTracking.addContactGroup(rightContactGroup)
        print(f"Added right foot contact group with {len(contact_forces['right_foot'])} forces")
    
    # Add left foot contact group if forces exist
    if contact_forces['left_foot']:
        leftContactGroup = osim.MocoContactTrackingGoalGroup(
            contact_forces['left_foot'], 'ExternalForce_1', ['/bodyset/toes_l'])
        leftContactGroup.setNormalizeTrackedForceMagnitude(True)
        contactTracking.addContactGroup(leftContactGroup)
        print(f"Added left foot contact group with {len(contact_forces['left_foot'])} forces")
    
    # Add the goal to the problem
    problem.addGoal(contactTracking)
    print("Contact tracking goal added successfully")
    
    return problem

# Calculate average walking speed from pelvis_tx motion data
def calculate_average_walking_speed(states_file, initial_time, final_time):
    """
    Calculate average walking speed from pelvis_tx translation data
    """
    import opensim as osim
    import numpy as np
    
    try:
        # Load the states/kinematics file
        kinematics_table = osim.TimeSeriesTable(states_file)
        
        # Get time column
        time_column = kinematics_table.getIndependentColumn()
        time_array = np.array([time_column.get(i) for i in range(time_column.size())])
        
        # Find pelvis_tx column (try different possible naming conventions)
        column_labels = kinematics_table.getColumnLabels()
        pelvis_tx_column = None
        pelvis_tx_name = None
        
        # Try different possible column names for pelvis_tx
        possible_names = [
            'pelvis_tx', 
            '/jointset/ground_pelvis/pelvis_tx/value',
            'ground_pelvis_pelvis_tx',
            'pelvis_tx/value'
        ]
        
        for name in possible_names:
            for i in range(column_labels.getSize()):
                if name in column_labels.get(i):
                    pelvis_tx_column = kinematics_table.getDependentColumn(column_labels.get(i))
                    pelvis_tx_name = column_labels.get(i)
                    break
            if pelvis_tx_column is not None:
                break
        
        if pelvis_tx_column is None:
            print("  WARNING: Could not find pelvis_tx column in states file")
            print(f"  Available columns: {[column_labels.get(i) for i in range(min(5, column_labels.getSize()))]}")
            return 1.25  # Default walking speed
        
        print(f"  Found pelvis_tx column: {pelvis_tx_name}")
        
        # Extract pelvis_tx position data
        pelvis_tx_data = np.array([pelvis_tx_column.get(i) for i in range(pelvis_tx_column.size())])
        
        # Find indices corresponding to the analysis time window
        start_idx = np.argmin(np.abs(time_array - initial_time))
        end_idx = np.argmin(np.abs(time_array - final_time))
        
        # Get time and position data for the analysis window
        analysis_time = time_array[start_idx:end_idx+1]
        analysis_position = pelvis_tx_data[start_idx:end_idx+1]
        
        # Calculate total displacement and time duration
        total_displacement = analysis_position[-1] - analysis_position[0]  # meters
        total_time = analysis_time[-1] - analysis_time[0]  # seconds
        
        # Calculate average speed
        average_speed = total_displacement / total_time
        
        print(f"  Analysis window: {initial_time:.3f} to {final_time:.3f} seconds")
        print(f"  Total displacement: {total_displacement:.3f} meters")
        print(f"  Total time: {total_time:.3f} seconds")
        print(f"  Calculated average speed: {average_speed:.3f} m/s")
        
        # Sanity check: typical walking speeds are 0.5-2.5 m/s
        if average_speed < 0.1 or average_speed > 3.0:
            print(f"  WARNING: Calculated speed ({average_speed:.3f} m/s) seems unrealistic")
            print(f"  Check your pelvis_tx data and time window")
            return max(0.5, min(2.5, abs(average_speed)))  # Clamp to reasonable range
        
        return abs(average_speed)  # Return absolute value
        
    except Exception as e:
        print(f"  ERROR calculating walking speed: {e}")
        print(f"  Using default speed of 1.25 m/s")
        return 1.25

# Calculate coordinate standard deviations and create intelligent weights
def calculate_coordinate_weights(states_file, model_path):
    """
    Calculate coordinate-specific weights based on motion variability and importance
    """
    import opensim as osim
    import numpy as np
    
    try:
        # Load the states/kinematics data
        kinematics_table = osim.TimeSeriesTable(states_file)
        
        # Load model to get coordinate information
        model = osim.Model(model_path)
        model.initSystem()
        coordSet = model.getCoordinateSet()
        
        # Get column labels from the data
        column_labels = kinematics_table.getColumnLabels()
        
        # Initialize weight storage
        coordinate_weights = {}
        coordinate_std_devs = {}
        
        print(f"  Analyzing {column_labels.getSize()} columns for weight calculation...")
        
        # =============================================================================
        # STEP 1: Calculate standard deviations for each coordinate
        # =============================================================================
        
        for i in range(column_labels.getSize()):
            col_name = column_labels.get(i)
            
            # Skip time column
            if 'time' in col_name.lower():
                continue
                
            # Get coordinate data
            try:
                column_data = kinematics_table.getDependentColumn(col_name)
                data_array = np.array([column_data.get(j) for j in range(column_data.size())])
                
                # Calculate standard deviation
                std_dev = np.std(data_array)
                coordinate_std_devs[col_name] = std_dev
                
            except Exception as e:
                print(f"    Warning: Could not process column {col_name}: {e}")
                continue
        
        # =============================================================================
        # STEP 2: Define coordinate importance categories
        # =============================================================================
        
        # Define importance weights based on coordinate type and body segment
        importance_weights = {
            # Pelvis coordinates (very important for overall motion)
            'pelvis_tx': {'value': 1.0, 'speed': 0.01},      # Forward progression
            'pelvis_ty': {'value': 2.0, 'speed': 0.02},      # Vertical motion (critical)
            'pelvis_tz': {'value': 1.5, 'speed': 0.015},     # Lateral motion
            'pelvis_tilt': {'value': 1.0, 'speed': 0.01},    # Pelvis rotations
            'pelvis_list': {'value': 1.0, 'speed': 0.01},
            'pelvis_rotation': {'value': 0.5, 'speed': 0.005},
            
            # Hip coordinates (important for gait)
            'hip_flexion': {'value': 1.0, 'speed': 0.01},
            'hip_adduction': {'value': 1.5, 'speed': 0.015},  # Important for stability
            'hip_rotation': {'value': 0.8, 'speed': 0.008},
            
            # Knee coordinates (critical for walking)
            'knee_angle': {'value': 2.0, 'speed': 0.02},     # Very important
            'knee_adduction': {'value': 1.0, 'speed': 0.01},
            'knee_rotation': {'value': 0.5, 'speed': 0.005},
            
            # Ankle coordinates (important for foot clearance and push-off)
            'ankle_angle': {'value': 2.0, 'speed': 0.02},    # Critical for walking
            # 'subtalar_angle': {'value': 1.0, 'speed': 0.01},
            'mtp_angle': {'value': 0.8, 'speed': 0.008},
            
            # Lumbar coordinates (often problematic, lower weights)
            'lumbar_extension': {'value': 0.1, 'speed': 0.001},  # Often causes issues
            'lumbar_bending': {'value': 0.1, 'speed': 0.001},
            'lumbar_rotation': {'value': 0.1, 'speed': 0.001},
            
            # Arm coordinates (if present, usually less important for walking)
            'arm_flex': {'value': 0.3, 'speed': 0.003},
            'arm_add': {'value': 0.2, 'speed': 0.002},
            'arm_rot': {'value': 0.1, 'speed': 0.001},
            'elbow_flex': {'value': 0.2, 'speed': 0.002},
            'pro_sup': {'value': 0.1, 'speed': 0.001},
            'wrist_flex': {'value': 0.1, 'speed': 0.001},
            'wrist_dev': {'value': 0.1, 'speed': 0.001},
        }
        
        # =============================================================================
        # STEP 3: Calculate final weights using both importance and variability
        # =============================================================================
        
        for col_name, std_dev in coordinate_std_devs.items():
            
            # Extract coordinate base name for lookup
            coord_base_name = None
            weight_type = None  # 'value' or 'speed'
            
            # Determine if this is a position or speed column
            if '/value' in col_name or col_name.endswith('_value'):
                weight_type = 'value'
            elif '/speed' in col_name or col_name.endswith('_speed'):
                weight_type = 'speed'
            else:
                # For raw coordinate names, assume position
                weight_type = 'value'
            
            # Extract base coordinate name
            for importance_key in importance_weights.keys():
                if importance_key in col_name:
                    coord_base_name = importance_key
                    break
            
            # Calculate weight based on importance and variability
            if coord_base_name and coord_base_name in importance_weights:
                importance_weight = importance_weights[coord_base_name][weight_type]
                
                # Weight calculation method similar to tracking_problem.py
                # Higher std_dev = lower weight (more variable = less strictly tracked)
                if std_dev > 0:
                    # Normalize by standard deviation (similar to tracking_problem.py approach)
                    variability_factor = 1.0 / (1.0 + std_dev)
                    final_weight = importance_weight * variability_factor
                else:
                    final_weight = importance_weight
                
                # Apply minimum weight threshold
                final_weight = max(final_weight, 0.001)
                
            else:
                # Default weight for unrecognized coordinates
                final_weight = 0.1 if weight_type == 'value' else 0.001
            
            coordinate_weights[col_name] = final_weight
            
            print(f"    {col_name}: std={std_dev:.4f}, weight={final_weight:.4f}")
        
        print(f"  Calculated weights for {len(coordinate_weights)} coordinates")
        return coordinate_weights
        
    except Exception as e:
        print(f"  ERROR calculating coordinate weights: {e}")
        return {}


def create_complete_states_reference(experimental_states_file, model_path, 
                                     output_file='complete_states_reference.sto',
                                     initial_time=0.81, final_time=1.65):
    """
    Create a complete states reference file with values for all state variables
    not already included in the experimental states file.
    
    Parameters:
    -----------
    experimental_states_file : str
        Path to the experimental states file (.sto) containing kinematic data
    model_path : str
        Path to the OpenSim model file (.osim)
    output_file : str
        Path for the output complete states reference file (default: 'complete_states_reference.sto')
    initial_time : float
        Start time for the reference (default: 0.81)
    final_time : float
        End time for the reference (default: 1.65)
        
    Returns:
    --------
    str
        Path to the created complete states reference file
    """
    import opensim as osim
    import numpy as np
    import os
    
    print(f"Creating complete states reference from experimental data...")
    print(f"  Experimental states: {experimental_states_file}")
    print(f"  Model: {model_path}")
    print(f"  Output: {output_file}")
    
    try:
        # Load the experimental states data
        experimentalStates = osim.TimeSeriesTable(experimental_states_file)
        
        # Load and initialize the model
        model = osim.Model(model_path)
        model.initSystem()
        
        # Get time data from experimental states
        experimentalTime = experimentalStates.getIndependentColumn()
        
        # Create time vector for the reference (can be different resolution)
        nrows = 100  # Number of time points for reference
        timeArray = np.linspace(initial_time, final_time, num=nrows)
        timeVector = osim.StdVectorDouble()
        for time_val in timeArray:
            timeVector.push_back(time_val)
        
        # Create the complete states reference table
        completeStatesRef = osim.TimeSeriesTable(timeVector)
        
        # Get all possible state variable names from the model
        allStateVars = model.getStateVariableNames()
        print(f"  Model has {allStateVars.getSize()} total state variables")
        
        # Get experimental state column labels
        expColumnLabels = experimentalStates.getColumnLabels()
        expStateNames = set()
        for i in range(expColumnLabels.getSize()):
            expStateNames.add(expColumnLabels.get(i))
        
        print(f"  Experimental data has {len(expStateNames)} state variables")
        
        # Track which states we're adding
        added_from_experimental = 0
        added_missing_coordinates = 0
        added_missing_activations = 0
        added_missing_tendon_forces = 0
        added_other_states = 0
        
        # Process each state variable from the model
        for i in range(allStateVars.getSize()):
            stateVar = allStateVars.get(i)
            
            if stateVar in expStateNames:
                # State exists in experimental data - interpolate from experimental data
                try:
                    # Get experimental data for this state
                    expColumn = experimentalStates.getDependentColumn(stateVar)
                    expTimeArray = np.array([experimentalTime.get(j) for j in range(experimentalTime.size())])
                    expDataArray = np.array([expColumn.get(j) for j in range(expColumn.size())])
                    
                    # Interpolate to new time points
                    interpolatedData = np.interp(timeArray, expTimeArray, expDataArray)
                    
                    # Add to reference table
                    dataVector = osim.Vector(nrows, 0.0)
                    for j, val in enumerate(interpolatedData):
                        dataVector.set(j, val)
                    
                    completeStatesRef.appendColumn(stateVar, dataVector)
                    added_from_experimental += 1
                    
                except Exception as e:
                    print(f"    Warning: Could not interpolate experimental data for {stateVar}: {e}")
                    # Fall back to default value
                    dataVector = osim.Vector(nrows, 0.0)
                    completeStatesRef.appendColumn(stateVar, dataVector)
                    added_other_states += 1
            
            else:
                # State missing from experimental data - create reference values
                reference_value = 0.0  # Default value
                
                # Determine appropriate reference value based on state type
                if 'activation' in stateVar.lower():
                    reference_value = 0.02  # Small baseline muscle activation
                    category = 'activation'
                    added_missing_activations += 1
                    
                elif 'tendon_force' in stateVar.lower() or 'normalized_tendon_force' in stateVar.lower():
                    reference_value = 0.02  # Small baseline tendon force
                    category = 'tendon_force'
                    added_missing_tendon_forces += 1
                    
                elif any(coord_type in stateVar.lower() for coord_type in ['/value', '/speed', 'coordinate']):
                    reference_value = 0.0  # Zero position/velocity for missing coordinates
                    category = 'coordinate'
                    added_missing_coordinates += 1
                    
                else:
                    reference_value = 0.0  # Default for other state types
                    category = 'other'
                    added_other_states += 1
                
                # Create constant reference data
                dataVector = osim.Vector(nrows, reference_value)
                completeStatesRef.appendColumn(stateVar, dataVector)
                
                print(f"    Added missing {category}: {stateVar} = {reference_value}")
        
        # Set table metadata
        completeStatesRef.addTableMetaDataString('inDegrees', 'no')  # All values in radians
        completeStatesRef.addTableMetaDataString('DataType', 'States')
        completeStatesRef.addTableMetaDataString('source', f'Generated from {experimental_states_file}')
        
        # Write the complete reference file
        osim.STOFileAdapter.write(completeStatesRef, output_file)
        
        # Print summary
        print(f"\n✓ Complete states reference created successfully:")
        print(f"    Total state variables: {allStateVars.getSize()}")
        print(f"    From experimental data: {added_from_experimental}")
        print(f"    Missing coordinates (set to 0.0): {added_missing_coordinates}")
        print(f"    Missing activations (set to 0.02): {added_missing_activations}")
        print(f"    Missing tendon forces (set to 0.02): {added_missing_tendon_forces}")
        print(f"    Other missing states (set to 0.0): {added_other_states}")
        print(f"    Output file: {output_file}")
        print(f"    Time range: {initial_time:.3f} to {final_time:.3f} seconds")
        print(f"    Time points: {nrows}")
        print(f"    Units: Radians (inDegrees=no)")
        
        return output_file
        
    except Exception as e:
        print(f"✗ Error creating complete states reference: {e}")
        return None

# def write_torso_tracking_reference()

def create_model_processor_base(model_path=path_model, 
                               external_loads_path=None,
                               # Muscle modeling options
                               use_degroote_muscles=True,
                               ignore_tendon_compliance=False,
                               ignore_passive_fiber_forces=True,
                               active_force_width_scale=1.5,
                               
                               # Reserve actuator options
                               reserve_strength=50.0,
                               reserve_coordinates=None,
                               
                               # Torque actuator options  
                               add_torque_actuators=False,
                               torque_actuator_coords=None,
                               torque_actuator_strengths=None,
                               
                               # Advanced options
                               implicit_tendon_dynamics=True,
                               lumbar_stiffness_scale=1.0,
                               
                               # Logging and debugging
                               verbose=True):
    """
    Create a ModelProcessor with comprehensive model modifications for optimization.
    This is a standalone version similar to create_model_processor_base from result.py
    
    Parameters:
    -----------
    model_path : str
        Path to the base OpenSim model file (.osim)
    external_loads_path : str, optional
        Path to external loads file (.xml) for ground reaction forces
    use_degroote_muscles : bool
        Replace muscles with DeGrooteFregly2016 muscles (default: True)
    ignore_tendon_compliance : bool
        Use rigid tendon assumption (default: True)
    ignore_passive_fiber_forces : bool
        Ignore passive muscle fiber forces (default: True)
    active_force_width_scale : float
        Scale factor for active force curve width (default: 1.5)
    reserve_strength : float
        Strength for reserve actuators (0 = no reserves) (default: 0.0)
    reserve_coordinates : list, optional
        List of coordinates to add reserve actuators to (default: None = all)
    add_torque_actuators : bool
        Add coordinate actuators for specified joints (default: False)
    torque_actuator_coords : list, optional
        List of coordinates for torque actuators (default: None)
    torque_actuator_strengths : dict, optional
        Dictionary of {coordinate: strength} for torque actuators (default: None)
    implicit_tendon_dynamics : bool
        Use implicit tendon compliance dynamics (default: False)
    lumbar_stiffness_scale : float
        Scale lumbar joint stiffness (default: 1.0)
    verbose : bool
        Print detailed processing information (default: True)
        
    Returns:
    --------
    osim.ModelProcessor
        Configured ModelProcessor ready for use in tracking/optimization
    """
    
    if verbose:
        print(f"\nCreating ModelProcessor from: {model_path}")
        osim.Logger.setLevelString('info')
    else:
        osim.Logger.setLevelString('error')
    
    # Initialize the ModelProcessor with the base model
    modelProcessor = osim.ModelProcessor(model_path)
    model = osim.Model(model_path)
    state = model.initSystem()
    mass = 60

    # # =============================================================================
    # # EXTERNAL LOADS (Ground Reaction Forces)
    # # =============================================================================
    
    if external_loads_path:
        if verbose:
            print(f"  Adding external loads: {external_loads_path}")
        modelProcessor.append(osim.ModOpAddExternalLoads(external_loads_path))
        
    
    # =============================================================================
    # TORQUE ACTUATORS
    # =============================================================================
    
    # Add in arm actuators
    if add_torque_actuators: #and torque_actuator_coords:
        if verbose:
            print("  Adding coordinate actuators (torque actuators)")
        
        coordNames = ['arm_flex', 'arm_add', 'arm_rot', 'elbow_flex', 'pro_sup']
        strengths = [0.5, 0.5, 0.5, 0.2, 0.2]
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

        if lumbar_stiffness_scale != 1.0:
            if verbose:
                print(f"  Scaling lumbar joint stiffness by {lumbar_stiffness_scale}")

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
        
        if reserve_strength > 0:
            if verbose:
                print(f"  Adding reserve actuators with strength: {reserve_strength}")
            
            coordNames = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz', 
                        'pelvis_list', 'pelvis_tilt', 'pelvis_rotation',
                        'hip_adduction_r', 'hip_rotation_r', 'hip_flexion_r',
                        'hip_adduction_l', 'hip_rotation_l', 'hip_flexion_l',
                        'knee_angle_r', 'ankle_angle_r', 
                        'knee_angle_l', 'ankle_angle_l', 
                        'mtp_angle_r', 'mtp_angle_l']
                        # 'subtalar_angle_r','subtalar_angle_l',
            for coordName in coordNames:
                actu = osim.CoordinateActuator()
                actu.set_coordinate(coordName)
                actu.setName(f'reserve_{coordName}')
                scale = 1.0
                if ((coordName == 'pelvis_tx') or
                        (coordName == 'pelvis_ty') or
                        (coordName == 'pelvis_tz')):
                    scale = 10.0 
                actu.setOptimalForce(scale * reserve_strength)
                actu.setMinControl(-1.0)
                actu.setMaxControl(1.0)
                model.addForce(actu)

            # Add reserve actuators to help with convergence
            # This requires processing the model first to get coordinate information
            # temp_model = modelProcessor.process()
            # temp_model.initSystem()
            
            # coordSet = temp_model.getCoordinateSet()
            
            # # Determine which coordinates get reserve actuators
            # if reserve_coordinates is None:
            #     # Add to all coordinates by default
            #     coords_for_reserves = [coordSet.get(i).getName() for i in range(coordSet.getSize())]
            # else:
            #     coords_for_reserves = reserve_coordinates
            
            # # Create reserve actuators for specified coordinates
            # for coord_name in coords_for_reserves:
            #     try:
            #         coord = coordSet.get(coord_name)
                    
            #         # Create a coordinate actuator (reserve)
            #         reserve = osim.CoordinateActuator()
            #         reserve.setName(f'reserve_{coord_name}')
            #         reserve.setCoordinate(coord)
            #         reserve.setOptimalForce(reserve_strength)
            #         reserve.setMinControl(-1.0)
            #         reserve.setMaxControl(1.0)
                    
            #         # Add to model via ModOp (this is simplified - actual implementation may vary)
            #         if verbose:
            #             print(f"    Added reserve actuator for {coord_name}")
                        
            #     except Exception as e:
            #         if verbose:
            #             print(f"    Warning: Could not add reserve for {coord_name}: {e}")
        
    # =============================================================================
    # FINALIZATION
    # =============================================================================

    model.finalizeConnections()

    modelProcessor = osim.ModelProcessor(model)

    # weld_joints = list()
    # for side in ['_l', '_r']:
    #     weld_joints.append(f"subtalar{side}")
    # modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weld_joints))
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))  

    
    # Print summary of modifications
    # temp_model = modelProcessor.process()
    # temp_model.initSystem()
    
    # print(f"  Model summary:")
    # print(f"    Coordinates: {modelProcessor.getNumCoordinates()}")
    # print(f"    Muscles: {modelProcessor.getMuscles().getSize()}")
    # print(f"    Forces: {modelProcessor.getForceSet().getSize()}")
    # print(f"    Bodies: {modelProcessor.getBodySet().getSize()}")



    # =============================================================================
    # MUSCLE MODELING MODIFICATIONS
    # =============================================================================

    # #Enable tendon compliance for ankle plantarflexors
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
    if implicit_tendon_dynamics:
        modelProcessor.append(
            osim.ModOpUseImplicitTendonComplianceDynamicsDGF())
        
    osim.Logger.setLevelString('info')

    return modelProcessor

    # if use_degroote_muscles:
    #     if verbose:
    #         print("  Replacing muscles with DeGrooteFregly2016 muscles")
        
    #     # Replace default muscles with optimization-friendly muscles
    #     modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
        
    #     # Configure DeGrooteFregly2016 muscle properties
    #     if ignore_passive_fiber_forces:
    #         if verbose:
    #             print("    Ignoring passive fiber forces")
    #         modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
        
    #     if active_force_width_scale != 1.0:
    #         if verbose:
    #             print(f"    Scaling active force curve width by {active_force_width_scale}")
    #         modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(active_force_width_scale))
    
    # # Tendon compliance settings
    # if ignore_tendon_compliance and not implicit_tendon_dynamics:
    #     if verbose:
    #         print("  Using rigid tendon assumption")
    #     modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    # elif implicit_tendon_dynamics:
    #     if verbose:
    #         print("  Using implicit tendon compliance dynamics")
    #     # Note: Implicit tendon dynamics may require additional setup
    
def set_state_bounds(state_bounds):

    # Define angular radian & velocity values for bounds (in radians for angles, meters for translations)
    pi = np.pi
    state_bounds = {
        # Pelvis translations (meters)
        '/jointset/ground_pelvis/pelvis_tx/value': ([-5, 10], [], []), 
        '/jointset/ground_pelvis/pelvis_ty/value': ([0, 1.5],[], []),
        '/jointset/ground_pelvis/pelvis_tz/value': ([-2, 2], [], []),
        
        # Pelvis rotations (radians)
        '/jointset/ground_pelvis/pelvis_tilt/value': ([-90*pi/180, 90*pi/180], [], []),
        '/jointset/ground_pelvis/pelvis_list/value': ([-90*pi/180, 90*pi/180], [], []),
        '/jointset/ground_pelvis/pelvis_rotation/value': ([-90*pi/180, 90*pi/180], [], []),
        
        # Hip joints (radians)
        '/jointset/hip_l/hip_flexion_l/value': ([-30*pi/180, 120*pi/180], [], []),
        '/jointset/hip_r/hip_flexion_r/value': ([-30*pi/180, 120*pi/180], [], []),
        '/jointset/hip_l/hip_adduction_l/value': ([-30*pi/180, 30*pi/180], [], []),
        '/jointset/hip_r/hip_adduction_r/value': ([-30*pi/180, 30*pi/180], [], []),
        '/jointset/hip_l/hip_rotation_l/value': ([-45*pi/180, 45*pi/180], [], []),
        '/jointset/hip_r/hip_rotation_r/value': ([-45*pi/180, 45*pi/180], [], []),
        
        # Knee joints (radians)
        '/jointset/walker_knee_l/knee_angle_l/value': ([0, 120*pi/180], [], []),
        '/jointset/walker_knee_r/knee_angle_r/value': ([0, 120*pi/180], [], []),
        
        # Ankle joints (radians)
        '/jointset/ankle_l/ankle_angle_l/value': ([-40*pi/180, 30*pi/180], [], []),
        '/jointset/ankle_r/ankle_angle_r/value': ([-40*pi/180, 30*pi/180], [], []),
        
        # MTP joints (radians) - if present
        '/jointset/mtp_l/mtp_angle_l/value': ([-30*pi/180, 60*pi/180], [], []),
        '/jointset/mtp_r/mtp_angle_r/value': ([-30*pi/180, 60*pi/180], [], []),
        
        # Lumbar joint (radians) - if present
        '/jointset/back/lumbar_extension/value': ([-30*pi/180, 30*pi/180], [], []),
        '/jointset/back/lumbar_bending/value': ([-30*pi/180, 30*pi/180], [], []),
        '/jointset/back/lumbar_rotation/value': ([-30*pi/180, 30*pi/180], [], []),
        
        # Pelvis translation speeds (m/s)
        '/jointset/ground_pelvis/pelvis_tx/speed': ([-10, 10], [], []),  # Forward walking speed
        '/jointset/ground_pelvis/pelvis_ty/speed': ([-1, 1], [], []),      # Vertical speed
        '/jointset/ground_pelvis/pelvis_tz/speed': ([-2, 2], [], []),     # Lateral speed
    }
    return state_bounds

def configure_objective_function_weights(weights):
    scale = 1.0
    weights = {
    'state_tracking_weight':    50 * scale,
    'control_weight':           25 * scale,
    'grf_tracking_weight':      7500 * scale,
    'torso_orientation_weight': 10 * scale,
    'feet_orientation_weight':  10 * scale,
    'control_tracking_weight':  0 * scale, 
    'aux_deriv_weight':         1000 * scale,
    'acceleration_weight':      1 * scale,
    }
    return weights
    
    # # Apply kinematic bounds
    # for state_path, bounds in kinematic_bounds.items():
    #     try:
    #         problem.setStateInfo(state_path, bounds)
    #         print(f"Set bounds for {state_path}: {bounds}")
    #     except Exception as e:
    #         print(f"Could not set bounds for {state_path}: {e}")
    
def muscleDrivenStateTracking(enable_contact_tracking=True, contact_weight=100, 
                              enable_advanced_solver=True, mesh_interval=0.01, 
                              enable_state_tracking=True,
                              target_walking_speed=2.0,
                              enable_speed_goal=False,
                              speed_goal_weight=10, 
                              enable_coordinate_weights=True, 
                              use_custom_weights=False,
                              custom_coordinate_weights=None,
                              inital_time=0.01,
                              final_time=0.752):


    log_file = setup_logging(log_file_path='test_log.txt',
                             also_print_to_console=True)


    # Construct a ModelProcessor and set it on the tool. The default
    # muscles in the model are replaced with optimization-friendly
    # DeGrooteFregly2016Muscles, and adjustments are made to the default muscle
    # parameters.
    modelProcessor = create_model_processor_base(path_model, 
                               external_loads_path=path_grf,
                               # Muscle modeling options
                               use_degroote_muscles=True,
                               ignore_tendon_compliance=False,
                               ignore_passive_fiber_forces=True,
                               active_force_width_scale=1.5,
                               
                               # Reserve actuator options
                               reserve_strength=50.0,
                               reserve_coordinates=None,
                               
                               # Torque actuator options  
                               add_torque_actuators=False,
                               torque_actuator_coords=None,
                               torque_actuator_strengths=None,
                               
                               # Advanced options
                               implicit_tendon_dynamics=True,
                               lumbar_stiffness_scale=1.0,
                               
                               # Logging and debugging
                               verbose=True)


    model = modelProcessor.process()
    model.initSystem()
    
    # modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))

    # modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    # modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    # Only valid for DeGrooteFregly2016Muscles.
    # modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    # Only valid for DeGrooteFregly2016Muscles.
    # modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    
    #^ After establishing object for model, count the amount of actuator/force objects present; this will be used to normalize control effort cost
    num_forces = 0
    for actu in model.getComponentsList():
        if (actu.getConcreteClassName().endswith('Muscle') or 
            actu.getConcreteClassName().endswith('Actuator')):
            num_forces += 1

    # Take this moment to convert the coordinates into radians and dump them into a temporary file
    # --------------------------------------------------------------------------------------------
    temp_states_file = 'temp_states_converted.sto'
    state_variances = 'state_standard_deviations.csv'

    kinematicsToStates(kinematicsFileName=path_states,
                       osimModelFileName=path_model,
                       outputFileName=temp_states_file,
                       inDegrees=True,
                       outDegrees=False,
                       excludeColumns=['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz'])

    compute_sto_stdevs(stoFileName=temp_states_file, outputCsvFileName=state_variances)


    # generate_guess_file(kinematicsFileName=path_states, osimModelFileName=path_model,
    #                 outputFileName='guess_file.sto',
    #                 muscle_activation_bounds=muscle_bounds,
    #                 normalized_tendon_force_bounds=tendon_bounds,
    #                 torque_actuator_bounds=torque_bounds,
    #                 inDegrees=True, filtFreq=None)
    
    expiremental_states = osim.TableProcessor('sub03_strifc_cond00000_speed02000_0001_0003_td_marks.sto')

    # Now create coordinates for tracking reference
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    exp_state_filter = expiremental_states.process()
    coordinatesTable = osim.TimeSeriesTable(exp_state_filter)
    # time = coordinatesTable.getIndependentColumn() #??'time' argument as default?
    # nrows = coordinatesTable.getNumRows()
    # pelvis_tx = coordinatesTable.getDependentColumn('/jointset/ground_pelvis/pelvis_tx/value')
    # pelvis_tx_new = osim.Vector(nrows, 0)
    # for i in np.arange(nrows):
    #     dx = (time[int(i)] - time[0]) * target_walking_speed
    #     pelvis_tx_new[int(i)] = pelvis_tx[int(i)] + dx

    # pelvis_tx_upd = coordinatesTable.updDependentColumn('/jointset/ground_pelvis/pelvis_tx/value')
    # for i in np.arange(nrows):
    #     pelvis_tx_upd[int(i)] = pelvis_tx_new[int(i)]
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    # for side in ['_l', 'r']:
    #     coordinatesTable.removeColumn(f'wrist_dev{side}')
    #     coordinatesTable.removeColumn(f'wrist_flex{side}')
    tableProcessor = osim.TableProcessor(coordinatesTable)
    tableProcessor.append(osim.TabOpLowPassFilter(6))
    tableProcessor.append(osim.TabOpUseAbsoluteStateNames())


    every_state_label = coordinatesTable.getColumnLabels()
    for label in every_state_label:
        if '/forceset' in label and '/activation' not in label:
            coordinatesTable.removeColumn(label)
            print(f"     Removed: {label}")
        if '/subtalar' in label:
            coordinatesTable.removeColumn(label)
            print(f" Removed: {label}")


    #     # Pelvis translations (meters)
    #     '/jointset/ground_pelvis/pelvis_tx/value': [-5, 10],
    #     '/jointset/ground_pelvis/pelvis_ty/value': [0, 1.5],
    #     '/jointset/ground_pelvis/pelvis_tz/value': [-2, 2],
        
    #     # Pelvis rotations (radians)
    #     '/jointset/ground_pelvis/pelvis_tilt/value': [-90*pi/180, 90*pi/180],
    #     '/jointset/ground_pelvis/pelvis_list/value': [-90*pi/180, 90*pi/180],
    #     '/jointset/ground_pelvis/pelvis_rotation/value': [-90*pi/180, 90*pi/180],
        
    #     # Hip joints (radians)
    #     '/jointset/hip_l/hip_flexion_l/value': [-30*pi/180, 120*pi/180],
    #     '/jointset/hip_r/hip_flexion_r/value': [-30*pi/180, 120*pi/180],
    #     '/jointset/hip_l/hip_adduction_l/value': [-30*pi/180, 30*pi/180],
    #     '/jointset/hip_r/hip_adduction_r/value': [-30*pi/180, 30*pi/180],
    #     '/jointset/hip_l/hip_rotation_l/value': [-45*pi/180, 45*pi/180],
    #     '/jointset/hip_r/hip_rotation_r/value': [-45*pi/180, 45*pi/180],
        
    #     # Knee joints (radians)
    #     '/jointset/walker_knee_l/knee_angle_l/value': [0, 120*pi/180],
    #     '/jointset/walker_knee_r/knee_angle_r/value': [0, 120*pi/180],
        
    #     # Ankle joints (radians)
    #     '/jointset/ankle_l/ankle_angle_l/value': [-40*pi/180, 30*pi/180],
    #     '/jointset/ankle_r/ankle_angle_r/value': [-40*pi/180, 30*pi/180],
        
    #     # MTP joints (radians) - if present
    #     '/jointset/mtp_l/mtp_angle_l/value': [-30*pi/180, 60*pi/180],
    #     '/jointset/mtp_r/mtp_angle_r/value': [-30*pi/180, 60*pi/180],
        
    #     # Lumbar joint (radians) - if present
    #     '/jointset/back/lumbar_extension/value': [-30*pi/180, 30*pi/180],
    #     '/jointset/back/lumbar_bending/value': [-30*pi/180, 30*pi/180],
    #     '/jointset/back/lumbar_rotation/value': [-30*pi/180, 30*pi/180],
        
    #     # Pelvis translation speeds (m/s)
    #     '/jointset/ground_pelvis/pelvis_tx/speed': [0, 3],  # Forward walking speed
    #     '/jointset/ground_pelvis/pelvis_ty/speed': [-1, 1],      # Vertical speed
    #     '/jointset/ground_pelvis/pelvis_tz/speed': [-2, 2],      # Lateral speed
    # }


    # Now get your predetermined state bounds
    state_bounds = {}
    state_borders = set_state_bounds(state_bounds)

    weight_set = {}
    weights = configure_objective_function_weights(weights=weight_set)
    # Create state tracking weights from scratch
    #-------------------------------------------
    model = modelProcessor.process()
    model.initSystem()
    model_coordinates = tableProcessor.process(model)
    paths = model_coordinates.getColumnLabels()

    ##  STILL HAVE TO GENERATE CSV STDEVS FILE IN COORDINATE CONVERSION
    coordinates_std = pd.read_csv(state_variances, index_col=0)
    state_weights = osim.MocoWeightSet()
    for valuePath in paths:
        speed_path = valuePath.replace('/value', '/speed')
        value_weight = 1.0
        speed_weight = 1.0

        for name in coordinates_std.index:
            stdev = coordinates_std.loc[name][0]
            denom = 10 * stdev

            if name in valuePath:
                if 'lumbar' in name:
                    if 'torso_orientation_weight' in weights:
                        value_weight = 0.0/denom
                        speed_weight = 0.00/denom
                    else:
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
                        value_weight = 1.0/denom
                        speed_weight = 0.1/denom
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

    
    
    # Continue onto initializing the base of the tracking problem
    # ===========================================================
    # Calculate coordinate-specific weights
    # coordinate_weights = calculate_coordinate_weights(temp_states_file, path_model)
    # tableProcessor = osim.TableProcessor(exp_states_with_references)
    # tableProcessor.append(osim.TabOpUseAbsoluteStateNames())

    # Create and name the basic instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName('single_run_track')
    track.setModel(modelProcessor)
    track.setStatesReference(osim.TableProcessor(coordinatesTable))
    track.set_states_global_tracking_weight(weights['state_tracking_weight'] / model.getNumCoordinates())
    track.set_allow_unused_references(True)
    track.set_apply_tracked_states_to_guess(True)
    track.set_states_weight_set(state_weights)

    debug_model_processor(modelProcessor, path_model)

    # This setting allows extra data columns contained in the states
    # reference that don't correspond to model coordinates.
    

    # Since there is only coordinate position data in the states references,
    # this setting is enabled to fill in the missing coordinate speed data using
    # the derivative of splined position data.
    track.set_track_reference_position_derivatives(True)
    track.set_control_effort_weight(weights['control_weight'] / num_forces)

    # Initial time, final time, and mesh interval.
    time = coordinatesTable.getIndependentColumn()
    if inital_time is not None:
        track.set_initial_time(inital_time)
    else:
        track.set_initial_time(time[0])
    if final_time is not None:
        track.set_final_time(final_time)
    else:
        track.set_final_time(time[-1])
    
    track.set_mesh_interval(mesh_interval)

    # Instead of calling solve(), call initialize() to receive a pre-configured
    # MocoStudy object based on the template settings. Use this to customize the
    # problem beyond the MocoTrack interface.
    study = track.initialize()
    # Get a reference to the MocoControlCost that is added to every MocoTrack
    # problem by default.
    problem = study.updProblem()

    # Now proceed with customizing the base tracking problem
    #**EXPERMINENTAL DATA NEEDS TO COME FROM ANOTHER GUESS SOLUTION??
    # expiremental_states = osim.TableProcessor('temp_states_converted.sto')

    # State bound setting
    # =============================================================================
    # Speed, muscle, torques, and tendons bounds
    # =============================================================================

    
    speed_bounds = [-20.0, 20.0]  # rad/s or m/s
    muscle_bounds = [0.01, 1.0]  # Muscle activations between 1% and 100%
    torque_bounds = [-1.0, 1.0]  # Normalized torque bounds
    tendon_bounds = [0.01, 1.8]

    # if '/value' in every_state:
    problem.setStateInfoPattern('/jointset/.*/speed', speed_bounds, [], [])
    problem.setControlInfoPattern('/forceset/.*', muscle_bounds, [], [])
    problem.setControlInfoPattern('/forceset/torque.*', torque_bounds, [], [])
    problem.setControlInfoPattern('/forceset/reserve.*', torque_bounds, [], [])
    problem.setStateInfoPattern('.*/activation', muscle_bounds, [], [])
    problem.setStateInfoPattern('.*torque.*/activation', torque_bounds, [], [])
    problem.setStateInfoPattern('.*reserve.*/activation', torque_bounds, [], [])
    problem.setStateInfoPattern('.*tendon_force', tendon_bounds, [], [])
    
    # Add these kinematic state bounds inot the problem
    for path in state_borders:
        bounds = state_borders[path][0]
        init_bounds = state_borders[path][1] if len(state_borders[path]) > 1 else []
        final_bounds = state_borders[path][2] if len(state_borders[path]) > 2 else []
        problem.setStateInfo(path, bounds, init_bounds, final_bounds)
    
    # Modify/update state tracking into problem
    state_tracking = osim.MocoStateTrackingGoal().safeDownCast(
        problem.updGoal('state_tracking'))
    state_tracking.setEnabled(enable_state_tracking and bool(weights['state_tracking_weight']))

    # Update control effort goal
    control_effort = osim.MocoControlGoal().safeDownCast(
        problem.updGoal('control_effort'))
    control_effort.setDivideByDisplacement(True)

    # # Average Speed Goal
    # average_speed = osim.MocoAverageSpeedGoal()
    # average_speed.setName('average_speed')
    # average_speed.set_desired_average_speed(target_walking_speed)
    
    # problem.addGoal(average_speed)

    # Torso Orientation Goal
    if 'torso_orientation_weight' in weights:

        # Set all activation values to .02 in the prescribed exp states file, and appending them if the actuator columns (almost all activation) don't exist then it just appends
        exp_state_table = osim.TableProcessor(temp_states_file).process()
        exp_state_labels = exp_state_table.getColumnLabels()
        torso_table = osim.TimeSeriesTable(exp_state_table)
        state_vars = model.getStateVariableNames()
        activation_vec = osim.Vector(torso_table.getNumRows(), 0.02)
        for isv in np.arange(state_vars.getSize()):
            state_var = state_vars.get(int(isv))
            if not state_var in exp_state_labels:
                torso_table.appendColumn(state_var, activation_vec)

        # torso_table.removeColumn('/forceset/addbrev_r')
        torso_orientation_goal = osim.MocoOrientationTrackingGoal('torso_orientation_goal', 
                                                                  weights['torso_orientation_weight'])
        torso_orientation_goal.setStatesReference(osim.TableProcessor(torso_table))
        paths = osim.StdVectorString()
        paths.push_back('/bodyset/torso')
        torso_orientation_goal.setFramePaths(paths)
        torso_orientation_goal.setEnabled(enable_state_tracking)
        
        problem.addGoal(torso_orientation_goal)

        torso_speed_goal = osim.MocoAngularVelocityTrackingGoal('torso_angular_velocity_goal', 
                                                                0.1 * weights['torso_orientation_weight'])
        torso_speed_goal.setStatesReference(osim.TableProcessor(torso_table))
        paths = osim.StdVectorString()
        paths.push_back('/bodyset/torso')
        torso_speed_goal.setFramePaths(paths)
        torso_speed_goal.setEnabled(enable_state_tracking)

        problem.addGoal(torso_speed_goal)

    coordinatesTable.removeColumn('lambda_cid35_p0')
    coordinatesTable.removeColumn('lambda_cid36_p0')
    coordinatesTable.removeColumn('gamma_cid35_p0')
    coordinatesTable.removeColumn('gamma_cid36_p0')
    
    # Foot Orientation Goal!
    if 'feet_orientation_weight' in weights:
        exp_state_table = expiremental_states.process()
        exp_state_labels = exp_state_table.getColumnLabels()
        foot_table = coordinatesTable
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
        foot_orientation_goal.setEnabled(enable_state_tracking)

        problem.addGoal(foot_orientation_goal)

        foot_speed_goal = osim.MocoAngularVelocityTrackingGoal('feet_speed_goal', 
                                                                 0.01 * weights['feet_orientation_weight'])
        foot_speed_goal.setStatesReference(osim.TableProcessor(foot_table))
        paths = osim.StdVectorString()
        paths.push_back('/bodyset/calcn_r')
        paths.push_back('/bodyset/calcn_l')
        foot_speed_goal.setFramePaths(paths)
        foot_speed_goal.setEnabled(enable_state_tracking)

        problem.addGoal(foot_speed_goal)        

    # # Distance constraint to prevent intersecting bodies
    # # --------------------------------------------------
    # distanceConstraint = osim.MocoFrameDistanceConstraint()
    # distanceConstraint.setName('distance_constraint')
    # foot_separation = 0.05
    
    # distanceConstraint.addFramePair(
    #         osim.MocoFrameDistanceConstraintPair(
    #         '/bodyset/calcn_l', '/bodyset/calcn_r', foot_separation, 1))
    # distanceConstraint.addFramePair(
    #         osim.MocoFrameDistanceConstraintPair(
    #         '/bodyset/toes_l', '/bodyset/toes_r', foot_separation, 1))
    # distanceConstraint.addFramePair(
    #         osim.MocoFrameDistanceConstraintPair(
    #         '/bodyset/calcn_l', '/bodyset/toes_r', foot_separation, 1))
    # distanceConstraint.addFramePair(
    #         osim.MocoFrameDistanceConstraintPair(
    #         '/bodyset/toes_l', '/bodyset/calcn_r', foot_separation, 1))

    # armDistance = 0.05
    # for body in ['humerus', 'ulna', 'radius', 'hand']:
    #     for side in ['_l', '_r']:
    #         distanceConstraint.addFramePair(
    #                 osim.MocoFrameDistanceConstraintPair(
    #                 '/bodyset/torso', f'/bodyset/{body}{side}', 
    #                 armDistance, np.inf))

    # distanceConstraint.setProjection('vector')
    # distanceConstraint.setProjectionVector(osim.Vec3(0, 0, 1))
    # problem.addPathConstraint(distanceConstraint)


    # # Define contact sphere force names for right foot
    # # Adjust these based on your specific model's contact sphere names
    # forceNamesRightFoot = [
    #     'forceset/contactHeel_r',
    #     'forceset/contactLateralRearfoot_r', 
    #     'forceset/contactLateralMidfoot_r',
    #     'forceset/contactLateralToe_r',
    #     'forceset/contactMedialToe_r',
    #     'forceset/contactMedialMidfoot_r'
    # ]
    
    # # Define contact sphere force names for left foot
    # forceNamesLeftFoot = [
    #     'forceset/contactHeel_l',
    #     'forceset/contactLateralRearfoot_l',
    #     'forceset/contactLateralMidfoot_l', 
    #     'forceset/contactLateralToe_l',
    #     'forceset/contactMedialToe_l',
    #     'forceset/contactMedialMidfoot_l'
    # ]

    # # Add contact tracking goal from function above
    # setup_contact_tracking(problem=problem, grf_file_path=path_grf, 
    #                        forceNamesRightFoot=forceNamesRightFoot, forceNamesLeftFoot=forceNamesLeftFoot,
    #                         contact_weight=weights['grf_tracking_weight'], 
    #                        normalize_tracking_err=True, tracking_enabled=enable_contact_tracking)

    # # exp_states_with_references = create_complete_states_reference(
    # #     experimental_states_file=temp_states_file,
    # #     model_path=path_model,
    # #     output_file='experimental_states_reference.sto',
    # #     initial_time=inital_time,
    # #     final_time=final_time
    # # )

    # SOLVER CONFIGURATION
    #====================
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    
    solver.set_optim_convergence_tolerance(1e-2)
    solver.set_optim_constraint_tolerance(1e-3)

    # solver.set_num_mesh_intervals(int(np.round((final_time - inital_time) / mesh_interval)))
    # solver.set_multibody_dynamics_mode('implicit')
    # solver.set_minimize_implicit_multibody_accelerations(True)
    # solver.set_implicit_multibody_accelerations_weight(config.acceleration_weight)
   
    solver.set_multibody_dynamics_mode('explicit')
    solver.set_minimize_implicit_auxiliary_derivatives(
        True)
    solver.set_implicit_auxiliary_derivatives_weight(
        weights['aux_deriv_weight'] / 6.0)
    solver.set_implicit_auxiliary_derivative_bounds(
        osim.MocoBounds(-100.0, 100.0))
    solver.set_parallel(28)
    solver.set_parameters_require_initsystem(False)
    solver.set_optim_max_iterations(10000)
    solver.set_scale_variables_using_bounds(True)
    solver.set_optim_finite_difference_scheme('forward')

    
    # # SETTING GUESS
    # # =============
    # # *For now the guess is being provided from a solution from an inverse problem
    # # BUT more work is needed still to determine where exactly guess files are produced
    # guess_file = 'sub03_strifc_cond00000_speed02000_0001_0003_td_marks.sto'   
    # input_guess = osim.MocoTrajectory(guess_file)
    # states_table = input_guess.exportToStatesTable()
    # control_table = input_guess.exportToControlsTable()
    # state_labels = states_table.getColumnLabels()
    # control_labels = control_table.getColumnLabels()

    # for label in state_labels:
    #     if 'reserve' in label:
    #         states_table.removeColumn(label)

    # for label in control_labels:
    #     if 'reserve' in label:
    #         control_table.removeColumn(label)

    # # for label in state_labels:
    # #     if '/forceset' in label and '/activation' not in label:
    # #         states_table.removeColumn(label)
    # #         print(f"     Removed: {label}")

    # # for label in control_labels:
    # #     if '/forceset' in label and '/activation' not in label:
    # #         control_labels.removeColumn(label)
    # #         print(f"     Removed: {label}")

    # sto_adapter = osim.STOFileAdapter()
    # sto_adapter.write(states_table, 'states_table.sto')
    # sto_adapter.write(control_table, 'controls_table.sto')
    
    # guess = solver.createGuess()
    # guess.insertControlsTrajectory(control_table, True)    
    # guess.insertStatesTrajectory(states_table, True)
    # # guess.randomizeAdd()
    # solver.setGuess(guess)

    # solver.resetProblem(problem)

    # Solve and visualize.
    solution_unsealed = study.solve()
    solution = solution_unsealed.unseal()
    solution.write(path_solution)

    # Make sure to close the log file
    if log_file:
        log_file.close()
    # Restore original stdout/stderr
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    # **Be sure to compute grf for each contact sphere
    # solution = osim.MocoTrajectory(path_solution)
    # contact_forces, copTable = create_contact_sphere_force_table(model=model, solution=solution)
    # osim.STOFileAdapter.write(contact_forces, os.path.join(f'{path_solution}', f'_contacts.sto'))

    # external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)

    # # Almost there, just archive the solution for the ground-contact model output
    # # ======================================
    # solution.write(path_solution)
    # osim.STOFileAdapter.write(external_loads, path_solution)
    # study.visualize(solution)

    # =============================================================================
    # ADD AVERAGE SPEED GOAL TO PROBLEM
    # =============================================================================

    # Calculate target walking speed from the input data
    # target_walking_speed = calculate_average_walking_speed(
    #     states_file=temp_states_file,  # Your input states/kinematics file
    #     initial_time=inital_time,        # Start time of analysis
    #     final_time=final_time           # End time of analysis
    # )
    # print("Adding average speed goal to optimization problem...")

    # try:
    #     # Create the average speed goal
    #     speed_goal_weight = 10.0  # Adjust this weight as needed
        
    #     average_speed_goal = osim.MocoAverageSpeedGoal()
    #     average_speed_goal.setName('average_speed_goal')
    #     average_speed_goal.setWeight(speed_goal_weight)
        
    #     # Set the desired average speed (calculated from data)
    #     average_speed_goal.set_desired_average_speed(target_walking_speed)
        
    #     # Add the goal to the problem
    #     problem.addGoal(average_speed_goal)
        
    #     print(f"✓ Added average speed goal:")
    #     print(f"  Target speed: {target_walking_speed:.3f} m/s")
    #     print(f"  Goal weight: {speed_goal_weight}")
        
    # except Exception as e:
    #     print(f"✗ Could not add average speed goal: {e}")

    # print(f"Target walking speed set to: {target_walking_speed:.3f} m/s")


    # Construct a TableProcessor of the coordinate data and pass it to the 
    # tracking tool. TableProcessors can be used in the same way as
    # ModelProcessors by appending TableOperators to modify the base table.
    # A TableProcessor with no operators, as we have here, simply returns the
    # base table.


#   # =============================================================================
#     # OPTIONAL: COORDINATE-SPECIFIC WEIGHTS SETUP
#     # =============================================================================
#     if enable_coordinate_weights:
#         print("\n" + "="*60)
#         print("SETTING UP COORDINATE-SPECIFIC WEIGHTS")
#         print("="*60)
        
#         if use_custom_weights and custom_coordinate_weights:
#             # =============================================================================
#             # OPTION 1: USE CUSTOM COORDINATE WEIGHTS
#             # =============================================================================
#             print("Using custom coordinate weights...")
            
#             try:
#                 # Create a MocoWeightSet with custom weights
#                 stateWeights = osim.MocoWeightSet()
                
#                 weights_applied = 0
#                 for coord_path, weight_value in custom_coordinate_weights.items():
#                     try:
#                         weight = osim.MocoWeight()
#                         weight.setName(coord_path)
#                         weight.setWeight(weight_value)
#                         stateWeights.cloneAndAppend(weight)
#                         weights_applied += 1
#                         print(f"  Custom weight: {coord_path} = {weight_value}")
#                     except Exception as e:
#                         print(f"  Warning: Could not add custom weight for {coord_path}: {e}")
                
#                 # Apply custom weights to tracking tool
#                 track.set_states_weight_set(stateWeights)
#                 print(f"✓ Applied {weights_applied} custom coordinate weights")
                
#             except Exception as e:
#                 print(f"✗ Could not apply custom weights: {e}")
#                 print("  Falling back to global tracking weight")
        
#         else:
#             # =============================================================================
#             # OPTION 2: CALCULATE WEIGHTS BASED ON STANDARD DEVIATION
#             # =============================================================================
#             print("Calculating coordinate weights based on standard deviation...")
            
#             try:
#                 # Load the states/kinematics data
#                 kinematics_table = osim.TimeSeriesTable(temp_states_file)
                
#                 # Load model to get coordinate information
#                 model = osim.Model(path_model)
#                 model.initSystem()
                
#                 # Get column labels from the data
#                 column_labels = kinematics_table.getColumnLabels()
                
#                 # Initialize weight storage
#                 coordinate_weights = {}
#                 coordinate_std_devs = {}
                
#                 print(f"  Analyzing {column_labels.getSize()} columns for weight calculation...")
                
#                 # Calculate standard deviations for each coordinate
#                 for i in range(column_labels.getSize()):
#                     col_name = column_labels.get(i)
                    
#                     # Skip time column
#                     if 'time' in col_name.lower():
#                         continue
                        
#                     # Get coordinate data
#                     try:
#                         column_data = kinematics_table.getDependentColumn(col_name)
#                         data_array = np.array([column_data.get(j) for j in range(column_data.size())])
                        
#                         # Calculate standard deviation
#                         std_dev = np.std(data_array)
#                         coordinate_std_devs[col_name] = std_dev
                        
#                     except Exception as e:
#                         print(f"    Warning: Could not process column {col_name}: {e}")
#                         continue
                
#                 # Define importance weights based on coordinate type and body segment
#                 importance_weights = {
#                     # Pelvis coordinates (very important for overall motion)
#                     'pelvis_tx': {'value': 1.0, 'speed': 0.01},
#                     'pelvis_ty': {'value': 2.0, 'speed': 0.02},
#                     'pelvis_tz': {'value': 1.5, 'speed': 0.015},
#                     'pelvis_tilt': {'value': 1.0, 'speed': 0.01},
#                     'pelvis_list': {'value': 1.0, 'speed': 0.01},
#                     'pelvis_rotation': {'value': 0.5, 'speed': 0.005},
                    
#                     # Hip coordinates (important for gait)
#                     'hip_flexion': {'value': 1.0, 'speed': 0.01},
#                     'hip_adduction': {'value': 1.5, 'speed': 0.015},
#                     'hip_rotation': {'value': 0.8, 'speed': 0.008},
                    
#                     # Knee coordinates (critical for walking)
#                     'knee_angle': {'value': 2.0, 'speed': 0.02},
#                     'knee_adduction': {'value': 1.0, 'speed': 0.01},
#                     'knee_rotation': {'value': 0.5, 'speed': 0.005},
                    
#                     # Ankle coordinates (important for foot clearance and push-off)
#                     'ankle_angle': {'value': 2.0, 'speed': 0.02},
#                     'subtalar_angle': {'value': 1.0, 'speed': 0.01},
#                     'mtp_angle': {'value': 0.8, 'speed': 0.008},
                    
#                     # Lumbar coordinates (often problematic, lower weights)
#                     'lumbar_extension': {'value': 0.1, 'speed': 0.001},
#                     'lumbar_bending': {'value': 0.1, 'speed': 0.001},
#                     'lumbar_rotation': {'value': 0.1, 'speed': 0.001},
                    
#                     # Arm coordinates (if present, usually less important for walking)
#                     'arm_flex': {'value': 0.3, 'speed': 0.003},
#                     'arm_add': {'value': 0.2, 'speed': 0.002},
#                     'arm_rot': {'value': 0.1, 'speed': 0.001},
#                     'elbow_flex': {'value': 0.2, 'speed': 0.002},
#                     'pro_sup': {'value': 0.1, 'speed': 0.001},
#                     'wrist_flex': {'value': 0.1, 'speed': 0.001},
#                     'wrist_dev': {'value': 0.1, 'speed': 0.001},
#                 }
                
#                 # Calculate final weights using both importance and variability
#                 for col_name, std_dev in coordinate_std_devs.items():
                    
#                     # Extract coordinate base name for lookup
#                     coord_base_name = None
#                     weight_type = None  # 'value' or 'speed'
                    
#                     # Determine if this is a position or speed column
#                     if '/value' in col_name or col_name.endswith('_value'):
#                         weight_type = 'value'
#                     elif '/speed' in col_name or col_name.endswith('_speed'):
#                         weight_type = 'speed'
#                     else:
#                         # For raw coordinate names, assume position
#                         weight_type = 'value'
                    
#                     # Extract base coordinate name
#                     for importance_key in importance_weights.keys():
#                         if importance_key in col_name:
#                             coord_base_name = importance_key
#                             break
                    
#                     # Calculate weight based on importance and variability
#                     if coord_base_name and coord_base_name in importance_weights:
#                         importance_weight = importance_weights[coord_base_name][weight_type]
                        
#                         # Weight calculation: higher std_dev = lower weight
#                         if std_dev > 0:
#                             variability_factor = 1.0 / (1.0 + std_dev)
#                             final_weight = importance_weight * variability_factor
#                         else:
#                             final_weight = importance_weight
                        
#                         # Apply minimum weight threshold
#                         final_weight = max(final_weight, 0.001)
                        
#                     else:
#                         # Default weight for unrecognized coordinates
#                         final_weight = 0.1 if weight_type == 'value' else 0.001
                    
#                     coordinate_weights[col_name] = final_weight
#                     print(f"    {col_name}: std={std_dev:.4f}, weight={final_weight:.4f}")
                
#                 # Create and apply MocoWeightSet
#                 stateWeights = osim.MocoWeightSet()
#                 weights_applied = 0
                
#                 for coord_path, weight_value in coordinate_weights.items():
#                     try:
#                         weight = osim.MocoWeight()
#                         weight.setName(coord_path)
#                         weight.setWeight(weight_value)
#                         stateWeights.cloneAndAppend(weight)
#                         weights_applied += 1
#                     except Exception as e:
#                         print(f"    Warning: Could not add weight for {coord_path}: {e}")
                
#                 # Apply calculated weights to tracking tool
#                 track.set_states_weight_set(stateWeights)
#                 print(f"✓ Applied {weights_applied} calculated coordinate weights")
                
#             except Exception as e:
#                 print(f"✗ Could not calculate coordinate weights: {e}")
#                 print("  Falling back to global tracking weight")
    
#     else:
#         print("Coordinate-specific weights disabled - using global tracking weight")



    # =============================================================================
    # ADVANCED WEIGHTING: Body-segment based adjustments
    # =============================================================================

    # print("Applying body-segment based weight adjustments...")

    # try:
    #     # Get the problem to apply additional weight adjustments
    #     # Note: This happens after track.initialize() has been called
        
    #     # Define segment-based weight multipliers (similar to tracking_problem.py)
    #     segment_multipliers = {
    #         'pelvis': 1.5,     # Critical for overall motion
    #         'hip': 1.2,        # Important for gait
    #         'knee': 1.3,       # Critical for walking
    #         'ankle': 1.3,      # Important for foot clearance
    #         'subtalar': 1.0,   # Normal importance
    #         'mtp': 0.8,        # Less critical
    #         'lumbar': 0.1,     # Often problematic
    #         'arm': 0.3,        # Less important for walking
    #         'elbow': 0.2,      # Less important
    #         'wrist': 0.1,      # Least important
    #     }
        
    #     # Apply segment-based multipliers to existing weights
    #     updated_weights = {}
        
    #     for coord_path, base_weight in coordinate_weights.items():
    #         # Find which segment this coordinate belongs to
    #         segment_multiplier = 1.0  # Default
            
    #         for segment, multiplier in segment_multipliers.items():
    #             if segment in coord_path.lower():
    #                 segment_multiplier = multiplier
    #                 break
            
    #         # Apply segment multiplier
    #         updated_weight = base_weight * segment_multiplier
    #         updated_weights[coord_path] = updated_weight
        
    #     print(f"  Applied segment-based multipliers to {len(updated_weights)} coordinates")
        
    #     # Update the weight set with new values
    #     stateWeights.clearAndDestroy()
        
    #     for coord_path, weight_value in updated_weights.items():
    #         try:
    #             weight = osim.MocoWeight()
    #             weight.setName(coord_path)
    #             weight.setWeight(weight_value)
    #             stateWeights.cloneAndAppend(weight)
    #         except:
    #             continue
        
    #     # Re-apply the updated weights
    #     track.set_states_weight_set(stateWeights)
        
    #     print("✓ Body-segment based weight adjustments applied")
        
    # except Exception as e:
    #     print(f"✗ Could not apply segment-based adjustments: {e}")
        

    ## =========================
    # Manual Weight override specific weights

    # print("Applying manual weight overrides...")

    # Define manual overrides for problematic or critical coordinates
    # manual_overrides = {
        # Increase weight for critical coordinates
        # '/jointset/ground_pelvis/pelvis_tx/value': 3.0
        # '/jointset/walker_knee_l/knee_angle_l/value': 3.0,
        # '/jointset/walker_knee_r/knee_angle_r/value': 3.0,
        # '/jointset/ankle_l/ankle_angle_l/value': 2.5,
        # '/jointset/ankle_r/ankle_angle_r/value': 2.5,
        
        # Decrease weight for problematic coordinates
        # '/jointset/back/lumbar_extension/value': 0.01,
        # '/jointset/back/lumbar_bending/value': 0.01,
        # '/jointset/back/lumbar_rotation/value': 0.01,
    # }

    # Apply manual overrides
    # for coord_path, override_weight in manual_overrides.items():
    #     try:
    #         # Find and update the weight in the weight set
    #         for i in range(stateWeights.getSize()):
    #             weight_obj = stateWeights.get(i)
    #             if coord_path in weight_obj.getName():
    #                 weight_obj.setWeight(override_weight)
    #                 print(f"  Override: {coord_path} = {override_weight}")
    #                 break
    #     except Exception as e:
    #         print(f"  Could not apply override for {coord_path}: {e}")

    # print("Manual weight overrides applied!")


    # =============================================================================
    # ADD AVERAGE WALKING SPEED GOAL
    # =============================================================================
    
    # if enable_speed_goal and target_walking_speed is not None:
    #     print(f"\n" + "="*50)
    #     print("ADDING AVERAGE WALKING SPEED GOAL")
    #     print("="*50)
        
    #     try:
    #         # Create average speed goal
    #         averageSpeedGoal = osim.MocoAverageSpeedGoal('average_speed_goal', speed_goal_weight)
    #         averageSpeedGoal.set_desired_average_speed(target_walking_speed)
            
    #         # Add the goal to the problem
    #         problem.addGoal(averageSpeedGoal)
            
    #         print(f"  ✓ Added average speed goal:")
    #         print(f"    Target speed: {target_walking_speed:.3f} m/s")
    #         print(f"    Goal weight: {speed_goal_weight}")
            
    #     except Exception as e:
    #         print(f"  ✗ Could not add average speed goal: {e}")
    # else:
    #     print("Average walking speed goal disabled")



    # =============================================================================
    # Continue with original problem setup
    # =============================================================================

    # effort = osim.MocoControlGoal.safeDownCast(problem.updGoal("control_effort"))
    # effort.setDivideByDisplacement(True)

    # =============================================================================
    # ADD WALKING SPEED GOAL (OPTIONAL)
    # =============================================================================

    # average_speed = osim.MocoAverageSpeedGoal()
    # average_speed.setName('average_speed')
    # average_speed.set_desired_average_speed(target_walking_speed)
    # problem.addGoal(average_speed)

    # Put a large weight on the pelvis CoordinateActuators, which act as the
    # residual, or 'hand-of-god', forces which we would like to keep as small
    # as possible.
    # model = modelProcessor.process()
    # model.initSystem()
    # forceSet = model.getForceSet()
    # for i in range(forceSet.getSize()):
    #     forcePath = forceSet.get(i).getAbsolutePathString()
    #     if 'pelvis' in str(forcePath):
    #         effort.setWeightForControl(forcePath, 10)

    

if __name__ == "__main__":
    # Solve the muscle-driven state tracking problem.
    # This problem could take an hour or more to solve, depending on the number of
    # processor cores available for parallelization. With 12 cores, it takes around
    # 25 minutes.
    muscleDrivenStateTracking(enable_contact_tracking=True, contact_weight=100, 
                              enable_advanced_solver=True, mesh_interval=0.01, 
                              enable_state_tracking=True,
                              target_walking_speed=2.0,
                              enable_speed_goal=False,
                              speed_goal_weight=10, 
                              enable_coordinate_weights=True, 
                              use_custom_weights=False,
                              custom_coordinate_weights=None,
                              inital_time=0.001,
                              final_time=0.746)

