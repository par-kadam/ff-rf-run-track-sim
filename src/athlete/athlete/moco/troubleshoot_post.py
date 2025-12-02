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

path_opensim_core = "_____"

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
path_solution       = os.path.join('./sub03_testing/', 'muscle_driven_state_tracking_sol_8.sto')


# CONSIDER PUTTING THE FOLLOWING FUNCTIONS IN ANOTHER SCRIPT AND IMPORTING FROM THERE FOR ORGANIZATION:
# ====================================================================================================


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
                               implicit_tendon_dynamics=False,
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
                        'knee_angle_r', 'ankle_angle_r', 'subtalar_angle_r',
                        'knee_angle_l', 'ankle_angle_l', 'subtalar_angle_l',
                        'mtp_angle_r', 'mtp_angle_l']
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

    # =============================================================================
    # FINALIZATION
    # =============================================================================

    model.finalizeConnections()

    modelProcessor = osim.ModelProcessor(model)

    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))  

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
                              inital_time=None,
                              final_time=None):


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
                               implicit_tendon_dynamics=False,
                               lumbar_stiffness_scale=1.0,
                               
                               # Logging and debugging
                               verbose=True)


    model = modelProcessor.process()
    model.initSystem()
    
    # modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))

    
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

    speed_bounds = [-20.0, 20.0]  # rad/s or m/s
    muscle_bounds = [0.01, 1.0]  # Muscle activations between 1% and 100%
    torque_bounds = [-1.0, 1.0]  # Normalized torque bounds
    tendon_bounds = [0.01, 1.8]


    expiremental_states = osim.TableProcessor('sub03_strifc_cond00000_speed02000_0001_0003_td_marks.sto')

    # Now create coordinates for tracking reference
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    exp_state_filter = expiremental_states.process()
    coordinatesTable = osim.TimeSeriesTable(exp_state_filter)

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
                        value_weight = 1.0/denom
                        speed_weight = 0.01/denom
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

    
    # Create and name the basic instance of the MocoTrack tool.
    track = osim.MocoTrack()
    track.setName('single_run_track')
    track.setModel(modelProcessor)
    track.setStatesReference(expiremental_states)
    track.set_states_global_tracking_weight(weights['state_tracking_weight'] / model.getNumCoordinates())
    track.set_apply_tracked_states_to_guess(True)
    track.set_states_weight_set(state_weights)

    debug_model_processor(modelProcessor, path_model)

    # This setting allows extra data columns contained in the states
    # reference that don't correspond to model coordinates.
    track.set_allow_unused_references(True)

    # Since there is only coordinate position data in the states references,
    # this setting is enabled to fill in the missing coordinate speed data using
    # the derivative of splined position data.
    track.set_track_reference_position_derivatives(True)
    track.set_control_effort_weight(weights['control_weight'] / num_forces)

    # Initial time, final time, and mesh interval.
    track.set_initial_time(inital_time)
    track.set_final_time(final_time)
    track.set_mesh_interval(0.01)

    # Instead of calling solve(), call initialize() to receive a pre-configured
    # MocoStudy object based on the template settings. Use this to customize the
    # problem beyond the MocoTrack interface.
    study = track.initialize()
    # Get a reference to the MocoControlCost that is added to every MocoTrack
    # problem by default.
    problem = study.updProblem()


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
    # control_effort = osim.MocoControlGoal().safeDownCast(
    #     problem.updGoal('control_effort'))
    # control_effort.setDivideByDisplacement(True)


    # SOLVER CONFIGURATION
    #====================
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    solver.resetProblem(problem)
    
    solver.set_optim_convergence_tolerance(1e-2)
    solver.set_optim_constraint_tolerance(1e-3)

    solver.set_num_mesh_intervals(int(np.round((final_time - inital_time) / mesh_interval)))

   
    solver.set_multibody_dynamics_mode('implicit')
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

    
    # SETTING GUESS
    # =============
    # *For now the guess is being provided from a solution from an inverse problem
    # BUT more work is needed still to determine where exactly guess files are produced
    guess_file = 'sub03_strifc_cond00000_speed02000_0001_0003_td_marks.sto'   
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

    # for label in state_labels:
    #     if '/forceset' in label and '/activation' not in label:
    #         states_table.removeColumn(label)
    #         print(f"     Removed: {label}")

    # for label in control_labels:
    #     if '/forceset' in label and '/activation' not in label:
    #         control_labels.removeColumn(label)
    #         print(f"     Removed: {label}")

    sto_adapter = osim.STOFileAdapter()
    sto_adapter.write(states_table, 'states_table.sto')
    sto_adapter.write(control_table, 'controls_table.sto')
    
    guess = solver.createGuess()
    guess.insertControlsTrajectory(control_table, True)    
    guess.insertStatesTrajectory(states_table, True)
    # guess.randomizeAdd()
    solver.setGuess(guess)


    # Solve and visualize.
    solution_unsealed = study.solve()
    solution = solution_unsealed.unseal()
    solution.write(path_solution)



    

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
                              inital_time=None,
                              final_time=None)

