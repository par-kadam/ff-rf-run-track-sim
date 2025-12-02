# -------------------------------------------------------------------------- #
# OpenSim Moco: exampleMocoInverse.py                                        #
# -------------------------------------------------------------------------- #
# Copyright (c) 2020 Stanford University and the Authors                     #
#                                                                            #
# Author(s): Christopher Dembia                                              #
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

# This example shows how to use the MocoInverse tool to exactly prescribe a
# motion and estimate muscle behavior for walking. The first example does not
# rely on electromyography data, while the second example penalizes deviation
# from electromyography data for a subset of muscles.
#
# See the README.txt next to this file for more information.

import argparse
import os
import sys
import numpy as np

# from build.opensim.core.Bindings.Python.simulation import TabOpUseAbsoluteStateNames
# from install.opensim.core.sdk.Python.opensim.simulation import TabOpLowPassFilter
import opensim as osim

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Geometry")

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/src/opensim-gui/opensim-models/Models/RajagopalModel/Geometry")

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

def solveMocoInverseKinematics(
        path_model          = None,
        path_musclepaths    = None,
        path_grf            = None,
        path_kinematics     = None,
        time_initial        = None,
        time_final          = None,
        trial_speed         = None,
        reserves            = 1.0,
        mesh_interval       = 0.01,
        path_solution_guess = None,
        path_solution       = None):


    # Construct the MocoInverse tool.
    inverse = osim.MocoInverse()

    # Construct a ModelProcessor and set it on the tool. The default
    # muscles in the model are replaced with optimization-friendly
    # DeGrooteFregly2016Muscles, and adjustments are made to the default muscle
    # parameters.
    if path_model is not None:
        model = osim.Model(path_model)
        state = model.initSystem()
        mass = model.getTotalMass(state)

        # # =============================================================================
        # # TORQUE ACTUATORS
        # # =============================================================================
        

        # coordNames = ['arm_flex', 'arm_add', 'arm_rot', 'elbow_flex', 'pro_sup']
        # strengths = [0.5, 0.5, 0.5, 0.2, 0.2]
        # for coordName, strength in zip(coordNames, strengths):
        #     for side in ['_l', '_r']:
        #         actu = osim.ActivationCoordinateActuator()
        #         actu.set_coordinate(f'{coordName}{side}')
        #         actu.setName(f'torque_{coordName}{side}')
        #         actu.setOptimalForce(strength*mass)
        #         actu.setMinControl(-1.0)
        #         actu.setMaxControl(1.0)
        #         model.addForce(actu)

        # # Then add in lumbar actuators
        # coordNames = ['lumbar_extension', 'lumbar_bending', 'lumbar_rotation']
        # for coordName in coordNames:
        #     actu = osim.ActivationCoordinateActuator()
        #     actu.set_coordinate(f'{coordName}')
        #     actu.setName(f'torque_{coordName}')
        #     actu.setOptimalForce(mass)
        #     actu.setMinControl(-1.0)
        #     actu.setMaxControl(1.0)
        #     model.addForce(actu)

        # =============================================================================
        # JOINT STIFFNESS MODIFICATIONS
        # =============================================================================

        # lumbar_stiffness_scale = 1.0
            
        # # Stiffness values for lumbar joint
        # stiffnesses = [1.0, 1.5, 0.5] # in N-m/rad*kg
        # for coordName, stiffness in zip(coordNames, stiffnesses):
        #     sgf = osim.SpringGeneralizedForce(coordName)
        #     sgf.setName(f'passive_stiffness_{coordName}')
        #     sgf.setStiffness(lumbar_stiffness_scale * stiffness * mass)
        #     sgf.setViscosity(2.0)
        #     model.addForce(sgf)

        
        # Custom ModOp to scale lumbar stiffness would go here
        # This is a placeholder - actual implementation would require custom ModOp
        # modelProcessor.append(osim.ModOpScaleLumbarStiffness(lumbar_stiffness_scale))

        # =============================================================================
        # RESERVE ACTUATORS
        # =============================================================================
        
        # reserve_strength = 50
            
        # coordNames = ['pelvis_tx', 'pelvis_ty', 'pelvis_tz', 
        #             'pelvis_list', 'pelvis_tilt', 'pelvis_rotation',
        #             'hip_adduction_r', 'hip_rotation_r', 'hip_flexion_r',
        #             'hip_adduction_l', 'hip_rotation_l', 'hip_flexion_l',
        #             'knee_angle_r', 'ankle_angle_r', 
        #             'knee_angle_l', 'ankle_angle_l', 
        #             'mtp_angle_r', 'mtp_angle_l']
        #             # 'subtalar_angle_r', 'subtalar_angle_l'
        # for coordName in coordNames:
        #     actu = osim.CoordinateActuator()
        #     actu.set_coordinate(coordName)
        #     actu.setName(f'reserve_{coordName}')
        #     scale = 1.0
        #     if ((coordName == 'pelvis_tx') or
        #             (coordName == 'pelvis_ty') or
        #             (coordName == 'pelvis_tz')):
        #         scale = 10.0 
        #     actu.setOptimalForce(scale * reserve_strength)
        #     actu.setMinControl(-1.0)
        #     actu.setMaxControl(1.0)
        #     model.addForce(actu)


        coord_set = model.updCoordinateSet()
        for coord_name in ['subtalar_angle']:
            for side in ['_l', '_r']:
                coord = coord_set.get(f'{coord_name}{side}')
                coord.set_locked(False)

        model.finalizeConnections()
        modelProcessor = osim.ModelProcessor(model)
        # weld_joints = list()
        # for side in ['_l', '_r']:
        #     weld_joints.append(f"subtalar{side}")
        # modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weld_joints))
        
        modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))
        modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
        modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
        # # Only valid for DeGrooteFregly2016Muscles.
        modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
        # # Only valid for DeGrooteFregly2016Muscles.
        modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
        # modelProcessor.append(osim.ModOpFiberDampingDGF(0.01))
        # # Use a function-based representation for the muscle paths. This is
        # # recommended to speed up convergence, but if you would like to use
        # # the original GeometryPath muscle wrapping instead, simply comment out
        # # this line. To learn how to create a set of function-based paths for
        # # your model, see the example 'examplePolynomialPathFitter.py'.
        modelProcessor.append(osim.ModOpReplacePathsWithFunctionBasedPaths(
                path_musclepaths))
        
        # modelProcessor.append(osim.ModOpAddReserves(reserves))
        

        
        # =============================================================================
        # MUSCLE MODELING MODIFICATIONS
        # =============================================================================

        #Enable tendon compliance for ankle plantarflexors
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
                    
                    muscle.set_active_force_width_scale(1.5)
                    muscle.set_max_contraction_velocity(10.0)

        model.finalizeConnections()
        modelProcessor = osim.ModelProcessor(model)
        
        modelProcessor.append(
            osim.ModOpUseImplicitTendonComplianceDynamicsDGF())

    # if path_grf is not None:
    #     modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))
        # modelProcessor.append(osim.ModOpDisableContactForces())

    # modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    # modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())

    # Only valid for DeGrooteFregly2016Muscles.
    # modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())

    # Only valid for DeGrooteFregly2016Muscles.
    # modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))

    # Use a function-based representation for the muscle paths. This is
    # recommended to speed up convergence, but if you would like to use
    # the original GeometryPath muscle wrapping instead, simply comment out
    # this line. To learn how to create a set of function-based paths for
    # your model, see the example 'examplePolynomialPathFitter.py'.
    # if path_musclepaths is not None:
    #     modelProcessor.append(
    #         osim.ModOpReplacePathsWithFunctionBasedPaths(path_musclepaths))

    if reserves is not None:
        modelProcessor.append(osim.ModOpAddReserves(reserves))
    


    inverse.setModel(modelProcessor)

    # Construct a TableProcessor of the coordinate data and pass it to the
    # inverse tool. TableProcessors can be used in the same way as
    # ModelProcessors by appending TableOperators to modify the base table.
    # A TableProcessor with no operators, as we have here, simply returns the
    # base table.
    if path_kinematics is not None:
        
        coords_table = osim.TimeSeriesTable(path_kinematics)
        time = coords_table.getIndependentColumn()
        nrows = coords_table.getNumRows()
        # pelvis_tx = coords_table.getDependentColumn('pelvis_tx')
        # pelvis_tx_new = osim.Vector(nrows, 0.0)
        # for i in np.arange(nrows):
        #     dx = (time[int(i)] - time[0]) * trial_speed
        #     pelvis_tx_new[int(i)] = pelvis_tx[int(i)] + dx

        # pelvis_tx_upd = coords_table.updDependentColumn('pelvis_tx')
        # for i in np.arange(nrows):
        #     pelvis_tx_upd[int(i)] = pelvis_tx_new[int(i)]
        
        TableProcessor = osim.TableProcessor(coords_table)
        TableProcessor.append(osim.TabOpLowPassFilter(6))
        TableProcessor.append(osim.TabOpUseAbsoluteStateNames())
        inverse.setKinematics(TableProcessor)

    if time_initial is not None:
        inverse.set_initial_time(time_initial)
    else :
        inverse.set_initial_time(time[0])

    if time_final is not None:
        inverse.set_final_time(time_final)
    else:
        inverse.set_final_time(time[-1])

    if mesh_interval is not None:
        inverse.set_mesh_interval(mesh_interval)

    # By default, Moco gives an error if the kinematics contains extra columns.
    # Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(True)

    # If a guess solution is available then initialize the solver with it.
    if path_solution_guess is not None:
        study = inverse.initialize()
        solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
        solver.setGuessFile(path_solution_guess)

        # Update the solver tolerances.
        # set_optim_convergence_tolerance(1e-6)
    # else:
    #     study = inverse.initialize()
    #     solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    #     solver.setGuessFile('example_muscle_drive_27.sto')
    #     solver.set_optim_convergence_tolerance(1e-2)
    #     solver.set_optim_constraint_tolerance(1e-2)
    #     solver.set_optim_max_iterations(10000)
    #     solver.set_parallel(28)
        # solver.set_multibody_dynamics_mode('implicit')
        # solver.set_optim_finite_difference_scheme('forward')

    # Solve the problem and write the solution to a Storage file.
    solution = inverse.solve()

    solution.getMocoSolution().write(path_solution)
    # If a solution is found, then Write it to a Storage file.
    # if not solution.isSealed():
    #     solution.write(path_solution)

#    # Generate a PDF with plots for the solution trajectory.
#    model = modelProcessor.process()
#    report = osim.report.Report(model, path_solution, bilateral=True)
#    # The PDF is saved to the working directory.
#    report.generate()

def solveMocoInverseKinematics_v44(
        path_model          = None,
        path_grf            = None,
        path_kinematics     = None,
        time_initial        = None,
        time_final          = None,
        force_reserve       = 5.0,
        mesh_interval       = 0.010,
        path_solution_guess = None,
        path_solution       = None):

    # Construct the MocoInverse tool.
    inverse = osim.MocoInverse()

    # Construct a ModelProcessor and set it on the tool. The default
    # muscles in the model are replaced with optimization-friendly
    # DeGrooteFregly2016Muscles, and adjustments are made to the default muscle
    # parameters.
    modelProcessor = osim.ModelProcessor(path_model)

    modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))


    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())

    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())

    # Only valid for DeGrooteFregly2016Muscles.
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    modelProcessor.append(osim.ModOpAddReserves(force_reserve))

    inverse.setModel(modelProcessor)

    # Construct a TableProcessor of the coordinate data and pass it to the
    # inverse tool. TableProcessors can be used in the same way as
    # ModelProcessors by appending TableOperators to modify the base table.
    # A TableProcessor with no operators, as we have here, simply returns the
    # base table.
    coordinates = osim.TableProcessor(path_kinematics)
    coordinates.append(TabOpLowPassFilter(6))
    coordinates.append(TabOpUseAbsoluteStateNames())
    inverse.setKinematics(
        coordinates)

#    # Initial time, final time, and mesh interval.
#    if time_initial is not None:
#        inverse.set_initial_time(time_initial)
#
#    if time_final is not None:
#        inverse.set_final_time(time_final)

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

    inverse.set_mesh_interval(mesh_interval)

    # By default, Moco gives an error if the kinematics contains extra columns.
    # Here, we tell Moco to allow (and ignore) those extra columns.
    inverse.set_kinematics_allow_extra_columns(True)

    study = inverse.initialize()

    # Update the solver tolerances.
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    if path_solution_guess is not None:
        solver.setGuessFile(path_solution_guess)
        inverse.set_optim_convergence_tolerance(1e-2)
        inverse.set_optim_constraint_tolerance(1e-2)

    # Solve the problem.
    solution = study.solve()

    model = modelProcessor.process()
    external_loads = osim.createExternalLoadsTableForGait(model, solution, forceNamesRightFoot, forceNamesLeftFoot)
    external_loads_file_path = path_solution.replace('.sto', '_external_loads.sto')
    osim.STOFileAdapter.write(external_loads, external_loads_file_path)

    # Write the solution to a Storage file.
    if not solution.getMocoSolution().isSealed():
        solution.getMocoSolution().write(path_solution)
        # Generate a PDF with plots for the solution trajectory.
        model = modelProcessor.process()
        report = osim.report.Report(model,
                                    path_solution,
                                    'bilateral', True)
        # # The PDF is saved to the working directory.
        report.generate()

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
        if 'sub' in subject and os.path.isdir(os.path.join(path_data, subject))]
    subjects.sort()

    # Initialize the MocoInverse parameters.
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

                if f.endswith('_SCALED_addbio_free_subt_mtp.osim'):
                    path_model = os.path.join(path_data, sub, ses, f)

                if f.endswith('_FunctionBasedPathSet.xml'):
                    path_musclepaths = os.path.join(path_data, sub, ses, f)

            file_cycles_kinematics = os.listdir(os.path.join(path_data, sub, ses))
            # Keep only the files having the sting '_ik.sto' in their names:
            # - filter out everything else.
            file_cycles_kinematics = [
                file_cycle_kinematics for file_cycle_kinematics in file_cycles_kinematics \
                if '_ik_move.sto' in file_cycle_kinematics and os.path.isfile(os.path.join(path_data, sub, ses, file_cycle_kinematics))]
            file_cycles_kinematics.sort()

            for file_cycle_kinematics in file_cycles_kinematics:
                file_cycle_stem = file_cycle_kinematics[:-12]
                trial_speed     = int(file_cycle_kinematics[-21])

                path_grf              = os.path.join(path_data, sub, ses, file_cycle_stem + '_grf_move.xml')
                path_cycle_kinematics = os.path.join(path_data, sub, ses, file_cycle_kinematics)
                path_solution         = os.path.join(path_data, sub, ses, file_cycle_stem + '_mi.sto')


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

                    solveMocoInverseKinematics(
                        path_model          = path_model,
                        path_musclepaths    = path_musclepaths,
                        path_grf            = path_grf,
                        path_kinematics     = path_cycle_kinematics,
                        time_initial        = None,
                        time_final          = None,
                        trial_speed         = 2,
                        reserves            = 20.0,
                        mesh_interval       = 0.05,
                        path_solution_guess = path_solution_guess,
                        path_solution       = path_solution)
