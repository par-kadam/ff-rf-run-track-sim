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
from osimFunctions import kinematicsToStates
import opensim as osim

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/sim/src/opensim-models/Geometry")

osim.ModelVisualizer.addDirToGeometrySearchPaths(
    "/home/lisca/biomec/sim/src/opensim-models/Models/RajagopalModel/Geometry")

def solveMocoInverseKinematicsEMG(
        path_model          = None,
        path_grf            = None,
        path_kinematics     = None,
        path_states         = None,
        path_musclepaths    = None,
        time_initial        = None,
        time_final          = None,
        force_reserve       = 1.0,
        mesh_interval       = 0.010,
        path_solution_guess = None,
        path_solution       = None,
        path_emg            = None):
    
    
    inverse = osim.MocoInverse()
    model = osim.Model(path_model)
    state = model.initSystem()
    mass = model.getTotalMass(state)
    initial_forceset = model.getForceSet()
    remove_actuators = list()

    coord_set = model.updCoordinateSet()
    for coord_name in ['subtalar_angle']:
        for side in ['_l', '_r']:
            coord = coord_set.get(f'{coord_name}{side}')
            coord.set_locked(False)

        model.finalizeConnections()

    weld_joints = list()
    for side in ['_l', '_r']:
        # weld_joints.append(f"subtalar{side}")
        weld_joints.append(f"radius_hand{side}")
        remove_actuators.append(f"wrist_flex{side}")
        remove_actuators.append(f"wrist_dev{side}")

    # Get & iterate through indices in reverse order to remove
    for i in range(initial_forceset.getSize() - 1, -1, -1):
        force = initial_forceset.get(i).getName()
        if force in remove_actuators:
            initial_forceset.remove(i)
            print(f"Removed: {force}")

    modelProcessor = osim.ModelProcessor(model)
    modelProcessor.append(osim.ModOpReplaceJointsWithWelds(weld_joints))

    modelProcessor.append(osim.ModOpAddExternalLoads(path_grf))
    modelProcessor.append(osim.ModOpIgnoreTendonCompliance())
    modelProcessor.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    modelProcessor.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    modelProcessor.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))

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


    modelProcessor.append(osim.ModOpAddReserves(force_reserve))
    inverse.setModel(modelProcessor)
    # inverse.setKinematics(osim.TableProcessor(path_kinematics))
    # inverse.set_initial_time(0.81)
    # inverse.set_final_time(1.79)
    inverse.set_mesh_interval(mesh_interval)
    
    path_states = path_states

    kinematicsToStates(kinematicsFileName=path_kinematics,
                    osimModelFileName=path_model,
                    outputFileName=path_states,
                    inDegrees=True,
                    outDegrees=False,
                    excludeColumns=['time', 'pelvis_tx', 'pelvis_ty', 'pelvis_tz'], 
                    include_accelerations=False,
                    include_speeds=False)

    coords_table = osim.TimeSeriesTable(path_kinematics)
    # for side in ['_l', '_r']:
    #     # coords_table.removeColumn(f'/jointset/radius_hand{side}/wrist_flex{side}/value')
    #     # coords_table.removeColumn(f'/jointset/radius_hand{side}/wrist_dev{side}/value')
    #     coords_table.removeColumn(f'wrist_flex{side}')
    #     coords_table.removeColumn(f'wrist_dev{side}')

    time = coords_table.getIndependentColumn()
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

    inverse.set_kinematics_allow_extra_columns(True)

    table_processor = osim.TableProcessor(coords_table)
    table_processor.append(osim.TabOpLowPassFilter(6))
    table_processor.append(osim.TabOpUseAbsoluteStateNames())
    inverse.setKinematics(table_processor)

    study = inverse.initialize()
    problem = study.updProblem()

    torque_bounds = [-1.0, 1.0]
    problem.setStateInfoPattern('.*reserve.*', torque_bounds, [], [])
    
    # Add electromyography tracking.
    emgTracking = osim.MocoControlTrackingGoal('emg_tracking')
    emgTracking.setWeight(1000.0)
    # Each column in electromyography.sto is normalized so the maximum value in
    # each column is 1.0.
    controlsRef = osim.TimeSeriesTable(path_emg)
    # Scale the tracked muscle activity based on peak levels from
    # "Gait Analysis: Normal and Pathological Function" by
    # Perry and Burnfield, 2010 (digitized by Carmichael Ong).
    soleus = controlsRef.updDependentColumn('right_soleus')
    gasmed = controlsRef.updDependentColumn('right_gastrocnemius_medialis')
    tibant = controlsRef.updDependentColumn('right_tibialis_anterior')
    peroneus = controlsRef.updDependentColumn('right_peroneus_longus')
    gaslat = controlsRef.updDependentColumn('right_gastrocnemius_lateralis')
    rectfem = controlsRef.updDependentColumn('right_rectus_femoris')
    bicfem = controlsRef.updDependentColumn('right_biceps_femoris')
    tfl = controlsRef.updDependentColumn('right_tensor_fascia_latae')
    for t in range(0, controlsRef.getNumRows()):
        soleus[t] = 1 * soleus[t]
        gasmed[t] = 1 * gasmed[t]
        tibant[t] = 1 * tibant[t]
        peroneus[t] = 1 * peroneus[t]
        gaslat[t] = 1 * gaslat[t]
        rectfem[t] = 1 * rectfem[t]
        bicfem[t] = 1 * bicfem[t]
        tfl[t]    = 1 * tfl[t]
    emgTracking.setReference(osim.TableProcessor(controlsRef))
    # Associate actuators in the model with columns in electromyography.sto.
    emgTracking.setReferenceLabel('/forceset/soleus_r', 'right_soleus')
    emgTracking.setReferenceLabel('/forceset/gasmed_r', 'right_gastrocnemius_medialis')
    emgTracking.setReferenceLabel('/forceset/gaslat_r', 'right_gastrocnemius_lateralis')
    emgTracking.setReferenceLabel('/forceset/tibant_r', 'right_tibialis_anterior')
    emgTracking.setReferenceLabel('/forceset/perlong_r', 'right_peroneus_longus')
    emgTracking.setReferenceLabel('/forceset/recfem_r', 'right_rectus_femoris')
    emgTracking.setReferenceLabel('/forceset/bflh_r', 'right_biceps_femoris')
    emgTracking.setReferenceLabel('/forceset/tfl_r', 'right_tensor_fascia_latae')
    problem.addGoal(emgTracking)

    # Update the solver tolerances.
    solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
    if path_solution_guess is not None:
        solver.setGuessFile(path_solution_guess)
        solver.set_optim_convergence_tolerance(1e-2)
        solver.set_optim_constraint_tolerance(1e-2)

    # Solve the problem and write the solution to a Storage file.
    solution = inverse.solve()
    
    # Write the solution to a Storage file.
    if not solution.getMocoSolution().isSealed():
        solution.getMocoSolution().write(path_solution)
        # Generate a PDF with plots for the solution trajectory.
        model = modelProcessor.process()
        # report = osim.report.Report(model,
        #                             path_solution,
        #                             'bilateral', True)
        # # # The PDF is saved to the working directory.
        # report.generate()

    # # Write the reference data in a way that's easy to compare to the solution.
    # controlsRef.removeColumn('medial_hamstrings')
    # controlsRef.removeColumn('biceps_femoris')
    # controlsRef.removeColumn('vastus_lateralis')
    # controlsRef.removeColumn('vastus_medius')
    # controlsRef.removeColumn('rectus_femoris')
    # controlsRef.removeColumn('gluteus_maximus')
    # controlsRef.removeColumn('gluteus_medius')
    # controlsRef.setColumnLabels(['/forceset/soleus_r', '/forceset/gasmed_r',
    #                              '/forceset/tibant_r'])
    # controlsRef.appendColumn('/forceset/gaslat_r', gasmed)
    # osim.STOFileAdapter.write(controlsRef, 'controls_reference.sto')

    # # Generate a report comparing MocoInverse solutions without and with EMG
    # # tracking.
    # model = modelProcessor.process()
    # output = 'example3DWalking_MocoInverseWithEMG_report.pdf'
    # ref_files = [
    #     'example3DWalking_MocoInverseWithEMG_solution.sto',
    #     'controls_reference.sto']
    # report = osim.report.Report(model,
    #                             'example3DWalking_MocoInverse_solution.sto',
    #                             output=output, bilateral=True,
    #                             ref_files=ref_files)
    # # The PDF is saved to the working directory.
    # report.generate()

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

    # Initialize the MocoInverse parameters.
    path_emg              = None
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

                # if f.endswith('grf_move.xml'):
                #     path_grf = os.path.join(path_data, sub, ses, f)

            file_cycles_kinematics = os.listdir(os.path.join(path_data, 'results_pool/'))
            # Keep only the files having the sting '_ik.sto', and '_emg.sto' in their names:
            # - filter out everything else.
            file_cycles_kinematics = [
                file_cycle_kinematics for file_cycle_kinematics in file_cycles_kinematics \
                if '_ik_move.sto' in file_cycle_kinematics and os.path.isfile(os.path.join(path_data, 'results_pool/', file_cycle_kinematics))]
            file_cycles_kinematics.sort()

            # file_cycles_grfs = os.listdir(os.path.join(path_data, sub, ses))
            # # Keep only the files having the sting '_ik.sto', and '_emg.sto' in their names:
            # # - filter out everything else.
            # file_cycles_grfs = [
            #     file_cycles_grfs for file_cycles_grfs in file_cycles_grfs \
            #     if '_grf_move.xml' in file_cycles_grfs and os.path.isfile(os.path.join(path_data, sub, ses, file_cycles_grfs))]
            # file_cycles_grfs.sort()

            file_count = 0

            for file_cycle_kinematics in file_cycles_kinematics:
                if file_count >=10:
                    break
                file_cycle_stem = file_cycle_kinematics[:-12]
                # file_grfs_stem = file_cycles_grfs[:-13]

                path_emg              = os.path.join(path_data, 'results_pool/', file_cycle_stem + '_emg.sto')
                path_states           = os.path.join(path_data, 'results_pool/', file_cycle_stem + '_rad.sto')
                path_grf              = os.path.join(path_data, sub, ses, file_cycle_stem + '_grf_move.xml')
                path_cycle_kinematics = os.path.join(path_data, 'results_pool/', file_cycle_kinematics)
                path_solution         = os.path.join(path_data, sub, ses, 'guesses/initial/', file_cycle_stem + '_mi_emg.sto')
                


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

                        # path_data = '/home/lisca/biomec/data/20250306_dataset/sub03/ses20250306/'
                        # path_inputs = '/home/lisca/biomec/data/20250306_dataset/results_pool'
                        # cycle_kinematics = os.path.join(path_inputs, 'sub03_strifc_cond00000_speed02000_0001_0003_ik_move.sto')
                        # cycle_grf = os.path.join(path_data, 'sub03_strifc_cond00000_speed02000_0001_0003_grf_move.xml')
                        # path_musclepaths = os.path.join(path_data, 'sub03_FunctionBasedPathSet.xml')
                        # path_model = os.path.join(path_data,'sub03_SCALED_addbio_free_subt_mtp.osim')
                        # file_cycle_stem = cycle_kinematics[:-12]
                        # setup_data_stem = cycle_grf[:-13]

                        # path_emg              = os.path.join(path_inputs, file_cycle_stem + '_emg.sto')
                        # path_states           = os.path.join(path_inputs, file_cycle_stem + '_rad.sto')
                        # path_grf              = os.path.join(path_data, setup_data_stem + '_grf_move.xml')
                        # path_cycle_kinematics = cycle_kinematics
                        # path_solution         = os.path.join(path_data, 'guesses/initial/', file_cycle_stem + '_mi_emg.sto')
                                    


                    solveMocoInverseKinematicsEMG(
                        path_model          = path_model,
                        path_musclepaths    = path_musclepaths,
                        path_grf            = path_grf,
                        path_kinematics     = path_cycle_kinematics,
                        path_states         = path_states,
                        time_initial        = None,
                        time_final          = None,
                        force_reserve       = 200.0,
                        mesh_interval       = 0.02,
                        path_solution_guess = None,
                        path_solution       = path_solution,
                        path_emg            = path_emg)

                    file_count += 1
