import sys
import os
import opensim as osim
import osimFunctions as helper
import test_moco as runSim

path_test = '/home/lisca/biomec/data/20240904_dataset/sub01/ses20240904'

helper.kinematicsToStates(kinematicsFileName = os.path.join(path_test, 'sub01_strifs_cond00000_speedp0300_0001_0001_ik.sto'),
    #kinematicsFileName = os.path.join(path_test, 'test_coordinates.sto'),
                       outputFileName = os.path.join(path_test, 'sub01_strifs_cond00000_speedp0300_0001_0001_rad.sto'),
                       osimModelFileName = os.path.join(path_test, 'sub01_SCALED.osim'),
                       inDegrees = True, outDegrees = False, filtFreq=None)

runSim.runMoco(
                            path_model=os.path.join(path_test, 'sub01_SCALED.osim'),
                            path_grf_xml=os.path.join(path_test, 'sub01_strifs_cond00000_speedp0300_0001_0001_grf.xml'),
                            path_kinematics=os.path.join(path_test, 'sub01_strifs_cond00000_speedp0300_0001_0001_ik.sto'),
                            path_states=os.path.join(path_test, 'sub01_strifs_cond00000_speedp0300_0001_0001_rad.sto'),
                            path_emg=os.path.join(path_test, 'sub01_strifs_cond00000_speedp0300_0001_0001_emg.sto'),
                            time_initial=None,
                            time_final=None,
                            path_solution_guess=None,
                            path_solution=os.path.join(path_test, 'test_states_track.sto'),
                            path_solution_kinematics=os.path.join(path_test, 'test_coordinates_track.sto'),
                            path_setup_file=os.path.join(path_test, 'test_coordinates_setup.xml'),
                            path_log=os.path.join(path_test, 'test_coordinates.log'),
                            force_reserve=1.0,
                            mesh_interval=0.10
)