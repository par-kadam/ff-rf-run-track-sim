function [ ...
    path_root, path_dataset, ...
    marker_list, calc_marker_list, ...
    path_model, file_model, ...
    name_emg_muscles, ...
    name_moco_reserves, name_moco_muscle_activation, name_moco_muscle_forces ...
    ] = data_pipeline_setup()

clc

% The 'install' directory of the 'biomec' directory structure.
addpath(genpath('./../../../install/btk/'));

% The path to the 'vicon2opensim/tools/' of the 'athlete' repository.
addpath(genpath('../lib/'));

path_root = '/home/lisca/biomec/sim/src/data2opensim/vicon/';
addpath(path_root);

% Path to MoCap data.
path_dataset = '/home/lisca/biomec/data/dataset_lutetium/vicon/';

% Path to the OpenSim model.
path_model = [pwd, '/model/'];

% define the model name - the task files etc must start with the same name
% file_model = 'lower_limb';
file_model = 'vicon2opensim2moco';

% GeometryPath for Matlab;
% Note: This is not enought for the external executables which the script
% calls. The inverse kinematics tool called in this way:
%  com = ['opensim-cmd run-tool ' data.Name '_Setup_InverseKinematics.xml'];
%  system(com);
% will stills spit warnings that it did not find files (meshes).
path_geometry_core =   '/home/lisca/biomec/sim/src/opensim-core/Bindings/Java/Matlab/OpenSenseExample/Geometry/';
path_geometry_models = '/home/lisca/biomec/sim/src/opensim-models/Geometry/';

org.opensim.modeling.ModelVisualizer.addDirToGeometrySearchPaths(path_geometry_core);
org.opensim.modeling.ModelVisualizer.addDirToGeometrySearchPaths(path_geometry_models);

% The markers which will be used.
marker_list = {
    'LASI'; 'RASI'; 'LPSI'; 'RPSI'; ...
    'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'; ...
    'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'};

calc_marker_list = {
    'LTIO'; 'LFEO'; 'LFEP'; ...
    'RFEP'; 'RFEO'; 'RTIO'};

name_emg_muscles = {
%     'Channel0',
    'RightBicepsFemoris',
    'RightGastrocnemius',
    'RightRectusFemoris',
    'RightTibialisAnterior',
%     'RightTensorFasciaLatae',
%     'RightGluteusMedius',
}

name_moco_reserves = {
%     'a_forceset_reserve_jointset_hip_r_hip_flexion_r',
%     'a_forceset_reserve_jointset_hip_r_hip_adduction_r',
%     'a_forceset_reserve_jointset_hip_r_hip_rotation_r',
    'a_forceset_reserve_jointset_walker_knee_r_knee_angle_r',
    'a_forceset_reserve_jointset_ankle_r_ankle_angle_r',
%     'a_forceset_reserve_jointset_hip_l_hip_flexion_l',
%     'a_forceset_reserve_jointset_hip_l_hip_adduction_l',
%     'a_forceset_reserve_jointset_hip_l_hip_rotation_l',
%     'a_forceset_reserve_jointset_walker_knee_l_knee_angle_l',
%     'a_forceset_reserve_jointset_ankle_l_ankle_angle_l',
    };

name_moco_muscle_activation = {
%     'a_forceset_addbrev_r_activation',
%     'a_forceset_addlong_r_activation',
%     'a_forceset_addmagDist_r_activation',
%     'a_forceset_addmagIsch_r_activation',
%     'a_forceset_addmagMid_r_activation',
%     'a_forceset_addmagProx_r_activation',
    'a_forceset_bflh_r_activation', % moco paper uses the short head
%     'a_forceset_bfsh_r_activation',
%     'a_forceset_edl_r_activation',
%     'a_forceset_ehl_r_activation',
%     'a_forceset_fdl_r_activation',
%     'a_forceset_fhl_r_activation',
%     'a_forceset_gaslat_r_activation',
    'a_forceset_gasmed_r_activation',
%     'a_forceset_glmax1_r_activation',
%     'a_forceset_glmax2_r_activation',
%     'a_forceset_glmax3_r_activation',
%     'a_forceset_glmed1_r_activation',
%     'a_forceset_glmed2_r_activation',
%     'a_forceset_glmed3_r_activation',
%     'a_forceset_glmin1_r_activation',
%     'a_forceset_glmin2_r_activation',
%     'a_forceset_glmin3_r_activation',
%     'a_forceset_grac_r_activation',
%     'a_forceset_iliacus_r_activation',
%     'a_forceset_perbrev_r_activation',
%     'a_forceset_perlong_r_activation',
%     'a_forceset_piri_r_activation',
%     'a_forceset_psoas_r_activation',
    'a_forceset_recfem_r_activation',
%     'a_forceset_sart_r_activation',
%     'a_forceset_semimem_r_activation',
%     'a_forceset_semiten_r_activation',
%     'a_forceset_soleus_r_activation',
%     'a_forceset_tfl_r_activation',
    'a_forceset_tibant_r_activation',
%     'a_forceset_tibpost_r_activation',
%     'a_forceset_vasint_r_activation',
%     'a_forceset_vaslat_r_activation',
%     'a_forceset_vasmed_r_activation',
%     'a_forceset_addbrev_l_activation',
%     'a_forceset_addlong_l_activation',
%     'a_forceset_addmagDist_l_activation',
%     'a_forceset_addmagIsch_l_activation',
%     'a_forceset_addmagMid_l_activation',
%     'a_forceset_addmagProx_l_activation',
%     'a_forceset_bflh_l_activation',
%     'a_forceset_bfsh_l_activation',
%     'a_forceset_edl_l_activation',
%     'a_forceset_ehl_l_activation',
%     'a_forceset_fdl_l_activation',
%     'a_forceset_fhl_l_activation',
%     'a_forceset_gaslat_l_activation',
%     'a_forceset_gasmed_l_activation',
%     'a_forceset_glmax1_l_activation',
%     'a_forceset_glmax2_l_activation',
%     'a_forceset_glmax3_l_activation',
%     'a_forceset_glmed1_l_activation',
%     'a_forceset_glmed2_l_activation',
%     'a_forceset_glmed3_l_activation',
%     'a_forceset_glmin1_l_activation',
%     'a_forceset_glmin2_l_activation',
%     'a_forceset_glmin3_l_activation',
%     'a_forceset_grac_l_activation',
%     'a_forceset_iliacus_l_activation',
%     'a_forceset_perbrev_l_activation',
%     'a_forceset_perlong_l_activation',
%     'a_forceset_piri_l_activation',
%     'a_forceset_psoas_l_activation',
%     'a_forceset_recfem_l_activation',
%     'a_forceset_sart_l_activation',
%     'a_forceset_semimem_l_activation',
%     'a_forceset_semiten_l_activation',
%     'a_forceset_soleus_l_activation',
%     'a_forceset_tfl_l_activation',
%     'a_forceset_tibant_l_activation',
%     'a_forceset_tibpost_l_activation',
%     'a_forceset_vasint_l_activation',
%     'a_forceset_vaslat_l_activation',
%     'a_forceset_vasmed_l_activation',
    };

name_moco_muscle_forces = {
%     'a_forceset_tau_pelvis_tilt',
%     'a_forceset_tau_pelvis_list',
%     'a_forceset_tau_pelvis_rotation',
%     'a_forceset_tau_pelvis_tx',
%     'a_forceset_tau_pelvis_ty',
%     'a_forceset_tau_pelvis_tz',

%     'a_forceset_addbrev_r',
%     'a_forceset_addlong_r',
%     'a_forceset_addmagDist_r',
%     'a_forceset_addmagIsch_r',
%     'a_forceset_addmagMid_r',
%     'a_forceset_addmagProx_r',
%     'a_forceset_bflh_r',
%     'a_forceset_bfsh_r',
%     'a_forceset_edl_r',
%     'a_forceset_ehl_r',
%     'a_forceset_fdl_r',
%     'a_forceset_fhl_r',
%     'a_forceset_gaslat_r',
%     'a_forceset_gasmed_r',
%     'a_forceset_glmax1_r',
%     'a_forceset_glmax2_r',
%     'a_forceset_glmax3_r',
%     'a_forceset_glmed1_r',
%     'a_forceset_glmed2_r',
%     'a_forceset_glmed3_r',
%     'a_forceset_glmin1_r',
%     'a_forceset_glmin2_r',
%     'a_forceset_glmin3_r',
%     'a_forceset_grac_r',
%     'a_forceset_iliacus_r',
%     'a_forceset_perbrev_r',
%     'a_forceset_perlong_r',
%     'a_forceset_piri_r',
%     'a_forceset_psoas_r',
%     'a_forceset_recfem_r',
%     'a_forceset_sart_r',
%     'a_forceset_semimem_r',
%     'a_forceset_semiten_r',
%     'a_forceset_soleus_r',
%     'a_forceset_tfl_r',
%     'a_forceset_tibant_r',
%     'a_forceset_tibpost_r',
%     'a_forceset_vasint_r',
%     'a_forceset_vaslat_r',
%     'a_forceset_vasmed_r',
%     'a_forceset_addbrev_l',
%     'a_forceset_addlong_l',
%     'a_forceset_addmagDist_l',
%     'a_forceset_addmagIsch_l',
%     'a_forceset_addmagMid_l',
%     'a_forceset_addmagProx_l',
%     'a_forceset_bflh_l',
%     'a_forceset_bfsh_l',
%     'a_forceset_edl_l',
%     'a_forceset_ehl_l',
%     'a_forceset_fdl_l',
%     'a_forceset_fhl_l',
%     'a_forceset_gaslat_l',
%     'a_forceset_gasmed_l',
%     'a_forceset_glmax1_l',
%     'a_forceset_glmax2_l',
%     'a_forceset_glmax3_l',
%     'a_forceset_glmed1_l',
%     'a_forceset_glmed2_l',
%     'a_forceset_glmed3_l',
%     'a_forceset_glmin1_l',
%     'a_forceset_glmin2_l',
%     'a_forceset_glmin3_l',
%     'a_forceset_grac_l',
%     'a_forceset_iliacus_l',
%     'a_forceset_perbrev_l',
%     'a_forceset_perlong_l',
%     'a_forceset_piri_l',
%     'a_forceset_psoas_l',
%     'a_forceset_recfem_l',
%     'a_forceset_sart_l',
%     'a_forceset_semimem_l',
%     'a_forceset_semiten_l',
%     'a_forceset_soleus_l',
%     'a_forceset_tfl_l',
%     'a_forceset_tibant_l',
%     'a_forceset_tibpost_l',
%     'a_forceset_vasint_l',
%     'a_forceset_vaslat_l',
%     'a_forceset_vasmed_l',
    };