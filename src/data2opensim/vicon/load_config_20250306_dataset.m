function config = load_config_20240904_dataset()

clc

config. visualize.pause                     = 0.0;
%
config.visualize.treadmill_data             = false;
config.visualize.force_plate_data           = false;
config.visualize.emg_data                   = false;
config.visualize.marker_data                = false;
%
config.visualize.subject                    = false;
config.visualize.session                    = false;
config.visualize.condition                  = false;
config.visualize.trial                      = false;
config.visualize.gait_cycle                 = false;
%
config.visualize.raw_data                   = false;
config.visualize.filter_and_resample_data   = false;
config.visualize.scale_model_to_data        = false;
config.visualize.extract_gait_cycles        = false;
config.visualize.compute_ik_and_ik          = false;
config.visualize.perform_all_normalizations = false;
config.visualize.compute_means_and_stds     = false;

% config. visualize.pause                     = 1.0;
% 
% config.visualize.treadmill_data             = true;
% config.visualize.force_plate_data           = true;
% config.visualize.emg_data                   = true;
% config.visualize.marker_data                = true;
% 
% config.visualize.subject                    = true;
% config.visualize.session                    = true;
% config.visualize.condition                  = true;
% config.visualize.trial                      = true;
% config.visualize.gait_cycle                 = true;
% 
% config.visualize.raw_data                   = true;
% config.visualize.filter_and_resample_data   = true;
% config.visualize.scale_model_to_data        = true;
% config.visualize.compute_ik_and_ik          = true;
% config.visualize.extract_gait_cycles        = true;
% config.visualize.perform_all_normalizations = true;
% config.visualize.compute_means_and_stds     = true;

% The 'install' directory of the 'biomec' work space.
addpath(genpath([getenv('WS_BIOMEC') '/install/btk/']));

% The path to the 'lib' directory of the this repository.
addpath(genpath('../lib/'));

config.path_root = [pwd '/../'];
addpath(config.path_root);

% Path to the root directory of the dataset.
% config.path_dataset = '/home/lisca/biomec/data/20240904_dataset/';
config.path_dataset = '/home/lisca/biomec/data/20250306_dataset/';
config.path_moco_inputs = '/home/lisca/biomec/data/20250306_dataset/results_pool/';

% Subject's metadata
% todo(lisca): Find a way to tell Vicon to export it automatically.

% 20240601_dataset
%name_subjects    = {'sub01',};
%process_subjects = [   true];
%masses           = [   78.0,];
%heights          = [ 1780.0,];
%trial_static     = { '0004',};

% 20250306_dataset
name_subjects      = {'sub03',};
process_subjects   = [    true];
masses             = [   58.0,];
heights            = [  178.0,];
trial_static       = { '0003',};

config.subjects.names        = name_subjects;
config.subjects.process      = process_subjects;
config.subjects.masses       = containers.Map(name_subjects, masses);
config.subjects.heights      = containers.Map(name_subjects, heights);
config.subjects.trial_static = containers.Map(name_subjects, trial_static);

% Filter the force plate data.
config.force_plate_data.cutoffFrequency  = 15; % Cutoff frequency in Hz
config.force_plate_data.order            =  4; % Filter order

% Keep only the following markers.
% config.marker_data.marker_list = {
%    'C_28'; 'C_29'; ...
%    'LASI'; 'LPSI'; 'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'; ...
%    'RASI'; 'RPSI'; 'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE';};

config.marker_data.marker_list = {
    'LASI'; 'LPSI'; 'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'; 'LMAnkle'; 'LSH1'; 'LSH2'; 'LSH3'; 'LTH1'; 'LTH2'; 'LTH3'; ...
    'RASI'; 'RPSI'; 'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'; 'RMT5'; 'RMAnkle'; 'RSH1'; 'RSH2'; 'RTH1'; 'RTH2'; 'RTH3'; ...
    'Pelvis_RightThigh_score'; 'Pelvis_LeftThigh_score'; ...
    'RD1M'; 'RMT5'; 'RSTL'; 'RLCA'; ...
    'LD1M'; 'LMT5'; 'LSTL'; 'LLCA';
    'LShoulder'; 'LClavicle'; 'C7'; 'LBiceps'; 'LElbow'; 'LMElbow'; 'LForearm'; 'LWrist'; 'RShoulder'; 'RBiceps'; 'RElbow'; 'RMElbow'; 'RForearm'; 'RWrist';}

config.marker_data.calc_marker_list = {
    'LTIO'; 'LFEO'; 'LFEP'; ...
    'RFEP'; 'RFEO'; 'RTIO'};

% Filter the marker data.
config.marker_data.moving_average_window_size = 10;

% EMG

% Rigth leg had EMG sensors attached on it.
name_muscles   = {                 'right_soleus', 'right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_peroneus_longus', 'right_tibialis_anterior', 'right_biceps_femoris', 'right_rectus_femoris', 'right_tensor_fascia_latae'};
channels.sub03 = {                     'R_Soleus',                'R_Med__Gastro',                 'R_Lat__Gastro',            'R_Peroneus',              'R_Tib_Ant_',        'R_Biceps_Fem_',        'R_Rectus_Fem_',              'R_Tensor_Fl_'};

maximum_contraction_trials.sub03 = ...
                 {                        '0901',                          '1001',                         '1101',                    '1201',                  '1301',                 '1401',                '1501',                       '1601'};

config.emg_data.name_muscles     = name_muscles;
config.emg_data.subjects.muscles = containers.Map( ...
    name_subjects, {...
        containers.Map(config.emg_data.name_muscles, channels.sub03),});

config.emg_data.maximum_contractions = containers.Map( ...
    name_subjects, {...
        containers.Map(config.emg_data.name_muscles, maximum_contraction_trials.sub03),});

% Gait cycle
% Left  foot stept on the forceplate 1
% Right food stept on the forceplate 2

% find when feet are on each forceplate and define the trial events as
% the first foot contact and the final foot contact
% -> when GRF's vertical component exceeds 5% of its maximum.
% -> when the GRF's vertical component exceeds 20Newtons.
%    - according to Adam this is the right criteria, not 5%.
%    - 20 Newtons does work properly on Adam's trials sub04.
%      - during swing there are values higher than 22.0 Newtons.
config.threshold_grf               =  30; % 5/100.0 * max(data.fp_data.GRF_data(i).F(:,3));
config.length_min_gait_cycle       = 700;
config.length_max_gait_cycle       = 900;
config.export_only_first_n_cycles  =  50;

% 

% Path to the OpenSim model.
config.path_model = [config.path_root, 'vicon/model/'];

% define the model name - the task files etc must start with the same name
% config.file_model = 'lower_limb';
config.file_model = 'vicon2opensim2moco';

config.path_model_physics_fitted = fullfile(pwd, "../../../data/20250306_dataset/sub03/sub03_addbiomechanics/Models/match_markers_and_physics.osim");
conifg.update_model_path = "../../../data/20250306_dataset/sub03/ses20250306/";
config.results_directory = [config.path_dataset 'sub03/ses20250306/'];


% GeometryPath for Matlab;
% Note: This is not enought for the external executables which the script
% calls. The inverse kinematics tool called in this way:
%  com = ['opensim-cmd run-tool ' data.Name '_Setup_InverseKinematics.xml'];
%  system(com);
% will stills spit warnings that it did not find files (meshes).
% path_geometry_core =   '/home/lisca/biomec/src/opensim-core/Bindings/Java/Matlab/OpenSenseExample/Geometry/';
% path_geometry_models = '/home/lisca/biomec/src/opensim-models/Geometry/';

% org.opensim.modeling.ModelVisualizer.addDirToGeometrySearchPaths(path_geometry_core);
% org.opensim.modeling.ModelVisualizer.addDirToGeometrySearchPaths(path_geometry_models);

config.name_emg_muscles = {
%     'Channel0',
    'R_Soleus',
    'R_Med__Gastro',
    'R_Lat__Gastro',
    'R_Peroneus',
    'R_Tib_Ant_',
    'R_Biceps_Fem_',
    'R_Rectus_Fem_',
    'R_Tensor_Fl_'
};

config.name_moco_reserves = {
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

config.name_moco_muscle_activation = {
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

config.name_moco_muscle_forces = {
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
