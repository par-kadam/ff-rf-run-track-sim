function config = load_config()

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
% addpath(genpath([getenv('WS_BIOMEC') '/install/btk/']));
addpath('C:\rl\btk\');

% The path to the 'lib' directory of the this repository.
addpath(genpath('../lib/'));

config.path_root = [pwd '/../'];
addpath(config.path_root);

% Path to the root directory of the dataset.
% 20240126_dataset_lutetium
% config.path_dataset = '/home/lisca/ssd/nvme0n1p9/lisca/datasets/20240126_dataset_lutetium/';
% 20240525_opencap_match_pilot
% config.path_dataset = '/home/lisca/ssd/nvme0n1p9/lisca/datasets/20240525_opencap_match_pilot/';
% config.path_dataset = 'C:\rl\20240525_opencap_match_pilot\';
% config.path_dataset = 'C:\rl\20240601_dataset\';
% config.path_dataset = '/home/lisca/ssd/nvme0n1p9/lisca/datasets/20240601_dataset/';
config.path_dataset = 'C:\rl\20240601_dataset\20240822_dataset\';
% Subject's metadata
% todo(lisca): Find a way to tell Vicon to export it automatically.

% 20240126_dataset_lutetium
% name_subjects    = {'sub01', 'sub02', 'sub03', 'sub04', 'sub05', 'sub06', 'sub07',};
% process_subjects = [   false,    true,    true,    true,   false,    true,   false,];
% process_subjects = [   false,   true,   false,   false,   false,   false,   false,];
% masses           = [   58.0,    58.0,    72.0,    78.0,    83.0,    85.0,    85.0,];
% heights          = [ 1820.0,  1820.0,  1760.0,  1780.0,  1740.0,  1820.0,  1820.0,];

% % 20240525_opencap_match_pilot
% name_subjects    = {'sub01',};
% process_subjects = [  false,];
% process_subjects = [   true,];
% masses           = [   78.0,];
% heights          = [ 1780.0,];
% static_trials    = {'_0001',};

% 20240601_dataset
% name_subjects    = {'sub01', 'sub02', 'sub03',};
% process_subjects = [  true,   true,      true,];
% masses           = [   78.0,    78.0,    85.0,];
% heights          = [ 1780.0,  1780.0,  1820.0,];
% trial_static     = {  '0001',  '0008',  '0002',};

% 20240822_dataset
name_subjects      = {'sub01'};
process_subjects   = [true];
masses             = [78.0];
heights            = [178.0];
trial_static       = { '0010' };

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
% %     'C_28';
%     'LASI'; 'LPSI'; 'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'; ...
%     'RASI'; 'RPSI'; 'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE';};

config.marker_data.marker_list = {
    'LASI'; 'LPSI'; 'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE'; 'LTOE'; 'LMT5'; 'LMAnkle'; 'LSH1'; 'LSH2'; 'LTH1'; 'LTH2'; ...
    'RASI'; 'RPSI'; 'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'; 'RMT5'; 'RMAnkle'; 'RSH1'; 'RSH2'; 'RTH1'; 'RTH2'; ...
    'RHJC'; 'LHJC'; ...
    'RD1M'; 'RMTB5'; 'RMTB1'; 'RSTL'; 'RLCA'; 'RPCA'; ...
    'LD1M'; 'LMTB5'; 'LMTB1'; 'LSTL'; 'LLCA'; 'LPCA';};

config.marker_data.calc_marker_list = {
    'LTIO'; 'LFEO'; 'LFEP'; ...
    'RFEP'; 'RFEO'; 'RTIO'};

% Filter the marker data.
config.marker_data.moving_average_window_size = 10;

% EMG

% Rigth leg had EMG sensors attached on it.

% 20240126_dataset_lutetium
% name_muscles   = {      'bicep_femoris',      'tensor_fascia_latae',      'tibialis_anterior', 'gastrocnemius_medial',     'glutenus_medius',      'rectus_femoris',};
% channels.sub01 = {           'Channel1',                'Channel21',              'Channel18',             'Channel4',           'Channel20',           'Channel22',};
% channels.sub02 = {           'Channel1',                'Channel21',              'Channel18',             'Channel4',           'Channel20',           'Channel22',};
% channels.sub03 = {           'Channel1',                'Channel21',              'Channel18',             'Channel4',           'Channel20',           'Channel22',};
% channels.sub04 = {'RightBiceps_Femoris', 'RightTensor_Fascia_Latae', 'RightTibialis_Anterior',   'RightGastrocnemius', 'RightGluteus_Medius', 'RightRectus_Femoris',};
% channels.sub05 = {           'Channel1',                'Channel21',              'Channel18',             'Channel4',           'Channel20',           'Channel22',};
% channels.sub06 = {'RightBiceps_Femoris', 'RightTensor_Fascia_Latae', 'RightTibialis_Anterior',   'RightGastrocnemius', 'RightGluteus_Medius', 'RightRectus_Femoris',};
% channels.sub07 = {'RightBiceps_Femoris', 'RightTensor_Fascia_Latae', 'RightTibialis_Anterior',   'RightGastrocnemius', 'RightGluteus_Medius', 'RightRectus_Femoris',};

% % 20240525_opencap_match_pilot
% name_muscles   = {      'bicep_femoris',      'tensor_fascia_latae',      'tibialis_anterior', 'gastrocnemius_medial',     'glutenus_medius',      'rectus_femoris',};
% channels.sub01 = {           'Channel1',                'Channel21',              'Channel18',             'Channel4',           'Channel20',           'Channel22',};

% 20240601_dataset
name_muscles   = {'right_gastrocnemius_medialis', 'right_gastrocnemius_lateralis', 'right_soleus', 'right_tibialis_anterior', 'right_vastus_medialis', 'right_vastus_lateralis', 'right_rectus_femoris',};
channels.sub01 = {               'R_Med__Gastro',                 'R_Lat__Gastro',     'R_Soleus',              'R_Tib_Ant_',                 'R_Vmo',                  'R_Vlo',        'R_Rectus_Fem_',};
channels.sub02 = {               'R_Med__Gastro',                 'R_Lat__Gastro',     'R_Soleus',              'R_Tib_Ant_',                 'R_Vmo',                  'R_Vlo',        'R_Rectus_Fem_',};
channels.sub03 = {               'R_Med__Gastro',                 'R_Lat__Gastro',     'R_Soleus',              'R_Tib_Ant_',                 'R_Vmo',                  'R_Vlo',        'R_Rectus_Fem_',};

maximum_contraction_trials.sub01 = ...
                 {                        '0201',                          '0301',         '0401',                    '0501',                  '0601',                    '0701',                '0801',};
maximum_contraction_trials.sub02 = ...
                 {                        '0201',                          '0301',         '0401',                    '0501',                  '0601',                    '0701',                '0801',};
maximum_contraction_trials.sub03 = ...
                 {                        '0201',                          '0301',         '0401',                    '0501',                  '0601',                    '0701',                '0801',};

config.emg_data.name_muscles     = name_muscles;
config.emg_data.subjects.muscles = containers.Map( ...
    name_subjects, {...
        containers.Map(config.emg_data.name_muscles, channels.sub01), ...
        containers.Map(config.emg_data.name_muscles, channels.sub02), ...
        containers.Map(config.emg_data.name_muscles, channels.sub03),});

config.emg_data.maximum_contractions = containers.Map( ...
    name_subjects, {...
        containers.Map(config.emg_data.name_muscles, maximum_contraction_trials.sub01), ...
        containers.Map(config.emg_data.name_muscles, maximum_contraction_trials.sub02), ...
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
config.threshold_grf               =  50; % 5/100.0 * max(data.fp_data.GRF_data(i).F(:,3));
config.length_min_gait_cycle       = 650;
config.export_only_first_n_cycles  =   2;

% 

% Path to the OpenSim model.
config.path_model = [config.path_root, 'vicon/model/'];

% define the model name - the task files etc must start with the same name
% config.file_model = 'lower_limb';
config.file_model = 'vicon2opensim2moco';

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
    'RightBicepsFemoris',
    'RightGastrocnemius',
    'RightRectusFemoris',
    'RightTibialisAnterior',
%     'RightTensorFasciaLatae',
%     'RightGluteusMedius',
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
