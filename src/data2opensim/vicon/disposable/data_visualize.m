function data_visualize(file_trial)

import org.opensim.modeling.*;

% Set the paths of the directories containg the resources to be used.
[path_root, path_dataset, marker_list, calc_marker_list, path_model, file_model] = data_pipeline_setup();


path_subjects = {
    'Adam/02.01/',
%     'Gheorghe/03.02/',
    };

conditions = {
    "bw_condn1000",
    "bw_condp0000",
    "fw_condp0000",
    "fw_condp1000"
    };

% Froce subject, trial and data.
path_subject = 'Adam/02.01/';
% path_subject = 'Gheorghe/03.02/';

condition = '_static';
condition = 'fw_condp0000';
% condition = 'bw_condp0000';
% condition = 'fw_condp1000';
% condition = 'bw_condn1000';

plot_grfs_raw           = false;
plot_markers_raw        = false;
plot_inverse_kinematics = false;
plot_moco_inverse       = false;

% plot_grfs_raw           = true;
% plot_markers_raw        = true;
% plot_inverse_kinematics = true;
plot_moco_inverse       = true;

%% Scan all subject directories and collect a list of all files.

for idx_subject = 1:length(path_subjects)
    disp([path_dataset path_subjects{idx_subject}]);

   % Get a list of all files matching the pattern
    files_all = dir(fullfile([path_dataset path_subjects{idx_subject}], '*.*'));
end

%% Visualize raw data: markers and forceplates

if plot_markers_raw
    % Create a pattern for dir function
    pattern_c3d = fullfile(path_dataset, path_subject, '*.c3d');

    % Get a list of all files matching the pattern
    files_c3d = dir(pattern_c3d);

    for i = 1:length(files_c3d)

        file_trial = files_c3d(i).name;
    %     disp(file_trial);

        if ~contains(file_trial, condition)
            continue
        end

        path_trial = [path_dataset path_subject file_trial];

        data_raw = btk_loadc3d(path_trial, 10);

        data_raw.marker_data = btk_sortc3d( ...
            data_raw.marker_data, marker_list, calc_marker_list);

        btk_animate_markers( ...
            data_raw.marker_data, data_raw.fp_data, 5);

%         btk_animate_markers(data_raw.marker_data);
    end
end

%% Visualzie the GRF data

if plot_grfs_raw

    figure_grfs = figure;
    set(figure_grfs, 'MenuBar', 'none');
    set(figure_grfs, 'ToolBar', 'none');

    % Iterate over all files

    pattern_grf = fullfile([path_dataset path_subject], '*_grf.sto');
    files_grf = dir(pattern_grf);
    for i = 1:length(files_grf)
        file_grf = files_grf(i).name;

        if ~contains(file_grf, condition)
            continue
        end

        % Consider only the OpenSim Moco solutions corresponding to the
        % selected trial and containing the '_mi' suffix into their name.
        path_grf = [path_dataset path_subject file_grf];
        disp(path_grf);

        info_path_grf = dir(path_grf);
        if info_path_grf.bytes == 0
            disp([file_grf ' file is empty']);
            continue;
        else
            disp(['Loading: ' file_grf]);
        end

        grf_data_table  = TimeSeriesTable(path_grf);

        figure(figure_grfs);
        hold on;
        title_grfs = strrep(file_grf, '_', '\_');
        title(title_grfs);
        legend(); % ('Location','northwest');

        grf_1_vy_label = string(grf_data_table.getColumnLabel(1));
        grf_1_vy_label = strrep(grf_1_vy_label, '_', '\_');
        grf_1_vy       = grf_data_table.getDependentColumnAtIndex(1).getAsMat();

        grf_2_vy_label = string(grf_data_table.getColumnLabel(10));
        grf_2_vy_label = strrep(grf_2_vy_label, '_', '\_');
        grf_2_vy       = grf_data_table.getDependentColumnAtIndex(10).getAsMat();

        plot(grf_1_vy, 'DisplayName', grf_1_vy_label);
        plot(grf_2_vy, 'DisplayName', grf_2_vy_label);
        drawnow;
    end

    disp('Click or press any key to plot the GRFs of the next gait cycle ...');
    w = waitforbuttonpress;
    close(figure_grfs);
end

%% Visualize OpenSim Inverse kinematics

if plot_inverse_kinematics
    % Create a pattern for dir function
    pattern_ik = fullfile([path_dataset path_subject], '*_ik.sto');

    % Get a list of all files matching the pattern
    files_ik = dir(pattern_ik);

    % Iterate over all files
    figure_iks = {};
    size_coordinates = 19;
    for idx_coordinate = 1:size_coordinates
        figure_iks{idx_coordinate} = figure;
        set(figure_iks{idx_coordinate}, 'MenuBar', 'none');
        set(figure_iks{idx_coordinate}, 'ToolBar', 'none');
    end

    for idx_file = 1:length(files_ik)

        file_ik = files_ik(idx_file).name;

        if ~contains(file_ik, condition)
            continue
        end

        % Consider only the OpenSim Moco solutions corresponding to the
        % selected trial and containing the '_mi' suffix into their name.
        path_ik = [path_dataset path_subject file_ik];
        disp(path_ik);

        info_path_ik = dir(path_ik);
        if info_path_ik.bytes == 0
            disp([file_ik ' file is empty!']);
            continue;
        end

        ik_data_table  = TimeSeriesTable(path_ik);
        ik_data        = osimTableToStruct(ik_data_table);
        ik_data_fields = fieldnames(ik_data);

        disp(ik_data_fields);

        for idx_coordinate =1:length(ik_data_fields)
            coordinate_name = ik_data_fields{idx_coordinate};

            figure(figure_iks{idx_coordinate});
            title_ik = [file_ik ' ' coordinate_name];
            title_ik = strrep(title_ik, '_', '\_');
            title(title_ik);
            hold on;
            plot(ik_data.(coordinate_name), 'DisplayName', 'rough');
            legend();
            drawnow;
        end

        disp('Press any space to process the next gait cycle ...')
        w = waitforbuttonpress;
        clf;
    end
end

%% Visualize OpenSim Moco Inverse

if plot_moco_inverse

    name_reserves = {
%         'a_forceset_reserve_jointset_hip_r_hip_flexion_r',
%         'a_forceset_reserve_jointset_hip_r_hip_adduction_r',
%         'a_forceset_reserve_jointset_hip_r_hip_rotation_r',
%         'a_forceset_reserve_jointset_walker_knee_r_knee_angle_r',
        'a_forceset_reserve_jointset_ankle_r_ankle_angle_r',
%         'a_forceset_reserve_jointset_hip_l_hip_flexion_l',
%         'a_forceset_reserve_jointset_hip_l_hip_adduction_l',
%         'a_forceset_reserve_jointset_hip_l_hip_rotation_l',
%         'a_forceset_reserve_jointset_walker_knee_l_knee_angle_l',
        'a_forceset_reserve_jointset_ankle_l_ankle_angle_l',
        };

    name_muscle_activation = {
%         'a_forceset_addbrev_r_activation',
%         'a_forceset_addlong_r_activation',
%         'a_forceset_addmagDist_r_activation',
%         'a_forceset_addmagIsch_r_activation',
%         'a_forceset_addmagMid_r_activation',
%         'a_forceset_addmagProx_r_activation',
%         'a_forceset_bflh_r_activation',
%         'a_forceset_bfsh_r_activation',
%         'a_forceset_edl_r_activation',
%         'a_forceset_ehl_r_activation',
%         'a_forceset_fdl_r_activation',
%         'a_forceset_fhl_r_activation',
%         'a_forceset_gaslat_r_activation',
%         'a_forceset_gasmed_r_activation',
%         'a_forceset_glmax1_r_activation',
%         'a_forceset_glmax2_r_activation',
%         'a_forceset_glmax3_r_activation',
%         'a_forceset_glmed1_r_activation',
%         'a_forceset_glmed2_r_activation',
        'a_forceset_glmed3_r_activation',
%         'a_forceset_glmin1_r_activation',
%         'a_forceset_glmin2_r_activation',
%         'a_forceset_glmin3_r_activation',
%         'a_forceset_grac_r_activation',
%         'a_forceset_iliacus_r_activation',
%         'a_forceset_perbrev_r_activation',
%         'a_forceset_perlong_r_activation',
%         'a_forceset_piri_r_activation',
%         'a_forceset_psoas_r_activation',
        'a_forceset_recfem_r_activation',
%         'a_forceset_sart_r_activation',
%         'a_forceset_semimem_r_activation',
%         'a_forceset_semiten_r_activation',
%         'a_forceset_soleus_r_activation',
%         'a_forceset_tfl_r_activation',
%         'a_forceset_tibant_r_activation',
%         'a_forceset_tibpost_r_activation',
%         'a_forceset_vasint_r_activation',
%         'a_forceset_vaslat_r_activation',
%         'a_forceset_vasmed_r_activation',
%         'a_forceset_addbrev_l_activation',
%         'a_forceset_addlong_l_activation',
%         'a_forceset_addmagDist_l_activation',
%         'a_forceset_addmagIsch_l_activation',
%         'a_forceset_addmagMid_l_activation',
%         'a_forceset_addmagProx_l_activation',
%         'a_forceset_bflh_l_activation',
%         'a_forceset_bfsh_l_activation',
%         'a_forceset_edl_l_activation',
%         'a_forceset_ehl_l_activation',
%         'a_forceset_fdl_l_activation',
%         'a_forceset_fhl_l_activation',
%         'a_forceset_gaslat_l_activation',
%         'a_forceset_gasmed_l_activation',
%         'a_forceset_glmax1_l_activation',
%         'a_forceset_glmax2_l_activation',
%         'a_forceset_glmax3_l_activation',
%         'a_forceset_glmed1_l_activation',
%         'a_forceset_glmed2_l_activation',
%         'a_forceset_glmed3_l_activation',
%         'a_forceset_glmin1_l_activation',
%         'a_forceset_glmin2_l_activation',
%         'a_forceset_glmin3_l_activation',
%         'a_forceset_grac_l_activation',
%         'a_forceset_iliacus_l_activation',
%         'a_forceset_perbrev_l_activation',
%         'a_forceset_perlong_l_activation',
%         'a_forceset_piri_l_activation',
%         'a_forceset_psoas_l_activation',
%         'a_forceset_recfem_l_activation',
%         'a_forceset_sart_l_activation',
%         'a_forceset_semimem_l_activation',
%         'a_forceset_semiten_l_activation',
%         'a_forceset_soleus_l_activation',
%         'a_forceset_tfl_l_activation',
%         'a_forceset_tibant_l_activation',
%         'a_forceset_tibpost_l_activation',
%         'a_forceset_vasint_l_activation',
%         'a_forceset_vaslat_l_activation',
%         'a_forceset_vasmed_l_activation',
        };

    name_muscle_forces = {
%         'a_forceset_tau_pelvis_tilt',
%         'a_forceset_tau_pelvis_list',
%         'a_forceset_tau_pelvis_rotation',
%         'a_forceset_tau_pelvis_tx',
%         'a_forceset_tau_pelvis_ty',
%         'a_forceset_tau_pelvis_tz',

%         'a_forceset_addbrev_r',
%         'a_forceset_addlong_r',
%         'a_forceset_addmagDist_r',
%         'a_forceset_addmagIsch_r',
%         'a_forceset_addmagMid_r',
%         'a_forceset_addmagProx_r',
%         'a_forceset_bflh_r',
%         'a_forceset_bfsh_r',
%         'a_forceset_edl_r',
%         'a_forceset_ehl_r',
%         'a_forceset_fdl_r',
%         'a_forceset_fhl_r',
%         'a_forceset_gaslat_r',
%         'a_forceset_gasmed_r',
%         'a_forceset_glmax1_r',
%         'a_forceset_glmax2_r',
%         'a_forceset_glmax3_r',
%         'a_forceset_glmed1_r',
%         'a_forceset_glmed2_r',
%         'a_forceset_glmed3_r',
%         'a_forceset_glmin1_r',
%         'a_forceset_glmin2_r',
%         'a_forceset_glmin3_r',
%         'a_forceset_grac_r',
%         'a_forceset_iliacus_r',
%         'a_forceset_perbrev_r',
%         'a_forceset_perlong_r',
%         'a_forceset_piri_r',
%         'a_forceset_psoas_r',
%         'a_forceset_recfem_r',
%         'a_forceset_sart_r',
%         'a_forceset_semimem_r',
%         'a_forceset_semiten_r',
%         'a_forceset_soleus_r',
%         'a_forceset_tfl_r',
%         'a_forceset_tibant_r',
%         'a_forceset_tibpost_r',
%         'a_forceset_vasint_r',
%         'a_forceset_vaslat_r',
%         'a_forceset_vasmed_r',
%         'a_forceset_addbrev_l',
%         'a_forceset_addlong_l',
%         'a_forceset_addmagDist_l',
%         'a_forceset_addmagIsch_l',
%         'a_forceset_addmagMid_l',
%         'a_forceset_addmagProx_l',
%         'a_forceset_bflh_l',
%         'a_forceset_bfsh_l',
%         'a_forceset_edl_l',
%         'a_forceset_ehl_l',
%         'a_forceset_fdl_l',
%         'a_forceset_fhl_l',
%         'a_forceset_gaslat_l',
%         'a_forceset_gasmed_l',
%         'a_forceset_glmax1_l',
%         'a_forceset_glmax2_l',
%         'a_forceset_glmax3_l',
%         'a_forceset_glmed1_l',
%         'a_forceset_glmed2_l',
%         'a_forceset_glmed3_l',
%         'a_forceset_glmin1_l',
%         'a_forceset_glmin2_l',
%         'a_forceset_glmin3_l',
%         'a_forceset_grac_l',
%         'a_forceset_iliacus_l',
%         'a_forceset_perbrev_l',
%         'a_forceset_perlong_l',
%         'a_forceset_piri_l',
%         'a_forceset_psoas_l',
%         'a_forceset_recfem_l',
%         'a_forceset_sart_l',
%         'a_forceset_semimem_l',
%         'a_forceset_semiten_l',
%         'a_forceset_soleus_l',
%         'a_forceset_tfl_l',
%         'a_forceset_tibant_l',
%         'a_forceset_tibpost_l',
%         'a_forceset_vasint_l',
%         'a_forceset_vaslat_l',
%         'a_forceset_vasmed_l',
        };

    figure_reserves = {};
    for idx_figure = 1:length(name_reserves)
        figure_reserves{idx_figure} = figure;
            set(figure_reserves{idx_figure}, 'MenuBar', 'none');
            set(figure_reserves{idx_figure}, 'ToolBar', 'none');
    end

    figure_muscle_activations = {};
    for idx_figure = 1:length(name_muscle_activation)
        figure_muscle_activations{idx_figure} = figure;
            set(figure_muscle_activations{idx_figure}, 'MenuBar', 'none');
            set(figure_muscle_activations{idx_figure}, 'ToolBar', 'none');
    end

    figure_muscle_forces = {};
    for idx_figure = 1:length(name_muscle_forces)
        figure_muscle_forces{idx_figure} = figure;
            set(figure_muscle_forces{idx_figure}, 'MenuBar', 'none');
            set(figure_muscle_forces{idx_figure}, 'ToolBar', 'none');
    end

    % Create a pattern for dir function
    pattern_mi = fullfile([path_dataset path_subject], '*_mi.sto');

    % Get a list of all files matching the pattern
    files_mi = dir(pattern_mi);

    % Iterate over all files
    for i = 1:length(files_mi)

        file_mi = files_mi(i).name;

        if ~contains(file_mi, condition)
            continue
        end

        % Consider only the OpenSim Moco solutions corresponding to the
        % selected trial and containing the '_mi' suffix into their name.
        path_mi = [path_dataset path_subject file_mi];
        disp(path_mi);

        info_path_mi = dir(path_mi);
        if info_path_mi.bytes == 0
            disp([file_mi ' file is empty!']);
            continue;
        end

        mi_data_table  = TimeSeriesTable(path_mi);
        mi_data        = osimTableToStruct(mi_data_table);
        mi_data_fields = fieldnames(mi_data);

        for idx_filed = 1:length(mi_data_fields)
            name_field = mi_data_fields{idx_filed};

            % Plot reservers
            for idx_reserve = 1:length(name_reserves)
                if contains(name_field, name_reserves{idx_reserve})
                    figure(figure_reserves{idx_reserve});
                    hold on;
                    title_reserves = ['Reserves ' file_mi];
                    title_reserves = strrep(title_reserves, '_', '\_');
                    title(title_reserves);
                    legend_reserves = strrep(name_field, 'a_forceset_reserve_jointset_', '');
                    legend_reserves = strrep(legend_reserves, 'ankle_', '');
                    legend_reserves = strrep(legend_reserves, '_', '\_');
                    legend_reserves = ['reserve\_' legend_reserves];
                    legend();
                    plot(mi_data.(name_field), 'DisplayName', legend_reserves);
                    drawnow;
                end
            end

            % Plot muscle activation
            for idx_muscle = 1:length(name_muscle_activation)
                if contains(name_field, name_muscle_activation{idx_muscle})
                    figure(figure_muscle_activations{idx_muscle});
                    hold on;
                    title_activation = ['Muscle activation ' file_mi];
                    title_activation = strrep(title_activation, '_', '\_');
                    title(title_activation);
                    legend_activation = strrep(name_field, 'a_forceset_', '');
                    legend_activation = strrep(legend_activation, '_activation', '');
                    legend_activation = strrep(legend_activation, '_', '\_');
                    legend_activation = ['activation\_' legend_activation];
                    legend();
                    plot(mi_data.(name_field), 'DisplayName', legend_activation);
                    drawnow;
                end
            end

            % Plot muscle forces
            for idx_muscle = 1:length(name_muscle_forces)
                if ~contains(name_field, "activation") && ...
                        contains(name_field, name_muscle_forces{idx_muscle})
                    figure(figure_muscle_forces{idx_muscle});
                    hold on;
                    title_force = ['Muscle forces ' file_mi];
                    title_force = strrep(title_force, '_', '\_');
                    title(title_force);
                    legend_force = strrep(name_field, 'a_forceset_', '');
                    legend_force = strrep(legend_force, '_', '\_');
                    legend_force = ['force\_' legend_force];
                    legend();
                    plot(mi_data.(name_field), 'DisplayName', legend_force);
                    drawnow;
                end
            end
        end

        disp('Click or press any space to plot the next gait cycle ...')
        w = waitforbuttonpress;

%         for idx_figure = 1:length(figure_reserves)
%             figure(figure_reserves{idx_figure});
%             clf;
%         end
% 
%         for idx_figure = 1:length(figure_muscle_activations)
%             figure(figure_muscle_activations{idx_figure});
%             clf;
%         end
% 
%         for idx_figure = 1:length(figure_muscle_forces)
%             figure(figure_muscle_forces{idx_figure});
%             clf;
%         end
    end

    disp('Click or press any key to close all figure and exit ...')
    w = waitforbuttonpress;

    for idx_figure = 1:length(figure_reserves)
        close(figure_reserves{idx_figure});
    end

    for idx_figure = 1:length(figure_muscle_activations)
        close(figure_muscle_activations{idx_figure});
    end

    for idx_figure = 1:length(figure_muscle_forces)
        close(figure_muscle_forces{idx_figure});
    end
end

