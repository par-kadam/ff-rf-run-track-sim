function data_process_subject()

%% Config the flags.
plot_emg_moco_comparison = false;

plot_emg_moco_comparison = true;

% Set the paths of the directories containg the resources to be used.
[path_root, path_dataset, ...
    marker_list, calc_marker_list, ...
    path_model, file_model ...
    name_emg_muscles, ...
    name_moco_reserves, name_moco_muscle_activation, name_moco_muscle_forces] = data_pipeline_setup();

% ---- OpenSim Matlab Tools example MoCap data. ------------------------- %
%
% data        = btk_loadc3d([mocap_path, 'AP_run_9.c3d'], 10);
% % data.Name   = 'AP';
% data.Mass   = 54.6;
% data.Height = 1650;
% 
% data_example = data
%
% ----------------------------------------------------------------------- %

% ---- Adamantium MoCap data. ------------------------------------------- %
%
% data        = btk_loadc3d([mocap_path, 'Robin_gait_3x_0.c3d'], 10);
% data.Name   = 'Robin';
% data.Mass   = 73.7;
% data.Height = 1780;
%
% ----------------------------------------------------------------------- %

% ---- M-Gait MoCap data. ----------------------------------------------- %
%
% data        = btk_loadc3d([mocap_path, 'Gheorghe_walking5.c3d'], 10);
% data.Name   = 'Gheorghe';
% data.Mass   = 84.8;
% data.Height = 1820;
%
% ----------------------------------------------------------------------- %

% ---- Semesterthesis Adam MoCap data. ---------------------------------- %
%
% data.Name   = 'Max';
% data.Mass   = 58;
% data.Height = 1820;
%
% data.Name   = 'Paul';
% data.Mass   = 58.0;
% data.Height = 1820;
% 
% data.Name   = 'Mudhit';
% data.Mass   = 72;
% data.Height = 1760;
% 
% data.Name   = 'Adam';
% data.Mass   = 78;
% data.Height = 1780;
% 
% data.Name   = 'JinYeong';
% data.Mass   = 83.0;
% data.Height = 1740;
%
% data.Name   = 'sub06';
% data.Mass   = 85.0;
% data.Height = 1820;
%
% data.Name   = 'sub07';
% data.Mass   = 85.0;
% data.Height = 1820;

% ----------------------------------------------------------------------- %

% Subject trials
path_subject = 'Adam/02.01/';

% Static trial
data.Name   = 'sub04';
data.Mass   = 78;
data.Height = 1780;

% Calibration
cd(path_root);
file_trial = 'sub04_static';
data_static_calibrate(path_subject, file_trial, data);

% Dynamic trials
file_trials = {
    'sub04_strifw_condp0000_speedp0100_0001', % 2-9, 10-24, 26, 27, 28
    % 'sub04_strifw_condp0000_speedp0100_0003',
    'sub04_stribw_condp0000_speedn0007_0003', % 2-6, 9-24, 26-29, 32,
%     'sub04_strifw_condp1000_speedp0100_0002', % 2-27    
%     'sub04_stribw_condn1000_speedn0007_0003', % ?24, ?25, 
};

% InverHse kinematics, inverse dynamics, etc.
for idx_trial = 1:length(file_trials)
    file_trial = file_trials{idx_trial};
    disp(file_trial);

    cd(path_root);
    [sub04_results_emg.(file_trial), ...
     sub04_results_moco.(file_trial)] = ...
        data_dynamic_process(path_subject, file_trial, data);
end

% ----------------------------------------------------------------------- %

% % Subject trials
% path_subject = 'Gheorghe/03.02/';
% 
% % Static trial
% data.Name   = 'sub06';
% data.Mass   = 85.0;
% data.Height = 1820;
% 
% % Calibrationh
% file_trial  = 'sub06_static';
% cd(path_root);
% data_static_calibrate(path_subject, file_trial, data);
% 
% % Dynamic trials
% file_trials = {
%     'sub06_strifw_condp0000_speedp0100_0001', % 2-12
% %     'sub06_strifw_condp0000_speedp0100_0002',
% %     'sub06_strifw_condp0000_speedp0100_0003',
%     'sub06_stribw_condp0000_speedn0060_0001', % 2-22
% 
% %     'sub06_strifw_condp1000_speedp0100_0001', % 2-25
% %     'sub06_stribw_condn1000_speedn0060_0004', % 2-21
% };
% 
% % Inverse kinematics, inverse dynamics, etc.
% for idx_trial = 1:length(file_trials)
%     file_trial = file_trials{idx_trial};
%     disp(file_trial);
% 
%     cd(path_root);
%     [sub06_results_emg.(file_trial), ...
%      sub06_results_moco.(file_trial)] = ...
%         data_dynamic_process(path_subject, file_trial, data);
% end

% ----------------------------------------------------------------------- %

if plot_emg_moco_comparison

    % Create a figure with four rows and four columns of subplots sharing the same axes
    figure;

    %
    position_title    = [0.5, 1.0, 0];
    size_font_title   = 22;
    size_font_xy      = 20;
    size_font_ticks   = 14;
    size_font_legend  = 18;
    weight_font       = 'Normal';
    
    % Define subplot sizes
    vertical_offset   = 0.06;
    horizontalSpacing = 0.045;  % Horizontal spacing between subplots
    verticalSpacing   = 0.065;  % Vertical spacing between subplots
    subplotWidth      = 0.20;  % Width of each subplot
    subplotHeight     = 0.125;  % Height of each subplot
    
    i = 1;
    for j = 1:4
        % Calculate position for each subplot
        position = [...
            horizontalSpacing + (j-1)*(subplotWidth+horizontalSpacing), ...
            vertical_offset   + (i-1)*(subplotHeight+verticalSpacing), ...
            subplotWidth, subplotHeight];
        % Create subplot
        emg_moco_axes{i, j} = axes('Position', position);

        % Remove frame around the subplot
        set(gca, 'box', 'off');
        % Hide unnecessary axis ticks and labels
        if j ~= 1
            set(gca, 'YTickLabel', []);
        end
        if i ~= 1
            set(gca, 'XTickLabel', []);
        end
    end
    
    % Define subplot sizes
    horizontal_offset = horizontalSpacing;
    vertical_offset   = verticalSpacing + subplotHeight + verticalSpacing;
    subplotWidth      = 0.2;  % Width of each subplot
    subplotHeight     = 0.3;  % Height of each subplot
    
    % Plot the subplots
    for i = 2:3
        for j = 1:4
            % Calculate position for each subplot
            position = [...
                horizontal_offset + (j-1)*(subplotWidth+horizontalSpacing), ...
                vertical_offset   + (i-2)*(subplotHeight+verticalSpacing), ...
                subplotWidth, subplotHeight];
            % Create subplot
            emg_moco_axes{i, j} = axes('Position', position);

            % Hide unnecessary axis ticks and labels
            if j ~= 1
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'XTickLabel', []);

        end
    end

    %% sub04

    file_trial_fw = 'sub04_strifw_condp0000_speedp0100_0001';
    file_trial_bw = 'sub04_stribw_condp0000_speedn0007_0003';
    
    %% gastrocnemius medial: RightGastrocnemius - a_forceset_gasmed_r_activation
    % forward flat set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{3, 1};    
    idx_muscle            = 2;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};    
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
    set(axis_muscle, 'fontsize', size_font_ticks);
%     xlabel(axis_muscle, 'percentage of gait cycle');
    ylabel(axis_muscle, 'activation', 'FontSize', size_font_xy);
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["gastrocnemius medial,", "forward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;
    % backward trial
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{3, 2};    
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["gastrocnemius medial,", "backward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;

    %% rectus femoris: RightRectusFemoris - a_forceset_recfem_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{3, 3};    
    idx_muscle            = 3;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["rectus femoris,", "forward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;
    % backward flat
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{3, 4};    
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;    
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill([1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);    
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["rectus femoris,", "backward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    legend(axis_muscle, ...
        'EMG mean', 'EMG activation', ...
        'Moco inverse mean', 'Moco inverse activation', 'FontSize', size_font_legend);
    drawnow;        

    %% tibialis anterior: RightTibialisAnterior - a_forceset_tibant_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{2, 1};    
    idx_muscle            = 4;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};    
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
    set(axis_muscle, 'fontsize', size_font_ticks);    
%     xlabel(axis_muscle, 'percentage of gait cycle');
    ylabel(axis_muscle, 'activation', 'FontSize', size_font_xy);
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["tibialis anterior,", "forward flat"];    
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;
    % backward trial
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{2, 2};
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];    
    title_emg = ["tibialis anterior,", "backward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;

    %% biceps femoris: RightBicepsFemoris - a_forceset_bflh_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{2, 3};
    idx_muscle            = 1;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
    set(axis_muscle, 'box', 'off');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["biceps femoris (long head),", "forward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;
    % backward flat
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{2, 4};    
    results_emg_mean      = sub04_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);    
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
    set(axis_muscle, 'box', 'off');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ["biceps femoris (long head),", "backward flat"];
    handle = title(axis_muscle, title_emg);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    drawnow;    

    %% reserve: a_forceset_reserve_jointset_ankle_r_ankle_angle_r
    % two left columns
    % forward flat
    file_trial            = file_trial_fw;
    axis                  = emg_moco_axes{1, 1};
    idx_reserve           = 2;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right ankle');    
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis, [0.0, 100.0]);
    ylim(axis, [-2.5,  +2.5]);
    set(axis, 'fontsize', size_font_ticks);
    xlabel(axis, 'percentage of gait cycle', 'FontSize', size_font_xy);
    ylabel(axis, 'Newton-meters', 'FontSize', size_font_xy);
%     title_reserve = [name_reserve_moco];
    title_reserve = ["reserve in the right ankle,", "forward flat"];
    handle = title(axis, title_reserve);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    legend(axis, 'Location', 'best', 'FontSize', size_font_legend);    
    drawnow;

    % backward trial
    file_trial            = file_trial_bw;
    axis                  = emg_moco_axes{1, 2};
    idx_reserve           = 2;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right ankle');
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis, [ 0.0, 100.0]);
    ylim(axis, [-2.5,  +2.5]);
    set(axis, 'fontsize', size_font_ticks);
    xlabel(axis, 'percentage of gait cycle', 'FontSize', size_font_xy);
%     ylabel(axis, 'activation');
%     title_reserve = [name_reserve_moco];
    title_reserve = ["reserve in the right ankle,", "backward flat"];
    handle = title(axis, title_reserve);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    legend(axis, 'Location', 'best', 'FontSize', size_font_legend);    
    drawnow;

    % two right columns
    % forward flat
    file_trial            = file_trial_fw;
    axis                  = emg_moco_axes{1, 3};
    idx_reserve           = 1;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right knee');    
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis, [0.0, 100.0]);
    ylim(axis, [-5.0,  +5.0]);
    set(axis, 'fontsize', size_font_ticks);
    xlabel(axis, 'percentage of gait cycle', 'FontSize', size_font_xy);
%     ylabel(axis, 'Newton-meters');
%     title_reserve = [name_reserve_moco];
    title_reserve = ["reserve in the right knee,", "forward flat"];
    handle = title(axis, title_reserve);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
    legend(axis, 'Location', 'best', 'FontSize', size_font_legend);    
    drawnow;

    % backward trial
    file_trial            = file_trial_bw;
    axis                  = emg_moco_axes{1, 4};
    idx_reserve           = 1;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub04_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub04_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right knee');
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis, [ 0.0, 100.0]);
    ylim(axis, [-2.5,  +2.5]);
    set(axis, 'fontsize', size_font_ticks);
    xlabel(axis, 'percentage of gait cycle', 'FontSize', size_font_xy);
%     ylabel(axis, 'activation');
%     title_reserve = [name_reserve_moco];
    title_reserve = ["reserve in the right knee,", "backward flat"];    
    handle = title(axis, title_reserve);
    set(handle, 'Units', 'normalized', 'Position', position_title);
    set(handle, 'FontSize', size_font_title, 'FontWeight', weight_font);
%     legend( ...
%         'EMG mean', 'EMG activation', ...
%         'Moco inverse mean', 'Moco inverse activation', ...
%         'reserve in right ankle')
    legend(axis, 'Location', 'best', 'FontSize', size_font_legend);
    drawnow;
end
    
    %% sub06

if false % plot_emg_moco_comparison

    % Create a figure with four rows and four columns of subplots sharing the same axes
    figure;
    
%     % Define subplot sizes
%     vertical_offset   = 0.04;
%     horizontalSpacing = 0.03;  % Horizontal spacing between subplots
%     verticalSpacing   = 0.03;  % Vertical spacing between subplots
%     subplotWidth      = 0.20;  % Width of each subplot
%     subplotHeight     = 0.125;  % Height of each subplot
    
    i = 1;
    for j = 1:4
        % Calculate position for each subplot
        position = [...
            horizontalSpacing + (j-1)*(subplotWidth+horizontalSpacing), ...
            vertical_offset   + (i-1)*(subplotHeight+verticalSpacing), ...
            subplotWidth, subplotHeight];
        % Create subplot
        emg_moco_axes{i, j} = axes('Position', position);

        % Remove frame around the subplot
        set(gca, 'box', 'off');
        % Hide unnecessary axis ticks and labels
        if j ~= 1
            set(gca, 'YTickLabel', []);
        end
        if i ~= 1
            set(gca, 'XTickLabel', []);
        end
    end
    
    % Define subplot sizes
    horizontal_offset = horizontalSpacing;
    vertical_offset   = verticalSpacing + subplotHeight + verticalSpacing;
    subplotWidth      = 0.2;  % Width of each subplot
    subplotHeight     = 0.3;  % Height of each subplot
    
    % Plot the subplots
    for i = 2:3
        for j = 1:4
            % Calculate position for each subplot
            position = [...
                horizontal_offset + (j-1)*(subplotWidth+horizontalSpacing), ...
                vertical_offset   + (i-2)*(subplotHeight+verticalSpacing), ...
                subplotWidth, subplotHeight];
            % Create subplot
            emg_moco_axes{i, j} = axes('Position', position);

            % Hide unnecessary axis ticks and labels
            if j ~= 1
                set(gca, 'YTickLabel', []);
            end
            set(gca, 'XTickLabel', []);

        end
    end

    file_trial_fw = 'sub06_strifw_condp0000_speedp0100_0001';
    file_trial_bw = 'sub06_stribw_condp0000_speedn0060_0001';

    %% gastrocnemius medial: RightGastrocnemius - a_forceset_gasmed_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{3, 1};    
    idx_muscle            = 2;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};    
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
    ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['gastrocnemius medial, forward flat'];
    title(axis_muscle, title_emg);
    drawnow;
    % backward trial
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{3, 2};    
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['gastrocnemius medial, backward flat'];
    title(axis_muscle, title_emg);    
    drawnow;

    %% rectus femoris: RightRectusFemoris - a_forceset_recfem_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{3, 3};    
    idx_muscle            = 3;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['rectus femoris, forward flat'];
    title(axis_muscle, title_emg);
    drawnow;
    % backward flat
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{3, 4};    
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;    
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill([1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);    
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['rectus femoris, backward flat'];
    title(axis_muscle, title_emg);
    legend(axis_muscle, ...
        'EMG mean', 'EMG activation', ...
        'Moco inverse mean', 'Moco inverse activation');    
    drawnow;        

    %% tibialis anterior: RightTibialisAnterior - a_forceset_tibant_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{2, 1};    
    idx_muscle            = 4;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};    
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
    ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['tibialis anterior, forward flat'];
    title(axis_muscle, title_emg);    
    drawnow;
    % backward trial
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{2, 2};
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);   
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];    
    title_emg = ['tibialis anterior, backward flat'];
    title(axis_muscle, title_emg);
    drawnow;

    %% biceps femoris: RightBicepsFemoris - a_forceset_bflh_r_activation
    % forward flat
    file_trial            = file_trial_fw;
    axis_muscle           = emg_moco_axes{2, 3};
    idx_muscle            = 1;
    name_muscle_emg       = name_emg_muscles{idx_muscle};
    name_muscle_moco      = name_moco_muscle_activation{idx_muscle};
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
    set(axis_muscle, 'box', 'off');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['biceps femoris (long head), forward flat'];
    title(axis_muscle, title_emg);
    drawnow;
    % backward flat
    file_trial            = file_trial_bw;
    axis_muscle           = emg_moco_axes{2, 4};    
    results_emg_mean      = sub06_results_emg.(file_trial).mean.(name_muscle_emg);
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_muscle_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    hold(axis_muscle, 'on');
    plot(axis_muscle, results_emg_mean, 'b-', 'LineWidth', 1.5);    
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_emg_mean, zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
    plot(axis_muscle, results_moco_mean, 'r-', 'LineWidth', 1.5);
    fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_muscle, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
    xlim(axis_muscle, [0.0, 100.0]);
    ylim(axis_muscle, [0.0,  +1.0]);
%     xlabel(axis_muscle, 'percentage of gait cycle');
%     ylabel(axis_muscle, 'activation');
    set(axis_muscle, 'box', 'off');
%     title_emg = [name_muscle_emg, ' ', name_muscle_moco];
    title_emg = ['biceps femoris (long head), backward flat'];
    title(axis_muscle, title_emg);
    drawnow;    

    %% reserve: a_forceset_reserve_jointset_ankle_r_ankle_angle_r
    % two left columns
    % forward flat
    file_trial            = file_trial_fw;
    axis_reserve          = emg_moco_axes{1, 1};
    idx_reserve           = 2;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis_reserve, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right ankle');    
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis_reserve, [0.0, 100.0]);
    ylim(axis_reserve, [-3.0,  +3.0]);
    xlabel(axis_reserve, 'percentage of gait cycle');
    ylabel(axis_reserve, 'Newton-meters');
%     title_reserve = [name_reserve_moco];
    title_reserve = ['reserve in the right ankle, forward flat'];
    title(axis_reserve, title_reserve);
    legend(axis_reserve, 'Location', 'best');    
    drawnow;

    % backward trial
    file_trial            = file_trial_bw;
    axis_reserve          = emg_moco_axes{1, 2};
    idx_reserve           = 2;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis_reserve, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right ankle');
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis_reserve, [ 0.0, 100.0]);
    ylim(axis_reserve, [-3.0,  +3.0]);
    xlabel(axis_reserve, 'percentage of gait cycle');
%     ylabel(axis_reserve, 'activation');
%     title_reserve = [name_reserve_moco];
    title_reserve = ['reserve in the right ankle, backward flat'];    
    title(axis_reserve, title_reserve);
    legend(axis_reserve, 'Location', 'best');    
    drawnow;

    % two right columns
    % forward flat
    file_trial            = file_trial_fw;
    axis_reserve          = emg_moco_axes{1, 3};
    idx_reserve           = 1;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis_reserve, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right knee');    
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis_reserve, [0.0, 100.0]);
    ylim(axis_reserve, [-10.0,  +10.0]);
    xlabel(axis_reserve, 'percentage of gait cycle');
%     ylabel(axis_reserve, 'Newton-meters');
%     title_reserve = [name_reserve_moco];
    title_reserve = ['reserve in the right knee, forward flat'];
    title(axis_reserve, title_reserve);
    legend(axis_reserve, 'Location', 'best');    
    drawnow;

    % backward trial
    file_trial            = file_trial_bw;
    axis_reserve          = emg_moco_axes{1, 4};
    idx_reserve           = 1;
    name_reserve_moco     = name_moco_reserves{idx_reserve};
    results_moco_mean     = sub06_results_moco.(file_trial).mean.(name_reserve_moco);
    results_moco_sd       = sub06_results_moco.(file_trial).sd.(name_muscle_moco);
    results_moco_sd_upper = results_moco_mean + results_moco_sd;
    results_moco_sd_lower = results_moco_mean - results_moco_sd;
    plot(axis_reserve, results_moco_mean, 'k-', 'LineWidth', 1.5, 'DisplayName', 'reserve in right knee');
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_mean, zeros(100, 1)'], 'r', 'FaceAlpha', 0.1);    
%     fill(axis_reserve, [1:100, fliplr(1:100)], [results_moco_sd_upper, results_moco_sd_lower], 'r', 'FaceAlpha', 0.1);
%     xlim(axis_reserve, [ 0.0, 100.0]);
    ylim(axis_reserve, [-3.0,  +3.0]);
    xlabel(axis_reserve, 'percentage of gait cycle');
%     ylabel(axis_reserve, 'activation');
%     title_reserve = [name_reserve_moco];
    title_reserve = ['reserve in the right knee, backward flat'];    
    title(axis_reserve, title_reserve);
%     legend( ...
%         'EMG mean', 'EMG activation', ...
%         'Moco inverse mean', 'Moco inverse activation', ...
%         'reserve in right ankle')
    legend(axis_reserve, 'Location', 'best');
    drawnow;

end
