function [results_emg, results_moco] = data_dynamic_process(path_subject, file_trial, data)

import org.opensim.modeling.*;
pause on;

%% Config the flags.

compute_inverse_kinematics       = false;
% compute_inverse_kinematics       = true;

compute_moco_inverse_statistics = false;
compute_moco_inverse_statistics = true;

plot_grf_raw_data            = false;
plot_marker_raw_data         = false;
plot_emg_raw_data            = false;
plot_grf_statistics          = false;
plot_emg_statistics          = false;
plot_inverse_kinematics      = false;
plot_moco_inverse            = false;
plot_moco_inverse_statistics = false;

% plot_grf_raw_data            = true;
% plot_emg_raw_data            = true;
% plot_marker_raw_data         = true;
% plot_grf_statistics          = true;
% plot_emg_statistics          = true;
% plot_inverse_kinematics      = true;
% plot_moco_inverse            = true;
% plot_moco_inverse_statistics = true;

mode_interactive             = false;
mode_interactive             = ...
    plot_grf_raw_data | plot_marker_raw_data | plot_emg_raw_data | ...
    plot_inverse_kinematics | plot_moco_inverse;

%% Load the *.c3d and *.csv files corresponding to the current trial.

% Set the paths of the directories containg the resources to be used.
[path_root, path_dataset, ...
    marker_list, calc_marker_list, ...
    path_model, file_model ...
    name_emg_muscles, ...
    name_moco_reserves, name_moco_muscle_activation, name_moco_muscle_forces] = data_pipeline_setup();

name_muscles           = name_emg_muscles;
name_reserves          = name_moco_reserves;
name_muscle_activation = name_moco_muscle_activation;
name_muscle_forces     = name_moco_muscle_forces;

cd([path_dataset, path_subject]);

% Load the trial's data.
data_c3d = btk_loadc3d([path_dataset, path_subject, file_trial, '.c3d'], 10);

% Copy the fields of c3d file into data received as parameter.
fields = fieldnames(data_c3d);
% Loop over each marker name
for idx_field = 1:length(fields)
    data.(fields{idx_field}) = data_c3d.(fields{idx_field});
end

emgTable = readtable([path_dataset, path_subject, file_trial, '.csv']);

%% Process the marker data.

% sort the C3D file so we know what is Marker data and what is calculated
% data -- THIS STEP ISN'T REALLY REQUIRED FOR THIS DATASET WHERE ALL
% THE MARKERS ARE REQUIRED
data.marker_data = btk_sortc3d(data.marker_data,marker_list,calc_marker_list);

% Filter marker data.
moving_average_window_size = 10;

% Get a cell array of marker names
markers = fieldnames(data.marker_data.Markers);

% Loop over each marker name
for i = 1:length(markers)
    % Get the current marker name
    marker = markers{i};
    
    % Copy the initial data into a new field.
    data.marker_data_initial.Markers.(marker) = data.marker_data.Markers.(marker);

    % Access the value of the current marker    
    marker_values = data.marker_data.Markers.(marker);
    
    % Apply movemean filter on marker's data.
    data.marker_data.Markers.(marker) = movmean(marker_values, moving_average_window_size);
end

% Upsample the marker data to the sampling rate of the force plate data.
factor_upsampling_marker = data.fp_data.Info(1).frequency / data.marker_data.Info.frequency;

% Get a cell array of field names
markers = fieldnames(data.marker_data.Markers);
markers = sort(markers);

% Loop over each marker name
for i = 1:length(markers)
    % Get the current marker name
    marker = markers{i};

    % Access the value of the current field
    marker_values = data.marker_data.Markers.(marker);

    % Apply interpolate each coordinate (x, y, z).
    x = interp(marker_values(:, 1), factor_upsampling_marker);
    y = interp(marker_values(:, 2), factor_upsampling_marker);
    z = interp(marker_values(:, 3), factor_upsampling_marker);

    data.marker_data.Markers.(marker) = [x, y, z];
end

data.marker_data.Info.frequency = data.marker_data.Info.frequency * factor_upsampling_marker;
data.marker_data.Info.NumFrames = length(data.marker_data.Markers.LANK);

% Update the timeline too.
data.marker_data.Time = interp(data.marker_data.Time, factor_upsampling_marker);

%% Process the GRF data.

% Filter GRF data. The design belongs to Mahaki, 2017
grf_SR           = data.fp_data.Info.frequency; %[Hz]
cutoffFrequency  = 15; % Cutoff frequency in Hz
order            =  4; % Filter order
normalizedCutoff = cutoffFrequency / (0.5 * grf_SR); % Normalize cutoff frequency
[b, a] = butter(order, normalizedCutoff, 'low');

for idx_forceplate = 1:length(data.fp_data.GRF_data)
    % Save the initial values.
    data.fp_data_initial.GRF_data(idx_forceplate).F = data.fp_data.GRF_data(idx_forceplate).F;

    % Rectify first
    grf_rect = abs(data.fp_data.GRF_data(idx_forceplate).F(:,3));

    % Apply zero-lag, low-pass Butterworth filter using filtfilt.
    grf_filt       = filtfilt(b, a, grf_rect);
    grf_moving_avg = movmean(grf_filt, 50);
    grf_smooth     = sgolayfilt(grf_moving_avg, 2, 101);

    % Save the smothed Z component of the GRF(idx_forceplate).
    data.fp_data.GRF_data(idx_forceplate).F(:,3) = grf_smooth;
end

%% Process the EMG data.

% muscle_names = ["Channel0" "RightBicepsFemoris" "RightTensorFasciaLatae" "RightTibialisAnterior" "RightGastrocnemius" "RightGluteusMedius" "RightRectusFemoris"];
emgChannel0Threshold = -2.0e+7;
emg_data = {};

% Filter out the segments for which the Channel0 is inactive:
% - the moving average is smaller than the threshold.
emg_moving_avg = movmean(emgTable.Channel0, 100);
idx_emg_active = emg_moving_avg <= emgChannel0Threshold;

% Extract the segment which is active, as signaled by Channel0.
emg_data.Channel0               = emgTable.Channel0(idx_emg_active);
emg_data.RightBicepsFemoris     = emgTable.RightBicepsFemoris(idx_emg_active);
emg_data.RightTensorFasciaLatae = emgTable.RightTensorFasciaLatae(idx_emg_active);
emg_data.RightTibialisAnterior  = emgTable.RightTibialisAnterior(idx_emg_active);
emg_data.RightGastrocnemius     = emgTable.RightGastrocnemius(idx_emg_active);
emg_data.RightGluteusMedius     = emgTable.RightGluteusMedius(idx_emg_active);
emg_data.RightRectusFemoris     = emgTable.RightRectusFemoris(idx_emg_active);

% Detrend, filter, rectify and normalize to the peak of the trial.
emg_SR           = 2000;  % Sampling rate of EMG data

% Loop over each muscle name
emg_data_muscles = fieldnames(emg_data);
for idx_muscle = 1:length(emg_data_muscles)
    name_muscle = emg_data_muscles{idx_muscle};
    trial_muscle = emg_data.(name_muscle);

    % Detrending
    trial_muscle_detrended = detrend(trial_muscle);

    % High-pass filter design (Butterworth, 8 Hz cutoff)
    [highPassB, highPassA] = butter(4, 8/(0.5*emg_SR), 'high');
    
    % Low-pass filter design for EMG (Butterworth, 500 Hz cutoff)
    [lowPassB, lowPassA] = butter(4, 500/(0.5*emg_SR), 'low');
    
    % Low-pass filter design for rectified EMG (Butterworth, 6 Hz cutoff)
    [low2b, low2a] = butter(4, 6 / (0.5*emg_SR), 'low');
    
    % Apply high-pass filter
    highPassEMG = filtfilt(highPassB, highPassA, trial_muscle_detrended);
    
    % Apply low-pass filter
    lowPassEMG = filtfilt(lowPassB, lowPassA, highPassEMG);
    
    % Rectify the signal
    rect_EMG = abs(lowPassEMG);
    
    % Apply low-pass filter to the rectified signal
    emg_filt_rect.(name_muscle) = filtfilt(low2b, low2a, rect_EMG);

    % Normalize relative to the peak of the trial.
    max_peak = max(abs(emg_filt_rect.(name_muscle)));
    emg_filt_rect_norm.(name_muscle) = emg_filt_rect.(name_muscle) / max_peak;

end

% Downsample the EMG data (2000Hz)to the sampling rate of the force plate
% data (1000Hz).
factor_downsampling = 2;
% Loop over each muscle name
emg_data_muscles = fieldnames(emg_data);
for idx_muscle = 1:length(emg_data_muscles)
    name_muscle = emg_data_muscles{idx_muscle};
    emg_filt_rect_norm.(name_muscle) = ...
        decimate(emg_filt_rect_norm.(name_muscle), factor_downsampling);
end

%% Extract the gait cycles based on force plates data.

% Left foot stept on the forceplate 1
% Right food stept on the forceplate 2
% Rigth leg had EMG sensors attached on it.

% find when feet are on each forceplate and define the trial events as
% the first foot contact and the final foot contact
% -> when GRF's vertical component exceeds 5% of its maximum.
% -> when the GRF's vertical component exceeds 20Newtons.
%    - according to Adam this is the right criteria, not 5%.
%    - 20 Newtons does work properly on Adam's trials sub04.
threshold_grf = 40; % 5/100.0 * max(data.fp_data.GRF_data(i).F(:,3));

% Extract the indexes which store force plate values higher than a
% threshold.
contact_idxs_grf_1 = find(data.fp_data.GRF_data(1).F(:,3) > threshold_grf);
contact_idxs_grf_2 = find(data.fp_data.GRF_data(2).F(:,3) > threshold_grf);

% Extract the indexes at which the strike starts.
lifts_grf_1 = find([diff(contact_idxs_grf_1.')~=1, 1]);
lifts_grf_2 = find([diff(contact_idxs_grf_2.')~=1, 1]);

% Extract the begining and the end of gait cycles: strike, strike, lift,
% lift.
gait_cycles_grf_1 = cell(1, numel(lifts_grf_1));
gait_cycles_grf_2 = cell(1, numel(lifts_grf_2));

startIdx = 1;
for i = 1:numel(lifts_grf_1)
    gait_cycles_grf_1{i} = contact_idxs_grf_1(startIdx:lifts_grf_1(i));
    startIdx = lifts_grf_1(i) + 1;
end

startIdx = 1;
for i = 1:numel(lifts_grf_2)
    gait_cycles_grf_2{i} = contact_idxs_grf_2(startIdx:lifts_grf_2(i));
    startIdx = lifts_grf_2(i) + 1;
end

%% Normalize the GRF data relative to the gait cycles extracted.

grfDataNormalized = cell(length(gait_cycles_grf_2), length(data.fp_data.GRF_data));
% Iterate the gait cycles.
idx_1st_gait_cycle = 2;
for idx_gait_cycle = idx_1st_gait_cycle:length(gait_cycles_grf_2)-1
    
    cycle_start = gait_cycles_grf_2{idx_gait_cycle}(1);
    cycle_end   = gait_cycles_grf_2{idx_gait_cycle + 1}(1);

    for idx_grf = 1:length(data.fp_data.GRF_data)
        data_grf = data.fp_data.GRF_data(idx_grf).F(:,3);
        data_grf = data_grf(cycle_start:cycle_end);

        % Create a linspace vector for the current gait cycle based on the
        % number of data points between strikes.
        length_data_grf = linspace(0, 100, length(data_grf));

        % Spline interpolation to normalize the data length to 100 points.
        data_grf_normamlized = spline(length_data_grf, data_grf, 0:99);

        grfDataNormalized{idx_gait_cycle, idx_grf} = data_grf_normamlized;
    end
end

% Compute the mean and standard deviation for the normalized GRF data.

grfDataMean    = cell(length(data.fp_data.GRF_data));
grfDataSD      = cell(length(data.fp_data.GRF_data));
grfDataLowerSD = cell(length(data.fp_data.GRF_data));
grfDataUpperSD = cell(length(data.fp_data.GRF_data));

for idx_force_plate = 1:length(data.fp_data.GRF_data)
    grfDataMean{idx_force_plate}    = mean(cell2mat(grfDataNormalized(:, idx_force_plate)), 1);
    grfDataSD{idx_force_plate}      =  std(cell2mat(grfDataNormalized(:, idx_force_plate)), 1);
    grfDataLowerSD{idx_force_plate} = grfDataMean{idx_force_plate} - grfDataSD{idx_force_plate};
    grfDataUpperSD{idx_force_plate} = grfDataMean{idx_force_plate} + grfDataSD{idx_force_plate};
end

% Plot the statistics on GRF data.
if plot_grf_statistics
    figure_grf_statistics = figure;
    set(figure_grf_statistics, 'MenuBar', 'none');
    set(figure_grf_statistics, 'ToolBar', 'none');

    hold on;
    title_grf = file_trial;
    title_grf = strrep(title_grf, '_', '\_');
    title_grf = [title_grf ' GRFs Mean and SD'];
    title(title_grf);
    plot(grfDataMean{1}, 'r-', 'LineWidth', 1.5, ...
        'DisplayName', ['GRF 1']);
    fill([1:100, fliplr(1:100)], [grfDataLowerSD{1}, fliplr(grfDataUpperSD{1})], 'g', 'FaceAlpha', 0.1);
    plot(grfDataLowerSD{1}, 'g--', 'LineWidth', 1);
    plot(grfDataUpperSD{1}, 'g--', 'LineWidth', 1);
    
    plot(grfDataMean{2}, 'r-', 'LineWidth', 1.5, ...
        'DisplayName', ['GRF 2']);
    fill([1:100, fliplr(1:100)], [grfDataLowerSD{2}, fliplr(grfDataUpperSD{2})], 'b', 'FaceAlpha', 0.1);
    plot(grfDataLowerSD{2}, 'b--', 'LineWidth', 1);
    plot(grfDataUpperSD{2}, 'b--', 'LineWidth', 1);

    xlabel('Percentage of gait cycle');
    % ylabel('Amplitude (V)');
    legend();
    drawnow;
end

%% Normalize the EMG data relative to the gait cycles extracted.

% emg_data_normalized_cycle = cell(length(gait_cycles_grf_2), length(fieldnames(emg_data)));
% Iterate the gait cycles.
idx_1st_gait_cycle = 2;
for idx_gait_cycle = idx_1st_gait_cycle:length(gait_cycles_grf_2)-1
    
    cycle_start = gait_cycles_grf_2{idx_gait_cycle}(1);
    cycle_end   = gait_cycles_grf_2{idx_gait_cycle + 1}(1);

    % Update the EMG figures.
    % Loop over each muscle name
    emg_data_muscles = fieldnames(emg_filt_rect_norm);
    for idx_muscle = 1:length(emg_data_muscles)
        name_muscle = emg_data_muscles{idx_muscle};
        emg_muscle  = emg_filt_rect_norm.(name_muscle)(cycle_start:cycle_end);
    
        % Create a linspace vector for the current gait cycle based on the
        % number of data points between strikes.
        emg_muscle_time = linspace(0, 99, length(emg_muscle));

        % Spline interpolation to normalize the data length to 100 points.
        emg_muscle_norm_cycle = spline(emg_muscle_time, emg_muscle, 0:99);

        emg_filt_rect_norm_norm.(name_muscle){idx_gait_cycle} = emg_muscle_norm_cycle;
    end
end

% Calculate mean and standard deviation for the normalized EMG data
% Loop over each muscle name listed in name_muscles.
emg_data_muscles = fieldnames(emg_data);
for idx_muscle = 1:length(emg_data_muscles)
    name_muscle = emg_data_muscles{idx_muscle};

    emg_data_mean.(name_muscle) = mean(cell2mat(emg_filt_rect_norm_norm.(name_muscle)'), 1);
    emg_data_sd.(name_muscle)   =  std(cell2mat(emg_filt_rect_norm_norm.(name_muscle)'), 1);

    emg_data_lower_sd.(name_muscle) = emg_data_mean.(name_muscle) - emg_data_sd.(name_muscle);
    emg_data_upper_sd.(name_muscle) = emg_data_mean.(name_muscle) + emg_data_sd.(name_muscle);

    % Update the value returned by this function.
    results_emg.mean.(name_muscle) = emg_data_mean.(name_muscle);
    results_emg.sd.(name_muscle) = emg_data_sd.(name_muscle);
end

% Plot the statistics on the EMG data.
if plot_emg_statistics
    % EMGs figures: one figure per muscle;
    figure_emgs_statistics = {};
    
    for idx_muscle=1:length(name_muscles)
        figure_emgs_statistics{idx_muscle} = figure;
        set(figure_emgs_statistics{idx_muscle}, 'MenuBar', 'none');
        set(figure_emgs_statistics{idx_muscle}, 'ToolBar', 'none');
    end

    % Loop over each muscle name listed in name_muscles.    
    for idx_muscle = 1:length(name_muscles)
        name_muscle = name_muscles{idx_muscle};
        figure(figure_emgs_statistics{idx_muscle});
        hold on;
        title_emg = [name_muscle];
        title(title_emg);
        fill([1:100, fliplr(1:100)], ...
            [emg_data_mean.(name_muscle), zeros(100, 1)'], 'b', 'FaceAlpha', 0.1);
%         xlim([0, 100]);
        ylim([0.0, +1.0]);
        xlabel('Percentage of gait cycle');
        ylabel('Amplitude relative to the peak');
        drawnow;
    end
end

%% Generate the files for OpenSim

if plot_grf_raw_data
    % GRFs figure.
    figure_grf_raw = figure;
    set(figure_grf_raw, 'MenuBar', 'none');
    set(figure_grf_raw, 'ToolBar', 'none');
end

if plot_marker_raw_data
    % Vicon marker figures: one figure per marker.
    figure_markers_raw = {};
    for idx_marker = 1:length(fieldnames(data.marker_data.Markers))
        figure_markers_raw{idx_marker} = figure;
        set(figure_markers_raw{idx_marker}, 'MenuBar', 'none');
        set(figure_markers_raw{idx_marker}, 'ToolBar', 'none');
    end
end

if plot_emg_raw_data
    % EMGs figures: one figure per muscle.
    figure_emgs_raw = {};
    for idx_muscle=1:length(name_muscles)
        figure_emgs_raw{idx_muscle} = figure;
        set(figure_emgs_raw{idx_muscle}, 'MenuBar', 'none');
        set(figure_emgs_raw{idx_muscle}, 'ToolBar', 'none');
    end
end

if plot_inverse_kinematics
    % IKs figures: one figure per coordinate.
    figure_iks = {};
    size_coordinates = 19;
    for idx_coordinate = 1:size_coordinates
        figure_iks{idx_coordinate} = figure;
        set(figure_iks{idx_coordinate}, 'MenuBar', 'none');
        set(figure_iks{idx_coordinate}, 'ToolBar', 'none');
    end
end

% Iterate the gait cycles.
idx_1st_gait_cycle = 2;
for idx_gait_cycle = idx_1st_gait_cycle:length(gait_cycles_grf_2)-1

    % Define start and end frame from the events to write the appropriate TRC
    % and MOT files for the OpenSim simulations
    data.Start_Frame = gait_cycles_grf_2{idx_gait_cycle}(1);
    data.End_Frame   = gait_cycles_grf_2{idx_gait_cycle + 1}(1);
    
    % Use the same color for data describing the right side.
    color_right_side = rand(1, 3);

    if plot_grf_raw_data
        % Update the GRFs figure.
        figure(figure_grf_raw);
        hold on;

        plot(data.fp_data_initial.GRF_data(1).F(data.Start_Frame:data.End_Frame, 3), ...
            'DisplayName', 'GRF1 - left initial');        
        plot(data.fp_data.GRF_data(1).F(data.Start_Frame:data.End_Frame, 3), ...
            'DisplayName', 'GRF1 - left filtered');

        plot(data.fp_data_initial.GRF_data(2).F(data.Start_Frame:data.End_Frame, 3), ...
            'DisplayName', 'GRF2 - right initial');
        plot(data.fp_data.GRF_data(2).F(data.Start_Frame:data.End_Frame, 3), ...
            'DisplayName', 'GRF2 - right filtered', 'Color', color_right_side);

        title_grf = [file_trial ' ' num2str(idx_gait_cycle)];
        title_grf = strrep(title_grf, '_', '\_');        
        title(title_grf);
        legend();
        drawnow;
    end

    % Plot marker data.
    if plot_marker_raw_data
        
        % Loop over each marker name
        for idx_marker = 1:length(markers)
            % Get the current marker name
            marker = markers{idx_marker};
        
            % Access the value of the current field
            marker_values = data.marker_data_initial.Markers.(marker);
    
            figure(figure_markers_raw{idx_marker});
            hold on;
            title_marker = [file_trial ' ' num2str(idx_gait_cycle) ' ' marker];
            title_marker = strrep(title_marker, '_', '\_');
            title(title_marker);
            frame_first = data.Start_Frame / factor_upsampling_marker;
            frame_last  = data.End_Frame / factor_upsampling_marker;
            legend();            
            plot(marker_values(frame_first:frame_last, 1), 'DisplayName', [marker ' x']);
            plot(marker_values(frame_first:frame_last, 2), 'DisplayName', [marker ' y']);
            plot(marker_values(frame_first:frame_last, 3), 'DisplayName', [marker ' z']);
            drawnow;
        end
    end

    if plot_emg_raw_data
        % Update the EMG figures.

        % Loop over each muscle name listed in name_muscles.
        for idx_muscle = 1:length(name_muscles)
            name_muscle = name_muscles{idx_muscle};
        
            figure(figure_emgs_raw{idx_muscle});
            hold on;
%             % Raw data
%             plot(emg_data.(name_muscle)(data.Start_Frame:data.End_Frame), ...
%                 'DisplayName', 'raw');
%             % Filtered data
%             plot(emg_filt_rect.(name_muscle)(data.Start_Frame:data.End_Frame), ...
%                 'DisplayName', 'filtered');
            % Normalized data
            plot(emg_filt_rect_norm.(name_muscle)(data.Start_Frame:data.End_Frame), ...
                'DisplayName', 'normalized');

%             ylim([-1.5, +1.5]);
            title_emg = ['cycle: ' num2str(idx_gait_cycle) ' ' name_muscle];
            title_emg = strrep(title_emg, '_', '\_');
            title(title_emg);
            xlabel('Percentage of gait cycle');
    %         ylabel('Amplitude (V)');            
            legend();
            drawnow;
        end
    end

    if compute_inverse_kinematics | plot_inverse_kinematics

        % now do conversions to TRC files using btk_c3d2trc
        TRCFileName = [file_trial '_' sprintf('%04d', idx_gait_cycle) '_trc.trc'];
        GRFFileName = [file_trial '_' sprintf('%04d', idx_gait_cycle) '_grf.sto'];
        GRFFileXML  = [file_trial '_' sprintf('%04d', idx_gait_cycle) '_grf.xml'];
        
        data_trc = btk_c3d2trc(data, 'off', TRCFileName, GRFFileName);
        
        % define the standard OpenSim files for IK
        ModelFile    = [data_trc.Name '_SCALED.osim'];
        IKTasksFile  = [path_model file_model '_IK_Tasks.xml'];
        IKOutputFile = [file_trial '_' sprintf('%04d', idx_gait_cycle) '_ik.sto'];
        
        setup_InverseKinematics( ...
            'data',data_trc, ...
            'ModelFile',ModelFile,...
            'IKTasksFile',IKTasksFile, ...
            'Accuracy',0.00005, ...
            'OutputFile',IKOutputFile);
        
        % run the ik tool from the command line
        com = ['opensim-cmd run-tool ' data_trc.Name '_Setup_InverseKinematics.xml'];
        system(com);
    
        % load the data from MOT file and tie kinematic data to data structure
        D = load_sto_file(IKOutputFile);
        fnames = fieldnames(D);
        
        clear a
        a = [strfind(lower(fnames),'tx') ...
             strfind(lower(fnames),'ty') ...
             strfind(lower(fnames),'tz') ...
             strfind(lower(fnames),'pelvis')];
        
        for i = 1:size(a,1)
            if isempty([a{i,1:3}])
                data_trc.Kinematics.Angles.(fnames{i}) = D.(fnames{i});
            elseif ~isempty([a{i,4}])
                data_trc.Kinematics.Angles.(fnames{i}) = D.(fnames{i});
            else
                data_trc.Kinematics.Markers.(fnames{i}) = D.(fnames{i});
            end
        end
            
        % assign forces to the specific segments using 'assign_force' function
        % the output creates a cell array with the name of the external force
        % (ExForce) and the body that is assigned to that force (ApBodies)
        data_trc = assign_forces(data_trc,{'RHEE','LHEE'},{'calcn_r','calcn_l'},0.25);
        
        % define the standard OpenSim files for ID   
        IDOutputFile = [file_trial '_' sprintf('%04d', idx_gait_cycle) '_id.sto'];
        
        % it is also necessary to generate an XML file containing the information
        % about which column is which and the bodies these forces are applied to
        % this is done with the grf2xml function
        grf2xml( ...
            data_trc, ...
            'ExternalLoadNames',data_trc.AssignForce.ExForce, ...
            'AppliedToBodies',data_trc.AssignForce.ApBodies,...
            'GRFFile',GRFFileName, ...
            'MOTFile',IKOutputFile, ...
            'LowPassFilterForKinematics',15,...
            'OutputFile',GRFFileXML);
    
        % make setup file for inverse dynamics
        setup_InverseDynamics(...
            'data',data_trc, ...
            'ModelFile',ModelFile, ...
            'MOTFile',IKOutputFile,...
            'GRFFile',GRFFileXML, ...
            'LowPassFilterForKinematics',15, ...
            'OutputFile',IDOutputFile); 
        
        % call ID tool from the command line using setup file
        com = ['opensim-cmd run-tool ' data_trc.Name '_Setup_InverseDynamics.xml'];
        system(com);
        
    %     % todo(lisca): why do we need this?
    %     clear D   
    %     % load the data from ID MOT file and tie kinetic data to data structure
    %     D = load_sto_file(IDOutputFile);
    %     data_trc.Kinetics = D;
    %     
    %     save([data_trc.TRC_Filename '.mat']');
    %     
    %     clear D
        
        %% Filter the IK computed by OpenSim
        disp("Filtering the IK which OpenSim computed ...");
        
        ik_data_table               = TimeSeriesTable(IKOutputFile);
        ik_data_table_metadata_keys = ik_data_table.getTableMetaDataKeys();
        
        ik_data       = osimTableToStruct(ik_data_table);
        
        grf_SR          = 1000; %[Hz]
        cutoffFrequency = 15; % Cutoff frequency in Hz
        order            = 4; % Filter order
        normalizedCutoff = cutoffFrequency / (0.5 * grf_SR); % Normalize cutoff frequency
        [b, a] = butter(order, normalizedCutoff, 'low');
        
        % Get a cell array of coordinate names.
        coordinates = fieldnames(ik_data);
        
        % Loop over each coordinate names.
        for idx_coordinate = 1:length(coordinates)

            % Get the current coordinate name.
            coordinate_name = coordinates{idx_coordinate};
    
            % Keep the time line as it is!
            if strcmp(coordinate_name, 'time') == 1
                continue;
            end
        
            % Access the value of the current coordinate
            coordinate_values = ik_data.(coordinate_name);
       
            % Apply zero-lag, low-pass Butterworth filter using filtfilt.
            filt            = filtfilt(b, a, coordinate_values);
            moving_avg      = movmean(filt, 50);
            smoothed_values = sgolayfilt(moving_avg, 2, 101);
            
            % todo(lisca): move the plot out of the filter.
            if plot_inverse_kinematics
                figure(figure_iks{idx_coordinate});
                title_ik = coordinate_name;
                title_ik = strrep(title_ik, '_', '\_');
                title(title_ik);
                hold on;
                plot(ik_data.(coordinate_name), 'DisplayName', 'rough');
                plot(smoothed_values, 'DisplayName', 'smoothed');
                legend();
                drawnow;
            end
    
            % Replace the coordinates with their smoothed values;
            ik_data.(coordinate_name) = smoothed_values;
        end
        
        ik_data_table_smoothed = osimTableFromStruct(ik_data);
        for i = 0:length(ik_data_table_metadata_keys)
            key   = ik_data_table_metadata_keys.get(i);
            value = ik_data_table.getTableMetaDataString(key);
            ik_data_table_smoothed.addTableMetaDataString(key, value);
        end
        
        disp(['Saving the filtered IK replacing ' IKOutputFile ' file ...']);
        STOFileAdapter.write(ik_data_table_smoothed, IKOutputFile);

    end

    % If some plot is active then wait to press any key before cleaning the
    % figures.
    if mode_interactive
        disp('Next gait cycle? ... just click or press the space key ...');
        w = waitforbuttonpress;
%         pause(0.5);
    end

    % Clean all open figures.
    if plot_grf_raw_data
        figure(figure_grf_raw);
        clf;
    end

    if plot_marker_raw_data
        for idx_marker = 1:length(figure_markers_raw)
            figure(figure_markers_raw{idx_marker});
            clf;
        end
    end

    if plot_emg_raw_data
        for idx_muscle=1:length(name_muscles)
            figure(figure_emgs_raw{idx_muscle});
            clf;
        end
    end

    if plot_inverse_kinematics
        for idx_figure =1:length(figure_iks)
            figure(figure_iks{idx_figure});
            clf;
        end
    end
end

%% Process the results of OpenSim Moco

if plot_moco_inverse
    figure_reserves = {};
    for idx_figure = 1:length(name_reserves)
        figure_reserves{idx_figure} = figure;
        set(figure_reserves{idx_figure}, 'MenuBar', 'none');
        set(figure_reserves{idx_figure}, 'ToolBar', 'none');
    end
    
    figure_muscle_activation = {};
    for idx_figure = 1:length(name_muscle_activation)
        figure_muscle_activation{idx_figure} = figure;
        set(figure_muscle_activation{idx_figure}, 'MenuBar', 'none');
        set(figure_muscle_activation{idx_figure}, 'ToolBar', 'none');
    end
    
    figure_muscle_forces = {};
    for idx_figure = 1:length(name_muscle_forces)
        figure_muscle_forces{idx_figure} = figure;
        set(figure_muscle_forces{idx_figure}, 'MenuBar', 'none');
        set(figure_muscle_forces{idx_figure}, 'ToolBar', 'none');
    end
end

if plot_moco_inverse_statistics
    figure_reserves_statistics = {};
    for idx_figure = 1:length(name_reserves)
        figure_reserves_statistics{idx_figure} = figure;
        set(figure_reserves_statistics{idx_figure}, 'MenuBar', 'none');
        set(figure_reserves_statistics{idx_figure}, 'ToolBar', 'none');
    end
    
    figure_muscle_activation_statistics = {};
    for idx_figure = 1:length(name_muscle_activation)
        figure_muscle_activation_statistics{idx_figure} = figure;
        set(figure_muscle_activation_statistics{idx_figure}, 'MenuBar', 'none');
        set(figure_muscle_activation_statistics{idx_figure}, 'ToolBar', 'none');
    end
    
    figure_muscle_forces_statistics = {};
    for idx_figure = 1:length(name_muscle_forces)
        figure_muscle_forces_statistics{idx_figure} = figure;
        set(figure_muscle_forces_statistics{idx_figure}, 'MenuBar', 'none');
        set(figure_muscle_forces_statistics{idx_figure}, 'ToolBar', 'none');
    end
end

if compute_moco_inverse_statistics | ...
    plot_moco_inverse_statistics | ...
    plot_moco_inverse

    % Collect all OpenSim Moco Inverse solutions corresponding to the
    % current file_trial parameter.
    files_mi = dir( ...
        fullfile([path_dataset path_subject], [file_trial '*_mi.sto']));
    
    % Iterate over all OpenSim Moco Inverse solutions corresponding to the
    % current file_trial parameter.
    for idx_file = 1:length(files_mi)

        % Extract the index of the gait cycle from the name the OpenSim
        % Moco *_mi.sto file.
        idx_gait_cycle = str2num(files_mi(idx_file).name(40:43));

        % Compose the full path of the OpenSim Moco solution.
        path_mi = [path_dataset path_subject files_mi(idx_file).name];
    
        % Skipp the *_mi.sto file if it is empy.
        info_path_mi = dir(path_mi);
        if info_path_mi.bytes == 0
            disp([files_mi(idx_file).name ' file is empty!']);
            continue;
        end
    
        mi_data_table  = TimeSeriesTable(path_mi);
        mi_data        = osimTableToStruct(mi_data_table);
        mi_data_fields = fieldnames(mi_data);
    
        for idx_filed = 1:length(mi_data_fields)
            name_field = mi_data_fields{idx_filed};
            mi_data_trial.(name_field){idx_file} = mi_data.(name_field);
    
            % Create a linspace vector for the current gait cycle based on the
            % number of data points between strikes.
            % todo(lisca): update this description in the context of OpenSim Moco
            length_mi = linspace(0, 100, length(mi_data.(name_field)));
    
            % Spline interpolation to normalize the data length to 100 points.
            mi_data_normalized.(name_field){idx_file} = spline(length_mi, mi_data.(name_field), 0:99);

        end
    end

    if plot_moco_inverse
        % Plot muscle activation
        mi_data_trial_fields = fieldnames(mi_data_trial);
        for idx_file = 1:length(files_mi)
            for idx_field = 1:length(mi_data_trial_fields)
                data_trial_field = mi_data_trial_fields{idx_field};
                for idx_muscle = 1:length(name_muscle_activation)
                    name_muscle = name_muscle_activation(idx_muscle);
                    if contains(data_trial_field, name_muscle)
                        figure(figure_muscle_activation{idx_muscle});
                        hold on;
                        plot(mi_data_trial.(data_trial_field){idx_file});
                        title_muscle_activation = [file_trial ' ' num2str(idx_file) ' ' name_muscle];
                        title_muscle_activation = strrep(title_muscle_activation, '_', '\_');
                        title(title_muscle_activation);
    %                     xlim([0, 100]);
%                         ylim([0.0, +1.0]);
    %                     xlabel('Percentage of gait cycle');
    %                     ylabel('Amplitude (V)');
                        drawnow;                    
                    end
                end
            end

            if mode_interactive
                disp('Next OpenSim Moco solution? Click or press any key.');
                w = waitforbuttonpress;
%                 pause(0.5);
            end

%             for idx_figure = 1:length(name_reserves)
%                 figure(figure_reserves{idx_figure});
%                 clf;
%             end
            
            for idx_figure = 1:length(name_muscle_activation)
                figure(figure_muscle_activation{idx_figure});
                clf;
            end
    
%             for idx_figure = 1:length(name_muscle_forces)
%                 figure(figure_muscle_forces{idx_figure});
%                 clf;
%             end

        end
    end

    % Compute the mean and standard devitation.
    for idx_filed = 1:length(mi_data_fields)
        name_field = mi_data_fields{idx_filed};
        % Compute the mean and standard deviation for the normalized Moco Inverse solutions.
        % Iterate the results of Moco Inverse solution.
        mi_mean.(name_field)     = mean(cell2mat(mi_data_normalized.(name_field)'), 1);
        mi_SD.(name_field)       =  std(cell2mat(mi_data_normalized.(name_field)'), 1);
        mi_lower_SD.(name_field) = mi_mean.(name_field) - mi_SD.(name_field);
        mi_upper_SD.(name_field) = mi_mean.(name_field) + mi_SD.(name_field);

        % Update the value returned by this function.
        results_moco.mean.(name_field) = mi_mean.(name_field);
        results_moco.sd.(name_field) = mi_SD.(name_field);
    end

    if plot_moco_inverse_statistics
    
        for idx_filed = 1:length(mi_data_fields)
            name_field = mi_data_fields{idx_filed};
    
            % Plot reservers
            for idx_reserve = 1:length(name_reserves)
                if contains(name_field, name_reserves{idx_reserve})
                    figure(figure_reserves_statistics{idx_reserve});
                    hold on;
                    plot(mi_mean.(name_field), 'r-', 'LineWidth', 1.5);
                    fill([1:100, fliplr(1:100)], ...
                        [mi_lower_SD.(name_field), fliplr(mi_upper_SD.(name_field))], 'b', 'FaceAlpha', 0.1);
                    plot(mi_lower_SD.(name_field), 'b--', 'LineWidth', 1);
                    plot(mi_upper_SD.(name_field), 'b--', 'LineWidth', 1);
                    title_reserve = [file_trial string(name_reserves(idx_reserve)) ' Mean and SD'];
                    title_reserve = strrep(title_reserve, '_', '\_');
                    title(title_reserve);                
%                     xlim([0, 100]);
%                     ylim([0.0, +1.0]);
                    xlabel('Percentage of gait cycle');
%                     ylabel('Amplitude (V)');
                    drawnow;
                end
            end
    
            % Plot muscle activation
            for idx_muscle = 1:length(name_muscle_activation)
                if contains(name_field, name_muscle_activation{idx_muscle})
                    figure(figure_muscle_activation_statistics{idx_muscle});
                    hold on;       
                    plot(mi_mean.(name_field), 'r-', 'LineWidth', 1.5);
%                     fill([1:100, fliplr(1:100)], ...
%                         [mi_lower_SD.(name_field), fliplr(mi_upper_SD.(name_field))], 'b', 'FaceAlpha', 0.1);
%                     plot(mi_lower_SD.(name_field), 'b--', 'LineWidth', 1);
%                     plot(mi_upper_SD.(name_field), 'b--', 'LineWidth', 1);
                    title_muscle_activation = [file_trial string(name_muscle_activation(idx_muscle)) ' Mean and SD'];
                    title_muscle_activation = strrep(title_muscle_activation, '_', '\_');
                    title(title_muscle_activation);
%                     xlim([0, 100]);
%                     ylim([0.0, +1.0]);
                    xlabel('Percentage of gait cycle');
%                     ylabel('Amplitude (V)');
                    drawnow;
                end
            end
    
            % Plot muscle forces
            for idx_muscle = 1:length(name_muscle_forces)
                if ~contains(name_field, "activation") && ...
                        contains(name_field, name_muscle_forces{idx_muscle})
    
                    figure(figure_muscle_forces_statistics{idx_muscle});
                    hold on;
                    plot(mi_mean.(name_field), 'r-', 'LineWidth', 1.5);
                    fill([1:100, fliplr(1:100)], ...
                        [mi_lower_SD.(name_field), fliplr(mi_upper_SD.(name_field))], 'b', 'FaceAlpha', 0.1);
                    plot(mi_lower_SD.(name_field), 'b--', 'LineWidth', 1);
                    plot(mi_upper_SD.(name_field), 'b--', 'LineWidth', 1);
                    title_muscle_force = [file_trial string(name_muscle_forces(idx_muscle)) ' Mean and SD'];
                    title_muscle_force = strrep(title_muscle_force, '_', '\_');
                    title(title_muscle_force);
%                     xlim([0, 100]);
%                     ylim([0.0, +1.0]);
                    xlabel('Percentage of gait cycle');
%                     ylabel('Amplitude (V)');
                    drawnow;                
                end
            end
        end
    end
end

% disp('Close all figures and process the next trial? Click or press any key.');
% w = waitforbuttonpress;

% Close all opene figures.
if plot_grf_raw_data
    close(figure_grf_raw);
end

if plot_marker_raw_data
    for idx_marker = 1:length(figure_markers_raw)
        close(figure_markers_raw{idx_marker});
    end
end

if plot_emg_raw_data
    for idx_muscle=1:length(figure_emgs_raw)
        close(figure_emgs_raw{idx_muscle});
    end
end

if plot_grf_statistics
    close(figure_grf_statistics);
end

if plot_emg_statistics
    for idx_muscle=1:length(figure_emgs_statistics)
        close(figure_emgs_statistics{idx_muscle});
    end
end

if plot_inverse_kinematics
    for idx_figure =1:length(figure_iks)
        close(figure_iks{idx_figure});
    end
end

if plot_moco_inverse
    for idx_figure = 1:length(figure_reserves)
        close(figure_reserves{idx_figure});
    end

    for idx_figure = 1:length(figure_muscle_activation)
        close(figure_muscle_activation{idx_figure});
    end

    for idx_figure = 1:length(figure_muscle_forces)
        close(figure_muscle_forces{idx_figure});
    end
end

if plot_moco_inverse_statistics
    for idx_figure = 1:length(figure_reserves_statistics)
        close(figure_reserves_statistics{idx_figure});
    end

    for idx_figure = 1:length(figure_muscle_activation_statistics)
        close(figure_muscle_activation_statistics{idx_figure});
    end

    for idx_figure = 1:length(figure_muscle_forces_statistics)
        close(figure_muscle_forces_statistics{idx_figure});
    end
end

