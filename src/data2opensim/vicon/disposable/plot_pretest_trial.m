function plot_pretest_trial(path_trial)
    % Ensure that all configuration variables are properly loaded.
    config = load_config();

    path_trial = '/home/lisca/ssd/nvme0n1p9/lisca/datasets/20231212_dataset_lutetium_pretest/Walking2.c3d';
    acq_trial = btkReadAcquisition(path_trial);
    [analogs, analogs_info] = btkGetAnalogs(acq_trial);

    name_muscles = { ...
        'L_GastrLat', 'R_GastrLat', 'L_GMed', 'R_BF',   'L_GastrMed', 'L_TFL',  'L_RF',       'L_VasMed', ...
        'L_BF',       'R_RF',       'R_TibA', 'L_TibA', 'R_TFL',      'R_GMed', 'R_GastrMed', 'R_VasMed'};

    frequency_emg = analogs_info.frequency;

    figure_emg = figure;

    for idx_channel = 1:16
        emg.(name_muscles{idx_channel}) = analogs.(['Emg_', char(string(idx_channel))]);

        % Detrending
        data_muscle_detrended = detrend(emg.(name_muscles{idx_channel}));

        % High-pass filter design (Butterworth, 8 Hz cutoff)
        [high_pass_b, high_pass_a] = butter(4,   8/(0.5*frequency_emg), 'high');
        % Low-pass filter design for EMG (Butterworth, 500 Hz cutoff)
        [low_pass_b, low_pass_a]   = butter(4, 500/(0.5*frequency_emg), 'low');
        % Low-pass filter design for rectified EMG (Butterworth, 6 Hz cutoff)
        [low_2b, low_2a]           = butter(4,   6/(0.5*frequency_emg), 'low');

        % Apply high-pass filter
        data_muscle_high_pass = filtfilt(high_pass_b, high_pass_a, data_muscle_detrended);
        % Apply low-pass filter
        data_muscle_low_pass  = filtfilt(low_pass_b, low_pass_a, data_muscle_high_pass);

        % Rectify the signal
        data_muscle_rect = abs(data_muscle_low_pass);

        % Apply low-pass filter to the rectified signal
        data_muscle_rect_filt = filtfilt(low_2b, low_2a, data_muscle_rect);

        emg.filtered.(name_muscles{idx_channel}) = data_muscle_rect_filt;

        figure(figure_emg);
        clf;
        hold on;
%         plot(emg.(name_muscles{idx_channel}), 'Display', 'raw');
        plot(emg.filtered.(name_muscles{idx_channel}), 'Display', 'filtered');
        drawnow;
        legend()
        waitforbuttonpress;
    end
end