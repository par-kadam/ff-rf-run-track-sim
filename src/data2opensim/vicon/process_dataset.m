function process_dataset()

    import org.opensim.modeling.*;

    % Ensure that all configuration variables are properly loaded.
    % config = load_config();
    % config = load_config_20240601_dataset();
    % config = load_config_20240822_dataset();
    % config = load_config_20240904_dataset();
    % config = load_config_20250306_dataset();
    config = load_config_20250524_dataset();
    

    % Note: the expected structure of the directory is:
    %       <subject directory>/<sesion directory>/<trials directory>

    %% Load all data from *.c3d files into RAM.
    % -> load all EMG and Vicon data into the memory.
    [dataset, figures] = load_trials(config);

    %% Bring all data at 1000Hz.
    % -> filter and upsample the treadmill data to 1000Hz.
    % -> filter and downsample the EMG data to 1000Hz.
    % -> filter the keep force plate data to 1000Hz.
    % -> filter and upsample the marker data to 1000Hz.
    [dataset, figures] = filter_and_resample_data(dataset, config, figures);

    %% Scale the OpenSim model to the data of each subject.
    [dataset, figures] = scale_model_to_data(dataset, config, figures);

    analyze_all_gait_cycles(dataset, config);
    %% Split each trial in gait cycles.
    [dataset, figures] = extract_gait_cycles(dataset, config, figures);

    %% Export each gait cycle extracted in OpenSim formats.
    [dataset, figures] = compute_ik_and_id(dataset, config, figures);

    %% Normalize all data
    % -> emg
    %    - relative to DMVC the maximum activation during a trial.
    %    - relative to MVIC the activation during maximum isokinetic contraction
    % -> force plate
    %    - relative to subject body weight.
    % -> emg + grf + marker
    %    - relative to the completeness of the gait cycle.
    % [dataset, figures] = perform_all_normalizations(dataset, config, figures);

    %% Compute the means and standard deviations on all data.
    % [dataset, figures] = compute_means_and_stds(dataset, config, figures);

end

%
function [figures] = create_figures(config)

    figures = [];

    % Treadmill figures
    if config.visualize.treadmill_data

        figures.treadmill.fig = figure;
        set(figures.treadmill.fig, 'MenuBar', 'none');
        set(figures.treadmill.fig, 'ToolBar', 'none');

        % Axes distances.
        figures.treadmill.distances = axes('Position', [0.05, 0.76, 0.93, 0.20]); % [left, bottom, width, height]
        set(figures.treadmill.distances, 'FontSize', 4);
        set(figures.treadmill.distances, 'XLabel', []);
        set(figures.treadmill.distances, 'XTickLabel', []);

        % xlim(figures.treadmill.distances, );
        % ylim(figures.treadmill.distances, );

        % label_x          = xlabel(figures.treadmill.distances, "Time (s)");
        label_y          = ylabel(figures.treadmill.distances, "Distance (m)");
        legend_distances = legend(figures.treadmill.distances);
        title_distances  =  title(figures.treadmill.distances, "Distances");

        % set(label_x,          'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(label_y,          'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(legend_distances, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(title_distances,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

        % Axes speeds.
        figures.treadmill.speeds = axes('Position', [0.05, 0.52, 0.93, 0.20]);
        set(figures.treadmill.speeds, 'FontSize', 4);
        set(figures.treadmill.speeds, 'XLabel', []);
        set(figures.treadmill.speeds, 'XTickLabel', []);

        % xlim(figures.treadmill.speeds, );
        % ylim(figures.treadmill.speeds, );

        % label_x       = xlabel(figures.treadmill.speeds, "Time (s)");
        label_y       = ylabel(figures.treadmill.speeds, "Speed (m/s)");
        legend_speeds = legend(figures.treadmill.speeds);
        title_speeds  =  title(figures.treadmill.speeds, "Speeds");

        % set(label_x,       'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(label_y,       'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(legend_speeds, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(title_speeds,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

        % Axes pitch.
        figures.treadmill.pitch = axes('Position', [0.05, 0.28, 0.93, 0.20]);
        set(figures.treadmill.pitch, 'FontSize', 4);
        set(figures.treadmill.pitch, 'XLabel', []);
        set(figures.treadmill.pitch, 'XTickLabel', []);

        % xlim(figures.treadmill.pitch, );
        % ylim(figures.treadmill.pitch, );

        title_pitch  = title(figures.treadmill.pitch, "PitchRotX");
        % label_x      = xlabel(figures.treadmill.pitch, "Time (s)");
        label_y      = ylabel(figures.treadmill.pitch, "PitchRotX (rad)");
        legend_pitch = legend(figures.treadmill.pitch);

        % set(label_x,      'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(label_y,      'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(legend_pitch, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(title_pitch,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

        % Axes sway.
        figures.treadmill.sway = axes('Position', [0.05, 0.04, 0.93, 0.20]);
        set(figures.treadmill.sway, 'FontSize', 4);
        % set(figures.treadmill.sway, 'XLabel', []);
        % set(figures.treadmill.sway, 'XTickLabel', []);

        % xlim(figures.treadmill.sway, );
        % ylim(figures.treadmill.sway, );

        label_x     = xlabel(figures.treadmill.sway, "Time (s)");
        label_y     = ylabel(figures.treadmill.sway, "SwayPosX (m)");
        legend_sway = legend(figures.treadmill.sway);
        title_sway  =  title(figures.treadmill.sway, "SwayPosX");

        set(label_x,     'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(label_y,     'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(legend_sway, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(title_sway,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

        drawnow;
    end

    % GRFs figure, both force plate data in the same figure.
    if config.visualize.force_plate_data

        figures.grf.fig = figure;
        set(figures.grf.fig, 'MenuBar', 'none');
        set(figures.grf.fig, 'ToolBar', 'none');

        figures.grf.forces = axes('Position', [0.05, 0.10, 0.93, 0.85]); % [left, bottom, width, height]
        set(figures.grf.forces, 'FontSize', 4);

        % xlim(figures.grf.forces, );
        ylim(figures.grf.forces, [-50.0, +2500.0]);

        label_x    = xlabel(figures.grf.forces, "Time (s)");
        label_y    = ylabel(figures.grf.forces, "Force (N)");
        legend_grf = legend(figures.grf.forces);
        title_grf  =  title(figures.grf.forces, "GRF");

        set(label_x,    'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(label_y,    'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(legend_grf, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
        set(title_grf,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

        drawnow;
    end

    % EMGs figures: one figure per muscle.
    if config.visualize.emg_data

        figures.emg.figs = cell(length(config.emg_data.name_muscles));
        for idx_mus=1:length(figures.emg.figs)

            figures.emg.figs{idx_mus}.fig = figure;
            set(figures.emg.figs{idx_mus}.fig, 'MenuBar', 'none');
            set(figures.emg.figs{idx_mus}.fig, 'ToolBar', 'none');

            figures.emg.figs{idx_mus}.activation = axes('Position', [0.05, 0.10, 0.93, 0.85]); % [left, bottom, width, height]
            set(figures.emg.figs{idx_mus}.activation, 'FontSize', 4);
            % xlim(fitures.emg.figs{idx_mus}.activation, );
            % ylim(fitures.emg.figs{idx_mus}.activation, );

            label_x    = xlabel(figures.emg.figs{idx_mus}.activation, 'Time (ms)');
            label_y    = ylabel(figures.emg.figs{idx_mus}.activation, 'Electric potential (mV)');
            legend_emg = legend(figures.emg.figs{idx_mus}.activation);
            title_emg  =  title(figures.emg.figs{idx_mus}.activation, "EMG");

            set(label_x,    'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
            set(label_y,    'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
            set(legend_emg, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
            set(title_emg,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

            drawnow;
        end
    end

    % Vicon marker figures: one figure per marker.
    if config.visualize.marker_data

        figures.mocap.figs = cell(length(config.marker_data.marker_list));
        for idx_mar = 1:length(figures.mocap.figs)

            figures.mocap.figs{idx_mar}.fig = figure;
            set(figures.mocap.figs{idx_mar}.fig, 'MenuBar', 'none');
            set(figures.mocap.figs{idx_mar}.fig, 'ToolBar', 'none');

            figures.mocap.figs{idx_mar}.coordinates = axes('Position', [0.05, 0.10, 0.93, 0.85]); % [left, bottom, width, height]
            set(figures.mocap.figs{idx_mar}.coordinates, 'FontSize', 4);

            % xlim(fitures.mocap.figs{idx_mar}.coordinates, [?, ?]);
            % ylim(fitures.mocap.figs{idx_mar}.coordinates, [?, ?]);

            label_x      = xlabel(figures.mocap.figs{idx_mar}.coordinates, "Time (s)");
            label_y      = ylabel(figures.mocap.figs{idx_mar}.coordinates, "Coordinate (mm)");
            legend_mocap = legend(figures.mocap.figs{idx_mar}.coordinates);
            title_mocap  =  title(figures.mocap.figs{idx_mar}.coordinates, "Marker");

            set(label_x,      'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
            set(label_y,      'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
            set(legend_mocap, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
            set(title_mocap,  'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');

            drawnow;
        end
    end
end

%
function [] = clear_figures(config, figures)

    if config.visualize.treadmill_data
        cla(figures.treadmill.distances);
        cla(figures.treadmill.speeds);
        cla(figures.treadmill.pitch);
        cla(figures.treadmill.sway);
    end

    if config.visualize.emg_data
        for idx_mus=1:length(figures.emg.figs)
            cla(figures.emg.figs{idx_mus}.activation);
        end
    end

    if config.visualize.force_plate_data
        cla(figures.grf.forces);
    end

    if config.visualize.marker_data
        for idx_mar = 1:length(figures.mocap.figs)
            cla(figures.mocap.figs{idx_mar}.coordinates);
        end
    end
end

%%
function [dataset, figures] = load_trials(config)

    fprintf('Loading subject data ...\n');

    path_dataset = config.path_dataset;
    subjects     = config.subjects.names(config.subjects.process);
    dataset = [];
    for idx_sub = 1:length(subjects)

        % Identify the directories storing the sessions.
        sessions      = dir([path_dataset subjects{idx_sub}]);
        sessions      = cellfun(@(x) x.name, num2cell(sessions), 'UniformOutput', false);
        sessions      = sort(sessions);
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);

        for idx_ses = 1:length(sessions)

            % Identify the *.c3d files holding EMG and Vicon data.
            trials = dir([path_dataset subjects{idx_sub} '/' sessions{idx_ses}]);
            trials = cellfun(@(x) x.name, num2cell(trials), 'UniformOutput', false);

            idxs_trials_enf = cellfun(@(x) contains(x, 'Trial.enf'), trials);
            trials_enf      = trials(idxs_trials_enf);
            trials_enf      = sort(trials_enf);

            for idx_tri = 1:length(trials_enf)

                file_trial_stem = trials_enf{idx_tri}(1:end - 10);

                fprintf("Loading %s trial ...", file_trial_stem);

                file_trial_treadmill = [file_trial_stem, '_treadmill.csv'];
                file_trial_emg       = [file_trial_stem, '_emg.c3d'];
                file_trial_mocap     = [file_trial_stem, '.c3d'];

                path_trial_treadmill = [ ...
                    path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/', file_trial_treadmill];
                path_trial_emg = [ ...
                    path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/', file_trial_emg];
                path_trial_mocap = [ ...
                    path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/', file_trial_mocap];

                if ~(exist(path_trial_treadmill, 'file') == 2) || ...
                   ~(exist(path_trial_emg,       'file') == 2) || ...
                   ~(exist(path_trial_mocap,     'file') == 2)
                    fprintf("skipped!\n");
                    continue
                end

                condition  = ['con', file_trial_stem(11:12), file_trial_stem(18:22), file_trial_stem(29:33)];
                trial      = ['tri', file_trial_stem(35:38)];

                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).name_trial = file_trial_stem;

                % Load treadmill data.
                treadmill_data = readtable(path_trial_treadmill);

                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).treadmill = treadmill_data;

                % Load the EMG data.
                [analogs, analogs_info] = btkGetAnalogs(btkReadAcquisition(path_trial_emg));

                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).emg.analogs       = analogs;
                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).emg.analogs_info  = analogs_info;

                % Load the motion capture data.
                trial_mocap = btk_loadc3d(path_trial_mocap, 10);

                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).mocap        = trial_mocap;
                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).mocap.Name   = subjects{idx_sub};
                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).mocap.Mass   = config.subjects.masses(subjects{idx_sub});
                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition).(trial).mocap.Height = config.subjects.heights(subjects{idx_sub});

                fprintf("done!\n");
            end
        end
    end

    % Print summary of the subjects, sessions and trials loaded.
    subjects = fieldnames(dataset);
    for idx_sub = 1:length(subjects)
        fprintf("%s\n", subjects{idx_sub});
        sessions = fieldnames(dataset.(subjects{idx_sub}));
        for idx_ses = 1:length(sessions)
            fprintf(" %s\n", sessions{idx_ses});
            conditions = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            for idx_con = 1:length(conditions)
                fprintf("  %s\n", conditions{idx_con});
                trials = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                for idx_tri = 1:length(trials)
                    fprintf("   %s\n", trials{idx_tri});
                end
            end
        end
    end

    figures = create_figures(config);

    % Wait of key pressed before plotting the next trial.
    if config.visualize.force_plate_data && ...
       config.visualize.raw_data
        fprintf([ ...
            '\nReady to visualize the raw data?' ...
            '\nJust click on the GRF figure and press the "n" key!\n']);
        figure(figures.grf.fig);
        while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
            waitforbuttonpress;
        end
        set(figures.grf.fig, 'CurrentCharacter', 'x');
    end

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) contains(x, 'con'), conditions);
            conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                for idx_tri = 1:length(trials)

                    name_trial      = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                    data_treadmill  = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).treadmill;
                    data_mocap      = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;
                    data_emg        = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;

                    % Plot raw data

                    % Use the same color for data describing the right side.
                    color_right_side = rand(1, 3);

                    if config.visualize.raw_data && ...
                       config.visualize.treadmill_data

                        hold(figures.treadmill.distances, "on");
                        plot(figures.treadmill.distances, ...
                            data_treadmill.Time, data_treadmill.Leftbelt_Distance, "DisplayName", "left");
                        plot(figures.treadmill.distances, ...
                            data_treadmill.Time, data_treadmill.Rightbelt_Distance, "DisplayName", "right", "Color", color_right_side);
                        hold(figures.treadmill.distances, "off");

                        hold(figures.treadmill.speeds, "on");
                        plot(figures.treadmill.speeds, ...
                            data_treadmill.Time, data_treadmill.Leftbelt_Speed, "DisplayName", "left");
                        plot(figures.treadmill.speeds, ...
                            data_treadmill.Time, data_treadmill.Rightbelt_Speed, "DisplayName", "right", "Color", color_right_side);
                        hold(figures.treadmill.speeds, "off");

                        hold(figures.treadmill.pitch, "on");
                        plot(figures.treadmill.pitch, ...
                            data_treadmill.Time, data_treadmill.PitchRotX, "DisplayName", "pitch");
                        hold(figures.treadmill.pitch, "off");

                        hold(figures.treadmill.sway, "on");
                        plot(figures.treadmill.sway, ...
                            data_treadmill.Time, data_treadmill.SwayPosX, "DisplayName", "sway");
                        hold(figures.treadmill.sway, "off");

                        drawnow;
                    end

                    if config.visualize.raw_data && ...
                       config.visualize.emg_data

                        channel_muscles = config.emg_data.subjects.muscles(subjects{idx_sub});
                        name_muscles    = config.emg_data.name_muscles;
                        for idx_mus = 1:length(name_muscles)
                            data_muscle_raw = data_emg.analogs.(channel_muscles(name_muscles{idx_mus}));

                            hold(figures.emg.figs{idx_mus}.activation, "on");

                            plot(figures.emg.figs{idx_mus}.activation, ...
                                data_muscle_raw, ...
                                'DisplayName', sprintf('%s raw', name_muscles{idx_mus}), 'Color', color_right_side);

                            hold(figures.emg.figs{idx_mus}.activation, "off");

                            title_emg = [name_trial ' ' name_muscles{idx_mus}];
                            title(figures.emg.figs{idx_mus}.activation, title_emg);

                            drawnow;
                        end
                    end

                    if config.visualize.raw_data && ...
                       config.visualize.force_plate_data

                        hold(figures.grf.forces, "on");
                        for idx_grf = 1:length(data_mocap.fp_data.GRF_data)
                            if idx_grf == 2
                                plot(figures.grf.forces, data_mocap.fp_data.GRF_data(idx_grf).F(:, 3), ...
                                    'DisplayName', 'GRF2 (right leg) raw', 'Color', color_right_side);
                            else
                                plot(figures.grf.forces, data_mocap.fp_data.GRF_data(idx_grf).F(:, 3), ...
                                    'DisplayName', 'GRF1 (left leg) raw');
                            end
                        end
                        hold(figures.grf.forces, "off");

                        title_grf = [name_trial, ' GRFs'];
                        title(figures.grf.forces, title_grf);

                        drawnow;
                    end

                    if config.visualize.raw_data && ...
                       config.visualize.marker_data

                        markers = fieldnames(data_mocap.marker_data.Markers);
                        for idx_mar = 1:min(length(figures.mocap.figs), length(markers))

                            hold(figures.mocap.figs{idx_mar}.coordinates, "on");
                            name_coords = {'x', 'y', 'z'};
                            for idx_coord = 1:length(name_coords)
                                legend_coord = sprintf('%s %s raw', ...
                                    markers{idx_mar}, name_coords{idx_coord});
                                plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                    data_mocap.marker_data.Markers.(markers{idx_mar})(:, idx_coord), ...
                                    'DisplayName', legend_coord);
                            end
                            hold(figures.mocap.figs{idx_mar}.coordinates, "off");

                            title_marker = [name_trial ' ' markers{idx_mar}];
                            title(figures.mocap.figs{idx_mar}.coordinates, title_marker);

                            drawnow;
                        end
                    end

%                     % Wait of key pressed before plotting the next trial.
%                     if config.visualize.force_plate_data && ...
%                        config.visualize.raw_data
%                         fprintf([ ...
%                             '\nReady to visualize the raw data of the next trial?' ...
%                             '\nJust click on the GRF figure and press the "n" key!\n']);
%                         figure(figures.grf.fig);
%                         while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
%                             waitforbuttonpress;
%                         end
%                         set(figures.grf.fig, 'CurrentCharacter', 'x');
%                     end

                    % Semi-interactive visualization
                    if config.visualize.raw_data&& ...
                       (config.visualize.subject || ...
                        config.visualize.session || ...
                        config.visualize.condition || ...
                        config.visualize.trial || ...
                        config.visualize.gait_cycle)

                        pause(config.visualize.pause)
                    end

                    if config.visualize.raw_data && ...
                       config.visualize.trial
                        clear_figures(config, figures);
                    end
                end
            end
        end
    end
end

%%
function [dataset, figures] = filter_and_resample_data(dataset, config, figures)

    if config.visualize.force_plate_data && ...
       config.visualize.filter_and_resample_data
        fprintf([ ...
            '\nReady to begin the filtering and resmpling of the data?' ...
            '\nJust click on the GRF figure and press the "n" key!\n']);
        figure(figures.grf.fig);
        while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
            waitforbuttonpress;
        end
        set(figures.grf.fig, 'CurrentCharacter', 'x');
    end

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            % Consider all conditions.
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) contains(x, 'con'), conditions);
            conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                for idx_tri = 1:length(trials)

                    % Note: data contains both: force plate signals and
                    % Vicon's marker trajectories.
                    name_trial = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                    % todo(lisca): filter the treadmill data too!
                    data_mocap = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;
                    data_emg   = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;

                    fprintf('%s filtering ...', name_trial);

                    frequency_force_plate = data_mocap.fp_data.Info.frequency;
                    frequency_emg         = data_emg.analogs_info.frequency;
                    factor_down           = frequency_emg / frequency_force_plate;

                    %
                    % 1. Filter and downsample the EMG data.
                    %


                    channel_muscles = config.emg_data.subjects.muscles(subjects{idx_sub});
                    name_muscles    = config.emg_data.name_muscles;
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.raw         = [];
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.rect_filt   = [];
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.downsampled = [];
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean       = [];
                    for idx_mus = 1:length(name_muscles)

                        channel_muscle = channel_muscles(name_muscles{idx_mus});
                        data_muscle    = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.analogs.(channel_muscle);

                        % Save a copy fo the original raw data.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.raw.(name_muscles{idx_mus}) = ...
                            data_muscle;

                        % Detrend
                        data_muscle_detrended = detrend(data_muscle);

                        % Check if this is a maximum voluntary contraction trial
                        is_mvc_trial = contains(conditions{idx_con}, 'conmc');

                        % Process regular trials with full filtering
                        if ~is_mvc_trial

                            % Band-pass filter first
                            low_cutoff = 50
                            high_cutoff = 500;
                            normalized_emg_cutoff = [low_cutoff, high_cutoff] / (frequency_emg / 2);
                            [b, a] = butter(config.force_plate_data.order, normalized_emg_cutoff, 'bandpass');
                            
                            % High-pass filter design (Butterworth, 8 Hz cutoff)
                            % [high_pass_b, high_pass_a] = butter(4, 8/(0.5*frequency_emg), 'high');
                            % Low-pass filter design for EMG (Butterworth, 500 Hz cutoff)
                            % [low_pass_b, low_pass_a] = butter(4, 500/(0.5*frequency_emg), 'low');
                            % Low-pass filter design for rectified EMG (Butterworth, 6 Hz cutoff)
                            [low_2b, low_2a] = butter(4, 15/(0.5*frequency_emg), 'low');

                            % Apply high-pass filter
                            % data_muscle_high_pass = filtfilt(high_pass_b, high_pass_a, data_muscle_detrended);
                            % Apply low-pass filter
                            % data_muscle_low_pass = filtfilt(low_pass_b, low_pass_a, data_muscle_high_pass);
                            
                            data_muscle_bandpassed = filtfilt(b, a, data_muscle_detrended);
                            % Rectify the signal
                            data_muscle_rect = abs(data_muscle_bandpassed);
                            
                            % Apply low-pass filter to the rectified signal
                            data_muscle_rect_filt = filtfilt(low_2b, low_2a, data_muscle_rect);

                            % else
                                % For MVC trials, skip filtering but still perform rectification
                            %     fprintf('Skipping filters for MVC trial: %s\n', trials{idx_tri});
                                
                                % Bandpass filter first

                                % Directly rectify the detrended signal without filtering
                            %     data_muscle_rect = abs(data_muscle_detrended);
                                
                                % Skip the final envelope filter as well
                            %     data_muscle_rect_filt = data_muscle_rect;
                        end

                       
                        % Save the rectified, then filtered data.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.rect_filt.(name_muscles{idx_mus}) = ...
                            data_muscle_rect_filt;

                        % Downsampling: 2000Hz -> 1000Hz.
                        data_muscle_rect_filt_downsampled = ...
                            decimate(data_muscle_rect_filt, factor_down);

                        % Save the filtered, detrended, and downsampled data, just for plotting it later.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.downsampled.(name_muscles{idx_mus}) = ...
                            data_muscle_rect_filt_downsampled;

                        % Save the filtered, detrended, and downsampled data.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus}) = ...
                            data_muscle_rect_filt_downsampled;

                    end

                    %
                    % 2. Filter the force plates data.
                    %

                    normalizedCutoff = config.force_plate_data.cutoffFrequency / (0.5 * frequency_force_plate);
                    [b, a]           = butter(config.force_plate_data.order, normalizedCutoff, 'low');

                    for idx_grf = 1:length(data_mocap.fp_data.GRF_data)

                        % Save a copy of the original raw data
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data_raw(idx_grf).F = ...
                        data_mocap.fp_data.GRF_data(idx_grf).F;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data_raw(idx_grf).P = ...
                        data_mocap.fp_data.GRF_data(idx_grf).P;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data_raw(idx_grf).M = ...
                        data_mocap.fp_data.GRF_data(idx_grf).M;

                        %***CONSIDER ADDING RECTIFICATION IN HERE
                        % Apply zero-lag, low-pass Butterworth filter using filtfilt.
                        grf_filt    = filtfilt(b, a, data_mocap.fp_data.GRF_data(idx_grf).F);
                        %grf_filt(:, 2)= abs(grf_filt(:, 2));
                        grf_movmean = movmean(grf_filt, 50);
                        grf_filt    = sgolayfilt(grf_movmean, 2, 101);

                        grf_filt_torque = filtfilt(b, a, data_mocap.fp_data.GRF_data(idx_grf).M);
                        grf_filt_position = filtfilt(b, a, data_mocap.fp_data.GRF_data(idx_grf).P);

                        % Save the filtered data.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).F = ...
                            grf_filt;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).M = ...
                            grf_filt_torque;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).P = ...
                            grf_filt_position;
                        
                        % Store back the smoothed data
                        % dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).F = combined_data;
                        
                    end

                    % ----------Threshold condition statement to zero out plate values during swing phase
                    % After both plates are filtered, handle threshold transitions and zeroing


                    %
                    % 3. Filter and upsample the marker data.
                    %

                    factor_upsampling_marker = ...
                        data_mocap.fp_data.Info(1).frequency / data_mocap.marker_data.Info.frequency;

                    markers = fieldnames(data_mocap.marker_data.Markers);
                    for idx_mar = 1:length(markers)

                        % Save a copy of the original raw data.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Markers_raw.(markers{idx_mar}) = ...
                            data_mocap.marker_data.Markers.(markers{idx_mar});

                        % Apply movemean filter on marker's data.
                        data_marker_filt = ...
                            movmean(data_mocap.marker_data.Markers.(markers{idx_mar}), config.marker_data.moving_average_window_size);

                        % Save the filtered data.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Markers_filt.(markers{idx_mar}) = ...
                            data_marker_filt;

                        % Apply interpolate each coordinate (x, y, z).
                        x = interp(data_marker_filt(:, 1), factor_upsampling_marker);
                        y = interp(data_marker_filt(:, 2), factor_upsampling_marker);
                        z = interp(data_marker_filt(:, 3), factor_upsampling_marker);

                        % Save the data ready to be processed further.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Markers_upsampled.(markers{idx_mar}) = ...
                            [x, y, z];

                        % ! ! ! Overwrite the Markers field.
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Markers.(markers{idx_mar}) = ...
                            [x, y, z];

                     end

                    % Upsample the timeline.
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Time = ...
                        interp(data_mocap.marker_data.Time, factor_upsampling_marker);

                    % Update the metadata.
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Info.frequency = ...
                        data_mocap.marker_data.Info.frequency * factor_upsampling_marker;
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Info.NumFrames = ...
                        length(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Time);

                    fprintf('done\n');

                    %
                    % Check the lengths of the arrays: EMG, GRF, Mocap and Time.
                    %

                    % EMG data
                    lengths = [];
                    for idx_mus = 1:length(name_muscles)
                        lengths = [lengths, ...
                            length(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus}))];
                    end

                    % GRF data
                    for idx_grf = 1:length(data_mocap.fp_data.GRF_data)
                        lengths = [lengths, ...
                            length( ...
                                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).F)];
                    end

                    % Mocap data
                    markers = fieldnames(data_mocap.marker_data.Markers);
                    for idx_mar = 1:length(markers)
                        lengths = [lengths, ...
                            length( ...
                                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Markers.(markers{idx_mar}))];
                    end

                    % Time data
                    lengths = [lengths, ...
                        length(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Time)];

                    lengths = lengths - min(lengths);

%                     if ~all(diff(lengths) == 0)
%                         disp(name_trial);
%                         disp(lengths);
%                         fprintf('%s: data arrays have different lengths ! ! !\n', name_trial);
%                     else
%                         fprintf('%s: data arrays have the same lengths -> all good.\n', name_trial);
%                     end

                    % Visualization
                    if config.visualize.filter_and_resample_data && ...
                       config.visualize.emg_data

                        for idx_mus = 1:length(name_muscles)

                            cla(figures.emg.figs{idx_mus}.activation);

                            hold(figures.emg.figs{idx_mus}.activation, "on");

                            plot(figures.emg.figs{idx_mus}.activation, ...
                                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus}), ...
                                'DisplayName', 'rectified, filtered, downsampled');

                            hold(figures.emg.figs{idx_mus}.activation, "off");

                            title_emg = sprintf("%s %s",name_trial, name_muscles{idx_mus});
                            title(figures.emg.figs{idx_mus}.activation, title_emg);

                            drawnow;
                         end
                    end

                    if config.visualize.filter_and_resample_data && ...
                       config.visualize.force_plate_data

                        % Use the same color for data describing the right side.
                        color_right_side = rand(1, 3);

                        % Update the GRFs figure.
                        hold(figures.grf.forces, "on");

                        for idx_grf = 1:length(data_mocap.fp_data.GRF_data)
                            if idx_grf == 2
                                plot(figures.grf.forces, ...
                                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).F(:, 3), ...
                                    'DisplayName', 'GRF2 - right filtered', 'Color', color_right_side);
                            else
                                plot(figures.grf.forces, ...
                                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).F(:, 3), ...
                                    'DisplayName', 'GRF1 - left filtered');
                            end
                        end

                        hold(figures.grf.forces, "off");

                        figure(figures.grf.fig);
                        title_grf = sprintf("%s GRF", name_trial);
                        title(figures.grf.forces, title_grf);

                        drawnow;
                    end

                    if config.visualize.filter_and_resample_data && ...
                       config.visualize.marker_data

                        markers = fieldnames(data_mocap.marker_data.Markers);
                        for idx_mar = 1:min(length(figures.mocap.figs), length(markers))

                            hold(figures.mocap.figs{idx_mar}.coordinates, "on");
                            name_coords = {'x', 'y', 'z'};
                            for idx_coord = 1:length(name_coords)
                                legend_coord = sprintf('%s %s filtered and upsampled', markers{idx_mar}, name_coords{idx_coord});
                                plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.marker_data.Markers.(markers{idx_mar})(:, idx_coord), ...
                                    'DisplayName', legend_coord);
                            end
                            hold(figures.mocap.figs{idx_mar}.coordinates, "off");

                            title_marker = sprintf("%s %s", name_trial, markers{idx_mar});
                            title(figures.mocap.figs{idx_mar}.coordinates, title_marker);

                            drawnow;
                        end
                    end

%                     % Plot one trial at the time.
%                     if config.visualize.force_plate_data && ...
%                        config.visualize.filter_and_resample_data
%                         fprintf([ ...
%                             '\nReady to visualize the next trial?' ...
%                             '\nJust click on the GRF figure and press the "n" key!\n']);
%                         figure(figures.grf.fig);
%                         while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
%                             waitforbuttonpress;
%                         end
%                         set(figures.grf.fig, 'CurrentCharacter', 'x');
%                     end

                    clear_figures(config, figures);

                end % loop through trials
            end % loop through conditions
        end % loop through sessions
    end % loop through subjects
end

%% Scale the OpenSim model to marker data.
function [dataset, figures] = scale_model_to_data(dataset, config, figures)

%     if config.visualize.force_plate_data && ...
%        config.visualize.filter_and_resample_data
%         fprintf([ ...
%             '\nReady to scale the OpenSim template model to the Mocap data of the subjects?' ...
%             '\nJust click on the GRF figure and press the "n" key!\n']);
%         figure(figures.grf.fig);
%         while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
%             waitforbuttonpress;
%         end
%         set(figures.grf.fig, 'CurrentCharacter', 'x');
%     end

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            % Consider only the static codition specified in the config file
            % -> filter out the conditions: maximum contraction and dynamic.
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) contains(x, 'const'), conditions);
            conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                % Filter out all static trials, except the trial specified in by config.subjects.trial_static.
                idxs_trials = cellfun(@(x) contains(x, 'tri') & contains(x, config.subjects.trial_static(subjects{idx_sub})) , trials);
                trials      = trials(idxs_trials);
                for idx_tri = 1:length(trials)

                    % Note: data contains both: force plate signals and
                    % Vicons marker trajectories.
                    name_trial = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                    data_mocap = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;

                    % Sort the C3D file so we know what is Marker data and
                    % what is calculated data.
                    data_mocap.marker_data = btk_sortc3d(...
                        data_mocap.marker_data, config.marker_data.marker_list, config.marker_data.calc_marker_list);

                    % Use a few frames of the static trial in order to scale the template model
                    % to the data.
                    data_mocap.Start_Frame = 11;
                    data_mocap.End_Frame   = 20;

                    % now do conversions to TRC files using btk_c3d2trc
                    TRCFileName = [name_trial '_trc.trc'];
                    GRFFileName = [name_trial '_grf.sto'];

                    % Switch to 'off' to turn animation off.
                    data = btk_c3d2trc(data_mocap, 'off', TRCFileName, GRFFileName);

                    subject_label = config.subjects.names;
                    current_subject = subject_label{1};

                    %model_port_path = ['../../../../biomec/data/20250306_dataset/sub03/ses20250306/sub03_SCALED_addbio_free_subt_mtp.osim'];
                    model_port_path = [config.results_directory current_subject '_SCALED_addbio_free_subt_mtp.osim'];
                    
                    %if ~exist(model_port_path, 'dir')
                    %        mkdir(model_port_path)
                    %end

                    % Execute the copy and rename operation
                    copyfile(config.path_model_physics_fitted, model_port_path);

                    % Define the standard config files OpenSims scaling tool.
                    ModelFile =          [config.path_model config.file_model '.osim'];   % ['sub03_SCALED_addbio_free_subt_mtp.osim']; 
                    ScaleTasksFile =     [config.path_model config.file_model '_Scale_Tasks.xml'];
                    MeasurementSetFile = [config.path_model config.file_model '_Scale_MeasurementSet.xml'];
                    MarkerSetFile =      [config.path_model config.file_model '_MarkerSetFile.xml'];

                    cd([config.path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/']);
                    setup_scale( ...
                        'data',data, ...
                        'ModelFile',ModelFile, ...
                        'ScaleTasksFile',ScaleTasksFile, ...
                        'MeasurementSetFile',MeasurementSetFile, ...
                        'MarkerSetFile',MarkerSetFile, ...
                        'PreserveMass','true');

                    com = ['opensim-cmd run-tool ' data_mocap.Name '_Setup_Scale.xml'];

                    % Get both stdout and stderr
                    [status, cmdout] = system(com);
                    fprintf('\n== Scale Tool Output ==\n');
                    fprintf('%s\n', cmdout);
                    fprintf('==End Scale Tool Output==\n\n');

                    fprintf('\n=== Running IK on Static Trial ===\n');

                    % Find the static trial TRC file created during scaling
                    subject_name = subjects{idx_sub}
                    static_trial_name = config.subjects.trial_static(subject_name);
                    static_trc_pattern = sprintf('%s_strist_cond00000_speed00000_%s_trc.trc', subject_name, static_trial_name);
                    static_files = dir(static_trc_pattern);

                    if ~isempty(static_files)
                        % Use the first matching static TRC file
                        static_trc_filename = static_files(1).name;
                        static_base_name = strrep(static_trc_filename, '_trc.trc', '');
                        
                        fprintf('Found static TRC file: %s\n', static_trc_filename);
                        
                        % Setup IK for static trial
                        static_ik_output = [static_base_name '_static_ik.sto'];
                        
                        % Create a temporary data structure for static IK setup
                        static_data.Name = static_base_name;
                        static_data.TRC_Filename = static_trc_filename;
                        
                        % Define the standard config files OpenSims IK tool for static trial
                        ModelFile = [subject_name '_SCALED.osim'];  % Use the scaled model
                        IKTasksFile = [config.path_model config.file_model '_IK_Tasks.xml'];
                        
                        % Setup static IK with higher accuracy for static pose
                        setup_InverseKinematics( ...
                            'data',        data, ...
                            'ModelFile',   ModelFile,...
                            'IKTasksFile', IKTasksFile, ...
                            'Accuracy',    0.000001, ...  % Higher accuracy for static trial
                            'OutputFile',  static_ik_output);
                        
                        % Run the static IK tool and capture output
                        static_com = ['opensim-cmd run-tool ' data_mocap.Name '_Setup_InverseKinematics.xml'];
                        [static_status, static_cmdout] = system(static_com);
                        
                        % Print the captured output
                        fprintf('\n=== Static IK Tool Output ===\n');
                        fprintf('%s\n', static_cmdout);
                        fprintf('=== End Static IK Tool Output ===\n\n');
                        
                        % Parse the static IK output for errors
                        [static_time_values, static_max_errors, static_max_error_markers, static_rms_errors] = parse_ik_output(static_cmdout);
                        
                        static_frame_numbers = 0:(length(static_time_values) - 1);
                        static_total_squared_errors = NaN(size(static_time_values));
                        % Write static IK results to files
                        if ~isempty(static_time_values)
                            write_ik_results_to_files(static_time_values, static_max_errors, static_max_error_markers, ...
                                                    static_rms_errors, static_frame_numbers, static_total_squared_errors, ...
                                                    static_base_name, 0, pwd); % Use 0 for static trial cycle number
                        end
                        
                        % Load and analyze the static IK results
                        if exist(static_ik_output, 'file')
                            fprintf('Static IK completed successfully. Output file: %s\n', static_ik_output);
                            
                            % Load the static IK data for analysis
                            static_ik_data = load_sto_file(static_ik_output);
                            
                            % Analyze static pose (optional - for quality assessment)
                            analyze_static_pose(static_ik_data, static_base_name, subjects{idx_sub});
                            
                            % Store static IK results in dataset for later reference
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).static_ik.filename = static_ik_output;
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).static_ik.data = static_ik_data;
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).static_ik.errors.time_values = static_time_values;
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).static_ik.errors.max_errors = static_max_errors;
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).static_ik.errors.rms_errors = static_rms_errors;
                            
                        else
                            fprintf('Warning: Static IK output file not found: %s\n', static_ik_output);
                        end
                        
                    else
                        fprintf('Warning: No static TRC file found matching pattern: %s\n', static_trc_pattern);
                    end


                    % Print captured output line by line
                    fprintf('\n== Scale Tool Output ==\n');
                    fprintf('%s\n', cmdout);
                    fprintf('==End Scale Tool Output==\n\n');

                    % Parse output for extracting scaling results
                    [scale_factors, measurements, marker_errors, generic_measurements, subject_measurements] = parse_scale_output(cmdout);
                    if ~isempty(scale_factors)
                        write_scale_results_to_files(scale_factors, measurements, marker_errors, generic_measurements, subject_measurements, ...
                            data_mocap.Name, subjects{idx_sub});
                    end
                    % system(com);
                end % loops throught all trials
            end % loops throught all conditions
        end % loops throught all sessions
    end % loops throught all subjects
end

% Function to parse Scale tool output and extract scaling information
function [scale_factors, measurements, marker_errors, generic_measurements, subject_measurements] = parse_scale_output(cmdout)
    % Initialize output structures
    scale_factors = struct();
    measurements = struct();
    marker_errors = struct();
    generic_measurements = [];
    subject_measurements = [];
    
    % Split output into lines
    if exist('splitlines', 'builtin') || exist('splitlines', 'file')    
        lines = splitlines(cmdout);
    else
        lines = regexp(cmdout, '\r?\n', 'split');
        if ischar(lines)
            lines = {lines}; % Ensure lines is a cell array
        end
    end
    
    current_section = '';
    
    % Parse each line looking for scaling information
    for i = 1:length(lines)
        if iscell(lines)
            line = strtrim(lines{i});
        else
            line = strim(lines(i));
        end
        
        % Skip empty lines
        if isempty(line)
            continue;
        end
        
        % Identify sections using case-insensitive matching
        line_lower = lower(line);
        if contains(line_lower, 'scale factors') || contains(line_lower, 'scale factor')
            current_section = 'scale_factors';
            continue;
        elseif contains(line_lower, 'model measurements') || contains(line_lower, 'generic measurements')
            current_section = 'model_measurements';
            continue;
        elseif contains(line_lower, 'subject measurements') || contains(line_lower, 'experimental measurements')
            current_section = 'subject_measurements';
            continue;
        elseif contains(line_lower, 'marker errors') || contains(line_lower, 'rms marker error')
            current_section = 'marker_errors';
            continue;
        elseif contains(line_lower, 'total squared error')
            current_section = 'total_error';
            continue;
        end
        
        % Parse based on current section
        switch current_section
            case 'scale_factors'
                % Look for scale factor lines: "body_name: X = 0.xxxx, Y = 0.xxxx, Z = 0.xxxx"
                scale_pattern = '([^:]+):\s*X\s*=\s*([\d\.-]+),?\s*Y\s*=\s*([\d\.-]+),?\s*Z\s*=\s*([\d\.-]+)';
                tokens = regexp(line, scale_pattern, 'tokens');
                if ~isempty(tokens)
                    body_name = strtrim(tokens{1}{1});
                    % Make valid field name and remove special characters
                    body_name = regexprep(body_name, '[^a-zA-Z0-9_]', '_');
                    if isempty(body_name) || ~isvarname(body_name)
                        body_name = sprintf('body_%d', i);
                    end
                    scale_factors.(body_name).X = str2double(tokens{1}{2});
                    scale_factors.(body_name).Y = str2double(tokens{1}{3});
                    scale_factors.(body_name).Z = str2double(tokens{1}{4});
                    
                    fprintf('Parsed scale factors for %s: X=%.4f, Y=%.4f, Z=%.4f\n', ...
                        body_name, scale_factors.(body_name).X, ...
                        scale_factors.(body_name).Y, scale_factors.(body_name).Z);
                end
                
            case 'model_measurements'
                % Look for measurement lines: "measurement_name = value" or "measurement_name: value"
                meas_pattern1 = '([^=]+)=\s*([\d\.-]+)';
                meas_pattern2 = '([^:]+):\s*([\d\.-]+)';
                tokens = regexp(line, meas_pattern1, 'tokens');
                if isempty(tokens)
                    tokens = regexp(line, meas_pattern2, 'tokens');
                end
                if ~isempty(tokens)
                    meas_name = strtrim(tokens{1}{1});
                    % Clean up measurement name
                    meas_name = regexprep(meas_name, '[^a-zA-Z0-9_]', '_');
                    if isempty(meas_name) || ~isvarname(meas_name)
                        meas_name = sprintf('measurement_%d', length(fieldnames(measurements)));
                    end
                    measurements.model.(meas_name) = str2double(tokens{1}{2});
                    generic_measurements(end+1) = str2double(tokens{1}{2});
                    
                    fprintf('Parsed model measurement %s: %.6f\n', meas_name, measurements.model.(meas_name));
                end
                
            case 'subject_measurements'
                % Look for measurement lines: "measurement_name = value" or "measurement_name: value"
                meas_pattern1 = '([^=]+)=\s*([\d\.-]+)';
                meas_pattern2 = '([^:]+):\s*([\d\.-]+)';
                tokens = regexp(line, meas_pattern1, 'tokens');
                if isempty(tokens)
                    tokens = regexp(line, meas_pattern2, 'tokens');
                end
                if ~isempty(tokens)
                    meas_name = strtrim(tokens{1}{1});
                    % Clean up measurement name
                    meas_name = regexprep(meas_name, '[^a-zA-Z0-9_]', '_');
                    if isempty(meas_name) || ~isvarname(meas_name)
                        meas_name = sprintf('measurement_%d', length(fieldnames(measurements)));
                    end
                    measurements.subject.(meas_name) = str2double(tokens{1}{2});
                    subject_measurements(end+1) = str2double(tokens{1}{2});
                    
                    fprintf('Parsed subject measurement %s: %.6f\n', meas_name, measurements.subject.(meas_name));
                end
                
            case 'marker_errors'
                % Look for marker error lines: "marker_name: error = value" or "marker_name = value"
                error_pattern1 = '([^:]+):\s*error\s*=\s*([\d\.-]+)';
                error_pattern2 = '([^:=]+)[:\s=]\s*([\d\.-]+)';
                tokens = regexp(line, error_pattern1, 'tokens');
                if isempty(tokens)
                    tokens = regexp(line, error_pattern2, 'tokens');
                end
                if ~isempty(tokens)
                    marker_name = strtrim(tokens{1}{1});
                    % Clean up marker name
                    marker_name = regexprep(marker_name, '[^a-zA-Z0-9_]', '_');
                    if isempty(marker_name) || ~isvarname(marker_name)
                        marker_name = sprintf('marker_%d', length(fieldnames(marker_errors)));
                    end
                    marker_errors.(marker_name) = str2double(tokens{1}{2});
                    
                    fprintf('Parsed marker error %s: %.6f\n', marker_name, marker_errors.(marker_name));
                end
                
                % Also look for RMS error summary
                rms_pattern = 'RMS.*error.*[=:]\s*([\d\.-]+)';
                rms_tokens = regexp(line, rms_pattern, 'tokens');
                if ~isempty(rms_tokens)
                    marker_errors.RMS_total = str2double(rms_tokens{1}{1});
                    fprintf('Parsed RMS total error: %.6f\n', marker_errors.RMS_total);
                end
                
            case 'total_error'
                % Look for total squared error
                total_pattern = 'Total.*squared.*error.*[=:]\s*([\d\.-]+)';
                tokens = regexp(line, total_pattern, 'tokens');
                if ~isempty(tokens)
                    marker_errors.total_squared_error = str2double(tokens{1}{1});
                    fprintf('Parsed total squared error: %.6f\n', marker_errors.total_squared_error);
                end
        end
    end
    
    fprintf('\nParsed Scale tool output:\n');
    if ~isempty(fieldnames(scale_factors))
        fprintf('  Scale factors for %d body segments\n', length(fieldnames(scale_factors)));
    end
    if isfield(measurements, 'model') && ~isempty(fieldnames(measurements.model))
        fprintf('  %d model measurements\n', length(fieldnames(measurements.model)));
    end
    if isfield(measurements, 'subject') && ~isempty(fieldnames(measurements.subject))
        fprintf('  %d subject measurements\n', length(fieldnames(measurements.subject)));
    end
    if ~isempty(fieldnames(marker_errors))
        fprintf('  Marker errors for %d markers\n', length(fieldnames(marker_errors)));
    end
end

% Function to write Scale results to XML and LOG files
function write_scale_results_to_files(scale_factors, measurements, marker_errors, ...
                                      generic_measurements, subject_measurements, ...
                                      trial_name, subject_name)
    
    % Create base filename
    base_filename = sprintf('%s_scale', trial_name);
    xml_filename = [base_filename '_results.xml'];
    log_filename = [base_filename '_results.log'];
    
    % Write XML file
    write_scale_xml(xml_filename, scale_factors, measurements, marker_errors, ...
                    generic_measurements, subject_measurements, trial_name, subject_name);
    
    % Write LOG file
    write_scale_log(log_filename, scale_factors, measurements, marker_errors, ...
                    generic_measurements, subject_measurements, trial_name, subject_name);
    
    fprintf('Scale tool results written to:\n');
    fprintf('  XML: %s\n', xml_filename);
    fprintf('  LOG: %s\n', log_filename);
end

% Function to write XML file with Scale data
function write_scale_xml(filename, scale_factors, measurements, marker_errors, ...
                         generic_measurements, subject_measurements, trial_name, subject_name)
    
    try
        % Create XML document
        docNode = com.mathworks.xml.XMLUtils.createDocument('ScaleAnalysis');
        rootElement = docNode.getDocumentElement;
        
        % Add metadata
        metadata = docNode.createElement('Metadata');
        rootElement.appendChild(metadata);
        
        % Subject/Trial information
        trialInfo = docNode.createElement('TrialInfo');
        trialInfo.setAttribute('trial_name', trial_name);
        trialInfo.setAttribute('subject_name', subject_name);
        trialInfo.setAttribute('timestamp', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        metadata.appendChild(trialInfo);
        
        % Scale factors section
        if ~isempty(fieldnames(scale_factors))
            scaleFactorsElement = docNode.createElement('ScaleFactors');
            rootElement.appendChild(scaleFactorsElement);
            
            body_names = fieldnames(scale_factors);
            for i = 1:length(body_names)
                bodyElement = docNode.createElement('Body');
                bodyElement.setAttribute('name', body_names{i});
                bodyElement.setAttribute('scale_X', sprintf('%.6f', scale_factors.(body_names{i}).X));
                bodyElement.setAttribute('scale_Y', sprintf('%.6f', scale_factors.(body_names{i}).Y));
                bodyElement.setAttribute('scale_Z', sprintf('%.6f', scale_factors.(body_names{i}).Z));
                scaleFactorsElement.appendChild(bodyElement);
            end
        end
        
        % Measurements section
        if isfield(measurements, 'model') || isfield(measurements, 'subject')
            measurementsElement = docNode.createElement('Measurements');
            rootElement.appendChild(measurementsElement);
            
            % Model measurements
            if isfield(measurements, 'model')
                modelElement = docNode.createElement('ModelMeasurements');
                measurementsElement.appendChild(modelElement);
                
                model_names = fieldnames(measurements.model);
                for i = 1:length(model_names)
                    measElement = docNode.createElement('Measurement');
                    measElement.setAttribute('name', model_names{i});
                    measElement.setAttribute('value', sprintf('%.6f', measurements.model.(model_names{i})));
                    modelElement.appendChild(measElement);
                end
            end
            
            % Subject measurements
            if isfield(measurements, 'subject')
                subjectElement = docNode.createElement('SubjectMeasurements');
                measurementsElement.appendChild(subjectElement);
                
                subject_names = fieldnames(measurements.subject);
                for i = 1:length(subject_names)
                    measElement = docNode.createElement('Measurement');
                    measElement.setAttribute('name', subject_names{i});
                    measElement.setAttribute('value', sprintf('%.6f', measurements.subject.(subject_names{i})));
                    subjectElement.appendChild(measElement);
                end
            end
        end
        
        % Marker errors section
        if ~isempty(fieldnames(marker_errors))
            errorsElement = docNode.createElement('MarkerErrors');
            rootElement.appendChild(errorsElement);
            
            error_names = fieldnames(marker_errors);
            for i = 1:length(error_names)
                errorElement = docNode.createElement('Error');
                errorElement.setAttribute('marker_or_type', error_names{i});
                errorElement.setAttribute('value', sprintf('%.6f', marker_errors.(error_names{i})));
                errorsElement.appendChild(errorElement);
            end
        end
        
        % Write XML to file
        xmlwrite(filename, docNode);
        
    catch ME
        fprintf('Error writing Scale XML file: %s\n', ME.message);
        % Fallback to simple XML writing
        write_simple_scale_xml(filename, scale_factors, measurements, marker_errors, ...
                               trial_name, subject_name);
    end
end

% Fallback function for simple XML writing
function write_simple_scale_xml(filename, scale_factors, measurements, marker_errors, ...
                                trial_name, subject_name)
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not create Scale XML file: %s', filename);
    end
    
    % Write XML header
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<ScaleAnalysis>\n');
    fprintf(fid, '  <Metadata>\n');
    fprintf(fid, '    <TrialInfo trial_name="%s" subject_name="%s" timestamp="%s"/>\n', ...
            trial_name, subject_name, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '  </Metadata>\n');
    
    % Scale factors
    if ~isempty(fieldnames(scale_factors))
        fprintf(fid, '  <ScaleFactors>\n');
        body_names = fieldnames(scale_factors);
        for i = 1:length(body_names)
            fprintf(fid, '    <Body name="%s" scale_X="%.6f" scale_Y="%.6f" scale_Z="%.6f"/>\n', ...
                    body_names{i}, scale_factors.(body_names{i}).X, ...
                    scale_factors.(body_names{i}).Y, scale_factors.(body_names{i}).Z);
        end
        fprintf(fid, '  </ScaleFactors>\n');
    end
    
    % Measurements
    if isfield(measurements, 'model') || isfield(measurements, 'subject')
        fprintf(fid, '  <Measurements>\n');
        
        if isfield(measurements, 'model')
            fprintf(fid, '    <ModelMeasurements>\n');
            model_names = fieldnames(measurements.model);
            for i = 1:length(model_names)
                fprintf(fid, '      <Measurement name="%s" value="%.6f"/>\n', ...
                        model_names{i}, measurements.model.(model_names{i}));
            end
            fprintf(fid, '    </ModelMeasurements>\n');
        end
        
        if isfield(measurements, 'subject')
            fprintf(fid, '    <SubjectMeasurements>\n');
            subject_names = fieldnames(measurements.subject);
            for i = 1:length(subject_names)
                fprintf(fid, '      <Measurement name="%s" value="%.6f"/>\n', ...
                        subject_names{i}, measurements.subject.(subject_names{i}));
            end
            fprintf(fid, '    </SubjectMeasurements>\n');
        end
        
        fprintf(fid, '  </Measurements>\n');
    end
    
    % Marker errors
    if ~isempty(fieldnames(marker_errors))
        fprintf(fid, '  <MarkerErrors>\n');
        error_names = fieldnames(marker_errors);
        for i = 1:length(error_names)
            fprintf(fid, '    <Error marker_or_type="%s" value="%.6f"/>\n', ...
                    error_names{i}, marker_errors.(error_names{i}));
        end
        fprintf(fid, '  </MarkerErrors>\n');
    end
    
    fprintf(fid, '</ScaleAnalysis>\n');
    fclose(fid);
end

% Function to write LOG file with Scale data
function write_scale_log(filename, scale_factors, measurements, marker_errors, ...
                         generic_measurements, subject_measurements, trial_name, subject_name)
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not create Scale LOG file: %s', filename);
    end
    
    % Write header
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'OPENSIM SCALE TOOL ANALYSIS LOG\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Trial Name: %s\n', trial_name);
    fprintf(fid, 'Subject Name: %s\n', subject_name);
    fprintf(fid, 'Analysis Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '================================================================================\n\n');
    
    % Scale factors section
    if ~isempty(fieldnames(scale_factors))
        fprintf(fid, 'SCALE FACTORS\n');
        fprintf(fid, '================================================================================\n');
        body_names = fieldnames(scale_factors);
        fprintf(fid, '%-20s %10s %10s %10s\n', 'Body Segment', 'X Factor', 'Y Factor', 'Z Factor');
        fprintf(fid, '%s\n', repmat('-', 1, 60));
        
        for i = 1:length(body_names)
            fprintf(fid, '%-20s %10.6f %10.6f %10.6f\n', ...
                    body_names{i}, scale_factors.(body_names{i}).X, ...
                    scale_factors.(body_names{i}).Y, scale_factors.(body_names{i}).Z);
        end
        fprintf(fid, '\n');
    end
    
    % Measurements section
    if isfield(measurements, 'model') || isfield(measurements, 'subject')
        fprintf(fid, 'MEASUREMENTS\n');
        fprintf(fid, '================================================================================\n');
        
        if isfield(measurements, 'model')
            fprintf(fid, 'Model (Generic) Measurements:\n');
            model_names = fieldnames(measurements.model);
            for i = 1:length(model_names)
                fprintf(fid, '  %-30s: %10.6f\n', model_names{i}, measurements.model.(model_names{i}));
            end
            fprintf(fid, '\n');
        end
        
        if isfield(measurements, 'subject')
            fprintf(fid, 'Subject (Experimental) Measurements:\n');
            subject_names = fieldnames(measurements.subject);
            for i = 1:length(subject_names)
                fprintf(fid, '  %-30s: %10.6f\n', subject_names{i}, measurements.subject.(subject_names{i}));
            end
            fprintf(fid, '\n');
        end
        
        % Measurement ratios (if both model and subject measurements exist)
        if isfield(measurements, 'model') && isfield(measurements, 'subject')
            fprintf(fid, 'Measurement Ratios (Subject/Model):\n');
            model_names = fieldnames(measurements.model);
            subject_names = fieldnames(measurements.subject);
            common_names = intersect(model_names, subject_names);
            
            for i = 1:length(common_names)
                ratio = measurements.subject.(common_names{i}) / measurements.model.(common_names{i});
                fprintf(fid, '  %-30s: %10.6f\n', common_names{i}, ratio);
            end
            fprintf(fid, '\n');
        end
    end
    
    % Marker errors section
    if ~isempty(fieldnames(marker_errors))
        fprintf(fid, 'MARKER ERRORS (meters)\n');
        fprintf(fid, '================================================================================\n');
        error_names = fieldnames(marker_errors);
        
        for i = 1:length(error_names)
            fprintf(fid, '%-20s: %12.8f\n', error_names{i}, marker_errors.(error_names{i}));
        end
        fprintf(fid, '\n');
    end
    
    % Summary statistics
    if ~isempty(fieldnames(scale_factors))
        fprintf(fid, 'SUMMARY STATISTICS\n');
        fprintf(fid, '================================================================================\n');
        
        % Calculate overall scale factor statistics
        all_x_factors = [];
        all_y_factors = [];
        all_z_factors = [];
        body_names = fieldnames(scale_factors);
        
        for i = 1:length(body_names)
            all_x_factors(end+1) = scale_factors.(body_names{i}).X;
            all_y_factors(end+1) = scale_factors.(body_names{i}).Y;
            all_z_factors(end+1) = scale_factors.(body_names{i}).Z;
        end
        
        fprintf(fid, 'Scale Factor Statistics:\n');
        fprintf(fid, '  X Direction - Mean: %.4f, Std: %.4f, Range: %.4f to %.4f\n', ...
                mean(all_x_factors), std(all_x_factors), min(all_x_factors), max(all_x_factors));
        fprintf(fid, '  Y Direction - Mean: %.4f, Std: %.4f, Range: %.4f to %.4f\n', ...
                mean(all_y_factors), std(all_y_factors), min(all_y_factors), max(all_y_factors));
        fprintf(fid, '  Z Direction - Mean: %.4f, Std: %.4f, Range: %.4f to %.4f\n', ...
                mean(all_z_factors), std(all_z_factors), min(all_z_factors), max(all_z_factors));
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'END OF SCALE ANALYSIS LOG\n');
    fprintf(fid, '================================================================================\n');
    
    fclose(fid);
end

%% Segment the trials into gait cycles.
function [dataset, figures] = extract_gait_cycles(dataset, config, figures)

    % Visualization
    if config.visualize.extract_gait_cycles && ...
       config.visualize.force_plate_data
         fprintf([ ...
            '\nReady for extracting the gait cycles?' ... 
            '\nJust select the GRF figure and press the "n" key\n']);       
        figure(figures.grf.fig);
        while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
            waitforbuttonpress;
        end
        set(figures.grf.fig, 'CurrentCharacter', 'x');
    end

    fprintf("\nExtracting all gait cycles ...");

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            % Filter out the conditions: maximum contraction and static
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'conmc') && ~contains(x, 'const'), conditions);
            conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                for idx_tri = 1:length(trials)

                    name_trial = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                    data_mocap = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;
                    data_emg   = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;

                    fprintf('%s\n', name_trial);

                    grf_filt = data_mocap.fp_data.GRF_data;
                    grf_raw = data_mocap.fp_data.GRF_data_raw;

                    % dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data_raw(idx_grf).F

                    % Extract the gait cycles based on force plates data.
                    % Extract the indexes which store force plate values higher than a
                    % threshold.
                    stance_grf_1_idxs        = find(grf_filt(1).F(:,3) > config.threshold_grf);
                    stance_grf_1_starts_ends = [1, diff(stance_grf_1_idxs.') == 1];
                    stance_grf_1_starts      = stance_grf_1_idxs(find([stance_grf_1_starts_ends, 0] & ~[0, stance_grf_1_starts_ends]));
                    stance_grf_1_ends        = stance_grf_1_idxs(find([0, stance_grf_1_starts_ends] & ~[stance_grf_1_starts_ends, 0]) - 1);

                    swing_grf_1_idxs          = find(grf_filt(1).F(:,3) < config.threshold_grf);
                    swing_grf_1_starts_ends   = [1, diff(  swing_grf_1_idxs.') == 1];
                    swing_grf_1_starts        = swing_grf_1_idxs(find([swing_grf_1_starts_ends, 0] & ~[0, swing_grf_1_starts_ends]));
                    swing_grf_1_ends          = swing_grf_1_idxs(find([0, swing_grf_1_starts_ends] & ~[swing_grf_1_starts_ends, 0]) - 1);

                    stance_grf_2_idxs        = find(grf_filt(2).F(:,3) > config.threshold_grf);
                    stance_grf_2_starts_ends = [1, diff(stance_grf_2_idxs.') == 1];
                    stance_grf_2_starts      = stance_grf_2_idxs(find([stance_grf_2_starts_ends, 0] & ~[0, stance_grf_2_starts_ends]));
                    stance_grf_2_ends        = stance_grf_2_idxs(find([0, stance_grf_2_starts_ends] & ~[stance_grf_2_starts_ends, 0]) - 1);

                    % Remove any incomplete stances at trial's start/end
                    if isempty(stance_grf_1_starts) || isempty(stance_grf_1_ends)
                        continue
                    end

                    lengths = [length(stance_grf_1_starts), length(stance_grf_1_ends)];
                    if ~all(diff(lengths) == 0)
                        disp(name_trial);
                        disp(lengths);
                        fprintf(' -> different numbers of stance starts and ends -> not good ! ! !\n');
                    else
                        fprintf(' -> same numbers of stance starts and ends -> all good.\n');
                    end

                    % Zero force plate data during swing phases
                    for idx_plate = 1:length(grf_raw)
                        if idx_plate == 1
                            % Zero left foot (plate 1) data between stance phases
                            if length(stance_grf_1_ends) > 1 && length(stance_grf_1_starts) > 1
                                for idx_cyc = 1:length(stance_grf_1_ends) - 1
                                    % Only zero data between end of current stance and start of next stance
                                    if stance_grf_1_starts(idx_cyc+1) > stance_grf_1_ends(idx_cyc)
                                        grf_raw(1).F(stance_grf_1_ends(idx_cyc):stance_grf_1_starts(idx_cyc+1), :) = 0;
                                        grf_raw(1).M(stance_grf_1_ends(idx_cyc):stance_grf_1_starts(idx_cyc+1), :) = 0;
                                        grf_raw(1).P(stance_grf_1_ends(idx_cyc):stance_grf_1_starts(idx_cyc+1), :) = 0;
                                        % fprintf('Zeroed force plate 1 data between frames %d and %d\n', stance_grf_1_ends(idx_cyc), stance_grf_1_starts(idx_cyc+1));
                                    end
                                end
                            end
                        else
                            % Zero right foot (plate 2) data between stance phases
                            if length(stance_grf_2_ends) > 1 && length(stance_grf_2_starts) > 1
                                for idx_cyc = 1:length(stance_grf_2_ends) - 1
                                    % Only zero data between end of current stance and start of next stance
                                    if stance_grf_2_starts(idx_cyc+1) > stance_grf_2_ends(idx_cyc)
                                        grf_raw(2).F(stance_grf_2_ends(idx_cyc):stance_grf_2_starts(idx_cyc+1), :) = 0;
                                        grf_raw(2).M(stance_grf_2_ends(idx_cyc):stance_grf_2_starts(idx_cyc+1), :) = 0;
                                        grf_raw(2).P(stance_grf_2_ends(idx_cyc):stance_grf_2_starts(idx_cyc+1), :) = 0;
                                        % fprintf('Zeroed force plate 2 data between frames %d and %d\n', stance_grf_2_ends(idx_cyc), stance_grf_2_starts(idx_cyc+1));
                                    end
                                end
                            end
                        end
                        
                        % Additionally zero all values below threshold
                        below_threshold = grf_raw(idx_plate).F(:,3) < config.threshold_grf;
                        grf_raw(idx_plate).F(below_threshold, :) = 0;
                        grf_raw(idx_plate).M(below_threshold, :) = 0;
                        grf_raw(idx_plate).P(below_threshold, :) = 0;
                        fprintf('Zeroed %d frames of force plate %d data below threshold\n', ...
                            sum(below_threshold), idx_plate);
                    end
                    
                    % Leave out stances shorter than minimum length
                    gait_cycles_short = [];
                    gait_cycles_long = [];
                    for idx_cyc = 1:length(stance_grf_1_starts)-1
                        % Skip cycles that were identified as too short or too long
                        length_stride = stance_grf_1_starts(idx_cyc+1) - stance_grf_1_starts(idx_cyc);

                        if length_stride < config.length_min_gait_cycle
                            fprintf('-> stride %d duration %d, shorter than %d -> excluded\n', ...
                                idx_cyc, length_stride, config.length_min_gait_cycle);
                            gait_cycles_short = [gait_cycles_short, idx_cyc];

                        elseif length_stride > config.length_max_gait_cycle
                            fprintf('-> stride duration too long\n', ...
                                idx_cyc, length_stride, config.length_max_gait_cycle);
                            gait_cycles_long = [gait_cycles_long, idx_cyc];
                        end
                     %end

                        if ismember(idx_cyc, gait_cycles_short) || ismember(idx_cyc, gait_cycles_long)
                            continue;
                        end

                        gait_cycle = ['cyc', num2str(idx_cyc)];
                        
                        % Start at the end of current left stance (left foot lift-off)
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).start = stance_grf_1_starts(idx_cyc);
                        
                        if idx_cyc+1 <= length(stance_grf_1_starts)
                            % End at the next left stance end (next left foot lift-off)
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).end = stance_grf_1_starts(idx_cyc+1);
                        else
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).end = length(grf_filt(1).F);
                        end
                        
                        cycle_start = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).start;
                        cycle_end = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).end;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).length = cycle_end - cycle_start + 1;
                    end
                    
                    % Ensure that we start with a (even incomplete) stance of
                    % the LEFT foot F(1).
                    % -> leave out the 1st stance.
                    if swing_grf_1_starts(1) < stance_grf_1_starts(1)
                        swing_grf_1_starts = swing_grf_1_starts(2:end);
                        swing_grf_1_ends   = swing_grf_1_ends(2:end);
                    end

                    % Ensure that each stance is followed by a (even incomplete)
                    % swing.
                    % -> leave out the last stance if it has not followed by
                    % (even an incomplete) swing.
                    if length(stance_grf_1_starts) > length(swing_grf_1_starts)
                        stance_grf_1_starts = stance_grf_1_starts(1:end - 1);
                        stance_grf_1_ends = stance_grf_1_ends(1:end - 1);
                    end

                    lengths = [ ...
                        length(stance_grf_1_starts), length(stance_grf_1_ends), ...
                        length(swing_grf_1_starts),  length(swing_grf_1_ends)];

                    if ~all(diff(lengths) == 0)
                        disp(name_trial);
                        disp(lengths);
                        fprintf(' -> different numbers of stances and swings -> not good ! ! !\n');
                    else
                        fprintf(' -> same numbers of stances and swings -> all good.\n');
                    end


                    stance_grf_1_starts(gait_cycles_short) = [];
                    stance_grf_1_ends(gait_cycles_short)   = [];

                    swing_grf_1_starts(gait_cycles_short)  = [];
                    swing_grf_1_ends(gait_cycles_short)    = [];

                    fprintf(" -> %3d total gait cycles\n", length(stance_grf_1_starts));

                    % Zero force plate 1 data between stance phases, including torques and COP
                    % for idx_cyc = 1:length(stance_grf_1_ends) - 1
                        % Zero all force plate 1 data between the end of stance and beginning of next stance
                        %grf_raw(1).F(stance_grf_1_ends(idx_cyc):stance_grf_1_starts(idx_cyc+1), :) = 0;
                        %grf_raw(1).M(stance_grf_1_ends(idx_cyc):stance_grf_1_starts(idx_cyc+1), :) = 0;
                        %grf_raw(1).P(stance_grf_1_ends(idx_cyc):stance_grf_1_starts(idx_cyc+1), :) = 0;
                        %fprintf('Zeroed force plate 1 data between frames %d and %d\n', stance_grf_1_ends(idx_cyc), stance_grf_1_starts(idx_cyc+1));
                    %end


                    % Zero all force plate values below threshold
                    for idx_plate = 1:length(grf_raw)

                        % Now apply the 4th order 20 Hz low pass filter to the zeroed raw data
                        frequency_force_plate = data_mocap.fp_data.Info.frequency;
                        normalizedCutoff = config.force_plate_data.cutoffFrequency / (0.5 * frequency_force_plate);
                        [b, a]           = butter(config.force_plate_data.order, normalizedCutoff, 'low');

                        below_threshold = grf_filt(idx_plate).F(:,3) < config.threshold_grf;
                        grf_raw(idx_plate).F(below_threshold, :) = 0;
                        grf_raw(idx_plate).M(below_threshold, :) = 0;
                        grf_raw(idx_plate).P(below_threshold, :) = 0;
                        fprintf('Zeroed %d frames of force plate %d data below threshold\n', ...
                            sum(below_threshold), idx_plate);

                        % Apply filter to zeroed raw data
                        grf_filt_F = filtfilt(b, a, grf_raw(idx_plate).F);
                        grf_filt_P = filtfilt(b, a, grf_raw(idx_plate).P);
                        grf_filt_M = filtfilt(b, a, grf_raw(idx_plate).M);
                        
                        % Update the data in the main data structure
                        data_mocap.fp_data.GRF_data(idx_plate).F = grf_filt_F;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_plate).F = grf_filt_F;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_plate).P = grf_filt_P;
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_plate).M = grf_filt_M;

                    end

                    % Validate data structure before calling btk_c3d2trc
                    for idx_grf = 1:length(grf_filt)
                        assert(size(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).F, 2) == 3, 'Force data should have 3 columns');
                        assert(size(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).P, 2) == 3, 'COP data should have 3 columns');
                        assert(size(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data(idx_grf).M, 2) == 3, 'Moment data should have 3 columns');
                    end

                    % Zero all force plate positional & torqu values below threshold
                    %for idx_plate = 1:length(data_mocap.fp_data.GRF_data)
                    %    below_threshold = data_mocap.fp_data.GRF_data(idx_plate).F(:,3) < config.threshold_grf;
                    %    data_mocap.fp_data.GRF_data(idx_plate).F(below_threshold, :) = 0;
                    %    fprintf('Zeroed %d frames of force plate %d data below threshold\n', ...
                    %        sum(below_threshold), idx_plate);
                    %end

                    % Now for plate 2 (right foot), ensure smooth transitions at zero boundaries
                    %for idx_cyc = 1:length(stance_grf_1_ends) - 1
                    %    % Only process cycles that meet length requirement
                    %    cycle_length = stance_grf_1_starts(idx_cyc+1) - stance_grf_1_ends(idx_cyc);
                    %    if cycle_length >= config.length_min_gait_cycle
                    %       % Process plate 2 (right foot)
                    %        plate_idx = 2;
                    %        
                    %       % Find non-zero regions for plate 2
                    %        non_zero_frames = find(data_mocap.fp_data.GRF_data(plate_idx).F(:,3) > 0);
                    %       
                    %       if ~isempty(non_zero_frames)
                    %            % Process transitions for each dimension
                    %            for dim = 1:3
                    %                % Find transitions from zero to non-zero and vice versa
                    %                transitions = diff([0; (data_mocap.fp_data.GRF_data(plate_idx).F(:,dim) > 0); 0]);
                    %                rise_frames = find(transitions == 1);    % Frame where force becomes non-zero
                    %                fall_frames = find(transitions == -1) - 1; % Frame where force becomes zero
                                    
                                    % Process each rise transition (zero to non-zero)
                               %     for i = 1:length(rise_frames)
                               %         start_smooth = max(1, rise_frames(i) - 10); % Start 10 frames before rise or at beginning
                               %         end_smooth = min(rise_frames(i) + 10, length(data_mocap.fp_data.GRF_data(plate_idx).F)); % End 10 frames after or at end
                                        
                                        % Extract region to smooth
                               %         region = start_smooth:end_smooth;
                               %         data_region = data_mocap.fp_data.GRF_data(plate_idx).F(region, dim);
                                        
                                        % Only smooth if we have enough points
                               %         if length(region) > 5
                               %             % Create window size (must be odd)
                               %             window_size = min(21, length(region));
                                %            if mod(window_size, 2) == 0
                                %                window_size = window_size - 1;
                                %            end
                                            
                                            % Apply Savitzky-Golay smoothing
                                %            smoothed_region = sgolayfilt(data_region, 3, window_size);
                                            
                                            % Update the data
                                %            data_mocap.fp_data.GRF_data(plate_idx).F(region, dim) = smoothed_region;
                                            
                                %            fprintf('Smoothed rise transition for plate 2, dimension %d at frame %d\n', ...
                                %                dim, rise_frames(i));
                                %        end
                                %    end
                                    
                                    % Process each fall transition (non-zero to zero)
                                  %  for i = 1:length(fall_frames)
                                 %       start_smooth = max(1, fall_frames(i) - 10); % Start 10 frames before fall or at beginning
                                %        end_smooth = min(fall_frames(i) + 10, length(data_mocap.fp_data.GRF_data(plate_idx).F)); % End 10 frames after or at end
                                        
                                        % Extract region to smooth
                               %         region = start_smooth:end_smooth;
                              %          data_region = data_mocap.fp_data.GRF_data(plate_idx).F(region, dim);
                                        
                                        % Only smooth if we have enough points
                             %           if length(region) > 5
                                            % Create window size (must be odd)
                            %                window_size = min(21, length(region));
                           %                 if mod(window_size, 2) == 0
                          %                      window_size = window_size - 1;
                         %                   end
                                            
                                            % Apply Savitzky-Golay smoothing
                        %                    smoothed_region = sgolayfilt(data_region, 3, window_size);
                                            
                                            % Update the data
                        %                    data_mocap.fp_data.GRF_data(plate_idx).F(region, dim) = smoothed_region;
                        %                    
                        %                    fprintf('Smoothed fall transition for plate 2, dimension %d at frame %d\n', ...
                        %                        dim, fall_frames(i));
                        %                end
                      %%              end
                     %           end
                     %       end
                    %    end
                    %end

                    % Update the dataset with the modified force plate data
                    %dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap.fp_data.GRF_data = data_mocap.fp_data.GRF_data;

                    % Find the maximum number of cycles we can extract
                    max_cycles = min(length(stance_grf_1_starts)-1, length(stance_grf_2_starts)-1);

                    for idx_cyc = 1:max_cycles
                        gait_cycle = ['cyc', num2str(idx_cyc)];
                        
                        % Skip cycles that were identified as too short or too long
                        if ismember(idx_cyc, gait_cycles_short) || ismember(idx_cyc, gait_cycles_long)
                            continue;
                        end
                        
                        % Start at the start/strike of current left stance
                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).start = stance_grf_1_starts(idx_cyc);
                        
                        % Check if we can safely access the next stance index
                        if idx_cyc+1 <= length(stance_grf_1_starts)
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).end = stance_grf_1_starts(idx_cyc+1);
                        else
                            % Use the last available frame if we cant access the next stance
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).end = length(grf_filt(1).F);
                            fprintf('Warning: Using last frame as end point for cycle %d\n', idx_cyc);
                        end
                        
                        cycle_start = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).start;
                        cycle_end = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(gait_cycle).end;
                
                        if config.visualize.extract_gait_cycles && ...
                           config.visualize.emg_data

                            name_muscles = config.emg_data.name_muscles;
                            for idx_mus = 1:length(name_muscles)

                                hold(figures.emg.figs{idx_mus}.activation, "on");
                                plot(figures.emg.figs{idx_mus}.activation, ...
                                    data_emg.clean.(name_muscles{idx_mus})(cycle_start:cycle_end), ...
                                    'DisplayName', sprintf("%s %2d rect.+filt.+down.", trials{idx_tri}, idx_cyc));
                                hold(figures.emg.figs{idx_mus}.activation, "off");

                                xlim(figures.emg.figs{idx_mus}.activation, [0., +900.0]);
                                % ylim(figures.emg.figs{idx_mus}.activation, ?);

                                title_emg = sprintf("%s, %s, cycle %2d of total %d cycles", ...
                                    name_trial, name_muscles{idx_mus}, idx_cyc, length(stance_grf_1_starts)-1);
                                title(figures.emg.figs{idx_mus}.activation, title_emg);

                                drawnow;
                            end
                        end

                        if config.visualize.extract_gait_cycles && ...
                           config.visualize.force_plate_data

                            % Use the same color for data describing the right side.
                            color_right_side = rand(1, 3);

                            hold(figures.grf.forces, "on");
                            for idx_grf = 1:length(data_mocap.fp_data.GRF_data)
                                if idx_grf == 2
                                    plot(figures.grf.forces, ...
                                        data_mocap.fp_data.GRF_data(2).F(cycle_start:cycle_end, 3), ...
                                        'LineWidth', 2,'LineStyle', '--', ...
                                        'DisplayName', sprintf("%s %2d GRF2 (right) filt.", trials{idx_tri}, idx_cyc), 'Color', color_right_side);
                                else
                                    plot(figures.grf.forces, ...
                                        data_mocap.fp_data.GRF_data(1).F(cycle_start:cycle_end, 3), ...
                                        'DisplayName', sprintf("%s %2d GRF1 (left) filt.", trials{idx_tri}, idx_cyc));
                                end
                            end

                            hold(figures.grf.forces, "off");

                            xlim(figures.grf.forces, [  0.0,  +900.0]);
                            ylim(figures.grf.forces, [-50.0, +2500.0]);

                            title_grf = sprintf( ...
                                "%s, GRFs, cycle %2d of %2d total cycles.", name_trial, idx_cyc, length(stance_grf_1_starts)-1);
                            title(figures.grf.forces, title_grf);

                            drawnow;
                        end

                        if config.visualize.extract_gait_cycles && ...
                           config.visualize.marker_data

                            markers = fieldnames(data_mocap.marker_data.Markers);
                            for idx_mar = 1:min(length(figures.mocap.figs), length(markers))

                                hold(figures.mocap.figs{idx_mar}.coordinates, "on");
                                name_coords = {'x', 'y', 'z'};
                                for idx_coord = 1:length(name_coords)
                                    legend_coord = sprintf("%s %2d %s %s filt.", ...
                                        trials{idx_tri}, idx_cyc, markers{idx_mar}, name_coords{idx_coord});
                                    plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                        data_mocap.marker_data.Markers.(markers{idx_mar})(cycle_start:cycle_end, idx_coord), ...
                                        'DisplayName', legend_coord);
                                end
                                hold(figures.mocap.figs{idx_mar}.coordinates, "off");

                                xlim(figures.mocap.figs{idx_mar}.coordinates, [0.0,  +900.0]);
                                % ylim(figures.mocap.figs{idx_mar}.coordinates, [?, ?]);

                                title_marker = sprintf("%s, %s, cycle %2d of total %d", ...
                                    name_trial, markers{idx_mar}, idx_cyc, length(stance_grf_1_starts)-1);
                                title(figures.mocap.figs{idx_mar}.coordinates, title_marker);

                                drawnow;
                            end
                        end

                        % Semi-interactive visualization
                        if config.visualize.extract_gait_cycles && ...
                           (config.visualize.subject || ...
                            config.visualize.session || ...
                            config.visualize.condition || ...
                            config.visualize.trial || ...
                            config.visualize.gait_cycle)

                            pause(config.visualize.pause)
                        end

                        % Visualize each gait cycles separately.
                        if config.visualize.gait_cycle
                            clear_figures(config, figures);
                        end

                    end % loop through gait cycles

                    % Visualize the gait cycles overlapped.
                    if config.visualize.trial
                        clear_figures(config, figures);
                    end

                end  % loops through trials

                % Visualize the gait conditions overlapped.
                if config.visualize.condition
                    clear_figures(config, figures);
                end

            end % loops through conditions
        end % loops through sessions
    end % loops through subjects
end

function analyze_all_gait_cycles(dataset, config)
    % Quick analysis of gait cycles across entire dataset
    
    fprintf('\n=== GAIT CYCLE ANALYSIS SUMMARY ===\n');
    
    subjects = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects = subjects(idxs_subjects);
    
    total_dataset_cycles = 0;
    
    for idx_sub = 1:length(subjects)
        fprintf('\nSubject: %s\n', subjects{idx_sub});
        subject_total = 0;
        
        sessions = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions = sessions(idxs_sessions);
        
        for idx_ses = 1:length(sessions)
            fprintf('  Session: %s\n', sessions{idx_ses});
            session_total = 0;
            
            conditions = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'conmc') && ~contains(x, 'const'), conditions);
            conditions = conditions(idxs_conditions);
            
            for idx_con = 1:length(conditions)
                fprintf('    Condition: %s\n', conditions{idx_con});
                condition_total = 0;
                
                trials = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials = trials(idxs_trials);
                
                for idx_tri = 1:length(trials)
                    data_mocap = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;
                    name_trial = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                    
                    % Quick count for this trial
                    trial_cycles = count_total_gait_cycles(data_mocap.fp_data.GRF_data, config.threshold_grf, config.length_min_gait_cycle);
                    
                    fprintf('      %s: %d cycles\n', trials{idx_tri}, trial_cycles);
                    
                    condition_total = condition_total + trial_cycles;
                end
                
                fprintf('    Condition Total: %d cycles\n', condition_total);
                session_total = session_total + condition_total;
            end
            
            fprintf('  Session Total: %d cycles\n', session_total);
            subject_total = subject_total + session_total;
        end
        
        fprintf('Subject Total: %d cycles\n', subject_total);
        total_dataset_cycles = total_dataset_cycles + subject_total;
    end
    
    fprintf('\n=== DATASET SUMMARY ===\n');
    fprintf('Total gait cycles in entire dataset: %d\n', total_dataset_cycles);
    fprintf('Cycles to be exported (first %d per trial): limited by config\n', config.export_only_first_n_cycles);
    fprintf('============================\n\n');
end

function total_cycles = count_total_gait_cycles(grf_data, threshold_grf, min_cycle_length)
    % Quick function to count total gait cycles in the entire trial
    
    % Extract stance phases for left foot (force plate 1)
    stance_grf_1_idxs = find(grf_data(1).F(:,3) > threshold_grf);
    
    if isempty(stance_grf_1_idxs)
        total_cycles = 0;
        return;
    end
    
    % Find stance start and end points
    stance_starts_ends = [1, diff(stance_grf_1_idxs.') == 1];
    stance_starts = stance_grf_1_idxs(find([stance_starts_ends, 0] & ~[0, stance_starts_ends]));
    stance_ends = stance_grf_1_idxs(find([0, stance_starts_ends] & ~[stance_starts_ends, 0]) - 1);
    
    % Ensure we have matching starts and ends'
    min_length = min(length(stance_starts), length(stance_ends));
    stance_starts = stance_starts(1:min_length);
    stance_ends = stance_ends(1:min_length);
    
    % Count cycles (each cycle goes from one stance start to the next)
    total_cycles = 0;
    for i = 1:length(stance_starts)-1
        cycle_length = stance_starts(i+1) - stance_starts(i);
        if cycle_length >= min_cycle_length
            total_cycles = total_cycles + 1;
        end
    end
    
    fprintf('  -> Total valid gait cycles found: %d\n', total_cycles);
end

%% Export the gayt cycles into the OpenSim data formats.
function [dataset, figures] = compute_ik_and_id(dataset, config, figures)

    fprintf("\nCalculating session-wide maximum EMG values from all trials within each subject's session...\n");

    % Loop through each subject and session to calculate session-wide maxima
    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);

    for idx_sub = 1:length(subjects)
        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        
        for idx_ses = 1:length(sessions)
            fprintf('Processing %s %s for session-wide maxima...\n', subjects{idx_sub}, sessions{idx_ses});
            
            % Initialize session maximum values for each muscle
            name_muscles = config.emg_data.name_muscles;
            session_muscle_max = struct();
            for idx_mus = 1:length(name_muscles)
                session_muscle_max.(name_muscles{idx_mus}) = 0;
            end
            
            % Get ALL conditions for this subject/session (including MVIC, static, and dynamic)
            all_conditions = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            
            for idx_con = 1:length(all_conditions)
                condition_name = all_conditions{idx_con};
                
                % Skip conditions that don`t contain trials
                if ~isstruct(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name))
                    continue;
                end
                
                trials = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials = trials(idxs_trials);
                
                for idx_tri = 1:length(trials)
                    fprintf('  Checking %s for session max...\n', trials{idx_tri});
                    
                    % Check each muscle in this trial
                    for idx_mus = 1:length(name_muscles)
                        muscle_name = name_muscles{idx_mus};
                        
                        % Get the rectified and filtered EMG data for the entire trial
                        if isfield(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name).(trials{idx_tri}).emg, 'rect_filt') && ...
                        isfield(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name).(trials{idx_tri}).emg.rect_filt, muscle_name)
                            
                            trial_data = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name).(trials{idx_tri}).emg.rect_filt.(muscle_name);
                            trial_max = max(trial_data);
                            
                            % Update session maximum if this trial`s max is higher
                            if trial_max > session_muscle_max.(muscle_name)
                                session_muscle_max.(muscle_name) = trial_max;
                                fprintf('    New session max for %s: %f (from %s %s)\n', ...
                                    muscle_name, trial_max, condition_name, trials{idx_tri});
                            end
                        end
                    end
                end
            end
            
            % Display final session maximum values
            fprintf('Final session maximum values for %s %s:\n', subjects{idx_sub}, sessions{idx_ses});
            for idx_mus = 1:length(name_muscles)
                fprintf('  %s: %f\n', name_muscles{idx_mus}, session_muscle_max.(name_muscles{idx_mus}));
            end
            
            % Store session maxima in ALL trials within this session (including MVIC, static, and dynamic)
            for idx_con = 1:length(all_conditions)
                condition_name = all_conditions{idx_con};
                
                if ~isstruct(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name))
                    continue;
                end
                
                trials = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials = trials(idxs_trials);
                
                for idx_tri = 1:length(trials)
                    % Store session maxima in each trial
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(condition_name).(trials{idx_tri}).emg.session_maximum = session_muscle_max;
                end
            end
        end
    end

    fprintf("Session-wide maximum calculation complete.\n\n");

    %% Add the global maximum calculation section here (from the first artifact)
    fprintf("\nCalculating global maximum EMG values from all dynamic trials...\n");

    fprintf("\nCalculating global maximum EMG values from all dynamic trials...\n");

    % Initialize global maximum values for each muscle
    name_muscles = config.emg_data.name_muscles;
    global_muscle_max = struct();
    for idx_mus = 1:length(name_muscles)
        global_muscle_max.(name_muscles{idx_mus}) = 0;
    end

    % Loop through all subjects, sessions, and dynamic trials to find global maxima
    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);

    for idx_sub = 1:length(subjects)
        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        
        for idx_ses = 1:length(sessions)
            % Filter to get only DYNAMIC conditions (exclude MVIC and static)
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'conmc') && ~contains(x, 'const'), conditions);
            conditions      = conditions(idxs_conditions);
            
            for idx_con = 1:length(conditions)
                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                
                for idx_tri = 1:length(trials)
                    fprintf('Processing %s %s %s %s for global max...\n', ...
                        subjects{idx_sub}, sessions{idx_ses}, conditions{idx_con}, trials{idx_tri});
                    
                    % Check each muscle in this trial
                    for idx_mus = 1:length(name_muscles)
                        muscle_name = name_muscles{idx_mus};
                        
                        % Get the rectified and filtered EMG data for the entire trial
                        if isfield(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg, 'rect_filt') && ...
                        isfield(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.rect_filt, muscle_name)
                            
                            trial_data = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.rect_filt.(muscle_name);
                            trial_max = max(trial_data);
                            
                            % Update global maximum if this trials max is higher
                            if trial_max > global_muscle_max.(muscle_name)
                                global_muscle_max.(muscle_name) = trial_max;
                                fprintf('  New global max for %s: %f\n', muscle_name, trial_max);
                            end
                        end
                    end
                end
            end
        end
    end

    % Display final global maximum values
    fprintf('\nFinal global maximum values from dynamic dataset:\n');
    for idx_mus = 1:length(name_muscles)
        fprintf('  %s: %f\n', name_muscles{idx_mus}, global_muscle_max.(name_muscles{idx_mus}));
    end

    % Store global maxima in all dynamic trials for later use
    for idx_sub = 1:length(subjects)
        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        
        for idx_ses = 1:length(sessions)
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'conmc') && ~contains(x, 'const'), conditions);
            conditions      = conditions(idxs_conditions);
            
            for idx_con = 1:length(conditions)
                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                
                for idx_tri = 1:length(trials)
                    % Store global maxima in each trial
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.global_dynamic_maximum = global_muscle_max;
                end
            end
        end
    end

    fprintf("Global maximum calculation complete.\n\n");

    %% Now proceed with the regular compute_ik_and_id processing
    fprintf("Starting IK and ID processing...\n");

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);

    for idx_sub = 1:length(subjects)
        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        
        for idx_ses = 1:length(sessions)
            % Filter out the conditions: maximum contraction and static.
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'conmc') && ~contains(x, 'const'), conditions);
            conditions      = conditions(idxs_conditions);
            
            for idx_con = 1:length(conditions)
                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                
                for idx_tri = 1:length(trials)
                    
                    % Re-define the data_emg variable to access cleaned EMG data
                    % NOTE: No need to update dataset structure here since global maxima are already stored
                    data_emg = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;

                    cycles      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}));
                    idxs_cycles = cellfun(@(x) contains(x, 'cyc'), cycles);
                    cycles      = cycles(idxs_cycles);
                    
                    for idx_cyc = 1:length(cycles)

                        % Compute IK and export it only for the first n gait cycles.
                        if config.export_only_first_n_cycles < idx_cyc
                            fprintf("Exporting only the first %2d cycles of each trial done!\n", ...
                                config.export_only_first_n_cycles);
                            continue
                        end

                        % Check if cycle meets minimum length requirement
                        % if isfield(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}), 'length')
                        %    cycle_length = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).length;
                        %    
                        %    if cycle_length < config.length_min_gait_cycle
                        %        fprintf('Skipping analysis for %s (length: %d < minimum: %d)\n', ...
                        %            cycles{idx_cyc}, cycle_length, config.length_min_gait_cycle);
                        %       continue;
                        %   end
                        %nd

                        name_trial = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                        data       = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;

                        % Define start and end frame from the events to write the appropriate TRC
                        % and MOT files for the OpenSim simulations
                        data.Start_Frame = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).start;
                        data.End_Frame   = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).end;

                        % Change to dataset directory for intermediate files
                        cd([config.path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/']);

                        % now do conversions to TRC files using btk_c3d2trc
                        TRCFileName = [name_trial '_' sprintf('%04d', idx_cyc) '_trc.trc'];
                        GRFFileName = [name_trial '_' sprintf('%04d', idx_cyc) '_grf.sto'];
                        GRFFileXML  = [name_trial '_' sprintf('%04d', idx_cyc) '_grf.xml'];
                        GRFFileXML_adjusted  = [name_trial '_' sprintf('%04d', idx_cyc) '_grf_move.xml'];
                        

                        data_trc = btk_c3d2trc(data, 'off', TRCFileName, GRFFileName);

                        % Define output files in moco_inputs directory
                        pre_moco_dir = config.path_moco_inputs;
                        if ~exist(pre_moco_dir, 'dir')
                            mkdir(pre_moco_dir)
                        end


                        % define the standard OpenSim files for IK and EMG
                        EMGFileName = fullfile(pre_moco_dir, [name_trial '_' sprintf('%04d', idx_cyc) '_emg.sto']);
                        ModelFile    = [data_trc.Name '_SCALED_addbio_free_subt_mtp.osim']; % [data_trc.Name '_SCALED.osim'];  
                        IKTasksFile  = [config.path_model config.file_model '_IK_Tasks.xml'];
                        IKOutputFile = fullfile(pre_moco_dir, [name_trial '_' sprintf('%04d', idx_cyc) '_ik.sto']);
                        

                        % Export EMG data to STO file using c3d2emg function before running the ik tool
                        data_sto_emg = c3d2emg(data_emg, data_trc, EMGFileName, 1000);
                        

                        % Uncommenting this causes all results to be dumped in same directory as c3ds
                        % cd([config.path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/']);
                        setup_InverseKinematics( ...
                            'data',        data_trc, ...
                            'ModelFile',   ModelFile,...
                            'IKTasksFile', IKTasksFile, ...
                            'Accuracy',    0.00005, ...
                            'OutputFile',  IKOutputFile);

                        % run the ik tool from the command line
                        com = ['opensim-cmd run-tool ' data_trc.Name '_Setup_InverseKinematics.xml'];
                        [status, cmdout] = system(com);
                        
                        fprintf('\n IK Tool Output \n');
                        fprintf('%s\n', cmdout);
                        fprintf('=End IK Tool Output=\n\n');

                        [time_values, max_errors, max_error_markers, rms_errors] = parse_ik_output(cmdout);

                        % Write parsed results to output files if data successfully parsed
                        if ~isempty(time_values)
                            frame_numbers = 0:(length(time_values) - 1);
                            total_squared_errors = NaN(size(time_values));
                            write_ik_results_to_files(time_values, max_errors, max_error_markers, rms_errors, ...
                                                    frame_numbers, total_squared_errors, name_trial, idx_cyc, pre_moco_dir);
                        end

                        % load the data from MOT file and tie kinematic data to data structure
                        D = load_sto_file(IKOutputFile);
                        
                        % Load the GRF file to get the original ground force position data
                        fprintf('Loading GRF file for position correction: %s\n', GRFFileName);
                        grf_data = load_sto_file(GRFFileName);

                        % Apply low-pass filter to all coordinates except time
                        fprintf('Applying 15 Hz low-pass filter to IK coordinates...\n');

                        % Filter parameters
                        sampling_freq = 1000; % Hz (assuming 1000 Hz sampling rate from your code)
                        cutoff_freq = 15; % Hz
                        filter_order = 4;

                        % Design the filter
                        nyquist_freq = sampling_freq / 2;
                        normalized_cutoff = cutoff_freq / nyquist_freq;
                        [b, a] = butter(filter_order, normalized_cutoff, 'low');

                        % Get all field names from the IK data
                        field_names = fieldnames(D);

                        % Filter each coordinate (skip the time column)
                        for i = 1:length(field_names)
                            field_name = field_names{i};
                            
                            % Skip the time column
                            if strcmp(field_name, 'time')
                                continue;
                            end
                            
                            % Apply zero-phase low-pass filter using filtfilt
                            D.(field_name) = filtfilt(b, a, D.(field_name));
                            
                            fprintf('Filtered coordinate: %s\n', field_name);
                        end

                        fprintf('Filtering complete. Now proceeding with pelvis_tx modification...\n');

                        % Store original pelvis_tx values AFTER filtering
                        original_pelvis_tx = D.pelvis_tx;

                        % Calculate differences between original (filtered) pelvis_tx and ground force positions
                        pelvis_grf1_px_diff = [];
                        pelvis_grf2_px_diff = [];

                        % Check if the required columns exist in the GRF file
                        if isfield(grf_data, 'ground_1_force_px')
                            pelvis_grf1_px_diff = original_pelvis_tx - grf_data.ground_1_force_px;
                            fprintf('Calculated differences between filtered pelvis_tx and ground_1_force_px\n');
                        else
                            fprintf('Warning: ground_1_force_px column not found in GRF file\n');
                        end

                        if isfield(grf_data, 'ground_2_force_px')
                            pelvis_grf2_px_diff = original_pelvis_tx - grf_data.ground_2_force_px;
                            fprintf('Calculated differences between filtered pelvis_tx and ground_2_force_px\n');
                        else
                            fprintf('Warning: ground_2_force_px column not found in GRF file\n');
                        end

                        % Modify pelvis_tx to represent the distance traveled by the subject
                        if isfield(D, 'pelvis_tx')
                            % Access treadmill data for this trial
                            data_treadmill = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).treadmill;
                            
                            % Get the first right belt speed value to use as a scaling factor
                            initial_speed = data_treadmill.Rightbelt_Speed(1);
                            
                            % Get the time column from the kinematics data
                            time_column = D.time;
                            
                            % Keep the first value as is
                            modified_pelvis_tx = zeros(size(original_pelvis_tx));
                            modified_pelvis_tx(1) = original_pelvis_tx(1);
                            
                            % For each subsequent time point, calculate the cumulative distance
                            for i = 2:length(time_column)
                                % Time elapsed since last frame
                                delta_t = time_column(i) - time_column(i-1);
                                
                                % Distance traveled in this time step (speed * time)
                                distance_increment = initial_speed * delta_t;
                                
                                % Calculate the change from the original data
                                original_increment = original_pelvis_tx(i) - original_pelvis_tx(i-1);
                                
                                % Subtract the original change and add the new distance increment
                                modified_pelvis_tx(i) = modified_pelvis_tx(i-1) + (distance_increment - original_increment);
                            end
                            
                            % Update the pelvis_tx values in both D and data_trc structures
                            D.pelvis_tx = modified_pelvis_tx;
                            data_trc.Kinematics.Angles.pelvis_tx = modified_pelvis_tx;
                            
                            % Create new filename for adjusted GRF data
                            [filepath, name, ext] = fileparts(GRFFileName);
                            if isempty(filepath)
                                new_grf_filename = [name '_move' ext];
                            else
                                new_grf_filename = fullfile(filepath, [name, '_move' ext]);
                            end
                            
                            % Create a new file for modified IK data and keep the initial 
                            [filepath,name, ext] = fileparts(IKOutputFile);
                            if isempty(filepath)
                                new_ik_filename = [name '_move' ext];
                            else
                                new_ik_filename = fullfile(filepath, [name '_move' ext]);
                            end

                            % After modifying pelvis_tx values, update the GRF file
                            if ~isempty(pelvis_grf1_px_diff) || ~isempty(pelvis_grf2_px_diff)
                                
                                fprintf('Updating GRF file with adjusted over-ground position values...\n');
                                
                                % Update the ground force position columns for ALL rows
                                if ~isempty(pelvis_grf1_px_diff)
                                    grf_data.ground_1_force_px = modified_pelvis_tx - pelvis_grf1_px_diff;
                                    fprintf('Updated ground_1_force_px column for %d data points\n', length(modified_pelvis_tx));
                                end
                                
                                if ~isempty(pelvis_grf2_px_diff)
                                    grf_data.ground_2_force_px = modified_pelvis_tx - pelvis_grf2_px_diff;
                                    fprintf('Updated ground_2_force_px column for %d data points\n', length(modified_pelvis_tx));
                                end                               
                            end
                            
                            % Write the updated GRF data back to file
                            write_updated_grf_file(new_grf_filename, grf_data);
                            
                            fprintf('Successfully updated GRF file: %s\n', GRFFileName);

                            % Write the modified (filtered + pelvis_tx corrected) data back to the IKOutputFile
                            try
                                % Write the STO file directly using MATLAB functions
                                fid = fopen(new_ik_filename, 'w');
                                if fid == -1
                                    error('Could not open file for writing: %s', IKOutputFile);
                                end
                                
                                % Write header lines (typical for OpenSim STO files)
                                fprintf(fid, 'Coordinates\n');
                                fprintf(fid, 'version=1\n');
                                fprintf(fid, 'nRows=%d\n', length(D.time));
                                fprintf(fid, 'nColumns=%d\n', length(fieldnames(D)));
                                fprintf(fid, 'inDegrees=yes\n');
                                fprintf(fid, 'endheader\n');
                                
                                % Write column headers
                                headers = fieldnames(D);
                                fprintf(fid, '%s', headers{1});
                                for i = 2:length(headers)
                                    fprintf(fid, '\t%s', headers{i});
                                end
                                fprintf(fid, '\n');
                                
                                % Write data rows with 6 decimal precision
                                for row = 1:length(D.time)
                                    fprintf(fid, '%.6f', D.(headers{1})(row));
                                    for i = 2:length(headers)
                                        fprintf(fid, '\t%.6f', D.(headers{i})(row));
                                    end
                                    fprintf(fid, '\n');
                                end
                                
                                fclose(fid);
                                fprintf('Successfully wrote filtered and modified data to %s\n', IKOutputFile);
                                
                            catch e
                                % If any error occurs, display it
                                fprintf('Error writing to file: %s\n', e.message);
                                if fid ~= -1 && fid ~= 0
                                fclose(fid);
                                end
                            end
                            
                            fprintf('Filtered and modified pelvis_tx values to represent cumulative distance based on treadmill speed.\n');
                        end
                        
                        % assign forces to the specific segments using 'assign_force' function
                        % the output creates a cell array with the name of the external force
                        % (ExForce) and the body that is assigned to that force (ApBodies)
                        data_trc = assign_forces(data_trc,{'LTOE','RTOE'},{'calcn_l','calcn_r'},0.25);

                        % define the standard OpenSim files for ID.
                        IDOutputFile = fullfile(pre_moco_dir, [name_trial '_' sprintf('%04d', idx_cyc) '_id.sto']);

                        cd([config.path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/']);

                        % it is also necessary to generate an XML file containing the information
                        % about which column is which and the bodies these forces are applied to
                        % this is done with the grf2xml function
                        grf2xml( ...
                            data_trc, ...
                            'ExternalLoadNames',          data_trc.AssignForce.ExForce, ...
                            'AppliedToBodies',            {'calcn_l','calcn_r'}, ... %data_trc.AssignForce.ApBodies,...
                            'GRFFile',                    GRFFileName, ...
                            'MOTFile',                    IKOutputFile, ...
                            'LowPassFilterForKinematics', -1,...
                            'OutputFile',                 GRFFileXML);

                        grf2xml( ...
                            data_trc, ...
                            'ExternalLoadNames',          data_trc.AssignForce.ExForce, ...
                            'AppliedToBodies',            {'calcn_l','calcn_r'}, ... %data_trc.AssignForce.ApBodies,...
                            'GRFFile',                    new_grf_filename, ...
                            'MOTFile',                    new_ik_filename, ...
                            'OutputFile',                 GRFFileXML_adjusted);

                        
                        GRFFileName_output = fullfile(pre_moco_dir, [name_trial '_' sprintf('%04d', idx_cyc) '_grf.sto']);
                        [~, grf_fname, grf_fext] = fileparts(new_grf_filename)
                        GRFFileName_output_adjusted = fullfile(pre_moco_dir, [grf_fname, grf_fext]);
                        % Also write over the GRFs that were initially extracted out from 
                        copyfile(GRFFileName, GRFFileName_output);
                        copyfile(new_grf_filename, GRFFileName_output_adjusted);


                        % make setup file for inverse dynamics
                        setup_InverseDynamics(...
                            'data',                       data_trc, ...
                            'ModelFile',                  ModelFile, ...
                            'MOTFile',                    IKOutputFile,...
                            'GRFFile',                    GRFFileXML, ...
                            'LowPassFilterForKinematics', 15, ...
                            'OutputFile',                 IDOutputFile); 

                        % Call the ID tool from the command line using setup file
                        com = ['opensim-cmd run-tool ' data_trc.Name '_Setup_InverseDynamics.xml'];
                        system(com);

                        %% Adding in first RRA iteration
                        % define the standard OpenSim files for RRA
                        RRATaskFile  = [config.path_model config.file_model '_RRA_Tasks.xml'];
                        RRAForceFile = [config.path_model config.file_model '_RRA_actuators.xml'];
                        RRAOutputModelFile = [name_trial '_' sprintf('%04d', idx_cyc) '_rra.osim'];

                        cd([config.path_dataset, subjects{idx_sub}, '/', sessions{idx_ses}, '/']);

                        setup_ReduceResiduals( ...
                            'data',        data_trc, ...
                            'ModelFile',   ModelFile,...
                            'MOTFile', IKOutputFile, ...
                            'GRFFile', GRFFileXML, ...
                            'RRATaskFile', RRATaskFile, ...
                            'RRAForceFile', RRAForceFile, ...
                            'AdjustedCOMBody', 'pelvis', ...
                            'LowPassFilterFreq', 6, ...
                            'Accuracy',    0.00005, ...
                            'OutputModelFile',  RRAOutputModelFile);

                        % run the rra tool from the command line
                        com = ['opensim-cmd run-tool ' data_trc.Name '_Setup_ReduceResiduals.xml'];
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

                        ik_data_table               = org.opensim.modeling.TimeSeriesTable(IKOutputFile);
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

                            % Adaptive Savitzky-Golay filter based on data length
                            data_length = length(moving_avg);

                            if data_length >= 101
                                % Use full frame length for normal-sized data
                                smoothed_values = sgolayfilt(moving_avg, 2, 101);
                            elseif data_length >= 21
                                % Use adaptive frame length for shorter data
                                frame_length = data_length;
                                if mod(frame_length, 2) == 0
                                    frame_length = frame_length - 1; % Frame length must be odd
                                end
                                % Ensure minimum reasonable frame length
                                frame_length = max(frame_length, 21);
                                if frame_length > data_length
                                    frame_length = data_length;
                                    if mod(frame_length, 2) == 0
                                        frame_length = frame_length - 1;
                                    end
                                end
                                smoothed_values = sgolayfilt(moving_avg, 2, frame_length);
                                fprintf('Using adaptive Savitzky-Golay filter (frame length: %d) for cycle with %d data points\n', ...
                                    frame_length, data_length);
                            elseif data_length >= 5
                                % For very short data, use simple moving average or basic smoothing
                                smoothed_values = smooth(moving_avg, min(5, data_length));
                                fprintf('Warning: Cycle too short (%d frames) for Savitzky-Golay filter, using simple smoothing\n', ...
                                    data_length);
                            else
                                % For extremely short data, use as-is
                                smoothed_values = moving_avg;
                                fprintf('Warning: Cycle extremely short (%d frames), no additional smoothing applied\n', ...
                                    data_length);
                            end

                            % todo(lisca): move the plot out of the filter.
                            if false % plot_inverse_kinematics
                                figure();
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
                        org.opensim.modeling.STOFileAdapter.write(ik_data_table_smoothed, IKOutputFile);
                    end % loop through gait cycles
                end 
            end % loop through trials
        end % loop through conditions
    end % loop through session
end % loop through subjects

% Function to write adjusted GRF file
function write_updated_grf_file(filename, grf_data)
    try
        % Create backup of original file
        backup_filename = [filename '.backup'];
        if exist(filename, 'file')
            copyfile(filename, backup_filename);
            fprintf('Created backup file: %s\n', backup_filename);
        end
        
        % Open file for writing
        fid = fopen(filename, 'w');
        if fid == -1
            error('Could not open file for writing: %s', filename);
        end
        
        % Get field names (column headers)
        field_names = fieldnames(grf_data);
        num_rows = length(grf_data.(field_names{1}));
        
        % Write header information (typical for OpenSim STO files)
        fprintf(fid, 'name %s\n', filename);
        fprintf(fid, 'datacolumns %d\n', length(field_names));
        fprintf(fid, 'datarows %d\n', num_rows);
        fprintf(fid, 'range %f %f\n', grf_data.time(1), grf_data.time(end));
        fprintf(fid, 'endheader\n');
        
        % Write column headers
        for i = 1:length(field_names)
            if i == 1
                fprintf(fid, '%s', field_names{i});
            else
                fprintf(fid, '\t%s', field_names{i});
            end
        end
        fprintf(fid, '\n');
        
        % Write data rows
        for row = 1:num_rows
            for col = 1:length(field_names)
                if col == 1
                    fprintf(fid, '%.8f', grf_data.(field_names{col})(row));
                else
                    fprintf(fid, '\t%.8f', grf_data.(field_names{col})(row));
                end
            end
            fprintf(fid, '\n');
        end
        
        fclose(fid);
        fprintf('Successfully wrote updated GRF data to: %s\n', filename);
        
    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        fprintf('Error writing GRF file: %s\n', ME.message);
        rethrow(ME);
    end
end

% Function to analyze static pose kinematics
function analyze_static_pose(static_ik_data, trial_name, subject_name)
    
    fprintf('\n=== Static Pose Analysis ===\n');
    
    % Get coordinate names (excluding time)
    coord_names = fieldnames(static_ik_data);
    coord_names = coord_names(~strcmp(coord_names, 'time'));
    
    % Create analysis filename
    analysis_filename = [trial_name '_static_pose_analysis.log'];
    
    fid = fopen(analysis_filename, 'w');
    if fid == -1
        fprintf('Warning: Could not create static pose analysis file\n');
        return;
    end
    
    % Write header
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'STATIC POSE KINEMATIC ANALYSIS\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Subject: %s\n', subject_name);
    fprintf(fid, 'Trial: %s\n', trial_name);
    fprintf(fid, 'Analysis Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, 'Total Frames: %d\n', length(static_ik_data.time));
    fprintf(fid, '================================================================================\n\n');
    
    % Analyze each coordinate
    fprintf(fid, 'COORDINATE ANALYSIS\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, '%-25s %10s %10s %10s %10s %10s\n', ...
            'Coordinate', 'Mean', 'Std Dev', 'Min', 'Max', 'Range');
    fprintf(fid, '%s\n', repmat('-', 1, 85));
    
    for i = 1:length(coord_names)
        coord_data = static_ik_data.(coord_names{i});
        
        coord_mean = mean(coord_data);
        coord_std = std(coord_data);
        coord_min = min(coord_data);
        coord_max = max(coord_data);
        coord_range = coord_max - coord_min;
        
        fprintf(fid, '%-25s %10.4f %10.4f %10.4f %10.4f %10.4f\n', ...
                coord_names{i}, coord_mean, coord_std, coord_min, coord_max, coord_range);
    end
    
    fprintf(fid, '\n');
    
    % Identify coordinates with high variability (potential issues)
    fprintf(fid, 'QUALITY ASSESSMENT\n');
    fprintf(fid, '================================================================================\n');
    
    high_variability_coords = {};
    variability_threshold = 5.0; % degrees for rotational coordinates
    
    for i = 1:length(coord_names)
        coord_data = static_ik_data.(coord_names{i});
        coord_range = max(coord_data) - min(coord_data);
        
        % Different thresholds for different coordinate types
        if contains(coord_names{i}, 'tx') || contains(coord_names{i}, 'ty') || contains(coord_names{i}, 'tz')
            % Translation coordinates - threshold in meters
            threshold = 0.01; % 1 cm
        else
            % Rotational coordinates - threshold in degrees
            threshold = variability_threshold;
        end
        
        if coord_range > threshold
            high_variability_coords{end+1} = coord_names{i};
            fprintf(fid, 'High variability detected in %s: range = %.4f\n', ...
                    coord_names{i}, coord_range);
        end
    end
    
    if isempty(high_variability_coords)
        fprintf(fid, 'All coordinates show acceptable stability for static pose.\n');
    else
        fprintf(fid, '\nRecommendation: Review marker placement and/or filtering for coordinates with high variability.\n');
    end
    
    fprintf(fid, '\n');
    
    % Summary statistics
    fprintf(fid, 'SUMMARY\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Total coordinates analyzed: %d\n', length(coord_names));
    fprintf(fid, 'Coordinates with high variability: %d\n', length(high_variability_coords));
    fprintf(fid, 'Time duration: %.3f seconds\n', max(static_ik_data.time) - min(static_ik_data.time));
    fprintf(fid, 'Sampling frequency: %.1f Hz\n', length(static_ik_data.time) / (max(static_ik_data.time) - min(static_ik_data.time)));
    
    fprintf(fid, '\n================================================================================\n');
    fprintf(fid, 'END OF STATIC POSE ANALYSIS\n');
    fprintf(fid, '================================================================================\n');
    
    fclose(fid);
    
    fprintf('Static pose analysis written to: %s\n', analysis_filename);
    
    % Print summary to console
    fprintf('Static pose summary:\n');
    fprintf('  - %d coordinates analyzed\n', length(coord_names));
    fprintf('  - %d coordinates with high variability\n', length(high_variability_coords));
    if ~isempty(high_variability_coords)
        fprintf('  - Problem coordinates: %s\n', strjoin(high_variability_coords, ', '));
    end
end

% Alternative function to run static IK with custom TRC file
function static_ik_output = run_static_ik_with_custom_trc(trc_filename, model_filename, ik_tasks_filename, output_name)
    
    fprintf('Running static IK with custom TRC file: %s\n', trc_filename);
    
    % Create temporary data structure
    temp_data.Name = output_name;
    temp_data.TRC_Filename = trc_filename;
    
    % Setup and run IK
    static_ik_output = [output_name '_static_ik.sto'];
    
    setup_InverseKinematics( ...
        'data',        temp_data, ...
        'ModelFile',   model_filename,...
        'IKTasksFile', ik_tasks_filename, ...
        'Accuracy',    0.000001, ...
        'OutputFile',  static_ik_output);
    
    % Run the IK tool
    com = ['opensim-cmd run-tool ' temp_data.Name '_Setup_InverseKinematics.xml'];
    [status, cmdout] = system(com);
    
    if status == 0
        fprintf('Static IK completed successfully: %s\n', static_ik_output);
    else
        fprintf('Static IK failed with status: %d\n', status);
        fprintf('Output: %s\n', cmdout);
    end
    
    % Note: static_ik_output is returned via the function signature above
end


% Added function to parse IK output and extract errors 
function [time_values, max_errors, max_error_markers, rms_errors] = parse_ik_output(cmdout)
    % Initialize output variables
    time_values = [];
    max_errors = [];
    max_error_markers = {};
    rms_errors = [];
    total_squared_errors = [];
    frame_numbers = [];

    % Split the command output into lines
    try
        lines = strsplitlines(cmdout);
    catch
    % using regexp as fallback
        lines = regexp(cmdout, '\r?\n', 'split');
        % Remove empty last element if existent
        if ~isempty(lines) && isempty(lines{end})
            lines = lines(1:end - 1);
        end
    end

    % Loop through each line to find relevant data
    for i = 1:length(lines)
        line = strtrim(lines{i});

        % Look exclusively for lines containing time and error info, which start with '[info] Frame'
        if startsWith(line, '[info] Frame')

            % Parse the line using regex to extract all components
            pattern = '\[info\]\s+Frame\s+(\d+)\s+\(t\s+=\s+([\d\.-]+)\):\s+total\s+squared\s+error\s+=\s+([\d\.-]+),\s+marker\s+error:\s+RMS\s+=\s+([\d\.-]+),\s+max\s+=\s+([\d\.-]+)\s+\(([^)]+)\)';
            tokens = regexp(line, pattern, 'tokens');

            if ~isempty(tokens)
                % Extract the values from the tokens
                frame_number = str2double(tokens{1}{1});
                time_value = str2double(tokens{1}{2});
                total_squared_error = str2double(tokens{1}{3});
                rms_error = str2double(tokens{1}{4});
                max_error = str2double(tokens{1}{5});
                marker_name = tokens{1}{6};

                % Append values to storre in variables
                frame_numbers(end+1) = frame_number;
                time_values(end+1) = time_value; 
                max_errors(end+1) = max_error; 
                max_error_markers{end+1} = marker_name;
                rms_errors(end+1) = rms_error; 

                fprintf('Frame %d: Time = %.3f, Total Squared Error = %.6f, RMS Error = %.6f, Max Error = %.6f, Marker = %s\n', ...
                    frame_number, time_value, total_squared_error, rms_error, max_error, marker_name);
            else
                % Try a more flexible pattern in case formatting varies slightly
                alt_pattern = '\[info\]\s+Frame\s+(\d+).*t\s+=\s+([\d\.-]+).*RMS\s+=\s+([\d\.-]+).*max\s+=\s+([\d\.-]+)\s+\(([^)]+)\)';
                alt_tokens = regexp(line, alt_pattern, 'tokens');
                
                if ~isempty(alt_tokens)
                    frame_num = str2double(alt_tokens{1}{1});
                    time_val = str2double(alt_tokens{1}{2});
                    rms_err = str2double(alt_tokens{1}{3});
                    max_err = str2double(alt_tokens{1}{4});
                    marker_name = alt_tokens{1}{5};
                    
                    % Store the values (without total squared error)
                    frame_numbers(end+1) = frame_num;
                    time_values(end+1) = time_val;
                    total_squared_errors(end+1) = NaN; % Not available in this pattern
                    rms_errors(end+1) = rms_err;
                    max_errors(end+1) = max_err;
                    max_error_markers{end+1} = marker_name;
                    
                    fprintf('Frame %d: t=%.6f, RMS=%.6f, Max=%.6f (%s) [TSE not parsed]\n', ...
                        frame_num, time_val, rms_err, max_err, marker_name);
                else
                    % If regex fails, print the line for debugging
                    fprintf('Warning: Could not parse IK frame line: %s\n', line);
                end
            end
        end 
    end

    % convert to column vectors
    time_values = time_values(:);
    max_errors = max_errors(:);
    rms_errors = rms_errors(:);
    frame_numbers = frame_numbers(:);
    total_squared_errors = total_squared_errors(:);

    if ~isempty(time_values)
        fprintf('Time range: %.6f to %.6f seconds\n', min(time_values), max(time_values));
        fprintf('Frame range: %d to %d\n', min(frame_numbers), max(frame_numbers));
        fprintf('RMS error range: %.6f to %.6f\n', min(rms_errors), max(rms_errors));
        fprintf('Max error range: %.6f to %.6f\n', min(max_errors), max(max_errors));
        
        % Find most problematic marker
        unique_markers = unique(max_error_markers);
        fprintf('Markers with highest errors: ');
        for j = 1:min(3, length(unique_markers))
            marker_indices = strcmp(max_error_markers, unique_markers{j});
            avg_error = mean(max_errors(marker_indices));
            fprintf('%s(%.4f) ', unique_markers{j}, avg_error);
        end
        fprintf('\n');
    end
end

% Function to write IK results to XML and LOG files
function write_ik_results_to_files(time_values, max_errors, max_error_markers, rms_errors, ...
                                   frame_numbers, total_squared_errors, trial_name, cycle_num, output_dir)
    
    % Create base filename with cycle number
    base_filename = sprintf('%s_%04d', trial_name, cycle_num);
    xml_filename = fullfile(output_dir, [base_filename '_ik_errors.xml']);
    log_filename = fullfile(output_dir, [base_filename '_ik_errors.log']);
    
    % Write XML file
    write_ik_xml(xml_filename, time_values, max_errors, max_error_markers, ...
                 rms_errors, frame_numbers, total_squared_errors, trial_name, cycle_num);
    
    % Write LOG file
    write_ik_log(log_filename, time_values, max_errors, max_error_markers, ...
                 rms_errors, frame_numbers, total_squared_errors, trial_name, cycle_num);
    
    fprintf('IK error results written to:\n');
    fprintf('  XML: %s\n', xml_filename);
    fprintf('  LOG: %s\n', log_filename);
end

% Function to write XML file with IK error data
function write_ik_xml(filename, time_values, max_errors, max_error_markers, ...
                      rms_errors, frame_numbers, total_squared_errors, trial_name, cycle_num)
    
    try
        % Create XML document
        docNode = com.mathworks.xml.XMLUtils.createDocument('IKErrorAnalysis');
        rootElement = docNode.getDocumentElement;
        
        % Add metadata
        metadata = docNode.createElement('Metadata');
        rootElement.appendChild(metadata);
        
        % Trial information
        trialInfo = docNode.createElement('TrialInfo');
        trialInfo.setAttribute('name', trial_name);
        trialInfo.setAttribute('cycle', num2str(cycle_num));
        trialInfo.setAttribute('timestamp', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        metadata.appendChild(trialInfo);
        
        % Summary statistics
        summary = docNode.createElement('Summary');
        summary.setAttribute('total_frames', num2str(length(time_values)));
        summary.setAttribute('time_start', sprintf('%.6f', min(time_values)));
        summary.setAttribute('time_end', sprintf('%.6f', max(time_values)));
        summary.setAttribute('rms_error_mean', sprintf('%.6f', mean(rms_errors)));
        summary.setAttribute('rms_error_std', sprintf('%.6f', std(rms_errors)));
        summary.setAttribute('rms_error_max', sprintf('%.6f', max(rms_errors)));
        summary.setAttribute('max_error_mean', sprintf('%.6f', mean(max_errors)));
        summary.setAttribute('max_error_std', sprintf('%.6f', std(max_errors)));
        summary.setAttribute('max_error_max', sprintf('%.6f', max(max_errors)));
        if ~all(isnan(total_squared_errors))
            summary.setAttribute('total_sq_error_mean', sprintf('%.6f', nanmean(total_squared_errors)));
            summary.setAttribute('total_sq_error_max', sprintf('%.6f', nanmax(total_squared_errors)));
        end
        metadata.appendChild(summary);
        
        % Frame data
        frameData = docNode.createElement('FrameData');
        rootElement.appendChild(frameData);
        
        for i = 1:length(time_values)
            frame = docNode.createElement('Frame');
            frame.setAttribute('number', num2str(frame_numbers(i)));
            frame.setAttribute('time', sprintf('%.8f', time_values(i)));
            frame.setAttribute('rms_error', sprintf('%.8f', rms_errors(i)));
            frame.setAttribute('max_error', sprintf('%.8f', max_errors(i)));
            frame.setAttribute('max_error_marker', max_error_markers{i});
            if ~isnan(total_squared_errors(i))
                frame.setAttribute('total_squared_error', sprintf('%.8f', total_squared_errors(i)));
            end
            frameData.appendChild(frame);
        end
        
        % Write XML to file
        xmlwrite(filename, docNode);
        
    catch ME
        fprintf('Error writing XML file: %s\n', ME.message);
        % Fallback to simple XML writing
        write_simple_xml(filename, time_values, max_errors, max_error_markers, ...
                         rms_errors, frame_numbers, total_squared_errors, trial_name, cycle_num);
    end
end

% Fallback function for simple XML writing (if Java XML fails)
function write_simple_xml(filename, time_values, max_errors, max_error_markers, ...
                          rms_errors, frame_numbers, total_squared_errors, trial_name, cycle_num)
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not create XML file: %s', filename);
    end
    
    % Write XML header and root element
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<IKErrorAnalysis>\n');
    
    % Metadata
    fprintf(fid, '  <Metadata>\n');
    fprintf(fid, '    <TrialInfo name="%s" cycle="%d" timestamp="%s"/>\n', ...
            trial_name, cycle_num, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '    <Summary total_frames="%d" time_start="%.6f" time_end="%.6f"', ...
            length(time_values), min(time_values), max(time_values));
    fprintf(fid, ' rms_error_mean="%.6f" rms_error_std="%.6f" rms_error_max="%.6f"', ...
            mean(rms_errors), std(rms_errors), max(rms_errors));
    fprintf(fid, ' max_error_mean="%.6f" max_error_std="%.6f" max_error_max="%.6f"', ...
            mean(max_errors), std(max_errors), max(max_errors));
    if ~all(isnan(total_squared_errors))
        fprintf(fid, ' total_sq_error_mean="%.6f" total_sq_error_max="%.6f"', ...
                nanmean(total_squared_errors), nanmax(total_squared_errors));
    end
    fprintf(fid, '/>\n');
    fprintf(fid, '  </Metadata>\n');
    
    % Frame data
    fprintf(fid, '  <FrameData>\n');
    for i = 1:length(time_values)
        fprintf(fid, '    <Frame number="%d" time="%.8f" rms_error="%.8f" max_error="%.8f" max_error_marker="%s"', ...
                frame_numbers(i), time_values(i), rms_errors(i), max_errors(i), max_error_markers{i});
        if ~isnan(total_squared_errors(i))
            fprintf(fid, ' total_squared_error="%.8f"', total_squared_errors(i));
        end
        fprintf(fid, '/>\n');
    end
    fprintf(fid, '  </FrameData>\n');
    fprintf(fid, '</IKErrorAnalysis>\n');
    
    fclose(fid);
end

% Function to write LOG file with IK error data
function write_ik_log(filename, time_values, max_errors, max_error_markers, ...
                      rms_errors, frame_numbers, total_squared_errors, trial_name, cycle_num)
    
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not create LOG file: %s', filename);
    end
    
    % Write header
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'IK ERROR ANALYSIS LOG\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Trial Name: %s\n', trial_name);
    fprintf(fid, 'Cycle Number: %d\n', cycle_num);
    fprintf(fid, 'Analysis Timestamp: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, 'Total Frames: %d\n', length(time_values));
    fprintf(fid, 'Time Range: %.6f to %.6f seconds\n', min(time_values), max(time_values));
    fprintf(fid, '================================================================================\n\n');
    
    % Summary statistics
    fprintf(fid, 'SUMMARY STATISTICS\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'RMS Error Statistics:\n');
    fprintf(fid, '  Mean:     %.8f\n', mean(rms_errors));
    fprintf(fid, '  Std Dev:  %.8f\n', std(rms_errors));
    fprintf(fid, '  Min:      %.8f\n', min(rms_errors));
    fprintf(fid, '  Max:      %.8f\n', max(rms_errors));
    fprintf(fid, '\n');
    
    fprintf(fid, 'Max Error Statistics:\n');
    fprintf(fid, '  Mean:     %.8f\n', mean(max_errors));
    fprintf(fid, '  Std Dev:  %.8f\n', std(max_errors));
    fprintf(fid, '  Min:      %.8f\n', min(max_errors));
    fprintf(fid, '  Max:      %.8f\n', max(max_errors));
    fprintf(fid, '\n');
    
    if ~all(isnan(total_squared_errors))
        fprintf(fid, 'Total Squared Error Statistics:\n');
        fprintf(fid, '  Mean:     %.8f\n', nanmean(total_squared_errors));
        fprintf(fid, '  Std Dev:  %.8f\n', nanstd(total_squared_errors));
        fprintf(fid, '  Min:      %.8f\n', nanmin(total_squared_errors));
        fprintf(fid, '  Max:      %.8f\n', nanmax(total_squared_errors));
        fprintf(fid, '\n');
    end
    
    % Marker analysis
    unique_markers = unique(max_error_markers);
    fprintf(fid, 'MARKER ERROR ANALYSIS\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'Markers causing maximum errors (frequency and average error):\n');
    for j = 1:length(unique_markers)
        marker_indices = strcmp(max_error_markers, unique_markers{j});
        frequency = sum(marker_indices);
        avg_error = mean(max_errors(marker_indices));
        max_error_for_marker = max(max_errors(marker_indices));
        fprintf(fid, '  %-10s: %3d occurrences, avg=%.6f, max=%.6f\n', ...
                unique_markers{j}, frequency, avg_error, max_error_for_marker);
    end
    fprintf(fid, '\n');
    
    % Detailed frame data
    fprintf(fid, 'DETAILED FRAME DATA\n');
    fprintf(fid, '================================================================================\n');
    if ~all(isnan(total_squared_errors))
        fprintf(fid, '%6s %12s %12s %12s %18s %10s\n', ...
                'Frame', 'Time', 'RMS_Error', 'Max_Error', 'Total_Sq_Error', 'Marker');
        fprintf(fid, '%s\n', repmat('-', 1, 82));
        for i = 1:length(time_values)
            fprintf(fid, '%6d %12.8f %12.8f %12.8f %18.8f %10s\n', ...
                    frame_numbers(i), time_values(i), rms_errors(i), max_errors(i), ...
                    total_squared_errors(i), max_error_markers{i});
        end
    else
        fprintf(fid, '%6s %12s %12s %12s %10s\n', ...
                'Frame', 'Time', 'RMS_Error', 'Max_Error', 'Marker');
        fprintf(fid, '%s\n', repmat('-', 1, 64));
        for i = 1:length(time_values)
            fprintf(fid, '%6d %12.8f %12.8f %12.8f %10s\n', ...
                    frame_numbers(i), time_values(i), rms_errors(i), max_errors(i), ...
                    max_error_markers{i});
        end
    end
    
    fprintf(fid, '\n');
    fprintf(fid, '================================================================================\n');
    fprintf(fid, 'END OF LOG\n');
    fprintf(fid, '================================================================================\n');
    
    fclose(fid);
end

%% Normalize the force plate, marker and EMG data relative to:
%  - the completeness of the gait cycle
%  - the maximum magnitude which the trial contains
function [dataset, figures] = perform_all_normalizations(dataset, config, figures)

    if config.visualize.perform_all_normalizations && ...
       config.visualize.force_plate_data
        fprintf([ ...
            '\nReady for performing all normalizations?' + ...
            '\nJust select the GRF figure and press the "n" key\n']);
        figure(figures.grf.fig);
        while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
            waitforbuttonpress;
        end
        set(figures.grf.fig, 'CurrentCharacter', 'x');
    end

    fprintf("\nNormalizing EMG, GRF and Mocap data ...\n");

    % Compute the normalization factor for the maximum voluntary isometric contraction (MVIC).
    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            % Filter out the conditions: static and dynamic.
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'conmc'), conditions);
            conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                for idx_tri = 1:length(trials)

                    name_muscles = config.emg_data.name_muscles;
                    for idx_mus = 1:length(name_muscles)

                        data_muscle = ...
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus});

                        data_muscle = max(data_muscle);

                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.trial_maximum.(name_muscles{idx_mus}) = ...
                            data_muscle;
                    end % loops through muscles

                end % loops throught trials

            end % loops throught conditions

        end % loops throught sessions

    end % loops throught subjects

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            % Consider all conditions: maximum contraction, static and dynamic.
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
%             idxs_conditions = cellfun(@(x) contains(x, 'const'), conditions);
%             conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);
                for idx_tri = 1:length(trials)

                    name_trial = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
                    data_mocap = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;
                    data_emg   = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;

                    fprintf('%s normalizing ...\n', name_trial);

                    % Prepare the normalization of the EMG data.

                    % 1. Compute / extract:
                    % - the maximum voluntary isometric contraction MVIC, based on the maximum contraction trials.
                    % - the .?. DMVC, based on the maximum activation during the whole trial.
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.mvic = [];
                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.dmvc = [];

                    name_muscles = config.emg_data.name_muscles;
                    % for idx_mus = 1:length(name_muscles)

                    %     data_muscle = ...
                    %         dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus});

                        % 1.1 Compute the maximum voluntary isometric contraction MVIC.
                    %     contraction_max_trials = config.emg_data.maximum_contractions(subjects{idx_sub});
                    %     contraction_max_trial  = contraction_max_trials(name_muscles{idx_mus});

                        % Extract the maximum activation from the maximum contraction trial which corresponds to the muscle which is processed in the current loop iteration.
                    %     contraction_max_activation_max = ...
                    %         dataset.(subjects{idx_sub}).(sessions{idx_ses}).('conmc0000000000').(['tri', contraction_max_trial]).emg.mvic_maximum.(name_muscles{idx_mus});

                        % Save the mvic_maximum which will be used in normalization.
                    %     dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.mvic_maximum.(name_muscles{idx_mus}) = ...
                    %         contraction_max_activation_max;

                        % 1.2 Compute the .?. DMVC.
                    %     data_muscle_max_peak = max(abs(data_muscle));

                        % Save the dmvc_maximum which will be used in normalization.
                    %     dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.dmvc_maximum.(name_muscles{idx_mus}) = ...
                    %         data_muscle_max_peak;

                    % end % loops throught muscles

                    % Update the data_emg variable, because (in the previous loop) we updated its source.
                    data_emg = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;

                    % 2. Perform the normalization of EMG, GRF and mocap data.
                    cycles      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}));
                    idxs_cycles = cellfun(@(x) contains(x, 'cyc'), cycles);
                    cycles      = cycles(idxs_cycles);
                    for idx_cyc = 1:length(cycles)

                        cycle_start = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).start;
                        cycle_end   = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).end;

                        % 2.1 Normalize the EMG data by computing:
                        %    - MVIC maximum voluntary isometric contraction and
                        %    - DMVC .?.
                        name_muscles = config.emg_data.name_muscles;
                        for idx_mus = 1:length(name_muscles)

                            data_muscle = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus})(cycle_start:cycle_end);
                            
                            % Get the trial peak value (instead of MVIC from separate trials)
                            if isfield(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg, 'trial_maximum')
                                trial_peak_max = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.trial_maximum.(name_muscles{idx_mus});
                            else
                                % Fallback: calculate trial peak if not already stored
                                trial_data = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.clean.(name_muscles{idx_mus});
                                trial_peak_max = max(trial_data);
                                
                                % Store it for future use
                                dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.trial_maximum.(name_muscles{idx_mus}) = trial_peak_max;
                                fprintf('Calculated trial peak for %s: %f\n', name_muscles{idx_mus}, trial_peak_max);
                            end
                            % mvic_max    = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.mvic_maximum.(name_muscles{idx_mus});
                            dmvc_max    =  trial_peak_max % dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg.dmvc_maximum.(name_muscles{idx_mus});

                            % 2.1.1 Perform trial peak normalization
                            data_muscle_trial_peak = data_muscle / trial_peak_max;

                            % 2.1.3.Perform the DMVC
                            data_muscle_dmvc = data_muscle / dmvc_max;

                            % 2.1.3 Normalize relative to the completeness of the gait cycle.
                            
                            % Create a linspace vector for the current gait cycle based on the number of data points between stances.
                            data_trial_peak_length = linspace(0, 100, length(data_muscle_trial_peak));
                            data_dmvc_length = linspace(0, 100, length(data_muscle_dmvc));

                            data_muscle_trial_peak_gait_cycle_completeness = spline(data_trial_peak_length, data_muscle_trial_peak, 0:99);
                            data_muscle_dmvc_gait_cycle_completeness = spline(data_dmvc_length, data_muscle_dmvc, 0:99);

                            % 2.1.4 Save the EMG data normalized relative to trial peak and completeness of the gait cycle.
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).emg.trial_peak_gcc.(name_muscles{idx_mus}) = ...
                                data_muscle_trial_peak_gait_cycle_completeness;

                            % 2.1.5 Save the EMG data normalized relative to DMVC and completeness of the gait cycle.
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).emg.dmvc_gcc.(name_muscles{idx_mus}) = ...
                                data_muscle_dmvc_gait_cycle_completeness;

                        end % loops through muscles

                        % After all muscles for a cycle are processed, export to .sto
                        c3d2emg(dataset, subjects{idx_sub}, sessions{idx_ses}, conditions{idx_con}, trials{idx_tri}, cycles{idx_cyc}, config, config.path_dataset);

                        % 2.2 Normalize the GRF data, relative to the body weight of the subject.
                        for idx_grf = 1:length(data_mocap.fp_data.GRF_data)

                            data_F = data_mocap.fp_data.GRF_data(idx_grf).F(cycle_start:cycle_end, :);

                            % 2.2.1 Normalize the GRF magnitude relative to subject's body weight.
                            %       - inclusive the lateral and forward components!
                            data_F_body_weight = data_F / config.subjects.masses(subjects{idx_sub}) / 9.8;

                            % 2.3.1 Normalize the GRF length relative to the completeness of the gait cycle.
                            data_F_body_weight_gait_cycle_completeness = zeros(100, 3);
                            for idx_coord = 1:3

                                % Spline interpolation to normalize the data length to 100 points.
                                data_F_body_weight_gait_cycle_completeness(:, idx_coord) = ...
                                    spline(linspace(0, 100, length(data_F_body_weight)), data_F_body_weight(:, idx_coord), 0:99);
                            end

                            % 2.1.3 Save the GRF data normalized relative to subject's body weight and completeness of the gait cycle.
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).mocap.GRF_data(idx_grf).F_bw_gcc = ...
                               data_F_body_weight_gait_cycle_completeness;

                        end % loops throught force plates

                        % 2.3 Normalize the mocap data relative to the completeness of the gait cycle.
                        markers = fieldnames(data_mocap.marker_data.Markers);
                        for idx_mar = 1:length(markers)

                            data_marker = data_mocap.marker_data.Markers.(markers{idx_mar})(cycle_start:cycle_end,:);

                            % 2.3.1 Normalize marker data to the completeness of the gait cycle.
                            data_marker_gait_cycle_completeness = zeros(100, 3);
                            for idx_coord = 1:3

                                % Spline interpolation to normalize the
                                % data length to 100 points.b
                                data_marker_gait_cycle_completeness(:, idx_coord) = ...
                                    spline(linspace(0, 100, length(data_marker)), data_marker(:, idx_coord), 0:99);

                            end % loops through 3D coordinates

                            % 2.3.2 Save the mocap data normalized relative to the completeness of the gait cycle.
                            dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).mocap.marker_data_gcc.Markers.(markers{idx_mar}) = ...
                                data_marker_gait_cycle_completeness;

                        end % loops through mocap markers

                        % Plotting

                        % Use the same color for data describing the right side.
                        color_right_side = rand(1, 3);

                        if config.visualize.perform_all_normalizations && ...
                        config.visualize.emg_data

                            name_muscles = config.emg_data.name_muscles;
                            for idx_mus = 1:length(name_muscles)

                                hold(figures.emg.figs{idx_mus}.activation, "on");

                                plot(figures.emg.figs{idx_mus}.activation, ...
                                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).emg.trial_peak_gcc.(name_muscles{idx_mus}), ...
                                    'DisplayName', sprintf("%s %2d normalized Trial Peak", trials{idx_tri}, idx_cyc));
                                plot(figures.emg.figs{idx_mus}.activation, ...
                                    dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).emg.dmvc_gcc.(name_muscles{idx_mus}), ...
                                    'DisplayName', sprintf("%s %2d normalized DMVC", trials{idx_tri}, idx_cyc));

                                hold(figures.emg.figs{idx_mus}.activation, "off");

                                xlim(figures.emg.figs{idx_mus}.activation, [-1.0, +101.0]);
                                % ylim(figures.emg.figs{idx_mus}.activation, [-0.1,   +1.1]);

                                xlabel(figures.emg.figs{idx_mus}.activation, 'Gait cycle (%)');
                                ylabel(figures.emg.figs{idx_mus}.activation, 'Activation (% Trial Peak / % DMVC)');

                                title_emg = sprintf("%s, %s, cycle %2d of total %d cycles", ...
                                    name_trial, name_muscles{idx_mus}, idx_cyc, length(cycles));
                                title(figures.emg.figs{idx_mus}.activation, title_emg);

                                drawnow;

                            end % loops through muscles

                        end

                        if config.visualize.perform_all_normalizations && ...
                           config.visualize.force_plate_data

                            hold(figures.grf.forces, "on");
                            for idx_grf = 1:length(data_mocap.fp_data.GRF_data)
                                if idx_grf == 2
                                    legend_grf = sprintf("%s %2d GRF%1d (right) normalized bw. gcc.", trials{idx_tri}, idx_cyc, idx_grf);
                                    plot(figures.grf.forces, ...
                                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).mocap.GRF_data(idx_grf).F_bw_gcc(:, 3), ...
                                        'DisplayName', legend_grf, 'Color', color_right_side, 'LineWidth', 2,'LineStyle', '--');
                                else
                                    legend_grf = sprintf("%s %2d GRF%1d (left) normalized bw. gcc.", trials{idx_tri}, idx_cyc, idx_grf);
                                    plot(figures.grf.forces, ...
                                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).mocap.GRF_data(idx_grf).F_bw_gcc(:, 3), ...
                                        'DisplayName', legend_grf);
                                end
                            end % loops throught force plates
                            hold(figures.grf.forces, "off");

                            xlim(figures.grf.forces, [-1.0,  +101.0]);
                            ylim(figures.grf.forces, [-0.1,    +3.5]);

                            xlabel(figures.grf.forces, 'Gait cycle (%)');
                            ylabel(figures.grf.forces, 'Body weight factor (BW)');

                            title_grf = sprintf( ...
                                "%s, GRFs, cycle %2d of %2d total cycles.", name_trial, idx_cyc, length(cycles));
                            title(figures.grf.forces, title_grf);

                            drawnow;
                        end

                        if config.visualize.perform_all_normalizations && ...
                           config.visualize.marker_data

                            markers = fieldnames(data_mocap.marker_data.Markers);
                            for idx_mar = 1:min(length(figures.mocap.figs), length(markers))

                                hold(figures.mocap.figs{idx_mar}.coordinates, "on");

                                name_coords = {'x', 'y', 'z'};
                                for idx_coord = 1:length(name_coords)
                                    legend_coord = sprintf("%s %2d %s %s normalized gcc.", ...
                                        trials{idx_tri}, idx_cyc, markers{idx_mar}, name_coords{idx_coord});
                                    plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                        dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).mocap.marker_data_gcc.Markers.(markers{idx_mar})(:, idx_coord), ...
                                        'DisplayName', legend_coord);
                                end % loops through 3D coordinates

                                hold(figures.mocap.figs{idx_mar}.coordinates, "off");

                                xlim(figures.mocap.figs{idx_mar}.coordinates, [-1.0,  +101.0]);
                                % ylim(figures.mocap.figs{idx_mar}.coordinates, [?, ?]);

                                xlabel(figures.mocap.figs{idx_mar}.coordinates, 'Gait cycle (%)');
                                ylabel(figures.mocap.figs{idx_mar}.coordinates, 'Marker coordinates (mm)');

                                title_marker = sprintf("%s, %s, cycle %2d of total %d", ...
                                    name_trial, markers{idx_mar}, idx_cyc, length(cycles));
                                title(figures.mocap.figs{idx_mar}.coordinates, title_marker);

                                drawnow;

                            end % loops through mocap markers
                        end

                        % Semi-interactive visualization.
                        if config.visualize.perform_all_normalizations && ...
                           (config.visualize.subject || ...
                            config.visualize.session || ...
                            config.visualize.condition || ...
                            config.visualize.trial || ...
                            config.visualize.gait_cycle)

                            pause(config.visualize.pause)
                        end

                        % Visualize each gait cycles separate.
                        if config.visualize.gait_cycle
                            clear_figures(config, figures);
                        end

                    end  % loops through gait cycles

                    if config.visualize.trial
                        clear_figures(config, figures);
                    end

                end  % loops through trials

                if config.visualize.condition
                    clear_figures(config, figures);
                end

            end  % loops through conditions

            if config.visualize.session
                clear_figures(config, figures);
            end

        end % loops through sessions

        if config.visualize.subject
            clear_figures(config, figures);
        end

    end % loops through subjects

end

%%
function [dataset, figures] = compute_means_and_stds(dataset, config, figures)

    if config.visualize.compute_means_and_stds && ...
       config.visualize.force_plate_data
        fprintf([ ...
            '\nReady for computing all means and standard deviations?' + ...
            '\nJust select the GRF figure and press the "n" key\n']);
        figure(figures.grf.fig);
        while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
            waitforbuttonpress;
        end
        set(figures.grf.fig, 'CurrentCharacter', 'x');
    end

    subjects      = fieldnames(dataset);
    idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
    subjects      = subjects(idxs_subjects);
    for idx_sub = 1:length(subjects)

        sessions      = fieldnames(dataset.(subjects{idx_sub}));
        idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
        sessions      = sessions(idxs_sessions);
        for idx_ses = 1:length(sessions)

            % Consider only the dynamic conditions.
            % - only the dynamic conditions contain gait cycles!
            conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
            idxs_conditions = cellfun(@(x) ~contains(x, 'const') && ~contains(x, 'conmc'), conditions);
            conditions      = conditions(idxs_conditions);
            for idx_con = 1:length(conditions)

                trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
                idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
                trials      = trials(idxs_trials);

                % Initialize the placeholder arrays for computing the means and the standard deviations.
                name_muscles = config.emg_data.name_muscles;
                for idx_mus = 1:length(name_muscles)
                    condition_emg_trial_peak_gcc.(name_muscles{idx_mus}) = [];
                    condition_emg_dmvc_gcc.(name_muscles{idx_mus}) = [];
                end

                for idx_grf = 1:length(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{1}).mocap.fp_data.GRF_data)
                    condition_F_bw_gcc.(['grf', char(num2str(idx_grf))]) = [];
                end

                markers = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{1}).mocap.marker_data.Markers);
                for idx_mar = 1:length(markers)
                    name_coords = {'x', 'y', 'z'};
                    for idx_coord = 1:length(name_coords)
                        condition_marker_gcc.(markers{idx_mar}).(name_coords{idx_coord}) = [];
                    end
                end

                for idx_tri = 1:length(trials)
                    cycles      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}));
                    idxs_cycles = cellfun(@(x) contains(x, 'cyc'), cycles);
                    cycles      = cycles(idxs_cycles);
                    for idx_cyc = 1:length(cycles)

                        cycle_emg   = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).emg;
                        cycle_mocap = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).(cycles{idx_cyc}).mocap;

                        name_muscles = config.emg_data.name_muscles;
                        for idx_mus = 1:length(name_muscles)

                            condition_emg_trial_peak_gcc.(name_muscles{idx_mus}) = [ ...
                                condition_emg_trial_peak_gcc.(name_muscles{idx_mus}); ...
                                cycle_emg.trial_peak_gcc.(name_muscles{idx_mus})];

                            condition_emg_dmvc_gcc.(name_muscles{idx_mus}) = [ ...
                                condition_emg_dmvc_gcc.(name_muscles{idx_mus}); ...
                                cycle_emg.dmvc_gcc.(name_muscles{idx_mus})];
                        end


                        %Processes GRF data for this cycle
                        %for idx_grf = 1:length(cycle_mocap.GRF_data)

                        %end

                        for idx_grf = 1:length(cycle_mocap.GRF_data)
                            condition_F_bw_gcc.(['grf', char(num2str(idx_grf))]) = [ ...
                                condition_F_bw_gcc.(['grf', char(num2str(idx_grf))]); ...
                                cycle_mocap.GRF_data(idx_grf).F_bw_gcc(:, 3).'];
                            condition_F_bw_gcc_mean.(['grf', char(num2str(idx_grf))]) = mean(condition_F_bw_gcc.(['grf', char(num2str(idx_grf))]), 1);
                            condition_F_bw_gcc_std.(['grf', char(num2str(idx_grf))])  =  std(condition_F_bw_gcc.(['grf', char(num2str(idx_grf))]), 1);

                            condition_F_bw_gcc_std_above.(['grf', char(num2str(idx_grf))]) = ...
                            condition_F_bw_gcc_mean.(['grf', char(num2str(idx_grf))]) + condition_F_bw_gcc_std.(['grf', char(num2str(idx_grf))]);
                            condition_F_bw_gcc_std_below.(['grf', char(num2str(idx_grf))]) = ...
                            condition_F_bw_gcc_mean.(['grf', char(num2str(idx_grf))]) - condition_F_bw_gcc_std.(['grf', char(num2str(idx_grf))]);
                        end

                        markers = fieldnames(cycle_mocap.marker_data_gcc.Markers);
                        for idx_mar = 1:length(markers)
                            name_coords = {'x', 'y', 'z'};
                            for idx_coord = 1:length(name_coords)
                                condition_marker_gcc.(markers{idx_mar}).(name_coords{idx_coord}) = [ ...
                                    condition_marker_gcc.(markers{idx_mar}).(name_coords{idx_coord});
                                    cycle_mocap.marker_data_gcc.Markers.(markers{idx_mar})(:, idx_coord).'];
                                condition_marker_gcc_mean.(markers{idx_mar}).(name_coords{idx_coord})       = mean(condition_marker_gcc.(markers{idx_mar}).(name_coords{idx_coord}));
                                condition_marker_gcc__std.(markers{idx_mar}).(name_coords{idx_coord})       =  std(condition_marker_gcc.(markers{idx_mar}).(name_coords{idx_coord}));
                                condition_marker_gcc__std_above.(markers{idx_mar}).(name_coords{idx_coord}) = condition_marker_gcc_mean.(markers{idx_mar}).(name_coords{idx_coord}) + condition_marker_gcc__std.(markers{idx_mar}).(name_coords{idx_coord});
                                condition_marker_gcc__std_below.(markers{idx_mar}).(name_coords{idx_coord}) = condition_marker_gcc_mean.(markers{idx_mar}).(name_coords{idx_coord}) - condition_marker_gcc__std.(markers{idx_mar}).(name_coords{idx_coord});
                            end
                        end

                        %Processes marker data for this cycle
                        %markers = fieldnames(cycle_mocap.marker_data_gcc.Markers);
                        %for idx_mar = 1:length(markers)

                         %   name_coords = {'x', 'y', 'z'};
                          %  for idx_coord = 1:length(name_coords)
                                
                           % end
                        %end

                    end
                end % loops through the trials

                name_muscles = config.emg_data.name_muscles;
                % for idx_mus = 1:length(name_muscles)
                %    condition_emg_trial_peak_gcc_mean.(name_muscles{idx_mus}) = mean(condition_emg_trial_peak_gcc.(name_muscles{idx_mus}), 1);
                %    condition_emg_trial_peak_gcc__std.(name_muscles{idx_mus}) =  std(condition_emg_trial_peak_gcc.(name_muscles{idx_mus}), 1);

                %    condition_emg_trial_peak_gcc_std_above.(name_muscles{idx_mus}) = ...
                %        condition_emg_trial_peak_gcc_mean.(name_muscles{idx_mus}) + condition_emg_trial_peak_gcc__std.(name_muscles{idx_mus});
                %    condition_emg_trial_peak_gcc_std_below.(name_muscles{idx_mus}) = ...
                %        condition_emg_trial_peak_gcc_mean.(name_muscles{idx_mus}) - condition_emg_trial_peak_gcc__std.(name_muscles{idx_mus});

                %    condition_emg_dmvc_gcc_mean.(name_muscles{idx_mus}) = mean(condition_emg_dmvc_gcc.(name_muscles{idx_mus}), 1);
                %    condition_emg_dmvc_gcc__std.(name_muscles{idx_mus}) =  std(condition_emg_dmvc_gcc.(name_muscles{idx_mus}), 1);

                %    condition_emg_dmvc_gcc_std_above.(name_muscles{idx_mus}) = ...
                %        condition_emg_dmvc_gcc_mean.(name_muscles{idx_mus}) + condition_emg_dmvc_gcc__std.(name_muscles{idx_mus});
                %    condition_emg_dmvc_gcc_std_below.(name_muscles{idx_mus}) = ...
                %        condition_emg_dmvc_gcc_mean.(name_muscles{idx_mus}) - condition_emg_dmvc_gcc__std.(name_muscles{idx_mus});
                % end

                
                if config.visualize.compute_means_and_stds && ...
                   config.visualize.emg_data

                    name_muscles = config.emg_data.name_muscles;
                    for idx_mus = 1:length(name_muscles)

                        hold(figures.emg.figs{idx_mus}.activation, "on");
                        plot(figures.emg.figs{idx_mus}.activation, ...
                            condition_emg_mvic_gcc_mean.(name_muscles{idx_mus}), ...
                            'DisplayName', sprintf("%s MVIC mean", name_muscles{idx_mus}));
                        fill(figures.emg.figs{idx_mus}.activation, ...
                            [1:100, fliplr(1:100)], ...
                            [condition_emg_mvic_gcc_std_below.(name_muscles{idx_mus}), fliplr(condition_emg_mvic_gcc_std_above.(name_muscles{idx_mus}))], 'b', 'FaceAlpha', 0.1);
                        plot(figures.emg.figs{idx_mus}.activation, ...
                            condition_emg_mvic_gcc_std_below.(name_muscles{idx_mus}), 'b--', 'LineWidth', 1);
                        plot(figures.emg.figs{idx_mus}.activation, ...
                            condition_emg_mvic_gcc_std_above.(name_muscles{idx_mus}), 'b--', 'LineWidth', 1);

                        plot(figures.emg.figs{idx_mus}.activation, ...
                            condition_emg_dmvc_gcc_mean.(name_muscles{idx_mus}), ...
                            'DisplayName', sprintf("%s DMVC mean", name_muscles{idx_mus}));
                        fill(figures.emg.figs{idx_mus}.activation, ...
                            [1:100, fliplr(1:100)], ...
                            [condition_emg_dmvc_gcc_std_below.(name_muscles{idx_mus}), fliplr(condition_emg_dmvc_gcc_std_above.(name_muscles{idx_mus}))], 'b', 'FaceAlpha', 0.1);
                        plot(figures.emg.figs{idx_mus}.activation, ...
                            condition_emg_dmvc_gcc_std_below.(name_muscles{idx_mus}), 'b--', 'LineWidth', 1);
                        plot(figures.emg.figs{idx_mus}.activation, ...
                            condition_emg_dmvc_gcc_std_above.(name_muscles{idx_mus}), 'b--', 'LineWidth', 1);
                        hold(figures.emg.figs{idx_mus}.activation, "off");

                        xlim(figures.emg.figs{idx_mus}.activation, [-1.0, +101.0]);
                        % ylim(figures.emg.figs{idx_mus}.activation, [-0.1,   +1.1]);

                        xlabel(figures.emg.figs{idx_mus}.activation, 'Gait cycle (%)');
                        ylabel(figures.emg.figs{idx_mus}.activation, 'Activation (?%)');

                        title_emg = sprintf("%s mean and standard deviation", conditions{idx_con});
                        title(figures.emg.figs{idx_mus}.activation, title_emg);

                        drawnow;
                    end
                end

                if config.visualize.compute_means_and_stds && ...
                   config.visualize.force_plate_data

                    hold(figures.grf.forces, "on");
                    for idx_grf = 1:length(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{1}).mocap.fp_data.GRF_data)
                        if idx_grf == 2
                            plot(figures.grf.forces, ...
                                condition_F_bw_gcc_mean.(['grf', char(num2str(idx_grf))]), ...
                                'DisplayName', sprintf('GRF %d mean', idx_grf));
                            fill(figures.grf.forces, ...
                                [1:100, fliplr(1:100)], ...
                                [condition_F_bw_gcc_std_below.(['grf', char(num2str(idx_grf))]), fliplr(condition_F_bw_gcc_std_above.(['grf', char(num2str(idx_grf))]))], 'b', 'FaceAlpha', 0.1);
                            plot(figures.grf.forces, ...
                                condition_F_bw_gcc_std_below.(['grf', char(num2str(idx_grf))]), 'b--', 'LineWidth', 1);
                            plot(figures.grf.forces, ...
                                condition_F_bw_gcc_std_above.(['grf', char(num2str(idx_grf))]), 'b--', 'LineWidth', 1);
                        else
                            plot(figures.grf.forces, ...
                                condition_F_bw_gcc_mean.(['grf', char(num2str(idx_grf))]), ...
                                'DisplayName', sprintf('GRF %d mean', idx_grf));
                            fill(figures.grf.forces, ...
                                [1:100, fliplr(1:100)], ...
                                [condition_F_bw_gcc_std_below.(['grf', char(num2str(idx_grf))]), fliplr(condition_F_bw_gcc_std_above.(['grf', char(num2str(idx_grf))]))], 'b', 'FaceAlpha', 0.1);
                            plot(figures.grf.forces, ...
                                condition_F_bw_gcc_std_below.(['grf', char(num2str(idx_grf))]), 'b--', 'LineWidth', 1);
                            plot(figures.grf.forces, ...
                                condition_F_bw_gcc_std_above.(['grf', char(num2str(idx_grf))]), 'b--', 'LineWidth', 1);
                        end
                    end
                    hold(figures.grf.forces, "off");

                    xlim(figures.grf.forces, [-1.0,  +101.0]);
                    ylim(figures.grf.forces, [-0.1,    +3.5]);

                    xlabel(figures.grf.forces, 'Gait cycle (%)');
                    ylabel(figures.grf.forces, 'Body weight factor (BW)');

                    title_grf = sprintf("%s mean and standard deviation", conditions{idx_con});
                    title(figures.grf.forces, title_grf);

                    drawnow;
                end

                if config.visualize.compute_means_and_stds && ...
                   config.visualize.marker_data

                    markers = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{1}).mocap.marker_data.Markers);
                    for idx_mar = 1:min(length(figures.mocap.figs), length(markers))

                        hold(figures.mocap.figs{idx_mar}.coordinates, "on");
                        name_coords = {'x', 'y', 'z'};
                        for idx_coord = 1:length(name_coords)
                            plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                condition_marker_gcc_mean.(markers{idx_mar}).(name_coords{idx_coord}),'DisplayName', sprintf("%s", conditions{idx_con}));
                            fill(figures.mocap.figs{idx_mar}.coordinates, ...
                                [1:100, fliplr(1:100)], ...
                                [      condition_marker_gcc__std_below.(markers{idx_mar}).(name_coords{idx_coord}), ...
                                fliplr(condition_marker_gcc__std_above.(markers{idx_mar}).(name_coords{idx_coord}))], 'b', 'FaceAlpha', 0.1);

                            plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                condition_marker_gcc__std_above.(markers{idx_mar}).(name_coords{idx_coord}),'DisplayName', sprintf("%s", name_coords{idx_coord}));
                            plot(figures.mocap.figs{idx_mar}.coordinates, ...
                                condition_marker_gcc__std_below.(markers{idx_mar}).(name_coords{idx_coord}),'DisplayName', sprintf("%s", name_coords{idx_coord}));
                        end
                        hold(figures.mocap.figs{idx_mar}.coordinates, "off");

                        xlim(figures.mocap.figs{idx_mar}.coordinates, [-1.0,  +101.0]);
                        % ylim(figures.mocap.figs{idx_mar}.coordinates, [?, ?]);

                        xlabel(figures.mocap.figs{idx_mar}.coordinates, 'Gait cycle (%)');
                        ylabel(figures.mocap.figs{idx_mar}.coordinates, 'Marker coordinates');

                        title_marker = sprintf("%s mean and standard deviation", conditions{idx_con});
                        title(figures.mocap.figs{idx_mar}.coordinates, title_marker);

                        drawnow;
                    end
                end

                if config.visualize.compute_means_and_stds && ...
                   (config.visualize.subject || config.visualize.session || config.visualize.condition) && ...
                   config.visualize.force_plate_data
                    fprintf(...
                        '\nReady for visualizing the means and standard deviations for the next condition?' + ...
                        '\nJust select the GRF figure and press the "n" key\n');
                    figure(figures.grf.fig);
                    while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
                        waitforbuttonpress;
                    end
                    set(figures.grf.fig, 'CurrentCharacter', 'x');
                end

                % Visualize the the means and standard deviation for each condition.
                if config.visualize.condition
                    clear_figures(config, figures);
                end

            end % loops through conditions

            % Visualize the gait cycles overlapped.
            if config.visualize.condition
                clear_figures(config, figures);
            end

        end % loops through sessions

        % Visualize the conditions overlapped.
        if config.visualize.session
            clear_figures(config, figures);
        end

    end % loops through subjects.

    % Visualize the the subjects overlapped.
    if config.visualize.subject
        clear_figures(config, figures);
    end

    fprintf("Processing done ... closing all figures.\n");
    close all;

end
