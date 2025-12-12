
function [dataset] = load_trials(config)

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

%     % Print summary of the subjects, sessions and trials loaded.
%     subjects = fieldnames(dataset);
%     for idx_sub = 1:length(subjects)
%         fprintf("%s\n", subjects{idx_sub});
%         sessions = fieldnames(dataset.(subjects{idx_sub}));
%         for idx_ses = 1:length(sessions)
%             fprintf(" %s\n", sessions{idx_ses});
%             conditions = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
%             for idx_con = 1:length(conditions)
%                 fprintf("  %s\n", conditions{idx_con});
%                 trials = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
%                 for idx_tri = 1:length(trials)
%                     fprintf("   %s\n", trials{idx_tri});
%                 end
%             end
%         end
%     end

%     figures = create_figures(config);
% 
%     % Wait of key pressed before plotting the next trial.
%     if config.visualize.force_plate_data && ...
%        config.visualize.raw_data
%         fprintf([ ...
%             '\nReady to visualize the raw data?' ...
%             '\nJust click on the GRF figure and press the "n" key!\n']);
%         figure(figures.grf.fig);
%         while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
%             waitforbuttonpress;
%         end
%         set(figures.grf.fig, 'CurrentCharacter', 'x');
%     end
% 
%     subjects      = fieldnames(dataset);
%     idxs_subjects = cellfun(@(x) contains(x, 'sub'), subjects);
%     subjects      = subjects(idxs_subjects);
%     for idx_sub = 1:length(subjects)
% 
%         sessions      = fieldnames(dataset.(subjects{idx_sub}));
%         idxs_sessions = cellfun(@(x) contains(x, 'ses'), sessions);
%         sessions      = sessions(idxs_sessions);
%         for idx_ses = 1:length(sessions)
% 
%             conditions      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}));
%             idxs_conditions = cellfun(@(x) contains(x, 'con'), conditions);
%             conditions      = conditions(idxs_conditions);
%             for idx_con = 1:length(conditions)
% 
%                 trials      = fieldnames(dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}));
%                 idxs_trials = cellfun(@(x) contains(x, 'tri'), trials);
%                 trials      = trials(idxs_trials);
%                 for idx_tri = 1:length(trials)
% 
%                     name_trial      = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).name_trial;
%                     data_treadmill  = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).treadmill;
%                     data_mocap      = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).mocap;
%                     data_emg        = dataset.(subjects{idx_sub}).(sessions{idx_ses}).(conditions{idx_con}).(trials{idx_tri}).emg;
% 
%                     % Plot raw data
% 
%                     % Use the same color for data describing the right side.
%                     color_right_side = rand(1, 3);
% 
%                     if config.visualize.raw_data && ...
%                        config.visualize.treadmill_data
% 
%                         hold(figures.treadmill.distances, "on");
%                         plot(figures.treadmill.distances, ...
%                             data_treadmill.Time, data_treadmill.Leftbelt_Distance, "DisplayName", "left");
%                         plot(figures.treadmill.distances, ...
%                             data_treadmill.Time, data_treadmill.Rightbelt_Distance, "DisplayName", "right", "Color", color_right_side);
%                         hold(figures.treadmill.distances, "off");
% 
%                         hold(figures.treadmill.speeds, "on");
%                         plot(figures.treadmill.speeds, ...
%                             data_treadmill.Time, data_treadmill.Leftbelt_Speed, "DisplayName", "left");
%                         plot(figures.treadmill.speeds, ...
%                             data_treadmill.Time, data_treadmill.Rightbelt_Speed, "DisplayName", "right", "Color", color_right_side);
%                         hold(figures.treadmill.speeds, "off");
% 
%                         hold(figures.treadmill.pitch, "on");
%                         plot(figures.treadmill.pitch, ...
%                             data_treadmill.Time, data_treadmill.PitchRotX, "DisplayName", "pitch");
%                         hold(figures.treadmill.pitch, "off");
% 
%                         hold(figures.treadmill.sway, "on");
%                         plot(figures.treadmill.sway, ...
%                             data_treadmill.Time, data_treadmill.SwayPosX, "DisplayName", "sway");
%                         hold(figures.treadmill.sway, "off");
% 
%                         drawnow;
%                     end
% 
%                     if config.visualize.raw_data && ...
%                        config.visualize.emg_data
% 
%                         channel_muscles = config.emg_data.subjects.muscles(subjects{idx_sub});
%                         name_muscles    = config.emg_data.name_muscles;
%                         for idx_mus = 1:length(name_muscles)
%                             data_muscle_raw = data_emg.analogs.(channel_muscles(name_muscles{idx_mus}));
% 
%                             hold(figures.emg.figs{idx_mus}.activation, "on");
% 
%                             plot(figures.emg.figs{idx_mus}.activation, ...
%                                 data_muscle_raw, ...
%                                 'DisplayName', sprintf('%s raw', name_muscles{idx_mus}), 'Color', color_right_side);
% 
%                             hold(figures.emg.figs{idx_mus}.activation, "off");
% 
%                             title_emg = [name_trial ' ' name_muscles{idx_mus}];
%                             title(figures.emg.figs{idx_mus}.activation, title_emg);
% 
%                             drawnow;
%                         end
%                     end
% 
%                     if config.visualize.raw_data && ...
%                        config.visualize.force_plate_data
% 
%                         hold(figures.grf.forces, "on");
%                         for idx_grf = 1:length(data_mocap.fp_data.GRF_data)
%                             if idx_grf == 2
%                                 plot(figures.grf.forces, data_mocap.fp_data.GRF_data(idx_grf).F(:, 3), ...
%                                     'DisplayName', 'GRF2 (right leg) raw', 'Color', color_right_side);
%                             else
%                                 plot(figures.grf.forces, data_mocap.fp_data.GRF_data(idx_grf).F(:, 3), ...
%                                     'DisplayName', 'GRF1 (left leg) raw');
%                             end
%                         end
%                         hold(figures.grf.forces, "off");
% 
%                         title_grf = [name_trial, ' GRFs'];
%                         title(figures.grf.forces, title_grf);
% 
%                         drawnow;
%                     end
% 
%                     if config.visualize.raw_data && ...
%                        config.visualize.marker_data
% 
%                         markers = fieldnames(data_mocap.marker_data.Markers);
%                         for idx_mar = 1:min(length(figures.mocap.figs), length(markers))
% 
%                             hold(figures.mocap.figs{idx_mar}.coordinates, "on");
%                             name_coords = {'x', 'y', 'z'};
%                             for idx_coord = 1:length(name_coords)
%                                 legend_coord = sprintf('%s %s raw', ...
%                                     markers{idx_mar}, name_coords{idx_coord});
%                                 plot(figures.mocap.figs{idx_mar}.coordinates, ...
%                                     data_mocap.marker_data.Markers.(markers{idx_mar})(:, idx_coord), ...
%                                     'DisplayName', legend_coord);
%                             end
%                             hold(figures.mocap.figs{idx_mar}.coordinates, "off");
% 
%                             title_marker = [name_trial ' ' markers{idx_mar}];
%                             title(figures.mocap.figs{idx_mar}.coordinates, title_marker);
% 
%                             drawnow;
%                         end
%                     end
% 
% %                     % Wait of key pressed before plotting the next trial.
% %                     if config.visualize.force_plate_data && ...
% %                        config.visualize.raw_data
% %                         fprintf([ ...
% %                             '\nReady to visualize the raw data of the next trial?' ...
% %                             '\nJust click on the GRF figure and press the "n" key!\n']);
% %                         figure(figures.grf.fig);
% %                         while ~strcmp(get(figures.grf.fig, 'CurrentCharacter'), 'n')
% %                             waitforbuttonpress;
% %                         end
% %                         set(figures.grf.fig, 'CurrentCharacter', 'x');
% %                     end
% 
%                     % Semi-interactive visualization
%                     if config.visualize.raw_data&& ...
%                        (config.visualize.subject || ...
%                         config.visualize.session || ...
%                         config.visualize.condition || ...
%                         config.visualize.trial || ...
%                         config.visualize.gait_cycle)
% 
%                         pause(config.visualize.pause)
%                     end
% 
%                     if config.visualize.raw_data && ...
%                        config.visualize.trial
%                         clear_figures(config, figures);
%                     end
%                 end
%             end
%         end
%     end
end
