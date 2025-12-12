function data_static_calibrate(path_subject, file_trial, data)

import org.opensim.modeling.*;

%% Config the plotting flags.
% plot_raw_data = true;
% plot_raw_data = false;

%% Load the *.c3d and *.csv files corresponding to the current trial.
[path_root, path_dataset, ...
    marker_list, calc_marker_list, ...
    path_model, file_model ...
    name_emg_muscles, ...
    name_moco_reserves, name_moco_muscle_activation, name_moco_muscle_forces] = data_pipeline_setup();

cd([path_dataset, path_subject]);

% Load the trial's data.
data_c3d = btk_loadc3d([path_dataset, path_subject, file_trial, '.c3d'], 10);

% Copy the fields of c3d file into data received as parameter.
fields = fieldnames(data_c3d);
% Loop over each marker name
for idx_field = 1:length(fields)
    data.(fields{idx_field}) = data_c3d.(fields{idx_field});
end

% ----------------------------------------------------------------------- %

% hold on;
% plot(data.fp_data.GRF_data(1).F(:,3));
% plot(data.fp_data.GRF_data(2).F(:,3));
% hold off;

% ----------------------------------------------------------------------- %

% sort the C3D file so we know what is Marker data and what is calculated
% data -- THIS STEP ISN'T REALLY REQUIRED FOR THIS DATASET WHERE ALL
% THE MARKERS ARE REQUIRED
data.marker_data = btk_sortc3d(data.marker_data, marker_list, calc_marker_list);

% Use a few frames of the static trial in order to scale the template model
% to the data.
data.Start_Frame = 11;
data.End_Frame =   20;

% now do conversions to TRC files using btk_c3d2trc
TRCFileName = [file_trial '_trc.trc'];
GRFFileName = [file_trial '_grf.sto'];

% Switch to 'off' to turn animation off.
data = btk_c3d2trc(data, 'off', TRCFileName, GRFFileName);

% Define the standard config files OpenSim's scaling tool.
ModelFile =          [path_model file_model '.osim'];
ScaleTasksFile =     [path_model file_model '_Scale_Tasks.xml'];
MeasurementSetFile = [path_model file_model '_Scale_MeasurementSet.xml'];

setup_scale( ...
    'data',data, ...
    'ModelFile',ModelFile, ...
    'ScaleTasksFile',ScaleTasksFile, ...
    'MeasurementSetFile',MeasurementSetFile, ...
    'PreserveMass','true');

com = ['opensim-cmd run-tool ' data.Name '_Setup_Scale.xml'];
system(com);
