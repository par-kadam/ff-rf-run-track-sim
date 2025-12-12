function data = c3d2emg(data_emg, data, EMGFileName, frequency)
% C3D2EMG - Converts EMG data from a structure to OpenSim .sto format
%
% Function to export EMG data to the STO file format for OpenSim, 
% similar to how btk_c3d2trc exports marker and force data
%
% This function normalizes all EMG values to each muscle's maximal 
% activation amplitude value stored in data_emg.mvic_maximum.(muscle_name)
%
% INPUTS:
%   data_emg      - Structure containing EMG data with both raw data in 'clean'
%                   and normalization values in 'mvic_maximum'
%   data          - Structure containing frame range information
%   EMGFileName   - Output filename for the .sto file
%
% OUTPUT:
%   data          - Updated structure (returned for consistency with btk_c3d2trc)
%                   Creates an EMG .sto file for OpenSim
%
% Based on btk_c3d2trc function by Glen Lichtwark

% Define the start and end frame for analysis
if ~isfield(data,'Start_Frame') || ~isfield(data,'End_Frame')
    error('Start_Frame and End_Frame must be defined in the data structure');
end

start_frame = data.Start_Frame;
end_frame = data.End_Frame;
nrows = end_frame - start_frame + 1;

%Just replacing the below if statement with a predetermined common sample rate denoted in filter_and_resample function
% Get frequency information for EMG data
%if isfield(data_emg, 'analogs_info') && isfield(data_emg.analogs_info, 'frequency')
%    frequency = data_emg.analogs_info.frequency;
%else
%    % Default to 1000Hz if frequency not found
%    frequency = 1000;
%    disp('Warning: EMG frequency not found, defaulting to 1000Hz');
%end

% Create time vector
emg_time = (0:nrows-1)' / frequency;

% Get muscle names (column labels)
% NOTE: This section extracts muscle names from the EMG data structure
% If you want to use custom column labels, replace this section
% ====================================================================
if isfield(data_emg, 'clean')
    name_muscles = fieldnames(data_emg.clean);
    % If you want to use custom column labels, replace the above line with:
    % name_muscles = {'your', 'custom', 'muscle', 'names', 'here'};
    % and ensure these match the keys in data_emg.clean
else
    error('No clean EMG data found in data_emg.clean');
end
% ====================================================================

% Initialize the EMG data matrix with time column
emg_data_out = emg_time;

% Create header for columns
emg_header = 'time\t';

% For each muscle, extract the data and add to the output matrix
for idx_mus = 1:length(name_muscles)
    muscle_name = name_muscles{idx_mus};
    
    % Get the raw EMG data for the specified frame range
    if isfield(data_emg.clean, muscle_name)
        raw_emg = data_emg.clean.(muscle_name)(start_frame:end_frame);
        
        % Get MVIC normalization value (maximal activation amplitude)
        if isfield(data_emg, 'mvic_maximum') && isfield(data_emg.mvic_maximum, muscle_name)
            mvic_max = data_emg.mvic_maximum.(muscle_name);
            
            % Check if the normalization value is valid
            if mvic_max > 0
                % Normalize the EMG data by dividing by the maximum value
                normalized_emg = raw_emg / mvic_max;
            else
                % If MVIC is zero or negative (invalid), use raw data
                normalized_emg = raw_emg;
                disp(['Warning: Invalid MVIC maximum (', num2str(mvic_max), ') for ' muscle_name ', using raw values']);
            end
        else
            % If no normalization value found, use raw data
            normalized_emg = raw_emg;
            disp(['Warning: No MVIC maximum found for ' muscle_name ', using raw values']);
        end
        
        % Add the EMG data to the output matrix
        emg_data_out = [emg_data_out, normalized_emg];
        
        % Add the muscle name to the header
        emg_header = [emg_header, muscle_name, '\t'];
    else
        disp(['Warning: Muscle ' muscle_name ' not found in EMG data']);
    end
end

% Remove the last tab from the header and add a newline
emg_header = [emg_header(1:end-2), '\n'];

% Format string for data output (tab-delimited with precision)
emg_format = '';
for i = 1:size(emg_data_out, 2)
    emg_format = [emg_format, '%20.6f\t'];
end
emg_format = [emg_format(1:end-2), '\n'];

disp('Writing emg.sto file...')

% Store the filename in the data structure 
data.EMGFileName = EMGFileName;

% Open file for writing
fid = fopen(EMGFileName, 'w');

if fid == -1
    error(['Could not open file for writing: ' EMGFileName]);
end

% Write the header information (5 lines as specified)
fprintf(fid, 'name %s\n', EMGFileName);
fprintf(fid, 'datacolumns %d\n', size(emg_data_out, 2));  % total # of datacolumns
fprintf(fid, 'datarows %d\n', length(emg_time)); % number of datarows
fprintf(fid, 'range %f %f\n', emg_time(1), emg_time(end)); % range of time data
fprintf(fid, 'endheader\n');
fprintf(fid, emg_header);

% Write the EMG data (write entire matrix at once, similar to btk_c3d2trc)
fprintf(fid, emg_format, emg_data_out');

% Close the file
fclose(fid);

disp(['EMG data written to ' EMGFileName]);
end