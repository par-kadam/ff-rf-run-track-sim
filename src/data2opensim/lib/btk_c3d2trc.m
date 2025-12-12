function data = btk_c3d2trc(data, anim, TRCFileName, GRFFileName)
% function btk_c3d2trc(file) OR
% function btk_c3d2trc(data)
%
% Function to convert data from a C3D file into the TRC and MOT file
% formats for OpenSim
%
% INPUT -   data - structure containing fields from from previously loaded
%               C3D file using btk_loadc3d.m
%           anim - animate 'on' or 'off' (default - 'on')
%           TRCFileName -
%           GRFFileName -
%
% OUTPUT -  data - structure containing the relevant data from the c3dfile
%                  Creates the TRC file and _grf.MOT file for OpenSim
%
% Written by Glen Lichtwark (University of Queensland)
% Updated September 2012

import org.opensim.modeling.*;

[pname, name, ext] = fileparts(data.marker_data.Filename);
if ispc
    pname = [pname '\'];
else 
    pname = [pname '/'];
end

%%
% if the mass, height and name aren't present then presribe - it is
% preferrable to have these defined in the data structure before running 
% this function - btk_loadc3d should try and do this for vicon data
if ~isfield(data,'Mass')
    data.Mass = 75;
end

if ~isfield(data,'Height')
    data.Height = 1750;
end

if ~isfield(data,'Name')
    data.Name = 'NoName';
end

%% define the start and end frame for analysis as first and last frame unless 
% this has already been done to change the analysed frames
if ~isfield(data,'Start_Frame')
    data.Start_Frame = 1;
    data.End_Frame = data.marker_data.Info.NumFrames;
end

%%
% define some parameters 
nrows = data.End_Frame-data.Start_Frame+1;
nmarkers = length(fieldnames(data.marker_data.Markers));

data.time = (1/data.marker_data.Info.frequency:1/data.marker_data.Info.frequency:(data.End_Frame-data.Start_Frame+1)/data.marker_data.Info.frequency)';

nframe = 1:nrows;

% anim the trial if animation = on
if strcmp(anim,'on')
    data.marker_data.First_Frame = data.Start_Frame;
    data.marker_data.Last_Frame = data.End_Frame;
    if isfield(data,'fp_data')
        btk_animate_markers(data.marker_data, data.fp_data, 5)
    else 
        btk_animate_markers(data.marker_data)
    end
end

%%
% we need to reorder the lab coordinate system to match that of the OpenSim
% system --> SKIP THIS STEP IF LAB COORDINATE SYSTEM IS SAME AS MODEL
% SYSTEM
markers = fieldnames(data.marker_data.Markers); % get markers names

if strcmp(data.marker_data.Info.units.ALLMARKERS,'mm')
    p_sc = 1000;
    data.marker_data.Info.units.ALLMARKERS = 'm';
else 
    p_sc = 1;
end

% go through each marker field and re-order from X Y Z to Y Z X
for i = 1:nmarkers
   data.marker_data.Markers.(markers{i}) =  [data.marker_data.Markers.(markers{i})(:,2)...
       data.marker_data.Markers.(markers{i})(:,3) data.marker_data.Markers.(markers{i})(:,1)]/p_sc;
end

%%
% now we need to make the headers for the column headings for the TRC file
% which are made up of the marker names and the XYZ for each marker

% first initialise the header with a column for the Frame # and the Time
% also initialise the format for the columns of data to be written to file
dataheader1 = 'Frame#\tTime\t';
dataheader2 = '\t\t';
format_text = '%i\t%2.4f\t';
% initialise the matrix that contains the data as a frame number and time row
data_out = [nframe; data.time'];

% now loop through each maker name and make marker name with 3 tabs for the
% first line and the X Y Z columns with the marker numnber on the second
% line all separated by tab delimeters
% each of the data columns (3 per marker) will be in floating format with a
% tab delimiter - also add to the data matrix
for i = 1:nmarkers
    dataheader1 = [dataheader1 markers{i} '\t\t\t'];    
    dataheader2 = [dataheader2 'X' num2str(i) '\t' 'Y' num2str(i) '\t'...
        'Z' num2str(i) '\t'];
    format_text = [format_text '%f\t%f\t%f\t'];
    % add 3 rows of data for the X Y Z coordinates of the current marker
    % first check for NaN's and fill with a linear interpolant - warn the
    % user of the gaps
    clear m
    m = find(isnan(data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),1)) >0);
    if ~isempty(m);
        clear t d
        disp(['Warning -' markers{i} ' data missing in parts. Frames ' num2str(m(1)) '-'  num2str(m(end))])
        t = time;
        t(m) = [];
        d = data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:);
        d(m,:) = [];
        data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:) = interp1(t,d,time,'linear','extrap');
    end
    data_out = [data_out; data.marker_data.Markers.(markers{i})((data.Start_Frame:data.End_Frame),:)'];
end
dataheader1 = [dataheader1 '\n'];
dataheader2 = [dataheader2 '\n'];
format_text = [format_text '\n'];

disp('Writing trc file...') 

%Output marker data to an OpenSim TRC file

data.TRC_Filename = [pname TRCFileName];

%open the file
fid_1 = fopen(data.TRC_Filename,'w');

% first write the header data
fprintf(fid_1,'PathFileType\t4\t(X/Y/Z)\t %s\n',data.TRC_Filename);
fprintf(fid_1,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid_1,'%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', data.marker_data.Info.frequency, data.marker_data.Info.frequency, nrows, nmarkers, data.marker_data.Info.units.ALLMARKERS, data.marker_data.Info.frequency,data.Start_Frame,data.End_Frame); 
fprintf(fid_1, dataheader1);
fprintf(fid_1, dataheader2);

% then write the output marker data
fprintf(fid_1, format_text,data_out);

% close the file
fclose(fid_1);

disp('Done.')

%%
% Write motion file containing GRFs

disp('Writing grf.mot file...')

if isfield(data,'fp_data')
    % assume that all force plates are collected at the same frequency!!!
    F = data.fp_data.Info(1).frequency/data.marker_data.Info.frequency;

    fp_time = 1/data.marker_data.Info.frequency:1/data.fp_data.Info(1).frequency:(F*(data.End_Frame-data.Start_Frame+1))/data.fp_data.Info(1).frequency;
    
    % initialise force data matrix with the time array and column header
    force_data_out = fp_time';
    force_header = 'time\t';
    force_format = '%20.6f\t';
    
    % go through each marker field and re-order from X Y Z to Y Z X and place
    % into data array and add data to the force data matrix --> also need to
    % divide by 1000 to convert to mm from m if necessary and Nmm to Nm 
    % these are the conversions usually used in most motion analysis systems
    % and if they are different, just change the scale factor below to p_sc 
    % value for the M data to 1. It should however get this from the file.
       
    for i = 1:length(data.fp_data.FP_data)
        
        % reoder data so lab coordinate system to match that of the OpenSim
        % system
       data.fp_data.GRF_data(i).P =  [data.fp_data.GRF_data(i).P(:,2)/p_sc ...
           data.fp_data.GRF_data(i).P(:,3)/p_sc data.fp_data.GRF_data(i).P(:,1)/p_sc];
       data.fp_data.GRF_data(i).F =  [data.fp_data.GRF_data(i).F(:,2) ...
           data.fp_data.GRF_data(i).F(:,3) data.fp_data.GRF_data(i).F(:,1)];
       data.fp_data.GRF_data(i).M =  [data.fp_data.GRF_data(i).M(:,2) ...
           data.fp_data.GRF_data(i).M(:,3) data.fp_data.GRF_data(i).M(:,1)]/p_sc;
       
       % do some cleaning of the COP before and after contact
       b = find(abs(diff(data.fp_data.GRF_data(i).P(:,3)))>0);
       if ~isempty(b)
           for j = 1:3
                data.fp_data.GRF_data(i).P(1:b(1),j) = data.fp_data.GRF_data(i).P(b(1)+1,j);
                data.fp_data.GRF_data(i).P(b(end):end,j) = data.fp_data.GRF_data(i).P(b(end)-1,j);
           end
       end
       
       % define the period which we are analysing
       K = (F*data.Start_Frame):1:(F*data.End_Frame);
          
       % add the force, COP and moment data for current plate to the force matrix 
       force_data_out = [force_data_out data.fp_data.GRF_data(i).F(K,:) data.fp_data.GRF_data(i).P(K,:) data.fp_data.GRF_data(i).M(K,:)];
       % define the header and formats
       force_header = [force_header num2str(i) '_ground_force_vx\t' num2str(i) '_ground_force_vy\t' num2str(i) '_ground_force_vz\t'...
           num2str(i) '_ground_force_px\t' num2str(i) '_ground_force_py\t' num2str(i) '_ground_force_pz\t' ...
           num2str(i) '_ground_torque_x\t' num2str(i) '_ground_torque_y\t' num2str(i) '_ground_torque_z\t'];
       force_format = [force_format '%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t'];
       
    end
    
    force_header = [force_header(1:end-2) '\n'];
    force_format = [force_format(1:end-2) '\n'];
    
    % assign a value of zero to any NaNs
    force_data_out(logical(isnan(force_data_out))) = 0;
    
    data.GRF_Filename = [pname GRFFileName];
    
    fid_2 = fopen(data.GRF_Filename,'w');
    
    % write the header information
    fprintf(fid_2,'name %s\n', data.GRF_Filename);
    fprintf(fid_2,'datacolumns %d\n', size(force_data_out,2));  % total # of datacolumns
    fprintf(fid_2,'datarows %d\n',length(fp_time)); % number of datarows
    fprintf(fid_2,'range %f %f\n',fp_time(1),fp_time(end)); % range of time data
    fprintf(fid_2,'endheader\n');
    fprintf(fid_2,force_header);
    
    % write the data
    fprintf(fid_2,force_format,force_data_out');
    
    fclose(fid_2);
    
    % % Load the GRF file and save it in a new format.
    % grf_data_table = TimeSeriesTable(data.GRF_Filename);
    % grf_data       = osimTableToStruct(grf_data_table);
    % 
    % grf_data_table_metadata_keys = grf_data_table.getTableMetaDataKeys();
    % 
    % % Create a new table with the same data;
    % grf_data_table_new = osimTableFromStruct(grf_data);
    % for i = 0:length(grf_data_table_metadata_keys) - 1
    %     key   = grf_data_table_metadata_keys.get(i);
    %     value = grf_data_table.getTableMetaDataString(key);
    % 
    %     grf_data_table_new.addTableMetaDataString(key, value);
    % end
    % 
    % STOFileAdapter.write(grf_data_table_new, data.GRF_Filename);
    
    disp('Done.');
else 
    disp('No force plate information available.')
end
