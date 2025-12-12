function data = opensim_walking_pipeline(file)

% This function will create the appropriate files to run an OpenSim
% simulation sequence for the walking data - use example model and data in
% Models\Plug_In_Gait_Lower_Limb

clc

% change this path to path where model and model files exist
% osim_path = '\\vmware-host\Shared Folders\Documents\MATLAB\Opensim\MatlabOpensimPipelineExample\Plug_In_Gait_Lower_Limb\';

% The 'install' directory of the 'biomec' directory structure.
addpath(genpath('../../../../install/btk/'));

% The path to the 'vicon2opensim/tools/' of the 'athlete' repository.
addpath(genpath('./tools/'));

% The path to the *.osim model.
osim_path = [pwd, '/model/'];

% define the model name - the task files etc must start with the same name
model = 'lower_limb';
    
if nargin == 1
if isempty(fileparts(file))
    pname = cd;
    pname = [pname '\'];
    fname = file;
else [pname, name, ext] = fileparts(file);
    fname = [name ext];
end
else [fname, pname] = uigetfile('*.c3d', 'Select C3D file');
end

cd (pname);

% load the c3d file using BTK
data = btk_loadc3d([pname, fname], 10);

marker_list = {'LASI'; 'RASI'; 'LPSI'; 'RPSI'; 'LTHI'; 'LKNE'; 'LTIB'; 'LANK'; 'LHEE';...
     'LTOE'; 'RTHI'; 'RKNE'; 'RTIB'; 'RANK'; 'RHEE'; 'RTOE'};
                     
calc_marker_list = {'LTIO'; 'LFEO'; 'LFEP'; 'RFEP'; 'RFEO'; 'RTIO'};
 
% sort the C3D file so we know what is Marker data and what is calculated
% data -- THIS STEP ISN'T REALLY REQUIRED FOR THIS DATASET WHERE ALL
% THE MARKERS ARE REQUIRED
data.marker_data = btk_sortc3d(data.marker_data,marker_list,calc_marker_list);

% IMPORTANT: A STATIC TRIAL MUST FIRST BE ANALYSED AND THIS FILENAME MUST
% CONTAIN THE STRING 'static' or 'cal'. IF THIS FILE EXISTS AND IS SELECTED
% THE STATIC TRIAL IS ANALYSED TO SCALE THE MODEL (ELSE STATEMENT BELOW).
% IF A DYNAMIC TRIAL IS SELECTED THE SCALED MODEL MUST ALREADY EXIST SO THE
% MODEL CAN RUN CORRECTLY (IF STATEMENT BELOW) - change to '01' or whatever
% else you require
if isempty(strfind(lower(fname),'static')) && isempty(strfind(lower(fname),'cal'))
    
    %find when feet are on each forceplate and define the trial events as
    %the first foot contact and the final foot contact
    E = [];
    for i = 1:length(data.fp_data.GRF_data)
        clear a
        a = find(data.fp_data.GRF_data(i).F(:,3) > 0.01*max(data.fp_data.GRF_data(i).F(:,3)));
        if ~isempty(a)
            P = round((a(1)*data.marker_data.Info.frequency/data.fp_data.Info(i).frequency)):round((a(end)*data.marker_data.Info.frequency/data.fp_data.Info(i).frequency));
            E(1) = min([min(E) P(1)]);
            E(2) = max([max(E) P(end)]);
        end
    end    
    
    % define start and end frame from the events to write the appropriate TRC
    % and MOT files for the OpenSim simulations
    data.Start_Frame = E(1);
    data.End_Frame = E(2);

    % now do conversions to TRC files using btk_c3d2trc
    data = btk_c3d2trc(data,'on');    % change to 'off' to turn animation off

    % define the standard OpenSim files for IK
    ModelFile = [data.Name '_SCALED.osim'];
    IKTasksFile = [osim_path model '_IK_Tasks.xml'];

    setup_InverseKinematics('data',data,'ModelFile',ModelFile,...
        'IKTasksFile',IKTasksFile,'Accuracy',0.00005);
    
    % run the ik tool from the command line
    com = ['opensim-cmd run-tool ' data.Name '_Setup_InverseKinematics.xml'];
    system(com);
    
    % load the data from MOT file and tie kinematic data to data structure
    D = load_sto_file([data.TRC_Filename(1:end-4) '_ik.mot']);
    fnames = fieldnames(D);

    clear a
    a = [strfind(lower(fnames),'tx') strfind(lower(fnames),'ty') strfind(lower(fnames),'tz') strfind(lower(fnames),'pelvis')];
    for i = 1:size(a,1)
        if isempty([a{i,1:3}])
            data.Kinematics.Angles.(fnames{i}) = D.(fnames{i});
        elseif ~isempty([a{i,4}])
            data.Kinematics.Angles.(fnames{i}) = D.(fnames{i});
        else data.Kinematics.Markers.(fnames{i}) = D.(fnames{i});
        end
    end
        
    % assign forces to the specific segments using 'assign_force' function
    % the output creates a cell array with the name of the external force
    % (ExForce) and the body that is assigned to that force (ApBodies)
    data = assign_forces(data,{'RHEE','LHEE'},{'calcn_r','calcn_l'},0.25);
    
    % define the standard OpenSim files for ID   
    MOTFile = [data.TRC_Filename(1:end-4) '_ik.mot'];
    GRFFile = [data.TRC_Filename(1:end-4) '_grf.mot'];
    GRFFileXML = [data.TRC_Filename(1:end-4) '_grf.xml'];
      
    % it is also necessary to generate an XML file containing the information
    % about which column is which and the bodies these forces are applied to
    % this is done with the grf2xml function
    grf2xml(data,'ExternalLoadNames',data.AssignForce.ExForce,'AppliedToBodies',data.AssignForce.ApBodies,...
        'GRFFile',GRFFile,'MOTFile',MOTFile,'LowPassFilterForKinematics',15,...
        'OutputFile',GRFFileXML);

    % make setup file for inverse dynamics
    setup_InverseDynamics('data',data, 'ModelFile',ModelFile,'MOTFile',MOTFile,...
        'GRFFile',GRFFileXML,'LowPassFilterForKinematics',15); 

    % call ID tool from the command line using setup file
    com = ['opensim-cmd run-tool ' data.Name '_Setup_InverseDynamics.xml'];
    system(com);
    
    clear D
    % load the data from ID MOT file and tie kinetic data to data structure
    D = load_sto_file([MOTFile(1:end-4) '_id.sto']);
    data.Kinetics = D;

    save([data.TRC_Filename(1:end-4) '.mat']');

    clear D
    
else
    
    data.Start_Frame = 2;
    data.End_Frame = 10;
    
    % the mass of the subject isn't in the c3d file, but can be calculated
    % from the static trial by using the force applied through the
    % forceplat that the subject stands on (forceplate 1)    
    data.Mass = mean(data.fp_data.GRF_data(1).F(:,3))/9.81;
    
    % convert the static data into a TRC file for the purposes of scaling
    data = btk_c3d2trc(data,'off'); % animation off
    
    % define the starndard OpenSim files for scaling
    ModelFile = [osim_path model '.osim'];
    MeasurementSetFile = [osim_path model '_Scale_MeasurementSet.xml'];
    ScaleTasksFile = [osim_path model '_Scale_Tasks.xml'];

    setup_scale('data',data,'ModelFile',ModelFile,'ScaleTasksFile',ScaleTasksFile,...
        'MeasurementSetFile',MeasurementSetFile,'PreserveMass','true');
    
    com = ['opensim-cmd run-tool ' data.Name '_Setup_Scale.xml'];
    system(com);
end





