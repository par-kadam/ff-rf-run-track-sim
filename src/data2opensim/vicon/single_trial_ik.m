function single_trial_ik()

    import org.opensim.modeling.*;

    % The path to the 'lib' directory of the this repository.
    addpath(genpath('../lib/'));


    % now do conversions to TRC files using btk_c3d2trc
    name_trial = 'sub03_strist_cond00000_speed00000_0003'
    TRCFileName = [name_trial '_trc.trc'];
    GRFFileName = [name_trial '_grf.sto'];
    % GRFFileXML  = [name_trial '_' sprintf('%04d', idx_cyc) '_grf.xml'];
    % EMGFileName = [name_trial '_' sprintf('%04d', idx_cyc) '_emg.sto'];


    data_trc = btk_c3d2trc(data, 'off', TRCFileName, GRFFileName);

    % Export EMG data to STO file using c3d2emg function
    % data_sto_emg = c3d2emg(data_emg, data_trc, EMGFileName, 1000);

    path_task_info = '/home/lisca/biomec/src/data2opensim/vicon/model/';
    path_dataset = '/home/lisca/biomec/data/20250306_dataset/';

    % define the standard OpenSim files for IK
    ModelFile    = ['sub03_SCALED.osim'];
    IKTasksFile  = ['vicon2opensim2moco_IK_Tasks.xml'];
    IKOutputFile = ['sub03_strist_cond00000_speed00000_0003_ik.sto'];

    cd([path_dataset, 'sub03', '/', 'ses20250306', '/']);

    setup_InverseKinematics( ...
        'data',        TRCFileName, ...
        'ModelFile',   ModelFile,...
        'IKTasksFile', IKTasksFile, ...
        'Accuracy',    0.00005, ...
        'OutputFile',  IKOutputFile);

    % run the ik tool from the command line
    com = ['opensim-cmd run-tool ' 'sub03_Setup_InverseKinematics.xml'];
    [status, cmdout] = system(com);
