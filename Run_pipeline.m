function Run_pipeline(txtpath, param)


    % Read in paths/ directories from summary_dirs.txt
    DirList = readtext(txtpath);
    DirList = DirList(~cellfun('isempty', DirList));
    nDir = length(DirList);
    errortxt = fullfile(pwd, 'Error_folders.txt');
    errorID = fopen(errortxt,'w');

    % Go over each folder to do the analysis
    
    for i = 1:nDir
        
        try

            cur_folder = DirList{i};

            %Set the current directory (data to analyze)
            cd(cur_folder)

            %Set parameters
            %For SVD analysis based on different initial dimensions, default = 0
            param.rechooseIniDim = 0;
            %Whether discard frames that contain large motions, recommend set to 0
            %When doing seed-based correlation analysis, those frames will be omitted,
            param.moveAssessFlag = 0; 
            %Flag = 0,1,2 (spontaneous, auditory stimuli, wavelength switching)
            param.flag = 0;
            %Spatial downsampling factor
            param.spacialFactor = 2;
            %Motion correction methods, 0:NoRMCorre, 1: dft
            param.motionCorrMethod = 0;
            %Lower limit for number of seeds
            param.total_seeds = 500;
            %Whether to use GPU, default = 0
            param.GPU_flag = 0;
            %Initial dimension to preserve during SVD analysis, default = 1
            param.iniDim = 1;
            %Partial correlation flag: 0 = regular, 1 = partial, 2 = partial correlation using lowest 1% std trace
            param.mean_flag = 0;
            %Timelag for correlation analysis, default = 0
            param.timelag = 0;
            %Threshold for generating connected components
            param.CCthreshold = 0.3;

            %Run the pipeline
            audPipe(param)

            %Run the seed-based correlation (regular)
            %byPassPreProcessing(3, param);

            %Run the seed-based correlation (partial)
            param.mean_flag = 1;
            byPassPreProcessing(3, param);

            %WAIP;
        catch
            warning('Unexpected error occurred, saving the folder path...')
            cur_folder = DirList{i};
            fprintf(errorID, [cur_folder, '\n']);
        end
            
    end
    
    fclose(errorID);

end





