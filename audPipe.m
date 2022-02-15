function audPipe(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auditory experiment analysis pipeline. Compatible with spontaneous
% activity preprocessing/seed-based correlation analysis. Auditory
% stimulation experiment. Wavelength switching experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All previous record is saved on Github
% Visit https://github.com/CrairLab/Yixiang_OOP_pipeline for more info
% Author: yixiang.wang@yale.edu
% Latest update:
% R10 09/11/19 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %If run on the HPC, use slurm to change the current directory
    try
        %currentFolder = 'E:\New folder';
        cd(currentFolder);
    catch
        disp('No such directory...Running on pwd...')
    end

    %Convert spike2 raw data to .mat files
    %Integration.Spike2Matlab(cd);
       
    %Detect movie/baphy/spike2/roi files 
    Integration.fileDetector();
    
    filelist = readtext('files.txt','\n');
    nmov = size(filelist,1);
    IdxInAll = cell(nmov,1);
    IdxInAll_1 = cell(nmov,1);
    IdxInAll_2 = cell(nmov,1);

    %Save current parameters
    save('parameter.mat','param');
    
    %Process each movie in parallel/sequentially
    try
        parfor f = 1:nmov
            disp('Load instance.mat file if existed...')
            [curLoad, ~, ~]  = Integration.readInSingleMatrix('instance',f);
            if ~isempty(curLoad)
                disp('Instance matrix detected, partially skip pre-processing...')
                curObj = curLoad.obj;
                curObj.prePipe(param);
            else
                disp('Did not detect instance.mat, try to generate obj from .tif files')
                [Idx,Idx1,Idx2] = Integration.processMovies(f,nmov,param);                   
                IdxInAll{f,1} = Idx;
                IdxInAll_1{f,1} = Idx1;
                IdxInAll_2{f,1} = Idx2;
            end
        end
    catch
        disp('Parallel running failed, try run the data sequentially!')
        parfor f = 1:nmov
            disp('Load instance.mat file if existed...')
            [curLoad, ~, ~]  = Integration.readInSingleMatrix('instance',f);
            if ~isempty(curLoad)
                disp('Instance matrix detected, partially skip pre-processing...')
                curObj = curLoad.obj;
                curObj.prePipe(param);
            else
                disp('Did not detect instance.mat, try to generate obj from .tif files')
                [Idx,Idx1,Idx2] = Integration.processMovies(f,nmov,param);                   
                IdxInAll{f,1} = Idx;
                IdxInAll_1{f,1} = Idx1;
                IdxInAll_2{f,1} = Idx2;
            end

        end
    end
        
    
    %Load exptparam(baphy) from saved object
    try
        [curLoad, ~, ~]  = Integration.readInSingleMatrix('instance',1);
        curObj = curLoad.obj;
        exptparam.PreStimSilence = curObj.PreStimSilence;
        exptparam.PostStimSilence = curObj.PostStimSilence;
        exptparam.Duration = curObj.Duration;
        exptparam.BlockDura = curObj.BlockDura;
    catch
        exptparam = [];
    end
        
    %Do averaging across all movies
    Integration.doAveragingAcrossMovies(flag,IdxInAll,IdxInAll_1,IdxInAll_2,exptparam);

end


