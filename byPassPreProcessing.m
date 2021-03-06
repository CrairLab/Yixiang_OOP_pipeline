function byPassPreProcessing(id,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function bypass pre-processing (assuming it has been done) and use
% filtered matrix to do further analysis: 1. re-calculate connected
% components. 2. Make pseudo-color movies 3. Seed-based correlation 
% 4-5 connectivity/dF over F -based K-means analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All previous record is saved on Github
% Visit https://github.com/CrairLab/Yixiang_OOP_pipeline for more info
% Author: yixiang.wang@yale.edu
% Latest update:
% R13 02/10/20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   id       determine which analysis to run
%   parSam    parameters for analysisS
    
    %Detect/create file list
    if isempty(dir(('files.txt')))
        Integration.fileDetector()
    end
    filelist = readtext('files.txt',' ');
    nmov = size(filelist,1);
    
    if ~exist('param','var')
        if isempty(dir(('parameter.mat')))
            param.rechooseIniDim = 0;
            param.moveAssessFlag = 0;
            param.flag = 0;
            param.spacialFactor = 2;
            param.total_seeds = 1000;
            param.GPU_flag = 0;
            param.iniDim = 1;
            param.mean_flag = 0;
            param.timelag = 0;
            param.motionCorrMethod = 0;
            param.CCthreshold = 0.3;
        else
            load('parameter.mat')
        end
    end
    
    switch id
        case 1 
        %Renew connected components
            for f = 1:nmov
                [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix('filtered',f);
                %Chop the matrix to contain only roi
                ppA_roi = movieData.focusOnroi(curLoad.A_dFoF);
                %Renew connected components
                Integration.renewCC(ppA_roi,param.CCthreshold, outputFolder,filename)
                disp(['Preprocessing done: ' filename]);
            end
        case 2
        %Make pseudo color movie
            for f = 1:nmov
                try
                    [curLoad,~,filename]  = Integration.readInSingleMatrix(['filtered' movTag], f);
                    if isempty(curLoad)
                        disp(['Movement tag = ' movTag]);
                        disp('Can not read in matrix with this tag, try a new tag...')
                        movTag = '';
                        [curLoad,~,filename]  = Integration.readInSingleMatrix(['filtered' movTag], f);
                    end
                catch
                    warning('Unable to read in filtered matrix!')
                end
                
                %Chop the matrix to contain only roi
                ppA_roi = movieData.focusOnroi(curLoad.A_dFoF);
                %make pseudo color movie
                movieData.makePseudoColorMovie(ppA_roi,filename(1:length(filename)-4))
                disp(['Preprocessing done: ' filename]);
            end
        case 3
        %Do Seed-based correlation 
            A_all = [];
            Avg_out_all = [];
            if ~isfield(param, 'rechooseIniDim')
                param.rechooseIniDim = 0;
            end
            if param.rechooseIniDim == 0  %Loading directly from filtered.mat         
                for f = 1:nmov
                    try
                        movTag = 'dsc';
                        %Detect if there is a filtered matrix with discarding process
                        [curLoad,~,~]  = Integration.readInSingleMatrix(['filtered' movTag], f);
                        
                        %Detect if there exists info of background dFoF 
                        Avg_out_dFoF = Integration.GenerateOutsideAverageFromScratch(f);
                        if isempty(Avg_out_dFoF)
                            warning('Can not detect/generate of background averaged trace!')
                        end                        
                        
                        if isempty(curLoad)
                            %If not, try to detect filtered matrix wo discarding process
                            disp(['Movement tag = ' movTag]);
                            disp('Can not detect matrix with this tag, try a new tag...')
                            movTag = '';
                            disp(['Movemetn tag = ' movTag]);
                            [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix(['filtered' movTag], f);
                            
                            if isempty(curLoad)
                                disp('Unable to detect pre-processed matrix!')
                            else
                                disp('Detect filtered matrix without discarding.')
                                disp('For seed-based correlation, apply frame discarding!')
                                %If detected filtered matrix wo discarding process, do frames discarding here
                                [moveAssessLoad,~,~]  = Integration.readInSingleMatrix('moveAssess', f);
                                if isempty(moveAssessLoad)
                                    disp('Unable to detect movAssess file!')
                                else
                                    [A_dFoF,~,~] = movieData.discardFrames(curLoad.A_dFoF, moveAssessLoad.NormTform_all);
                                    if size(A_dFoF,3) ~= size(curLoad.A_dFoF,3) 
                                        curLoad.A_dFoF = A_dFoF;
                                        checkname = [filename(1:length(filename)-4) '_filtereddsc.mat'];
                                        %save the frame-discarded filtered matrix
                                        save(fullfile(outputFolder,checkname), 'A_dFoF', '-v7.3')
                                    end
                                end
                            end
                        end                     
                    catch
                        warning('Unable to read in pre-processed matrix!')
                    end
                    
                    try
                        A_dFoF = curLoad.A_dFoF;
                        %If is not a spontaneous trail, do gross dF over F
                        if param.flag
                            disp('Do gross dF over F here!')
                            A_dFoF = movieData.grossDFoverF(A_dFoF);                       
                        end
                        clear curLoad
                        A_all = cat(3, A_all, A_dFoF);
                        Avg_out_all = cat(3, Avg_out_all, Avg_out_dFoF);
                    catch
                        disp('Unable to read A_dFoF!')
                        disp('')
                    end

                end
            else  %Reconstruct from SVD components
                for f = 1:nmov
                    [curLoad,~,~]  = Integration.readInSingleMatrix('SVD',f);
                    U = curLoad.U;
                    S = curLoad.S;
                    V = curLoad.V;
                    %Rechoose initial dimension
                    iniDim = param.rechooseIniDim;
                    sz = size(V);
                    V = reshape(V, [sz(1)*sz(2) sz(3)]);
                    A_rcs = V(:,iniDim:end)*S(iniDim:end,iniDim:end)*U(:,iniDim:end)';
                    A_rcs = reshape(A_rcs,[sz(1),sz(2),size(A_rcs,2)]);
                    A_z = zscore(A_rcs,1,3);
                    disp('Z-scored reconstructed matrix')
                    %Gaussian smoothing
                    A_dFoF = Integration.GauSmoo(A_z,1); %set sigma = 1
                    disp('Gaussian smoothing is done');
                    disp(' ')
                    %If detected filtered matrix wo discarding process, do frames discarding here
                    [moveAssessLoad,outputFolder,filename]  = Integration.readInSingleMatrix('moveAssess', f);
                    if isempty(moveAssessLoad)
                        disp('Unable to detect movAssess file!')
                    else
                        [A_dFoF,~,~] = movieData.discardFrames(A_dFoF, moveAssessLoad.NormTform_all);
                        checkname = [filename(1:length(filename)-4) '_filtereddsc.mat'];
                        %save the frame-discarded filtered matrix
                        save(fullfile(outputFolder,checkname), 'A_dFoF', '-v7.3')
                    end
                  
                    %prepare for seed-based correlation analysis
                    A_all = cat(3, A_all, A_dFoF);
                    Avg_out_all = [];
                    clear curLoad
                end
            end
            
            %See if length of background trace is consistent with the movie
            if size(Avg_out_all,3) ~= size(A_all)
                warning('Movie length and background trace length do not agree!')
                Avg_out_all = [];
            else
                Avg_out_all = Avg_out_all(:);
            end
            
            movieData.SeedBasedCorr_GPU(A_all,param.spacialFactor,param.total_seeds,param.GPU_flag,0,...
            param.mean_flag, Avg_out_all, param.timelag);
        case 4
        %Do connectivity K-means analysis 
            A_all = [];
            for f = 1:nmov
                    [curLoad,~,~]  = Integration.readInSingleMatrix('filtered',f);
                    A_all = cat(3, A_all, curLoad.A_dFoF);
            end
            curLoad = load('Correlation_Matrix.mat');
            [A_all,corrMatrix,~,~] = connectivityKmeans(A_all,curLoad.corrMatrix,0);
            save('AllMatrix.mat','A_all','-v7.3')
            save('CorrelationMatrix_all.mat','corrMatrix','-v7.3');
        %Do dFoF K-means analysis
        case 5
            A_all = [];
            for f = 1:nmov
                    [curLoad,~,~]  = Integration.readInSingleMatrix('filtered',f);
                    A_all = cat(3, A_all, curLoad.A_dFoF);
            end
            save('AllMatrix.mat','A_all','-v7.3')
            dFoFKmeans(A_all);
          
    end

    disp(['Processing done at:' pwd]);
    
    clearvars;  

end
%%

