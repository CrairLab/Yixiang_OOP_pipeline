function byPassPreProcessing(id,param)
%This function bypass pre-processing (assuming it has been done) and use
%filtered matrix to do further analysis
%R1 07/11/18 can do renewCC, seed-based correlation, or generate pseudo color
%movie based on the filtered matrix.
%R2 07/13/18 Require methods from Integration and Movie classes 
%R3 10/02/18 Move the renewCC functino to Integration class, compatiable
%with Integration class R11 or higher
%R4 10/24/18 Feed in param to be compatiable with R15 movieData class,
%which allow GPU computing when doing seed-based correlation maps
%R5 01/21/19 Changed function structure, allow reconstruction of matrix
%based on previously saved SVD.mat file and rechoose initial dimension when
%doing reconstruction

%Inputs:
%id       determine which function to run
%(id == 1, renewCC; id == 2, pseudo-clolor movie; id == 3 seed-based corr)
%

    if ~exist('files.txt','file')
        Integration.fileDetector()
    end
    filelist = readtext('files.txt',' ');
    nmov = size(filelist,1);
    A_all = [];
    
    if ~exist('param','var')
        param.spacialFactor = 2;
        param.total_seeds = 500;
        param.GPU_flag = 0;
        param.rechooseIniDim = 1;
    end


    for f = 1:nmov
        cur_Names = Names(f,0);
        filename = cur_Names.filename;
        currentFolder = pwd;
        outputFolder = fullfile(currentFolder,cur_Names.outputFolder);
        %If param.rechooseIniDim == 0 don't change the initial dimension
        if param.rechooseIniDim == 0
            checkname = [filename(1:length(filename)-4) '_filtered.mat'];               
                if exist(fullfile(outputFolder,checkname),'file')
                    %Check whether pre-processing has been done before
                    disp('Filtered matrix detected, loading .mat file...')
                    curLoad = load(fullfile(outputFolder,checkname));
                    if id == 3
                        A_all = cat(3, A_all, curLoad.Ga_TH_A);
                    else
                        chooseAction(id,Ga_TH_A,outputFolder,filename)
                    end
                    disp('')
                else
                    disp('')
                    disp('No filetered matrix detected!!!')
                    disp('')
                end
        else
            checkname = [filename(1:length(filename)-4) '_SVD.mat'];
            if exist(fullfile(outputFolder,checkname),'file')
                %Check whether pre-processing has been done before
                disp('SVD matrix detected, loading .mat file...')
                curLoad = load(fullfile(outputFolder,checkname),'U','S','V');
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
                Ga_TH_A = Integration.GauSmoo(A_z,2); %set sigma = 2
                disp('Gaussian smoothing is done');
                disp(' ')
                %prepare for seed-based correlation analysis
                A_all = cat(3, A_all, Ga_TH_A);
            else
                disp('')
                disp('No svd matrix detected!!!')
                disp('')
            end
        end

    end
    
    if id == 3
        movieData.SeedBasedCorr_GPU(A_all,param.spacialFactor,param.total_seeds,param.GPU_flag);
    end
    disp(['Processing done at:' pwd]);
    
    clearvars;  
end


%%
function filelist = fileDetector()
%Detect the filtered matrix with pre-fix 'AveragedMatrix_'
    temp_info = dir;
    filelist = {};
    n = 0;
    for i = 1:size(temp_info,1)  
        if regexp(temp_info(i,1).name,'AveragedMatrix_')
            n = n+1;
            filelist{n,1} = temp_info(i,1).name;
        end
    end
end

%%
function chooseAction(id,A_input,outputFolder,filename)
%Choose different actions specified by id number
    switch id
            case 1
                %Chop the matrix to contain only roi
                ppA_roi = movieData.focusOnroi(A_input);
                %Renew connected components
                Integration.renewCC(ppA_roi,outputFolder,filename)
                disp(['Preprocessing done: ' filename]);
            case 2
                %Chop the matrix to contain only roi
                ppA_roi = movieData.focusOnroi(A_input);
                %make pseudo color movie
                movieData.makePseudoColorMovie(ppA_roi,filename(1:length(filename)-4))
                disp(['Preprocessing done: ' filename]);
    end   
end