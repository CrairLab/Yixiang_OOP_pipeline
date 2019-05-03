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
%R6 01/24/19 Changed function structure to optimize running efficiency and
%readability only compatible with Integration class R16+
%R7 01/31/19 Changed Ga_TH_A to A_dFoF to be compatible with new
%Integration class. Only compatible with Integration R19+
%R8 05/03/19 Modify seed-based correlation procedures. Only compatible with
%Integration R23 +.

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
    
    switch id
        case 1 
        %Renew connected components
            for f = 1:nmov
                [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix('filtered',f);
                %Chop the matrix to contain only roi
                ppA_roi = movieData.focusOnroi(curLoad.A_dFoF);
                %Renew connected components
                Integration.renewCC(ppA_roi,outputFolder,filename)
                disp(['Preprocessing done: ' filename]);
            end
        case 2
        %Make pseudo color movie
            for f = 1:nmov
                [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix('filtered',f);
                %Chop the matrix to contain only roi
                ppA_roi = movieData.focusOnroi(curLoad.A_dFoF);
                %make pseudo color movie
                movieData.makePseudoColorMovie(ppA_roi,filename(1:length(filename)-4))
                disp(['Preprocessing done: ' filename]);
            end
        case 3
        %Do Seed-based correlation 
            A_all = [];
            if param.moveAssessFlag
                movTag = 'dsc';
            else 
                movTag = '';
            end             
            
            if param.rechooseIniDim == 0           
                for f = 1:nmov
                    try
                        [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix(['filtered' movTag], f);
                    catch
                        disp(['Movement tag = ' movTag]);
                        disp('Can not read in matrix with this tag, try a new tag...')
                        movTag = ''
                        [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix(['filtered' movTag], f);
                    end
                    A_all = cat(3, A_all, curLoad.A_dFoF);
                end
            else
                for f = 1:nmov
                    [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix('SVD',f);
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
                    %prepare for seed-based correlation analysis
                    A_all = cat(3, A_all, A_dFoF);
                end
            end
            movieData.SeedBasedCorr_GPU(A_all,param.spacialFactor,param.total_seeds,param.GPU_flag,0);
        case 4
        %Do connectivity K-means analysis 
            A_all = [];
            for f = 1:nmov
                    [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix('filtered',f);
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
                    [curLoad,outputFolder,filename]  = Integration.readInSingleMatrix('filtered',f);
                    A_all = cat(3, A_all, curLoad.A_dFoF);
            end
            save('AllMatrix.mat','A_all','-v7.3')
            dFoFKmeans(A_all);
          
    end

    disp(['Processing done at:' pwd]);
    
    clearvars;  
end

%%

