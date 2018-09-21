function byPassPreProcessing(id)
%This function bypass pre-processing (assuming it has been done) and use
%filtered matrix to do further analysis
%R1 07/11/18 can do renewCC, seed-based correlation, or generate pseudo color
%movie based on the filtered matrix.
%R2 07/13/18 Require methods from Integration and Movie classes 
%
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


    for f = 1:nmov
        cur_Names = Names(f,0);
        filename = cur_Names.filename;
        currentFolder = pwd;
        outputFolder = fullfile(currentFolder,cur_Names.outputFolder);
        checkname = [filename(1:length(filename)-4) '_filtered.mat'];               
            if exist(fullfile(outputFolder,checkname),'file')
                %Check whether pre-processing has been done before
                disp('Filtered matrix detected, loading .mat file...')
                load(fullfile(outputFolder,checkname));
                
                switch id
                    case 1
                        %Chop the matrix to contain only roi
                         ppA_roi = movieData.focusOnroi(Ga_TH_A);
                        %Renew connected components
                        renewCC(ppA_roi,outputFolder,filename)
                        disp(['Preprocessing done: ' filename]);
                    case 2
                        %Chop the matrix to contain only roi
                         ppA_roi = movieData.focusOnroi(Ga_TH_A);
                        %make pseudo color movie
                        movieData.makePseudoColorMovie(ppA_roi,filename(1:length(filename)-4))
                        disp(['Preprocessing done: ' filename]);
                    case 3
                        %prepare for seed-based correlation analysis
                        A_all = cat(3, A_all, Ga_TH_A);
                end             
                
                disp('')
            else
                disp('')
                disp('No filetered matrix detected!!!')
                disp('')
            end

    end
    
    if id == 3
        movieData.SeedBasedCorr(A_all)
    end
    disp(['Processing done at:' pwd]);
    
    clearvars;  
end
%%
function renewCC(ppA_roi,outputFolder,filename)
%renew connected components

        %Black-white thresholding of pre-processed A
        [BW_ppA,~] = Integration.bwThresholding_10prctPixels(ppA_roi); %ppA is short for pre-processed A
        clear Ga_TH_A;

        %Generate connected component
        [region,BW_ppA] = Integration.GenerateCC(ppA_roi,BW_ppA);
        checkname = ['CC_' filename(1:length(filename)-4) '_region.mat'];
        save(fullfile(outputFolder,checkname),'region');
        clear TH_A region

        %Save binary movie
        checkname = ['Binary_' filename(1:length(filename)-4) '.mat'];
        save(fullfile(outputFolder,checkname),'BW_ppA');
        clear BW_ppA
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