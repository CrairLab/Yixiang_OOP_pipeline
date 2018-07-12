function byPassPreProcessing(id)
%This function bypass pre-processing (assuming it has been done) and use
%filtered matrix to do further analysis
%07/11/18 can do renewCC, seed-based correlation, or generate pseudo color
%movie based on the filtered matrix.
%
%
%Inputs:
%id       determine which function to run
%(id == 1, renewCC; id == 2, pseudo-clolor movie; id == 3 seed-based corr)
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
                         ppA_roi = focusOnroi(Ga_TH_A);
                        %Renew connected components
                        renewCC(ppA_roi,outputFolder,filename)
                        disp(['Preprocessing done: ' filename]);
                    case 2
                        %Chop the matrix to contain only roi
                         ppA_roi = focusOnroi(Ga_TH_A);
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
        SeedBasedCorr(A_all)
    end
    disp(['Processing done at:' pwd]);
    
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
function A = focusOnroi(A)
%Chop the size to the roi part
    [dim1_lower,dim1_upper,dim2_lower,dim2_upper] = movieData.getROIBoundsFromImage(A(:,:,1)); 
    A = A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:); %ppA: pre-processed 

end

%%
function SeedBasedCorr(A,spacialFactor)
% Generate seed-based correlation maps based on seeds and filtered matrix
% Read in seeds(rois) from 'Seeds.zip'
%
% Inputs:
%   A                input matrix 
%   spacialFactor    factor previously used for downsampling
%
% Outpus:
%   seed-based correlation maps

    %Detect seeds
    temp = ROI('Seeds.zip');
    ROI_all = temp.ROIData;
    if isempty(ROI_all)
        disp('Seeds not provided! Cannot generate correlation maps!')
        return
    end

    %downSampleRatio = 1/spacialFactor
    if nargin == 2
        downSampleRatio = 1/spacialFactor;
    else 
        downSampleRatio = 0.5;
    end

        sz = size(A);
        imgall = reshape(A, sz(1)*sz(2), sz(3));

        disp('Note: ROI coordinates should be based on the downsampled movie size!')
        [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz);

        %Generate correlation map for each seed
        for r = 1:length(roi)     
            if size(roi{r}, 1) ~= sz(1)
                mask = imresize(roi{r}, downSampleRatio, 'bilinear');
            else
                mask = roi{r};
            end
            maskId = find(mask > 0);

            %Get seed time series 
            seedTrace(r, :) = squeeze(nanmean(imgall(maskId, :), 1));     

            %Generate correlation matrix
            corrM = corr(imgall', seedTrace(r, :)');

            corrMatrix(:, :, r) = reshape(corrM, sz(1), sz(2));

            %Plot correlation map
            h = figure; 
            cur_img = corrMatrix(:, :, r);
            imagesc(cur_img); colormap jet; colorbar; axis image
            caxis([-0.3, 1]); title(['roi:', num2str(r)]);
            hold on

            %Label seed position
            fill(roiPolygon{r}(:, 1), roiPolygon{r}(:, 2), 'y')

            %Save the plot           
            saveas(h, ['roi', num2str(r), '.png'])     
        end
    
end

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