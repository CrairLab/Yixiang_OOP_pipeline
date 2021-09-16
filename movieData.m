classdef movieData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movieData stores the input matrix 
% It has dozens of static and normal functions to process the matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All previous record is saved on Github
% Visit https://github.com/CrairLab/Yixiang_OOP_pipeline for more info
% Author: yixiang.wang@yale.edu
% Latest update:
% R36 08/27/20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    properties
        A;   %Input matrix        
    end
    
    methods
        
        function obj = movieData(filename,flag)
        %movieData object constructor
            if ~(flag == 2)
                obj.A = single(movieData.inputMovie(filename));
            end
            %obj.dfOverFA = movieData.fetchAverageImage(obj.A,filename);
  
        end
        
    end
    
    
    
    methods(Static)
        
         function A = fetchAverageImage(A,filename)
        
         %   This function computes and outputs the F/F0 - 1, where F0
         %   represents the matrix of averaged intensities over time at
         %   individual pixels
         %   
         %   Input:
         %       filename
         %   
         %   Output:
         %       Normalized movie (dF/F0, F0 is the mean intensity matrix
         %       averaged over time)
    
        
            %%%Make deltaF/F movie%%%
            sz = size(A); szZ=sz(3);
            
            %Transform 2D images to 1D arrays
            npix = prod(sz(1:2));           
            A = reshape(A, npix, szZ); 
            %avg at each pixel location in the image over time
            Amean = mean(A,2); 
            % F/F0 - 1 == ((F-F0)/F0);
            A = A ./ (Amean * ones(1,szZ)) - 1;               
            Amean = reshape(Amean,sz(1),sz(2));
            A = reshape(A, sz(1), sz(2), szZ);

            %%%Write average image%%%
            switch nargin
                case 2
                    fnm = [filename(1:length(filename)-4) '_AVG.tif'];
                    Amean=mat2gray(Amean);
                    [AmeanInd, ~] = gray2ind(Amean,65536);
                    imwrite(AmeanInd, fnm);
                    disp(filename)
            end

        end
        
        function A = inputMovie(filename)
        
        %    This function takes advantage of the openMovie function
        %    from the wholeBrainDX-master and uses it to input frames in
        %    a specific movie. It allows concatenation of several movies
        %    with a common prefix.
        %
        %    input: filename
        %
        %    output: 3D-matrix A containing data from a specific movie
        
            A = openMovie(filename);
            %tmp = dir([filename(1:length(filename)-4) '@00*.tif']);
            %Make a combined matrix for each recording
            %if ~isempty(tmp)
            %    for j = 1:numel(tmp)
            %        fn = tmp(j).name;
            %        B = openMovie(fn);
            %        A = cat(3, A, B);
            %        clear B
            %    end
            %end             
        end
             
        
        function timeProjectMaps(A,filename)     
        
        %    This function takes advantage of the timeColorMapProj functino
        %    from the wholeBrainDX-master and outputs the time-projection
        %    avi movies
        
            
            frStart=1;
            frEnd=size(A,3);
            [~, Iarr] = timeColorMapProj(A,frStart, frEnd, filename);
            %clear A
            %Iarr2montage(Iarr, frStart, frEnd, 10, filename);
            
            %Write avi movie
            Iarr2avi(Iarr, frStart, frEnd, filename);
        end
        
        
        function A = TopHatFiltering(A,hat)
        
        %    Doing top-hat filtering here, default hat value equals to 150
        %    
        %    Inputs:
        %        hat   parameter for top hat filtering  
        %    
        %    Outputs:
        %        A     filtered matrix
            
        
            if nargin == 1
                hat = 150;
            end
            
            sz = size(A);
            if length(sz) == 2
                sz(3) = 1;
            end
            
            A = reshape(A,sz(1)*sz(2),sz(3));
            se = strel('line', hat, 0);
            flpA = horzcat(fliplr(A(:,1:hat)),A(:,:));  
            % A_ln = zeros(sz(1)*sz(2), (sz(3))+hat);
            
            % Tophat filtering here, the filter is acting on the time axis
            % http://utam.gg.utah.edu/tomo03/03_mid/HTML/node120.html                 
            A_ln = imtophat(flpA, se);  %opening:https://en.wikipedia.org/wiki/Opening_(morphology)
          
            A = A_ln(:,(hat+1):end); 
            A=reshape(A,sz(1),sz(2),sz(3));    
        end
        
        
        function A = GauSmoo(A,sigma)
        
        %    Do Gaussian smoothing here
        %    
        %    Inputs:
        %        sigma
        %    
        %    Outputs:
        %        A     Smoothed matrix
        
            if nargin == 1
                sigma = 1;
            end
                                   
            sz = size(A);
            %A = zscore(reshape(A,sz(1)*sz(2)*sz(3),1)); %covert to standardized zscores so that mean=0, and sd = 1;
            A = reshape(A,sz);   
            Iarr = zeros(size(A));

            parfor fr = 1:sz(3) %option:for, do NOT use for it is 3x slower
               I = A(:,:,fr);
                   I2 = gaussSmooth(I,sigma,'same');
               Iarr(:,:,fr) = I2;
            end
            A = Iarr;
        end
        
        %{
        %% Substituted by ROIMask function in the ROI class 11/24/17 R3
            function [index, bothMasks3D] = FetchROI(sz)

                 ROI = ReadImageJROI('RoiSet.zip');       
                 bothMasks = poly2mask(ROI{1,1}.mnCoordinates(:,1), ROI{1,1}.mnCoordinates(:,2),sz(1),sz(2))+poly2mask(ROI{1,2}.mnCoordinates(:,1), ROI{1,2}.mnCoordinates(:,2),sz(1),sz(2));
                 index = find(bothMasks);
                 bothMasks3D = repmat(bothMasks,[1 1 sz(3)]);

            end
        %}
        
        
        function MeanFrame = MeanDuringSti(FreqVideo, exptparam)
        
        %    Calculate the mean of frames during stimulation.
        %   (during-stimulaton period)
        %    
        %    Inputs:
        %       FreqVideo     Averaged movie
        %       exptparam     experiment parameters from baphy files
        %
        %    Outputs:
        %       MeanFrame       mean intensities averaged over certain
        %       period
            
            ini_frame = round(exptparam.PreStimSilence) * 10;
            last_frame = round(exptparam.PreStimSilence + exptparam.Duration*2) * 10;
            MeanFrame = mean(FreqVideo(:,:,ini_frame+1:last_frame),3);
        end
        
        
        %Intra-trails dF/F
        function outputVideo = InsideVideodFoverF(curVideo, firstFrame)
        
        %    Calculate the df/f for each frame in a block video
        %    grouped by either frequencies or frequency/volume
        %    combination. Take the average of the pre-stimulation
        %    frames as background. Using this averaged matrix to
        %    calculate df/f for each frame.
        %
        %    Inputs:
        %        curVideo     currently working video
        %        firstFrame   first frame when stimulation starts
        %    
        %    Outputs:
        %        outputVideo     df/f video
        
        if nargin<2
            firstFrame = 10;
        end
        
            meanVideo = mean(curVideo(:,:,1:firstFrame),3);
            sz= size(curVideo);
            meanVideo = repmat(meanVideo,[1,1,sz(3)]);
            outputVideo = curVideo./meanVideo - 1;
         
            outputVideo(isnan(outputVideo)) = 0;
            %outputVideo = mat2gray(outputVideo);
        end
        
        
         function FreqColorMap(A,filename)
         
         %    This function output the time-color projeciting maps based on frames
         %    sorted by different frequencies. Say 5 different frequencies
         %    were presented to the animal, Frames would be sorted and
         %    stored in _FramesByFreq.mat. This function read in this .mat
         %    data and use these indices to extract frames with the same
         %    frequencies. These frames will be averaged across the whole
         %    movie and made into 40-frame para-stimulation video clips.
         %    
         %    Inputs:
         %       A          input movie
         %       filename   corresponding filename
         %    
         %    Outputs:
         %       Averaged df over f images over frames 20-29
         %       Averaged 40-frame para-stimulation video clips
            
            LoadName = [filename(1:length(filename)-4) '_FramesByFreq.mat'];
            load(LoadName);
            FreqNumber = size(FramesByFreq,2);
            sz = size(A);
            temp_length = [];

            % The intrinsic pseudorandom property of Baphy could render slightly 
            % different trail numbers for differnt frequencies.
            % This step get the max trail number across differen frequencies for
            % further constrcutin of a unified matrix
            for i = 1:FreqNumber
                temp_length = [temp_length size(FramesByFreq{2,i},1)];
            end
            MaxFrameLength = ceil(max(temp_length)/40)*40;

            for i = 1:FreqNumber
                temp = FramesByFreq(2,i);
                temp_Idx = temp{1};
                temp_Idx = [temp_Idx;zeros(MaxFrameLength - length(temp_Idx),1)];
                temp_Idx = reshape(temp_Idx,40,[]);
                FreqPic = zeros(sz(1),sz(2));
                FreqVideo = zeros(sz(1),sz(2),40);
                
                for j = 1:40
                    temp_Idx_j = temp_Idx(j,:);
                    temp_Idx_j(temp_Idx_j ==  0) = [];
                    NonZero = temp_Idx_j;
                    FreqVideo(:,:,j) = mean(A(:,:,NonZero),3);       
                end
                
                cur_dFoverF = movieData.InsideVideodFoverF(FreqVideo); 
                FreqVideo = cur_dFoverF;
                FreqPic = movieData.MeanDuringSti(FreqVideo);
                pictype = mat2gray(squeeze(FreqPic));
                [pictypeInd, ~] = gray2ind(pictype,65536);
                tmp_filename = strcat([filename(1:length(filename)-4)],['_',num2str(FramesByFreq{1,i}),'Hz_AVG-dFoF.png']);
                imwrite(pictypeInd, tmp_filename);
                tmp_filename = strcat([filename(1:length(filename)-4)],['_',num2str(FramesByFreq{1,i}),'Hz_dFoF.avi']);
                [~, Iarr] = timeColorMapProj(squeeze(FreqVideo),1,40,tmp_filename);            
                Iarr2avi(Iarr,1,40,tmp_filename);
            end

        end
        
        
        
        function FreqVoluColorMap(A,filename,nmov,exptparam)
            
        %     This function output the time-color projeciting maps based on frames
        %     sorted by different frequencies and volume combinations. 
        %     Say 5 different frequencies and 5 different sound levels
        %     were presented to the animal, Frames would be sorted and
        %     stored in _FramesByFreqVolu.mat. This function read in this .mat
        %     data and use these indices to extract frames with the same
        %     frequency and volume combination. These frames will be averaged across the whole
        %     movie and made into 40-frame para-stimulation video clips.
        %     
        %     Inputs:
        %        A          input movie
        %        filename   corresponding filename
        %     
        %     Outputs:
        %        Averaged df over f images over frames during stimulation
        %        Averaged para-stimulation video clips

        if isempty(exptparam)
            exptparam.PreStimSilence = 1;
            exptparam.PostStimSilence = 3.5;
            exptparam.Duration = 0.5;
            exptparam.BlockDura = 60;
        end
                 
            ini_frame = round(exptparam.PreStimSilence) * 10;
            last_frame = round(exptparam.PreStimSilence + exptparam.Duration*2) * 10 - 1;
            BlockDura = exptparam.BlockDura;
       
            LoadName = [filename(1:length(filename)-4) '_FramesByFreqVolu.mat'];
            load(LoadName);
            FreqNumber = size(FramesByFreqVolu,2);
            VoluNumber = size(FramesByFreqVolu,1);
            sz = size(A);
            temp_length = [];

            % The intrinsic pseudorandom property of Baphy could render slightly 
            % different trail numbers for differnt frequencies.
            % This step get the max trail number across differen frequencies for
            % further constrcutin of a unified matrix
            for i = 1:FreqNumber
                for j = 1:VoluNumber            
                    temp_length = [temp_length size(FramesByFreqVolu{j,i},1)];
                end
            end
            MaxFrameLength = ceil(max(temp_length)/BlockDura)*BlockDura;
            
            %FreqVoluPic = zeros(sz(1),sz(2));
            FreqVoluVideo = zeros(sz(1),sz(2));
            
            for i = 1:FreqNumber
                for j = 1:VoluNumber          
                    temp_Idx = FramesByFreqVolu{j,i};
                    ThisSignature = temp_Idx(1:2,1);
                    temp_Idx = temp_Idx(3:size(temp_Idx,1),1);
                    temp_Idx = [temp_Idx;zeros(MaxFrameLength - length(temp_Idx),1)];
                    temp_Idx = reshape(temp_Idx,BlockDura,[]);

                    for k = 1:BlockDura
                        temp_Idx_k = temp_Idx(k,:);
                        temp_Idx_k(temp_Idx_k ==  0) = [];
                        NonZero = temp_Idx_k;
                        FreqVoluVideo(:,:,k) = mean(A(:,:,NonZero),3);
                    end
                    
                    cur_dFoverF = movieData.InsideVideodFoverF(FreqVoluVideo, ini_frame);
                    FreqVoluVideo = cur_dFoverF;
                    
                    %FreqVoluPic(:,:,j,i) = movieData.MeanDuringSti(cur_dFoverF);
                    %CurrentMatrix = FreqVoluVideo(:,:,:,j,i);
                    movieData.AverageAcrossMovies(ThisSignature,FreqVoluVideo,nmov);
                    %{
                    pictype = mat2gray(squeeze(FreqVoluPic(:,:,j,i)));
                    [pictypeInd, cmap] = gray2ind(pictype,65536);
                    tmp_filename = [filename(1:length(filename)-4) '_' num2str(FramesByFreqVolu{j,i}(1,1)) 'Hz_' num2str(FramesByFreqVolu{j,i}(2,1)) 'dB_AVG-dFoF.png'];
                    imwrite(pictypeInd, tmp_filename);

                    tmp_filename = [filename(1:length(filename)-4) '_' num2str(FramesByFreqVolu{j,i}(1,1)) 'Hz_' num2str(FramesByFreqVolu{j,i}(2,1)) 'dB_AVG-dFoF.avi'];
                    [maxProj, Iarr] = timeColorMapProj(squeeze(FreqVoluVideo(:,:,:,j,i)),1,40,tmp_filename);            
                    Iarr2avi(Iarr,1,40,tmp_filename);
                    %}


                end
            end
        end


        function AverageAcrossMovies(ThisSignature,CurrentMatrix,nmov)
        
        %    Averaged across different movies. Grouped by specific frequencies and
        %    volumes
        %
        %    Inputs:
        %        ThisSignature      A matrix that stores all possible combination of
        %                           frequences and voluems for a specific movie
        %        CurrentMatrix      Currently working matrix (movie)
        %        nmov               Amount of movies 
        %
        %    Outputs
        %        AveragedMatrix.mat      The ongoing add-up averaged matrix.

       
            current_name = ['AveragedMatrix_' num2str(ThisSignature(1,1)) 'Hz_' num2str(ThisSignature(2,1)) 'dB.mat'];
            
            if nargin == 2
                filelist = readtext('files.txt',' ');
                nmov = size(filelist,1);
            end
            
            %Reasign 0 values to nan
            CurrentMatrix(CurrentMatrix == 0) = nan;
            
            %If the .mat file existed, load and update
            if exist(current_name,'file')
                load(current_name);
                AveragedMatrix = AveragedMatrix + CurrentMatrix./nmov;
                save(current_name,'AveragedMatrix');
            else
                AveragedMatrix = CurrentMatrix./nmov;
                save(current_name,'AveragedMatrix');
            end

        end
        
        
        function AveragedMapsAcrossMovies(exptparam)
            
            %    Calling this function will read in all
            %    AveragedMatrix_freq_volu.mat data (should be updated with
            %    all movies in advance, see the AverageAcrossMovies
            %    function) and output maps and videos averaged across movies
            %    It should be called at the end of the main pipeline. 
            %    Input: exptparam from baphy file
            
            BlockDura = exptparam.BlockDura; %Duraion of a single block
            
            temp_info = dir;
            filelist = {};
            n = 0;
            for i = 1:size(temp_info,1)  
                if regexp(temp_info(i,1).name,'AveragedMatrix_')
                    n = n+1;
                    filelist{n,1} = temp_info(i,1).name;
                end
            end

            for i = 1:n
                load(filelist{i,1});
                filename = filelist{i,1};
                AveragedMatrix(isnan(AveragedMatrix)) = 0;
                AveragedMap = movieData.MeanDuringSti(AveragedMatrix, exptparam);
                pictype = mat2gray(squeeze(AveragedMap));
                [pictypeInd, ~] = gray2ind(pictype,65536);
                tmp_filename = [filename(1:length(filename)-4) '.png'];
                imwrite(pictypeInd, tmp_filename);
                tmp_filename = [filename(1:length(filename)-4) '.avi'];
                %VideoDFoverF = movieData.InsideVideodFoverF(AveragedMatrix); %modified in R3
                [~, Iarr] = timeColorMapProj(squeeze(AveragedMatrix),1,BlockDura,tmp_filename);
                Iarr2avi(Iarr,1,BlockDura,tmp_filename);

            end

        end
        
                
        %%Downsample
        function outputA = downSampleMovie(A,s,t)
        
        %    Dowansample the input matrix A
        %    Temporal downsampling: choose frames spanning by a specified factor
        %    spatial downsampling: averaged pixels in a s-by-s window
        %
        %    Inputs:
        %        A   Input movies(3D matrix)
        %        s   spatial downsampling factor 
        %        t   temporal downsampling factor
        %
        %    Outputs:
        %        outputA    downsampled matrix

            switch nargin
                case 1
                    disp(['Spatial Factor: 2' ])
                    outputA = movieData.spatialDown(A,2);                    
                case 2
                    disp('Only do spatial downsampling here!')
                    disp(['Spatial Factor:' num2str(s)])
                    outputA = movieData.spatialDown(A,s);
                case 3
                    disp(['Spatial Factor:' num2str(s)])
                    disp(['Temporal Factor:' num2str(t)])
                    outputA = movieData.temporalDown(A,t);
                    outputA = movieData.spatialDown(outputA,s);
                    
            end
        end


        function downA = spatialDown(A,s)
        
        %    Downsample the input matrix A by factor s. The idea is assign a new pixel
        %    by averaging pixels within a s-by-s whindow.
        %
        %    Inputs: 
        %        A   input matrix
        %        s   spatial downsampling factor
        %
        %    Output:
        %        downA    spatially downsampled matrix
            
            downA = imresize(A, 1/s, 'bilinear');
        %    sz = size(A);
        %    newsz1 = mod(-sz(1),s)+sz(1);
        %    newsz2 = mod(-sz(2),s)+sz(2);
        %    IdxEnd1 = ceil(sz(1)/s)*s;
        %    IdxEnd2 = ceil(sz(2)/s)*s;
            
            %downA = nan(newsz1/s,newsz2/s,sz(3));
        %    downA = [];
            
        %    if length(sz)<3
        %        sz(3) = 1;
        %    end

        %    parfor k = 1:sz(3)  
        %        A_newk_avg = nan(newsz1/s,newsz2/s,s^2);
        %        A_k = A(:,:,k);
        %        A_newk = nan(newsz1,newsz2);
        %        A_newk(1:sz(1),1:sz(2)) = A_k;
        %        count = 0;
        %        for i = 1:s
        %            Idx1 = [1+i-1:s:IdxEnd1];
        %            for j = 1:s
        %                count = count + 1;
        %                Idx2 = [1+j-1:s:IdxEnd2];
        %                A_newk_avg(:,:,count) = A_newk(Idx1,Idx2);
        %            end
        %        end
        %        downA_k = nanmean(A_newk_avg,3);
        %        downA(:,:,k) = downA_k;
        %    end
        end


        function A_new = temporalDown(A,t)
        
        %    This function downsamples the original movie simply by selecting a
        %    subset of the input A (save every t-th frame)
        %
        %    Inputs:
        %        A   input matrix
        %        t   temporal downsampling factor
        %
        %    Outputs:
        %       A_new   new subset of the original movie
        
            sz = size(A);
            A_ln = reshape(A,sz(1)*sz(2),sz(3));
            A_new = [];
            downIdx = downsample((1:sz(3)),t);
            A_ln = A_ln(:,downIdx);
            A_new = reshape(A_ln,sz(1),sz(2),[]);
            disp(['Temporally downsampling by the factor of ' num2str(t)]);

        end
        
        
      
        function ImageThresholding(thresh)
        
        %    Thresholding all averaged images started with 'AveragedMatrix_'
        %    zscoring the matrix than using default thresholding factor
        %    thrsh = 1
        %    
        %    Inputs:
        %        Read in images with prefix 'AveragedMatrix'
        %        thresh     thresholding factor
        %    
        %    Ouputs:
        %        Thresholded bw images
        
        
        if nargin == 0
            thresh = 1;
        end

        temp_info = dir;
            filelist = {};
            n = 0;
            for i = 1:size(temp_info,1)  
                if ~isempty(findstr(temp_info(i,1).name,'AveragedMatrix_'))...
                        && isempty(findstr(temp_info(i,1).name,'threshfactor')) %#ok<FSTR>
                    n = n+1;
                    filelist{n,1} = temp_info(i,1).name;
                end
            end

            for i = 1:n
                filename = filelist{i,1};
                A = imread(filename);
                sz = size(A);
                pixN = sz(1) * sz(2);
                doubleA = reshape(double(A),pixN,1);
                bwA = movieData.simpleBwThresholding(doubleA,thresh);
                bwA = reshape(bwA,sz(1),sz(2));
                pictype = mat2gray(bwA);
                savename = [filename '_threshfactor' num2str(thresh) '.png'];
                imwrite(pictype,savename,'png');

            end
        end
        
        
        
        function A = simpleBwThresholding(A,thresh)
            
        %    Do simple black-white thresholding for 3-D matrix or 2-D matrix
        %    Note that the mean values are calculated from each frame
        %    
        %    Inputs:
        %        A    input matrix, either 3D or 2D
        %        thresh    thresholding factor
        % 
        %    Outputs:
        %        A    output thresholded matrix
           
            if nargin == 1
                thresh = 1;
            end
            sz = size(A);            
            if length(sz) == 2
                sz(3) = 1;
            end

            %singleThresholding = @(A, thresh) reshape(zscore(A(:)),size(A))>thresh;
            parfor k = 1:sz(3)
                curImg = A(:,:,k);
                curImg(~isnan(curImg)) = zscore(curImg(~isnan(curImg)));
                tmp = curImg > thresh;
                tmp(tmp == 0) = nan;
                A(:,:,k) = tmp;
                
                %A(:,:,k) = singleThresholding(A(:,:,k),thresh);                            
            end

        end
        
        
        
        function [bwA,postiveRatio] = bwThresholding(A,factor)
        
        %    Filter out pixels <= mean + factor * std. Note that the mean
        %    and std values here are calculated from the whole movie.
        %
        %    Inputs:
        %        A     Input movie (3D matrix). Should be a matrix after 
        %              bleaching correction, top-hat filtering etc
        %        factor    Threshold = factor*std above mean
        %
        %    Outputs:
        %        bwA             Thresholded movie
        %        postiveRatio    Fraction of pixels that are above threshold


            sz = size(A);
            bwA = zeros(sz);
            postiveRatio = zeros(1,sz(3));

            %Acuqire current matrix
            bgMatrix = A;
            bgMatrix(bgMatrix == 0) = nan;
            bgMean = nanmean(bgMatrix(:));
            bgStd = nanstd(bgMatrix(:));
            if nargin == 1
                factor = 1;
            end
            
            disp(['Thresholding factor = ' num2str(factor)]);
            
            %filter out pixels <= mean + factor * std
            parfor i = 1:sz(3)
                currMatrix = A(:,:,i);
                pixelsN = sz(1) * sz(2) - sum(isnan(currMatrix(:)));
                pixelsOverT = sum(currMatrix(:) > (bgMean + factor * bgStd));
                overThreRatio = pixelsOverT / pixelsN;
                postiveRatio(1,i) = overThreRatio;

                bwCurrMatrix = currMatrix > (bgMean + factor * bgStd);
                bwA(:,:,i) = bwCurrMatrix;              
            end

        end
        
        
        function [bwA,factors_all] = bwThresholding_10prctPixels(A)
        
        %    Implementing binary thresholding here for bwconncomp by filtering 
        %    out pixels with intensity < mean + factor * std. Mean and std 
        %    are calculated based on each single frame. factors are incrementally
        %    decreased to preserve at least 10% of overall pixels (within in ROIs)
        %    
        %    Inputs:
        %        A     Input movie (3D matrix). Should be a matrix processed
        %              after applying ROIs, Gaussian smoothing and top-hat filtering
        %    Outputs:
        %       bwA          Thresholded movie
        %       factors_all  Factors applied for each frame               
        
            sz = size(A);
            bwA = zeros(sz);
            factors_all = zeros(sz(3),1);
            parfor i = 1:sz(3)
                %Acuqire current matrix
                currMatrix = A(:,:,i);
                currMatrix(currMatrix == 0) = nan;
                currMean = nanmean(currMatrix(:));
                currStd = nanstd(currMatrix(:));
                factor = 1.6;
                
                pixelsN = sz(1) * sz(2) - sum(isnan(currMatrix(:)));
                pixelsOverT = 0;
                overThreRatio = pixelsOverT / pixelsN;
                
                %while proportion of over-thresholding ratio pixels is
                %smaller than 10%. Reduce thresholding factor.
                while (overThreRatio < 0.1) && (factor >= 0)
                    pixelsOverT = sum(currMatrix(:) > (currMean + factor * currStd));
                    overThreRatio = pixelsOverT / pixelsN;
                    factor = factor - 0.2;                    
                end
                factors_all(i) = factor;
                
                bwCurrMatrix = currMatrix > (currMean + factor * currStd);
                bwA(:,:,i) = bwCurrMatrix;              
               
            end
            
        end
        
        function A = grossDFoverF(A, flag)
        %    Doing gross dFoverF calculation.
        if nargin == 1
            flag = 1;
        end
        
            sz = size(A);
            A_re = reshape(A,[sz(1)*sz(2),sz(3)]);
                
            if flag == 1
            %Use bottom 5th percentile as F0
                A_F0 = prctile(A_re,5,2);
                A_F0 = repmat(A_F0,[1,sz(3)]);
            else
            % Use mean values as F0
                A_F0 = repmat(mean(A_re,2),[1,sz(3)]);
            end
            A = reshape(A_re./A_F0 - 1,sz);
            
        end
        
        
        function A = movieRigReg(A, A_fixed, movIdx)
        
        %    Do rigid registration to input matrix A, actually only do
        %    translation to save time
        %    
        %    Inputs:
        %    A        Original movie without rigid registration (3D)
        %    A_fixed  fixed frame as reference
        %    movIdx   Indices indicating frames to be registered
        %    
        %    Outputs:
        %    A        Movie after rigid registration
            A_fixed(isnan(A_fixed)) = 0;
            A(isnan(A)) = 0;
            
            [optimizer, metric] = imregconfig('monomodal');
            if nargin == 1
                A_fixed = nanmean(A,3);
            end
            if nargin == 2
                tic;
                parfor i = 1:size(A,3)
                    if ~mod(i,400)
                        disp(['Finish rigid registration at frame #' num2str(i)]);
                    end
                    A(:,:,i) = imregister(A(:,:,i), A_fixed, 'translation', optimizer, metric);
                end
                toc;
            elseif nargin == 3
                Idx = find(movIdx);
                A_toRegister = A(:,:,Idx);
                tic;
                A_Registered = [];
                parfor i = 1:length(Idx)
                    if ~mod(i,200)
                        disp(['Finish rigid registration at frame #' num2str(i)]);
                    end
                    A_Registered(:,:,i) = imregister(A_toRegister(:,:,i),...
                        A_fixed, 'translation', optimizer, metric);
                end
                A(:,:,Idx) = A_Registered;
                toc;
            end
            
            A(A == 0) = nan;
            
        end
        
        
        function A = bleachCorrection(A)
        
        %    Correcting fluorescence bleaching, assuming exponential decrease
        %   
        %   Inputs:
        %        A          input matrix
        %   
        %    Outputs:
        %        A          bleaching corrected matrix
            
            sz= size(A);
            A_re = reshape(A,[sz(1)*sz(2),sz(3)]);
            x = 1:sz(3);
            
            mean_series = mean(A_re,1);
            
            f = fit(x',mean_series','exp1');
            trend = f.a.*exp(f.b.*x);
            trend = trend./min(trend);
            
            
            trend = repmat(trend,[sz(1)*sz(2),1]);
            A_corrct = A_re./trend;
            A = reshape(A_corrct,sz);
                
        end
        
        
        
        function [A,U,S,V,iniDim,PC_exp] = roiSVD(A, iniDim)

        %    This function automatically identify the minimum rectangle containing roi
        %    and then do SVD denosing to only the roi part of the original matrix. It 
        %    preseves the 4th to the psv_dim-th principle components and recover the
        %    reconstructed roi part to original matrix size
        %
        %    Inputs:
        %    A        3D matrix
        %    iniDim   First preserved component
        %   
        %    Outputs:
        %    A        processed matrix
            
            %iniDim defines from which dimension shall the function starts
            %to preserve
            if nargin<2
                iniDim = 1;
            end
            disp(['Initial Dimention = ' num2str(iniDim)]);
            
            %Identify the vertex of the minimum rectangle containing roi
            cur_img = A(:,:,1);
            [dim1_lower,dim1_upper,dim2_lower,dim2_upper] = movieData.getROIBoundsFromImage(cur_img);

            %Do SVD to the roi only
            A_roi = A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:);
            sz = size(A_roi);
            A_roi = reshape(A_roi,[sz(1)*sz(2),sz(3)]);
            A_roi(isnan(A_roi)) = 0;
            psv_dim = round(size(cur_img,1)*size(cur_img,2)/800); %empirically suffcient to have a good-quality reconstruction
            psv_dim = max(psv_dim,120); %At least preserve first 120 dimensions
            disp(['Preserved dimensions = ' num2str(psv_dim)]);
            %tic; [U,S,V] = svds(A_roi,psv_dim); toc;
            tic; [U,S,V] = svds(A_roi',psv_dim); toc;
            
            %Reconstruct the roi pairt using the the 4th to the last component
            A_roi_rcs = V(:,iniDim:end)*S(iniDim:end,iniDim:end)*U(:,iniDim:end)';
            %A_roi_rcs = reshape(A_roi_rcs,[sz(1),sz(2),sz(3)]);
            A_roi_rcs = reshape(A_roi_rcs,[sz(1),sz(2),sz(3)]);

            %Recover the roi to the original image size
            A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:) = A_roi_rcs;
            %A(A == 0) = nan;
            
            %Save U,S,V
            %U = reshape(U,[sz(1), sz(2), size(U,2)]);
            V = reshape(V,[sz(1), sz(2), size(V,2)]);
            %save('SVD_decomp.mat','U','S','V');
            
            % Calculate the propotion of variance explained by each PCs
            eigs = diag(S).^2;
            Sum_var = sum(diag(A_roi'*A_roi));
            PC_exp = eigs./Sum_var;
            
            %save('Var_explained_by_PCs.mat', 'PC_exp');
            
            disp('');

        end
        
        
        
        function [nZ_1_lower,nZ_1_upper,nZ_2_lower,nZ_2_upper] = getROIBoundsFromImage(cur_img)
        
        %    This function identify the coordinates of vertex of the minimum
        %    rectangle containing the roi
        %
        %    Inputs:
        %        cur_img          A 2D image containg roi
        %
        %    Outputs:
        %        nZ_1_lower      lower left vertex of the minimum rectangle
        %        nZ_1_upper      upper left vertex of the minimum rectangle
        %        nZ_2_lower      lower right vertex of the minimum rectangle
        %        nZ_2_upper      upper right vertex of the minimum rectangle
        
            
            
            if (~any(isnan(cur_img(:)))) && (cur_img(1) == 0) && (cur_img(end) == 0)
                cur_img = ~(cur_img == 0);
            else
                cur_img = ~isnan(cur_img);
            end

            nZ_2 = find(mean(cur_img,1));
            nZ_2_upper = max(nZ_2);
            nZ_2_lower = min(nZ_2);
            nZ_1 = find(mean(cur_img,2));
            nZ_1_upper = max(nZ_1);
            nZ_1_lower = min(nZ_1);

        end
        
        function makePseudoColorMovie(A,movieName,FrameRate)
            
        % Make pseudo-color (jet 256) movie from zscored matrix A
        % Inputs:
        %   A    3D matrx (grayscale)
        %   movieName   Name of the movie
        %
            if ~exist('movieName','var')
                movieName = 'outputPseudocolorMovie';
            end
            
            if ~exist('FrameRate','var')
                FrameRate = 50;
            end
            
            %Doing zscoring here
            [X,~] = gray2ind(zscore(A,1,3));
            sz = size(X);
            new_X = zeros(sz(1),sz(2),3,sz(3));
            for i = 1:sz(3)
                %disp(num2str(i))
                new_X(:,:,:,i) = ind2rgb(X(:,:,i), jet(256));
            end
            mov = immovie(new_X);
            v = VideoWriter(movieName,'Motion JPEG AVI');
            v.FrameRate = FrameRate;
            open(v)
            writeVideo(v,mov);
            close(v)
        end

        
        function corrMatrix = SeedBasedCorr_GPU(A,spatialFactor,...
                total_seeds,GPU_flag,plot_flag,mean_flag, Avg_out_all, timelag)
        % Generate seed-based correlation maps based on seeds and filtered matrix
        % Read in seeds(rois) from 'Seeds.zip'. If manually defined seeds
        % are not available, try to automatically generate seeds that evenly
        % cover given rois or the whole image. Compatiable with GPU
        % computing if the GPU_flg == 1.
        %
        % Inputs:
        %   A                input matrix 
        %   spatialFactor    factor previously used for downsampling
        %   total_seeds      number of seeds to be generated
        %   GPU_flag         whether run on GPU
        %   plot_flag        whether plot the correlation maps or not
        %   mean_flag        whether regress out background correlation (1 yes)
        %   Avg_out_all      background trace
        %   timelag          number of frames for time-lag correlation.
        %                    Can be positive or negative.
        %   
        % Outpus:
        %   seed-based correlation maps


            %downSampleRatio = 1/spatialFactor
            if ~exist('spatialFactor','var')
                spatialFactor = 2;
            end

            if ~exist('total_seeds','var')
                total_seeds = 400;
            end

            if ~exist('GPU_flag', 'var')
                GPU_flag = 0;
            end
            
            if ~exist('plot_flag','var')
                plot_flag = 0;
            end
            
            if ~exist('mean_flag','var')
                mean_flag = 0;
            end
            
            if ~exist('Avg_out_all','var')
                Avg_out_all = [];
            end
            
            if ~exist('timelag','var')
                timelag = 0;
            end

            disp(['Note: Defalut spatial factor == ' num2str(spatialFactor)]);

            sz = size(A);
            imgall = reshape(A, sz(1)*sz(2), sz(3));

            %Detect seeds
            temp = ROI('Seeds.zip');
            ROI_all = temp.ROIData;
            if isempty(ROI_all)
                disp('Seeds not provided (no manually defined seeds!)')
                disp('Generating seeds evenly covering current rois...')
                roi = ROI.genSeedingROIs(total_seeds,sz(1:2));
                sflag = 1; % for automatically generated seeds
            else

                %Generate cell array of rois
                try
                    [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz);

                    %When ROIs are generated based on original movie size, it
                    %might not cause error. These lines try to detect whether
                    %roi size agrees with sz. If not, it indicates that roi
                    %need to be generated with original movie size.
                    flag = 0;
                    for r = 1:length(roi)     
                        if size(roi{r}, 1) ~= sz(1)
                            flag = 1;
                            break
                        end
                    end
                    if flag == 1
                        [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz.*spatialFactor);
                        disp('Program thinks that rois are generated based on ORIGINAL movie size')
                    else
                        disp('Program thinks that rois are generated based on DOWNSAMPLED movie size')
                    end

                catch
                    [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz.*spatialFactor);
                    disp('Program thinks that rois are generated based on ORIGINAL movie size')
                end
                disp(['Recommending generate ROIs based on downsampled movie/frames!'])
                disp(['Warning: If ROIs are generated based on original sized movies, correlation maps could be wrong!'])
                sflag = 2; %for manually defined seeds
            end

            %Generate correlation map for each seed
            disp(['Detected ' num2str(length(roi)) ' seeds...'])

            for r = 1:length(roi)

                %Report progress
                if mod(r,100) == 0
                    disp(['Generating seeds#' num2str(r)])
                end

                %sflag == 1 automatically generated seeds; sflag == 2
                %manually generated seeds
                if sflag == 2
                    if size(roi{r}, 1) ~= sz(1)
                        mask = imresize(roi{r}, 1/spatialFactor, 'bilinear');
                    else
                        mask = roi{r};
                    end
                    maskId = find(mask > 0);
                elseif sflag == 1
                    %All automatically generated rois are based on
                    %downsampled image size
                    maskId = roi(r,:);
                end
                %Get seed time series 
                seedTrace(r, :) = squeeze(nanmean(imgall(maskId, :), 1));     
            end

            %Generate correlation matrix
            [corrMatrix, roi] = movieData.generateCorrMatrix(sz, roi, seedTrace,...
                imgall, GPU_flag, mean_flag, Avg_out_all, timelag);
            
            if plot_flag
                if sflag == 1
                    movieData.plotCorrM(roi, corrMatrix)
                elseif sflag == 2
                    movieData.plotCorrM(roi, corrMatrix, 'roipolygon',...
                        roiPolygon, 'downsample', spatialFactor)
                end
            end
            
        end
        
        function plotCorrM(roi, corrMatrix, varargin)       
        % Plot correlation maps based on input seeds & correlation matrix
        % If you have manually defined seeds in .zip file, please do follow:
        % 1. Get ROI_all from your .zip roi file using ROI class
        % 2. [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz);
        % 3. plotCorrM(roi, corrMatrix, 'roipolygon'', roiPolygon value,
        % 'downsample', downsample ratio)
        %
        % Inputs:
        %       roi           seeds indices 
        %       corrMatrix    correlation matrix
        %       spatialFactor    spatialFactor of the original matrix
        %       roipolygon    polygon region of the input seeds
        %
        % Outputs:
        %      correlation maps based on input information
           
            disp('If you have manually defined seeds, please do follow:')
            disp('1. Get ROI_all from your .zip roi file using ROI class')
            disp('2. [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz);')
            disp('plotCorrM(roi, corrMatrix, ''roipolygon'', roiPolygon value, ''spatialFactor'', spatial factor)')
            
            sflag = 1;
            if nargin < 2
               error ('You must supply at least rois and correlation matrix');
            end
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            elseif nargin > 3
                sflag = 2;
                for i=1:2:(length(varargin)-1)
                    if ~ischar (varargin{i})
                     error (['Unknown type of optional parameter name (parameter' ...
                     ' names must be strings).']);
                    end
                    switch lower (varargin{i})
                        case 'roipolygon'
                            roiPolygon = varargin{i+1};
                        case 'downsample'
                            spatialFactor = varargin{i+1};
                    end
                end
            end
                     
                          
           %Plot correlation map
           disp(['Actual number of correlation maps = ' num2str(length(roi))])
           
           parfor r = 1:length(roi)
                h = figure; 
                set(gcf,'Visible', 'off');
                cur_img = corrMatrix(:, :, r);
                imagesc(cur_img); colormap jet; colorbar; axis image
                caxis([-0.2, 1]); title(['roi:', num2str(r)]);
                hold on                   

                %Label seed position
                if sflag == 2
                    if size(roi{r}, 1) ~= sz(1)
                        %fill(roiPolygon{r}(:, 1)./spatialFactor, roiPolygon{r}(:, 2)./spatialFactor, 'y')
                    else
                        %fill(roiPolygon{r}(:, 1), roiPolygon{r}(:, 2), 'y')
                    end
                elseif sflag == 1
                    %create a polygon region to label the location of
                    %current roi
                    maskId = roi(r,:);
                    x1 = ceil(maskId(1)/size(cur_img,1));
                    y1 = maskId(1) - floor(maskId(1)/size(cur_img,1))*size(cur_img,1);
                    try
                        fill([x1-2,x1-2,x1+2,x1+2],[y1-2,y1+2,y1+2,y1-2], 'y')
                    catch
                        disp('Index exceed matrix limit')
                    end
                end

                %Save the plot
                num_str = num2str(10000 + r);
                num_str(1) = '0';
                saveas(h, ['roi', num_str, '.png'])     
            end 
            
        end
      
        
        
        function corrMatrix = batchSeedsGPU(sz,seedTrace,imgall)
        %Break a big matrix into smaller batch matrix due to limited GPU
        %memory limitation
        %
        %Inputs:
        %   sz             size of the original matrix
        %   seedsTrace     traces of seeds
        %   imgall         all images combined
        %
        %Outputs:
        %   corrMatrix     correlation matrix
        %
            reset(gpuDevice);
            pixelNum = sz(1)*sz(2);
            %Determine batch size based on the length of the time series
            batchSize = floor(10000./sz(3)*6000); %10000 is empirical, 6000 is from typical 10min movie
            batchIdx = [0:batchSize:pixelNum];
            seedTrace_GPU = gpuArray(seedTrace);
            corrM = [];
            for i = 1:length(batchIdx)
                if i < length(batchIdx)
                    cur_batch = imgall(batchIdx(i)+1:batchIdx(i+1),:);
                else
                    cur_batch = imgall(batchIdx(i)+1:end,:);
                end

                cur_batch_GPU = gpuArray(cur_batch);
                cur_corrM = corr(cur_batch_GPU', seedTrace_GPU');
                corrM = [corrM; gather(cur_corrM)];
            end
            sz_3 = size(corrM,2);
            corrMatrix = reshape(corrM,sz(1), sz(2),sz_3);
            reset(gpuDevice);
        end

        
        function A = focusOnroi(A)
        
        %Chop the matirx to the roi part
        %Inputs/Outputs:
        %A     3D Matrix
        
            [dim1_lower,dim1_upper,dim2_lower,dim2_upper] = movieData.getROIBoundsFromImage(A(:,:,1)); 
            A = A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:); %ppA: pre-processed 
        end
        
        function A_pooled = maxPooling2D(A,r)

        %   This function ouput the maximum-pooled matrix. The raidus of sliding
        %   window for maximum pooling is defined by r (or in default == 2).
        %
        %   Input:
        %       A   Input matrix (2D)
        %       r   radius of the sliding window
        %    
        %   Output:
        %       A_pooled = maximum pooled matrix
        %
            if nargin == 1
                r = 2;
            end
            sz = size(A);
            %Create a padding matrix
            width_new = sz(1) + 2*r;
            length_new = sz(2) + 2*r;
            sz_new = [width_new, length_new];
            A_pd = NaN(sz_new);
            A_pd(1+r:width_new-r,1+r:length_new-r) = A;

            A_pooled = zeros(sz);

            %Slide the window through every each pixels of the original matrix
            for i = 1+r:width_new-r
                for j = 1+r:length_new-r
                    %Calculate the maximum value in the current window
                    curWindow = A_pd(i-r:i+r,j-r:j+r);
                    max_curWindow = max(curWindow(:));

                    %Asign this value to the corresponding pixel in the pooled matrix
                    A_pooled(i-r,j-r) = max_curWindow;
                end
            end
        end
        
        
        
        function [icasig, M, W, corr_map] = getICA(A)
        %   Do fastICA and generate correlation maps with computed independent
        %   components
        %
        %   Inputs:
        %       A          3D matrix
        %
        %   Outputs:
        %       Save correlation maps with each independent components
        %       Save the computed ICA results
        
            %further downsampled the movie
            disp('for ICA analysis, further downsample the data by the factor of 2(spatial+temporal)')
            A_DS = Integration.downSampleMovie(A,2,2);
            %Restrict the area to only roi region
            [dim1_lower,dim1_upper,dim2_lower,dim2_upper] = movieData.getROIBoundsFromImage(A_DS(:,:,1));
            A_roi = A_DS(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:);
            sz = size(A_roi);
            A_re = reshape(A_roi,[sz(1)*sz(2),sz(3)]);

            if ~exist('ICA_results.mat','file')    
                %Get rid of nan in the matrix
                A_re_zv = A_re;
                A_re_zv(isnan(A_re_zv)) = 0;
                lastEig = max(round(sz(1)*sz(2)/800),50);
                tic;
                [icasig, M, W] = fastica(A_re_zv , 'verbose','off','lastEig', lastEig, 'numOfIC', 50);
                toc;
                save('ICA_results.mat','icasig','M','W')
            else
                load('ICA_results.mat');
            end

            %Generating correlation map with each independent component
            corr_map = corr(icasig',A_re');     
            
            disp('ICA processing is done')
        end
        
        
        function all_features = getFeatures2D(A)
        %Get 2D connected components features from binary matrix A
        %   Inputs:
        %       A          Input binary matrix
        %   Outputs:
        %       all_feature      Features from regionprops defined in the function

            sz = size(A);
            all_features = [];
            parfor i = 1:sz(3)
                cur_img = A(:,:,i);
                cur_CC =  bwconncomp(cur_img);
                cur_STATS = regionprops(cur_CC,'Centroid','Extrema','Eccentricity','Orientation'...
                    ,'FilledArea','MajorAxisLength','MinorAxisLength'); 
                for j = 1:length(cur_STATS)
                    curr_features = [cur_STATS(j).Centroid...
                        ,cur_STATS(j).Extrema(:)'...
                        ,cur_STATS(j).Orientation...
                        ,cur_STATS(j).Eccentricity ...
                        ,cur_STATS(j).FilledArea, cur_STATS(j).MajorAxisLength...
                        ,cur_STATS(j).MinorAxisLength];
                    all_features = [all_features; curr_features];
                end
            end
        end
        
        
        
        function [A, output_all] = dftReg(A, loadtag)
        %Use discrete fourier transform algo to register (rigid) input
        %matrix A. Output the registered matrix
        %
        %   Inputs:
        %     A                   3D matrix
        %     loadftag            tag for loading specific .mat file     
        %
        %   Outputs:
        %     A                   Registered matrix
        %     output_all          output from dftregistration
            
            if nargin < 2
                loadflag = 0;
            else
                loadflag = 1;
                disp(['loadtag detected: ' loadtag])
            end
        
            mask = isnan(A);
            A(mask) = 0;
            %Use the median frame as reference

            if loadflag
                if isfile(['Reference_frame_' loadtag '.mat']) 
                    load(['Reference_frame_' loadtag '.mat'])
                else
                    A_median = nanmedian(A,3);
                    save(['Reference_frame_' loadtag '.mat'],'A_median')
                end
            else
                A_median = nanmedian(A,3);
            end
            
            %Do 2d fft to A_median
            A_mf = fft2(A_median);
            parfor i = 1:size(A,3)
                A_c = A(:,:,i);
                %Do 2d fft to current frame
                A_cf = fft2(A_c);
                [output, Greg] = dftregistration(A_mf,A_cf,2);
                output_all(i,:) = output;
                Greg_all(:,:,i) = abs(ifft2(Greg));
            end
            A = Greg_all.*~mask;
        
            disp('Rigid registration finished...')
        end
        
                
                
        
        function [A, movTag, output_all, NormTform_all, movIdx_saved] = ...
                movAssess(A, param, loadtag)
        %   Assess movie and get rid of frames with large movements 
        %   Use discrete fourier transform algorithm to compare phase
        %   correlation. When filtering frames, use a step function as the
        %   convolution kernel (width == 5).
        %
        %   Inputs:
        %     A                   3D matrix
        %     param               parameters, should contain moveAssessFlag
        %                         outputFolder, filename, spacialFactor
        %
        %   Outputs:
        %     A                   Registered matrix
        %     movTag              Tag for generating filename
        %     output_all           The array that stores all transformation matrices
        %     NormTform_all       Norm of each matrices (minus I)
        %     movIdx_saved        Indices of matrix that will be saved
                          
            %Do dftregistration
            if nargin < 3
                %Whether provided loadtag or not
                [A_ori, output_all] = movieData.dftReg(A);
                if nargin < 2
                    moveAssessFlag = 0;
                    outputFolder = cd;
                    filename = 'temporaryFile.mat';
                    spacialFactor = 1;
                else
                    moveAssessFlag = param.moveAssessFlag;
                    outputFolder = param.outputFolder;
                    filename = param.filename;  
                    spacialFactor = param.spacialFactor;
                end
            else
                [A_ori, output_all] = movieData.dftReg(A, loadtag);
                moveAssessFlag = param.moveAssessFlag;
                outputFolder = param.outputFolder;
                filename = param.filename;  
                spacialFactor = param.spacialFactor;
            end
                           
            NormTform_all = sqrt(output_all(:,3).^2 + output_all(:,4).^2);
            
            if moveAssessFlag
                [A, movIdx_saved, ~] = movieData.discardFrames(A_ori, NormTform_all);
                %If there are frames being discarded set movTag to dsc
                disp('Trying to discard moving frames...')
                if size(A,3)== size(A_ori,3)
                    movTag = '';
                    disp('No frame discarded')
                else
                    movTag = 'dsc';
                    A_ori_DS = movieData.downSampleMovie(A_ori,spacialFactor);
                    A_ori_DS = reshape(A_ori_DS, [size(A_ori_DS,1)*size(A_ori_DS,2),...
                                size(A_ori_DS,3)]);
                    checkname = [filename(1:length(filename)-4) '_ori_DS_registered.mat'];
                    save(fullfile(outputFolder,checkname),'A_ori_DS','-v7.3');
                end
            else
                A = A_ori;
                movIdx_saved = ones(size(A,3),1);
                movTag = '';
            end
                   
            %if flag == 1 replace moving frames with mean-intensity frame
            %if flag == 1 discard the neighbouring 5 frames
            %A_mean = reshape(A_mean,[sz(1) sz(2)]);

            disp(['Mean tform magnitude (minus I) = ' num2str(mean(NormTform_all))]);
        end
        
      
        
        
        function [A, movIdx_saved, saveRatio] = discardFrames(A, NormTform_all)
        %Discard frames in the input movie A based on NormTform_all
        
        flag = 0; %default flag: do not discard frames
        sz = size(A);
        if length(NormTform_all) == sz(3)
            %If the norm is larger than 0.49 (0.5 pixel at either direction)
            %label as large-movement frame. Save frames that do not
            %move that much as movIdx_saved
            movIdx_saved = NormTform_all < 0.49;
            saveRatio = sum(movIdx_saved)/sz(3);
            disp(['saveRatio at 0.49 threshold: ' num2str(saveRatio)])
            filter = [1,1,1,1,1,1,1,1,1,1,1]; %Default filter: discard the neighbouring 11 frames 
            
            if saveRatio < 0.8
                warning('More than half of the movie will be discarded given current threshold!')
                disp('Increase motion detection threshold to 0.7!')
                disp('Double check movie quality!')
                %If the norm is larger than 1 (>1 pixel at both directions)
                %label as large-movement frame. Save frames that do not
                %move that much as movIdx_saved
                movIdx_saved = NormTform_all < 1.5;
                saveRatio = sum(movIdx_saved)/sz(3);
                disp(['saveRatio at 2.0 threshold: ' num2str(saveRatio)])
                filter = [1,1,1,1,1,1,1,1,1,1,1,1,1]; %Default filter: discard the neighbouring 13 frames
            end
            
            %If more than 5% of the movie have substantial movements, warn the user
            if (saveRatio < 0.95) || (max(NormTform_all) > 0.51)
                warning('This movie contains more than 5% moving frames/ movement with big jitters!')
                disp('Turn on discarding processing...')
                flag = 1;
            end
            
            if flag
                warning('Discarding moving frames here...')
                movIdx_replace =  ~movIdx_saved;
                movIdx_replace = conv(movIdx_replace, filter, 'same');
                movIdx_replace = movIdx_replace > 0;
                %A(:,:,movIdx_replace) = repmat(A_mean, [1,1,sum(movIdx_replace)]);
                A(:,:,movIdx_replace) = [];
                movIdx_saved = ~movIdx_replace;
            end   
            
            disp(['Relatively stable frames ratio = ' num2str(saveRatio)]);
            
        else
            warning('Movie size does not agree with movAsess file!')
            saveRatio = 1;
            movIdx_saved = NormTform_all > -1;
        end
     
        end
                
        
        
        function [A, movTag,output_all, NormTform_all, movIdx_saved] = ...
                movAssess_NoRMCorre(Yf, param)
        %   Assess movie and get rid of frames with large movements 
        %   Use NoRMCorre algo. See https://github.com/flatironinstitute/NoRMCorre
        %
        %   Inputs:
        %     Yf                  3D matrix
        %     param               parameters, should contain moveAssessFlag
        %                         outputFolder, filename, spacialFactor
        %
        %   Outputs:
        %     A                   Registered matrix
        %     movTag              Tag for generating filename
        %     output_all           The array that stores all transformation matrices
        %     NormTform_all       Norm of each matrices (minus I)
        %     movIdx_saved        Indices of matrix that will be saved 
        
            if nargin < 2
                moveAssessFlag = 0;
                outputFolder = cd;
                filename = 'temporaryFile.mat';
                spacialFactor = 1;
            else
                moveAssessFlag = param.moveAssessFlag;
                outputFolder = param.outputFolder;
                filename = param.filename;
                spacialFactor = param.spacialFactor;
            end

            
            try %In case rigid registration failed
                %Change from double to single
                Yf = single(Yf);
                [d1,d2,T] = size(Yf);

                % perform some sort of deblurring/high pass filtering
                if (0)    
                    hLarge = fspecial('average', 40);
                    hSmall = fspecial('average', 2); 
                    for t = 1:T
                        Y(:,:,t) = filter2(hSmall,Yf(:,:,t)) - filter2(hLarge, Yf(:,:,t));
                    end
                    %Ypc = Yf - Y;
                    bound = size(hLarge,1);
                else
                    gSig = 7; 
                    gSiz = 3*gSig; 
                    psf = fspecial('gaussian', round(2*gSiz), gSig);
                    ind_nonzero = (psf(:)>=max(psf(:,1)));
                    psf = psf-mean(psf(ind_nonzero));
                    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
                    %Y = imfilter(Yf,psf,'same');
                    %bound = 2*ceil(gSiz/2);
                    Y = imfilter(Yf,psf,'symmetric');
                    bound = 0;
                end
                % first try out rigid motion correction
                    % exclude boundaries due to high pass filtering effects
                options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);

                % register using the high pass filtered data and apply shifts to original data
                tic; [M1,shifts1,template1] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r); toc % register filtered data
                    % exclude boundaries due to high pass filtering effects
                tic; Mr = apply_shifts(Yf,shifts1,options_r,bound/2,bound/2); toc % apply shifts to full dataset
                    % apply shifts on the whole movie
                % compute metrics 
                [cY,mY,vY] = motion_metrics(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_r.max_shift);
                [cYf,mYf,vYf] = motion_metrics(Yf,options_r.max_shift);

                [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
                [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);

                %shifts_r = squeeze(cat(3,shifts1(:).shifts));
                try % In case non-rigid registration cannot run through
                    % now apply non-rigid motion correction
                    % non-rigid motion correction is likely to produce very similar results
                    % since there is no raster scanning effect in wide field imaging

                    options_nr = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',50, ...
                        'grid_size',[128,128]*2,'mot_uf',4,'correct_bidir',false, ...
                        'overlap_pre',32,'overlap_post',32,'max_shift',20);

                    tic; [M2,shifts2,template2] = normcorre_batch(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_nr,template1); toc % register filtered data
                    tic; Mpr = apply_shifts(Yf,shifts2,options_nr,bound/2,bound/2); toc % apply the shifts to the removed percentile

                    % compute metrics

                    [cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
                    [cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);

                    % plot shifts        

                    shifts_r = squeeze(cat(3,shifts1(:).shifts));
                    shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
                    shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
                    shifts_x = squeeze(shifts_nr(:,2,:))';
                    shifts_y = squeeze(shifts_nr(:,1,:))';
                    A_ori = Mpr;
                    disp('Non-rigid registration succeeded!')
                catch
                    disp('Non-rigid registration failed but rigid registration succeeded!')
                    A_ori = Mr;
                end
                %Save movement assessment result
                shifts_r = squeeze(cat(3,shifts1(:).shifts));
                NormTform_all = sqrt(shifts_r(:,1).^2 + shifts_r(:,2).^2);
                movIdx_saved = NormTform_all < 0.49;
                output_all = shifts_r;
                
                %Whether discard frames or not
                if moveAssessFlag
                    [A, movIdx_saved, ~] = movieData.discardFrames(A_ori, NormTform_all);
                    disp('Trying to discard moving frames...')
                    %If there are frames being discarded set movTag to dsc
                    if size(A,3)== size(A_ori,3)
                        movTag = '';
                        disp('No frame discarded')
                    else
                        movTag = 'dsc';
                        A_ori_DS = movieData.downSampleMovie(A_ori,spacialFactor);
                        A_ori_DS = reshape(A_ori_DS, [size(A_ori_DS,1)*size(A_ori_DS,2),...
                                    size(A_ori_DS,3)]);
                        checkname = [filename(1:length(filename)-4) '_ori_DS_registered.mat'];
                        save(fullfile(outputFolder,checkname),'A_ori_DS','-v7.3');
                    end
                else
                    A = A_ori;
                    movTag = '';
                    movIdx_saved = ones(size(A,3),1);
                end
            catch
                %If there is a glich...
                A_ori = [];
                A = [];
                NormTform_all = [];
                output_all = [];
                movIdx_saved = [];
                movTag = '';
                disp('NoRMCorre rigid registration failed!')
            end
            

            
            
        end
        
        
        function [A,tform_all,NormTform_all] = movAssessUsingEdge(A, flag, factor)
        %   Movement assessment using edge detection. Detect edge of each
        %   frames and the reference frame (median). Do registration using
        %   edges and detect movement spikes based on magnitute of
        %   transforming matrix. Discard neighbouring frames around the 
        %   spikes if flag == 1.
        %
        %   Inputs:
        %      A       input matrix, 3D
        %      flag    whether discard frames
        %      factor  fludge factor to control edge detection sensibility
        %
        %   Outputs:
        %      A       filtered matrix
        %      tform_all           Array stored all transformation matrices
        %      NormTform_all       Norm of each matrices (minus I)


        %fludge factor to control edge detection sensibility
        if nargin == 1
            factor = 0.5;
            flag = 0;
        end

        if nargin == 2
            factor = 0.5;
        end

        %Further downsample movie if it's too large
        %if (round(size(A,1)/100)>=2)||(round(size(A,2)/100)>=2)
        %    A = movieData.downSampleMovie(A, 2, 1);
        %    disp('Further downsample the movie for movement detection...')
        %end

        %Filtering out top 20% pixels in the median image
        sz = size(A);
        A_median = nanmedian(A,3);
        if ceil(sz(2)/sz(1)) >= 2
            threshold = 50;
        else
            threshold = 80;
        end
        Mask = A_median<prctile(A_median(:),threshold);

        %Only do edge detection in the region excluding the top 20% pixels
        %This will exclude pixels with high activity going on
        A_median_mask = A_median.*Mask;   
        A_mask = A.*Mask;

        %Generate threshold from the median image 
        [~,threshold] = edge(A_median_mask,'Canny');

        %For image dilation
        %se90 = strel('line',3,90);
        %se0 = strel('line',3,0);

        parfor i = 1:sz(3)
            A_cur = A_mask(:,:,i);
            BW = edge(A_cur,'Canny',threshold * factor);
            %BWdil = imdilate(BW,[se90 se0]);
            %BWs(:,:,i) = BWdil;
            BWs(:,:,i) = BW;
            if mod(i,1000) == 0
                disp(['Edge detection finished at frame ' num2str(i)])
            end
        end

        %Generate reference frame
        BW_ref = nanmedian(BWs,3);

        %Filtering connected componets with inconsistent duration
        [~,BWs] = Integration.GenerateCC(BWs,BWs,0, 1, round(size(BWs,3)./10));

        %Prepare for rigid registration
        [optimizer, metric] = imregconfig('monomodal');
        NormTform_all = zeros(1,sz(3));
        tform_all = cell(1,sz(3));

        parfor i = 1:sz(3)
            %Get transformation matrix using BW matrices
            curI = BWs(:,:,i);
            tform = imregtform(double(curI), BW_ref ,'translation', optimizer, metric);
            A_cur = A(:,:,i);
            A_cur(isnan(A_cur)) = 0;
            %Get translation vector
            T = tform.T;
            A_cur = imwarp(A_cur,tform,'OutputView',imref2d(size(BW_ref)));
            A_cur(A_cur == 0) = nan;
            A(:,:,i) = A_cur;
            %Calculate norm of the translation vector
            NormTform_all(i) = norm(T(3,1:2),'fro');
            tform_all{i} = tform;
        end

        %Detect movment spikes from background
        filter_sp = [-1,2,-1];
        NormTform_conv = conv(NormTform_all, filter_sp, 'same');

        %Top 10% (zscore > 1.5) are considered movement spikes 
        NTz = zscore(NormTform_conv);
        movSpikes  = abs(NTz) > 1.5;

        %Neibouring 7 frames of the spikes are considered moving frames
        filter_mov = [1,1,1,1,1,1,1];
        movFrames = conv(movSpikes, filter_mov, 'same');
        movFrames = movFrames > 0;
        saveRatio = 1 - sum(movFrames)/sz(3);

        %If more than 5% of movie have substantial movements, change
        %flag to 1 so moving frames will be replaced at the next step.
        if saveRatio < 0.95
            disp('This movie contains more than 5% moving frames!')
            disp('Consider to discard moving frames')
        end

        %if flag == 1 discard the neighbouring 7 frames
        if flag
            A(:,:,movFrames) = [];
            disp('DISCARD MOVING FRAMES!')
        end

        disp(['Mean tform magnitude (minus I) = ' num2str(mean(NormTform_all))]);
        disp(['Relatively stable frames ratio = ' num2str(saveRatio)]);

    end
        
        
        
        
        function corrM = corrMapThresholding(corrM,threshpct)
        %Thresholding correlation matrix, preserve top 20% correlated region in
        %each frame in default
        %
        %Inputs:
        %   corrM          correlation matrix
        %   threshpct      thresholding percentile
        %
        %Outputs:
        %   corrM          threholded matrix
        %

        if nargin == 1
            threshpct = 80;
        end

            parfor i = 1:size(corrM,3)
                %normalizing correlation based on distance between pixels and seeds
                %nomarlized correlation = correlation * exp((xi -x_seed)^2/2*delta)
                %in which delta = max(square(distances btw pixels)) 
                currMap = corrM(:,:,i);

                [seed_x, seed_y] = ind2sub(size(currMap),find(currMap == max(currMap(:))));
                [all_x, all_y] = ind2sub(size(currMap),1:length(currMap(:)));
                dis_sq = (all_x - seed_x).^2 + (all_y - seed_y).^2;
                dis_sq = reshape(dis_sq, size(currMap));
                delta_sq = max(dis_sq);
                e_term = exp(dis_sq./(2*delta_sq));
                nmMap = currMap.*e_term;

                %Preserve top 20% correlated region in default
                tholdedM = nmMap>prctile(nmMap(:),threshpct);
                corrM(:,:,i) = tholdedM;
                %imagesc([currMap;nmMap;tholdedM]); colormap jet; colorbar; axis image
                %caxis([-0.5, 1]); 
            end
        end
        
        function [mean_CCprops,allCCprops] = getAllCCfromCorrM(corrM)
        %Crop and reshape input correlation matrix and get all connectd components
        %properties from each frame of the input matrix
        %
        %Inputs:
        %   corrM          correlation matrix
        %
        %Outputs:
        %   allCCprops     all connected components properties
        %

            sz = size(corrM);
            %Get rid of maps that are all NaN (due to imperfect seeds sampling)
            corrM_re = reshape(corrM,[sz(1)*sz(2),sz(3)]);
            corrM_re = corrM_re(:,any(~isnan(corrM_re)));
            corrM = reshape(corrM_re,[sz(1),sz(2),size(corrM_re,2)]);
            %Confine the analysis only in non-nan (roi) region
            %corrM = movieData.focusOnroi(corrM);
            %Thresholding the matrix 
            thCorrM = movieData.corrMapThresholding(corrM,90);
            %CCfiltered_thCorrrM = Integration.filterCC_byFrame(thCorrM);
            CCfiltered_thCorrrM = thCorrM;
            allCCprops = movieData.getFeatures2D(CCfiltered_thCorrrM);
            mean_CCprops = mean(allCCprops);
            allCCprops = allCCprops - mean_CCprops;
        end
        
        
        
        function [trainSet,testSet] = generateBatch(A)
        %Generate training set and test set from permutated input matrix A
        %A can be 2D or 3D matrix
        %Inputs:
        %   A        2D or 3D matrix
        %
        %Outputs:
        %   trainSet    Subset of input A (75% of the data)
        %   testSet     Subset of input A (25% of the data)

            if length(size(A)) == 3
                n = size(A,3);
                %Permutation
                Idx = randperm(n);
                %training set contains 75% of input data
                trainSize = round(n/4*3);
                permA = A(:,:,Idx);
                trainSet = permA(:,:,1:trainSize);
                testSet = permA(:,:,trainSize+1:end);
            elseif length(size(A)) == 2
                n = size(A,2);
                Idx = randperm(n);
                trainSize = round(n/4*3);
                permA = A(:,Idx);
                trainSet = permA(:,1:trainSize);
                testSet = permA(:,trainSize+1:end);
            end

        end

        
        
        function [corrMatrix,seedsList] = seedsCorrelation(A, seedsList)
        %Generate seedsbased correaltion matrix from A given seedsList
        %
        %Inputs:
        %   A            3D matrix
        %   seedsList    list of given seeds (2D subscripts or 1D indices)
        %
        %Outputs:
        %   seedsList    renewed seedsList
        %   corrMatrix   correlation matrix

            try
                seedsTrace = [];
                sz = size(A);
                if length(size(seedsList)) == 2
                    seedsList = sub2ind(sz(1:2),seedsList);
                end
                imgall = reshape(A, sz(1)*sz(2), sz(3));
                for r = 1:length(seedsList)
                    seedTrace(r, :) = squeeze(nanmean(imgall(seedsList(r), :), 1));
                end
                corrM = corr(imgall',seedTrace');
                realSeeds = any(~isnan(corrM));
                corrM = corrM(:,realSeeds);
                seedsList = seedsList(realSeeds);
                sz_3 = size(corrM,2);
                corrMatrix = reshape(corrM,sz(1),sz(2),sz_3);
                save('CorrelationMatrix_provided_list.mat','corrMatrix');
            catch
                warning('SeedsList provids seeds that exceed current matrix dimensions!');
            end
        end
        
        
        function [corrMatrix, roi] = generateCorrMatrix(sz, roi, ...
                seedTrace, imgall, GPU_flag, mean_flag, Avg_out_all, timelag)
        %Generate correlation matrix using parameters specified by user
        %Inputs:
        %   sz          size of the matrix
        %   roi         roi defined by the user
        %   seedTrace   trace of the seed
        %   imgall      all pixels reshaped in 2D
        %   GPU_flag    whether use GPU for computing
        %   mean_flag   whether regress out background
        %   Avg_out_all background trace (if not provided put [])
        %   timelag     for timelag analysis
        %
            
       %Truncate traces using timelag
        [seedTrace, imgall] = movieData.timelagTruncate(seedTrace, imgall, timelag);          

       %Control on mean-activity-trace (lower 1% std) using partial correlation if mean_flag == 1
        switch mean_flag
            case 0 %Do not do partial correlation
                %Use GPU if GPU_flag == 1
                if GPU_flag
                    try
                        corrMatrix = movieData.batchSeedsGPU(sz,seedTrace,imgall);
                        disp('Run seeds based correlation on GPU')
                        warning('Can not do partial correlation when using GPU!')
                    catch
                        [corrM, pvalM] = corr(imgall',seedTrace');
                        [corrMatrix, roi] = filterNaNCorrMap(corrM, roi, sz);
                        disp('Run on GPU failed, run seeds based correlation on CPU')
                    end
                else
                    [corrM, pvalM] = corr(imgall',seedTrace');
                    [corrMatrix, roi] = filterNaNCorrMap(corrM, roi, sz);
                    disp('Run seeds based correlation on CPU')
                end
            case 1 %Do partial correlation using background trace or lowest 1% trace
                disp('Doing partial correlation on CPU...')
                if isempty(Avg_out_all) 
                   disp('Background trace not provided!')
                   disp('Define background trace as mean of the lowest 1% std pixels!')
                   std_all = nanstd(imgall,0,2);
                   lowPixels = std_all <= prctile(std_all,1);
                   %avg_trace = nanmean(imgall,1);
                   avg_trace = nanmean(imgall(lowPixels,:),1);
                   [corrM, pvalM] = partialcorr(imgall',seedTrace', avg_trace');
                else
                   [Avg_out_all, ~] = movieData.timelagTruncate(Avg_out_all', Avg_out_all', timelag);
                   [corrM, pvalM] = partialcorr(imgall',seedTrace', Avg_out_all');
                   disp('Generated correlation matrix with provided background trace!')
                end
                [corrMatrix, roi] = filterNaNCorrMap(corrM, roi, sz);
            case 2 %Do partial correlation using lowest 1% trace
                disp('Doing partial correlation using mean trace...')
                %avg_trace = nanmean(imgall,1);
                std_all = nanstd(imgall,0,2);
                lowPixels = std_all <= prctile(std_all,1);
                avg_trace = nanmean(imgall(lowPixels,:),1);
                [corrM, pvalM] = partialcorr(imgall',seedTrace', avg_trace');
                [corrMatrix, roi] = filterNaNCorrMap(corrM, roi, sz);
            end
                      
           %Create new savenames
           %c = clock;
           %timetag = [num2str(c(1)) num2str(c(2)) num2str(c(3)) num2str(c(4)) num2str(c(5))];
           %nametag = [num2str(GPU_flag) num2str(mean_flag) '_' num2str(timelag) '_' timetag];
           if mean_flag && ~isempty(Avg_out_all)
               extra_tag = '_withAvgOut';
           else
               extra_tag = '';
           end
           
           nametag = [num2str(GPU_flag) num2str(mean_flag) '_' num2str(timelag) '_' num2str(size(roi,1))  extra_tag];
           savename1 = ['Correlation_Matrix_' nametag '.mat'];
           savename2 = ['Seeds_' nametag '.mat'];
           
           try
               pvalM = reshape(pvalM, size(corrMatrix));
           catch
               disp('Can not reshape pval matrix!')
           end
 
           save(savename1,'corrMatrix','pvalM');
           save(savename2,'roi');
           
           %Build-in function to filter nan correlation map
           function [corrMatrix, roi]= filterNaNCorrMap(corrM, roi, sz)
            %Filter out correlation maps that are essentially nan
                realSeeds = any(~isnan(corrM));
                corrM = corrM(:,realSeeds);
                roi = roi(realSeeds);
                sz_3 = size(corrM,2);
                corrMatrix = reshape(corrM,sz(1),sz(2),sz_3);
           end
                     
        end
        

        
        function [seedTrace, imgall] = timelagTruncate(seedTrace, imgall, timelag)
        %Truncate the calcium traces for time-lag correaltion
       
            if nargin == 2
                %zero lag in default
                timelag = 0;
            end
            
            duration = size(seedTrace,2);
            
            %Allow postive or negative time-lag (lagging forward or backward) 
            if timelag >= 0
                %forward timelap
                seedTrace = seedTrace(:, 1 : duration - timelag);
                imgall = imgall(:, 1 + timelag : duration);
            else
                %backward timelap
                timelag = abs(timelag);
                seedTrace = seedTrace(:, 1 + timelag : duration);
                imgall = imgall(:, 1 : duration - timelag);
            end
        end
        
        
        
        function OutputMovie(A, filename, rate)
        % Write the input matrix (3D) into a .avi movie
        % Inputs:
        %       A          3D matrix
        %       filename   name of the output movie
        %       rate       frame rate

            if nargin < 2
                filename = 'Output';
            end

            if nargin < 3
                rate = 50;%default rate 50 frames/sec
            end

            writerObj = VideoWriter([filename '.avi']);
            writerObj.FrameRate = rate;
            % set the seconds per image

            % open the video writer
            open(writerObj);
            % write the frames to the video

            parfor i=1:size(A,3)
                % convert the image to a frame
                frame = A(:,:,i);
                figure('visible','off');
                imshow(mat2gray(frame))
                frame = getframe;
                F(i) = frame;
            end

            for i=1:length(F)
                % convert the image to a frame
                frame = F(i);    
                writeVideo(writerObj, frame);
            end

            % close the writer object
            close(writerObj);

        end
        
        
    end         
end
       