classdef movieData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   movieData stores the input matrix 
%   It has several static and normal functions to process the matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%R1 11/23/17 fetchAverageImage is a normal function 
%R2 11/24/17 dfOverFA is no longer a property of class movieData 
%R2 11/24/17 fetchAverageImage is a static function 
%R2 11/24/17 More built-in static functions for class movieData 
%R3 11/24/17 Major change in the way calculating df/f. Now using intra-video
%df/f. see InsideVideodFOverF and MeanDuringSti 
%R3 11/24/17 In gaussSmooth function, the zscoring step is commented out 
%R3 11/24/17 More built_in static functions for class movieData 
%R3 11/24/17 Rename classname from movie to movieData due to conflict with
%matlab default function movie 
%R3 11/24/17 Remove function MoviedFoverF 
%R3 11/24/17 Modify tophat function, remove the ROI checking line 
%R4 11/25/17 Modify tophat function, remove index variable in the function 
%R5 11/25/17 Modify insideVideodFoverF function by adding 'outputVideo(isnan(outputVideo)) = 0;' 
%R6 12/01/17 Change function bwthresholding's name to simpleBwThresholding 
%R6 12/01/17 Add new function bwThresholding 
%R6 12/01/17 Modify the way to aquire variable nmov in function AverageAcrossMovies 
%!!!THE PREVIOUS MODIFICATION MAKING THE CLASS COMPATIBLE WITH VERSIONS after
%audpipe R4 and Integration R2!!!
%R7 12/01/17 Add grossDFoverF 
%R7 12/11/17 Complement comments for several important function 
%R7 01/11/18 Revise the bwThresholding function, add the condition (factor
%>=0) to stay in the while loop 
%R8 05/24/18 Conditioning on flag in the class constructor, in order to
%prevent double-call in the base class. Compatible with Integration R5+/ROI R1+/
%audPipe R5+ Incorporate wavelength switching. This version is R8 
%R8 Directly input nmov when using the function AverageAcrossMovies 
%R9 New bwThresholding function (fixed threshold, defalut = 2 std above mean)
%Compatiable with Integration R7          
%R10 06/13/18 New static functions allowing rigid registration and 
%bleaching correction        
%R11 06/19/18 SVD deposition (to roi part, faster)
%R11 07/05/18 new roiSVD function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    properties
        A;   %Input matrix        
    end
    
    methods
        
        function obj = movieData(filename,flag)
        %movieData object constructor
            if ~(flag == 2)
                obj.A = movieData.inputMovie(filename);
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
            tmp = dir([filename(1:length(filename)-4) '@00*.tif']);
            %Make a combined matrix for each recording
            if ~isempty(tmp)
                for j = 1:numel(tmp)
                    fn = tmp(j).name;
                    B = openMovie(fn);
                    A = cat(3, A, B);
                    clear B
                end
            end             
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
                sigma = 3;
            end
                                   
            sz = size(A);
            %A = zscore(reshape(A,sz(1)*sz(2)*sz(3),1)); %covert to standardized zscores so that mean=0, and sd = 1;
            A = reshape(A,sz);   
            Iarr = zeros(size(A));

            for fr = 1:sz(3) %option:parfor, do NOT use parfor it is 3x slower
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
        
        
        function MeanFrame = MeanDuringSti(FreqVideo)
        
        %    Calculate the mean of 20th-29th frames in a 40-frame movie.
        %   (during-stimulaton period)
        %    
        %    Inputs:
        %       FreqVideo     Averaged movie containing 40 frames
        %
        %    Outputs:
        %       MeanFrame       mean intensities averaged over frames 20-29
        
            MeanFrame = mean(FreqVideo(:,:,20:29),3);
        end
        
        
        %Intra-trails dF/F
        function outputVideo = InsideVideodFoverF(curVideo)
        
        %    Calculate the df/f for each frame in a 40-frame video
        %    grouped by either frequencies or frequency/volume
        %    combination. Take the average of the pre-stimulation 20
        %    frames as background. Using this averaged matrix to
        %    calculate df/f for each frame.
        %
        %    Inputs:
        %        curVideo     currently working video
        %    
        %    Outputs:
        %        outputVideo     df/f video
        
            meanVideo = mean(curVideo(:,:,1:19),3);
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
        
        
        
        function FreqVoluColorMap(A,filename,nmov)
            
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
        %        Averaged df over f images over frames 20-29
        %        Averaged 40-frame para-stimulation video clips
        
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
            MaxFrameLength = ceil(max(temp_length)/40)*40;
            
            %FreqVoluPic = zeros(sz(1),sz(2));
            FreqVoluVideo = zeros(sz(1),sz(2));
            
            for i = 1:FreqNumber
                for j = 1:VoluNumber          
                    temp_Idx = FramesByFreqVolu{j,i};
                    ThisSignature = temp_Idx(1:2,1);
                    temp_Idx = temp_Idx(3:size(temp_Idx,1),1);
                    temp_Idx = [temp_Idx;zeros(MaxFrameLength - length(temp_Idx),1)];
                    temp_Idx = reshape(temp_Idx,40,[]);

                    for k = 1:40
                        temp_Idx_k = temp_Idx(k,:);
                        temp_Idx_k(temp_Idx_k ==  0) = [];
                        NonZero = temp_Idx_k;
                        FreqVoluVideo(:,:,k) = mean(A(:,:,NonZero),3);
                    end
                    
                    cur_dFoverF = movieData.InsideVideodFoverF(FreqVoluVideo);
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
            
            if exist(current_name,'file')
                load(current_name);
                AveragedMatrix = AveragedMatrix + CurrentMatrix./nmov;
                save(current_name,'AveragedMatrix');
            else
                AveragedMatrix = CurrentMatrix./nmov;
                save(current_name,'AveragedMatrix');
            end

        end
        
        
        function AveragedMapsAcrossMovies()
            
            %    Calling this function will read in all
            %    AveragedMatrix_freq_volu.mat data (should be updated with
            %    all movies in advance, see the AverageAcrossMovies
            %    function) and output maps and videos averaged across movies
            %    It should be called at the end of the main pipeline. 
            

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
                AveragedMap = movieData.MeanDuringSti(AveragedMatrix);
                pictype = mat2gray(squeeze(AveragedMap));
                [pictypeInd, ~] = gray2ind(pictype,65536);
                tmp_filename = [filename(1:length(filename)-4) '.png'];
                imwrite(pictypeInd, tmp_filename);
                tmp_filename = [filename(1:length(filename)-4) '.avi'];
                %VideoDFoverF = movieData.InsideVideodFoverF(AveragedMatrix); %modified in R3
                [~, Iarr] = timeColorMapProj(squeeze(AveragedMatrix),1,40,tmp_filename);
                Iarr2avi(Iarr,1,40,tmp_filename);

            end

        end
        
                
        %%Downsample
        function outputA = downSampleMovie(A,s,t)
        
        %    Dowansample the input matrix A
        %    Temporal downsampling: choose frames spanning by a specified factor
        %    Spacial downsampling: averaged pixels in a s-by-s window
        %
        %    Inputs:
        %        A   Input movies(3D matrix)
        %        s   spacial downsampling factor 
        %        t   temporal downsampling factor
        %
        %    Outputs:
        %        outputA    downsampled matrix

            switch nargin
                case 1
                    disp(['Spaceil Factor: 2' ])
                    outputA = movieData.spacialDown(A,2);                    
                case 2
                    disp('Only do spacial downsampling here!')
                    disp(['Spaceil Factor:' num2str(s)])
                    outputA = movieData.spacialDown(A,s);
                case 3
                    disp(['Spaceil Factor:' num2str(s)])
                    disp(['Temporal Factor:' num2str(t)])
                    outputA = movieData.temporalDown(A,t);
                    outputA = movieData.spacialDown(outputA,s);
                    
            end
        end


        function downA = spacialDown(A,s)
        
        %    Downsample the input matrix A by factor s. The idea is assign a new pixel
        %    by averaging pixels within a s-by-s whindow.
        %
        %    Inputs: 
        %        A   input matrix
        %        s   spacial downsampling factor
        %
        %    Output:
        %        downA    spacially downsampled matrix

            sz = size(A);
            newsz1 = mod(-sz(1),s)+sz(1);
            newsz2 = mod(-sz(2),s)+sz(2);
            IdxEnd1 = ceil(sz(1)/s)*s;
            IdxEnd2 = ceil(sz(2)/s)*s;
            A_newk = nan(newsz1,newsz2);
            A_newk_avg = nan(newsz1/s,newsz2/s,s^2);
            downA = nan(newsz1/s,newsz2/s,sz(3));
            downA = [];

            for k = 1:sz(3)    
                A_k = A(:,:,k);
                A_newk(1:sz(1),1:sz(2)) = A_k;
                count = 0;
                for i = 1:s
                    Idx1 = [1+i-1:s:IdxEnd1];
                    for j = 1:s
                        count = count + 1;
                        Idx2 = [1+j-1:s:IdxEnd2];
                        A_newk_avg(:,:,count) = A_newk(Idx1,Idx2);
                    end
                end
                A_newk_avg = nanmean(A_newk_avg,3);
                downA(:,:,k) = A_newk_avg;
            end

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
                if ~isempty(findstr(temp_info(i,1).name,'AveragedMatrix_')) && isempty(findstr(temp_info(i,1).name,'threshfactor')) %#ok<FSTR>
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
            
            singleThresholding = @(A, thresh) reshape(zscore(A(:)),size(A))>thresh;
            for k = 1:sz(3)               
                A(:,:,k) = singleThresholding(A(:,:,k),thresh);                            
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
            bgMatrix = A(:,:,:);
            bgMatrix(bgMatrix == 0) = nan;
            bgMean = nanmean(bgMatrix(:));
            bgStd = nanstd(bgMatrix(:));
            if nargin == 1
                factor = 1;
            end
            
            disp(['Thresholding factor = ' num2str(factor)]);
            
            %filter out pixels <= mean + factor * std
            for i = 1:sz(3)
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
            for i = 1:sz(3)
                %Acuqire current matrix
                currMatrix = A(:,:,i);
                currMatrix(currMatrix == 0) = nan;
                currMean = nanmean(currMatrix(:));
                currStd = nanstd(currMatrix(:));
                factor = 2;
                
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
        
        function A = grossDFoverF(A)
        %    Doing gross dFoverF calculation.
        
            A_mean = mean(A,3);
            for i = 1:size(A,3)
            	A(:,:,i) = (A(:,:,i) - A_mean)./A_mean;
            end
        end
        
        
        function A = movieRigReg(A_fixed, A)
        
        %    Do rigid registration to input matrix A, actually only do
        %    translation to save time
        %    
        %    Inputs:
        %    A        Original movie without rigid registration
        %    
        %    Outputs:
        %    A        Movie after rigid registration
 
            [optimizer, metric] = imregconfig('monomodal');
            tic;
            for i = 1:size(A,3)
                if ~mod(i,400)
                    disp(['Finish rigid registration at frame #' num2str(i)]);
                end
                A(:,:,i) = imregister(A(:,:,i), A_fixed, 'translation', optimizer, metric);
            end
            toc;
        end
        
        
        function A = bleachCorrection(A)
        
        %    Correcting fluorescence bleaching, assuming exponential decrease
        %   
        %   Inputs:
        %        A          input matrix
        %   
        %    Outputs:
        %        A          bleaching corrected matrix
            
            mean_series = zeros(1,size(A,3));
            x = [1:size(A,3)];
            
            for i = 1:size(A,3)
                cur_A = A(:,:,i);
                mean_series(i) = nanmean(cur_A(:));
            end
            
            f = fit(x',mean_series','exp1');
            trend = f.a.*exp(f.b.*x);
            trend = trend./min(trend);
            
            for i = 1:size(A,3)
                A(:,:,i) = A(:,:,i)./trend(i);
            end
                
        end
        
        
        
        function A = roiSVD(A, iniDim)

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
                iniDim = 4;
            end
            
            %Identify the vertex of the minimum rectangle containing roi
            cur_img = A(:,:,1);
            [dim1_lower,dim1_upper,dim2_lower,dim2_upper] = movieData.getROIBoundsFromImage(cur_img);

            %Do SVD to the roi only
            A_roi = A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:);
            sz = size(A_roi);
            A_roi = reshape(A_roi,[sz(1)*sz(2),sz(3)]);
            A_roi(isnan(A_roi)) = 0;
            psv_dim = round(size(cur_img,1)*size(cur_img,2)/800); %empirically suffcient to have a good-quality reconstruction
            disp(['Preserved dimensions = ' num2str(psv_dim)]);
            tic; [U,S,V] = svds(A_roi,psv_dim); toc;

            %Reconstruct the roi pairt using the the 4th to the last component
            A_roi_rcs = U(:,iniDim:end)*S(iniDim:end,iniDim:end)*V(:,iniDim:end)';
            A_roi_rcs = reshape(A_roi_rcs,[sz(1),sz(2),sz(3)]);

            %Recover the roi to the original image size
            A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:) = A_roi_rcs;
            %A(A == 0) = nan;

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

        
        
    end    
end
       