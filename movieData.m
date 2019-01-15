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
%R11 07/05/18 new roiSVD function, new grossDFoverF (mean->nanmean)
%R11 07/07/18 relax decremental factor starting from 1.6 in function bwThresholding_10prctPixels
%R11 07/11/18 new function makePseudoColorMovie
%R12 07/13/18 new function SeedBasedCorr and focusOnroi
%R12 07/14/18 Modify function SeedBasedCorr
%R12 08/30/18 Modify function SeedBasedCorr
%R12 08/31/18 Modify function SeedBasedCorr
%R13 09/15/18 Major improvement of SeedBased Corr, now can automatically
%generate seeds if there is not manually defined ones. Compatiable with ROI
%class R3 or higher 
%R13 09/19/18 Modify function SeedBasedCorr
%R14 09/19/18 Previous algo to calculate coordinates of downsampled seeds is
%wrong. Corrected in related functions in both ROI and movieData Class
%compatible with ROI R4 or higher
%R15 10/24/18 New function to allow seed based correlation maps generated on GPU
%Compatiable with byPassPreProcessing R4 or higher
%R16 12/12/18 New function maxPooling2D to allow maximum pooling
%R16 12/28/18 Renew the function makePseudoColorMovie by adding zscoring
%R16 01/03/19 Modify the roiSVD function
%R17 01/08/19 Add ICA analysis function getICA 
%R18 01/15/19 modify several functino to allow better parellel computing.
%Add new function plotCorrM 
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
                    disp(['Spaceil Factor: 2' ])
                    outputA = movieData.spatialDown(A,2);                    
                case 2
                    disp('Only do spatial downsampling here!')
                    disp(['Spaceil Factor:' num2str(s)])
                    outputA = movieData.spatialDown(A,s);
                case 3
                    disp(['Spaceil Factor:' num2str(s)])
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

            sz = size(A);
            newsz1 = mod(-sz(1),s)+sz(1);
            newsz2 = mod(-sz(2),s)+sz(2);
            IdxEnd1 = ceil(sz(1)/s)*s;
            IdxEnd2 = ceil(sz(2)/s)*s;
            A_newk = nan(newsz1,newsz2);
            A_newk_avg = nan(newsz1/s,newsz2/s,s^2);
            %downA = nan(newsz1/s,newsz2/s,sz(3));
            downA = [];
            
            if length(sz)<3
                sz(3) = 1;
            end

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
            parfor k = 1:sz(3)               
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
        
        function A = grossDFoverF(A)
        %    Doing gross dFoverF calculation.
        
            A_mean = nanmean(A,3);
            sz = size(A);
            A_re = reshape(A,[sz(1)*sz(2),sz(3)]);
            A_mean = repmat(reshape(A_mean,[sz(1)*sz(2),1]),[1,sz(3)]);
            A = reshape(A_re./A_mean - 1,sz);
            
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
            parfor i = 1:size(A,3)
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
        
        
        
        function [A,U,S,V] = roiSVD(A, iniDim)

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
            disp(['Preserved dimensions = ' num2str(psv_dim)]);
            %tic; [U,S,V] = svds(A_roi,psv_dim); toc;
            tic; [U,S,V] = svds(A_roi',psv_dim); toc;
            
            %Reconstruct the roi pairt using the the 4th to the last component
            A_roi_rcs = U(:,iniDim:end)*S(iniDim:end,iniDim:end)*V(:,iniDim:end)';
            %A_roi_rcs = reshape(A_roi_rcs,[sz(1),sz(2),sz(3)]);
            A_roi_rcs = reshape(A_roi_rcs',[sz(1),sz(2),sz(3)]);

            %Recover the roi to the original image size
            A(dim1_lower:dim1_upper,dim2_lower:dim2_upper,:) = A_roi_rcs;
            %A(A == 0) = nan;
            
            %Save U,S,V
            %U = reshape(U,[sz(1), sz(2), size(U,2)]);
            V = reshape(V,[sz(1), sz(2), size(V,2)]);
            %save('SVD_decomp.mat','U','S','V');
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

        
        function SeedBasedCorr_GPU(A,spatialFactor,total_seeds,GPU_flag)
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
        %
        % Outpus:
        %   seed-based correlation maps


            %downSampleRatio = 1/spatialFactor
            if exist('spatialFactor','var')
                downSampleRatio = 1/spatialFactor;
            else 
                downSampleRatio = 0.5;
            end

            if ~exist('total_seeds','var')
                total_seeds = 400;
            end

            if ~exist('GPU_flag', 'var')
                GPU_flag = 0;
            end

            disp(['Note: Defalut downsampled ratio == ' num2str(downSampleRatio)]);

            sz = size(A);
            imgall = reshape(A, sz(1)*sz(2), sz(3));

            %Detect seeds
            temp = ROI('Seeds.zip');
            ROI_all = temp.ROIData;
            if isempty(ROI_all)
                disp('Seeds not provided (no manually defined seeds!)')
                disp('Generating seeds evenly covering current rois...')
                roi = ROI.genSeedingROIs(total_seeds,downSampleRatio);
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
                        [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz./downSampleRatio);
                        disp('Program thinks that rois are generated based on ORIGINAL movie size')
                    else
                        disp('Program thinks that rois are generated based on DOWNSAMPLED movie size')
                    end

                catch
                    [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz./downSampleRatio);
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
                        mask = imresize(roi{r}, downSampleRatio, 'bilinear');
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
            if GPU_flag
                try
                    corrMatrix = movieData.batchSeedsGPU(sz,seedTrace,imgall);
                    disp('Run seeds based correlation on GPU')
                catch
                    corrM = corr(imgall',seedTrace');
                    realSeeds = any(~isnan(corrM));
                    corrM = corrM(:,realSeeds);
                    roi = roi(realSeeds);
                    sz_3 = size(corrM,2);
                    corrMatrix = reshape(corrM,sz(1),sz(2),sz_3);
                    disp('Run on GPU failed, run seeds based correlation on CPU')
                end
            else
                corrM = corr(imgall',seedTrace');
                realSeeds = any(~isnan(corrM));
                corrM = corrM(:,realSeeds);
                roi = roi(realSeeds);
                sz_3 = size(corrM,2);
                corrMatrix = reshape(corrM,sz(1),sz(2),sz_3);
                disp('Run seeds based correlation on CPU')
            end
            save('Correlation_Matrix.mat','corrMatrix');  
            
            if sflag == 1
                movieData.plotCorrM(roi, corrMatrix)
            elseif sflag == 2
                movieData.plotCorrM(roi, corrMatrix, 'roipolygon', roiPolygon, 'downsample', downSampleRatio)
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
        %       downsample    downsample ratio of the original matrix
        %       roipolygon    polygon region of the input seeds
        %
        % Outputs:
        %      correlation maps based on input information
           
            disp('If you have manually defined seeds, please do follow:')
            disp('1. Get ROI_all from your .zip roi file using ROI class')
            disp('2. [roi,roiPolygon] = ROI.generateROIArray(ROI_all,sz);')
            disp('plotCorrM(roi, corrMatrix, ''roipolygon'', roiPolygon value, ''downsample'', downsample ratio)')
            
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
                            downSampleRatio = varargin{i+1};
                    end
                end
            end
                     
                          
           %Plot correlation map
           disp(['Actual number of correlation maps = ' num2str(length(roi))])
            for r = 1:length(roi)
                h = figure; 
                set(gcf,'Visible', 'off');
                cur_img = corrMatrix(:, :, r);
                imagesc(cur_img); colormap jet; colorbar; axis image
                caxis([-0.5, 1]); title(['roi:', num2str(r)]);
                hold on                   

                %Label seed position
                if sflag == 2
                    if size(roi{r}, 1) ~= sz(1)
                        fill(roiPolygon{r}(:, 1).*downSampleRatio, roiPolygon{r}(:, 2).*downSampleRatio, 'y')
                    else
                        fill(roiPolygon{r}(:, 1), roiPolygon{r}(:, 2), 'y')
                    end
                elseif sflag == 1
                    %create a polygon region to label the location of
                    %current roi
                    maskId = roi(r,:);
                    x1 = ceil(maskId(1)/size(cur_img,1));
                    y1 = maskId(1) - floor(maskId(1)/size(cur_img,1))*size(cur_img,1);
                    fill([x1,x1,x1+2,x1+2],[y1,y1+2,y1+2,y1], 'y')
                end

                %Save the plot
                num_str = num2str(1000 + r);
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
                [icasig, M, W] = fastica(A_re_zv ,'lastEig', lastEig, 'numOfIC', 50);
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
                cur_STATS = regionprops(cur_CC,'Eccentricity','Orientation','Centroid'...
                    ,'Extrema','FilledArea','MajorAxisLength','MinorAxisLength'); 
                for j = 1:length(cur_STATS)
                    curr_features = [cur_STATS(j).Centroid,cur_STATS(j).Orientation...
                        ,cur_STATS(j).Eccentricity, cur_STATS(j).Extrema(:)'...
                        ,cur_STATS(j).FilledArea, cur_STATS(j).MajorAxisLength...
                        ,cur_STATS(j).MinorAxisLength];
                    all_features = [all_features; curr_features];
                end
            end
        end
  
    end
           
end
       