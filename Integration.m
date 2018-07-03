classdef Integration < spike2 & baphy & movieData & Names & ROI & wlSwitching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Integration is a subclass that inherits from all superclasses to
%analyze auditory experiment data. Specifically, Integration inherits
%from spike2\baphy\movidData\Names\ROI\wlSwitching classes
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Warning:
%Change backslashes in any path to forwardslashes when running on HPC!!! 
%Warning lifted:
%05/28/18 now the program can handle both cases (UNIX or WINDOWS syntax)    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R2 12/01/2017 Add a pipeline processing function to differentially
%process different data sets (for future spontaneous experimental set-up
%and stimuli-evoked set-up) 
%!!!NOTE THAT INTEGRATION R2 AND HIGHER ARE ONLY COMPATIABLE WITH 
%movieData R6 AND HIGHER!!!
%R3 12/08/2017 Add fileDetector function to allow automatically txt files
%generation 
%R3 12/08/2017 Add try catch module to tackle idx mismatch 
%R4 12/08/2017 Impose dFOverF to tophat filtered matrix 
%R4 12/12/2017 Add function minSpotSize to compute the minimum size of CC
%to be preserved. 
%R4 12/12/2017 Move makeCCMovie function to static 
%R4 05/23/2018 Tolerant to inperfect movie recording (when less than 6000
%frames were recorded) 
%R5 05/24/2018 New superclass wlSwitching. New static functions
%processMovies, doAveragingAcrossMovies. Compatible with audiPipe
%R5+/movieData R8+/ROI R1+/wlSwitching R1+.   This version is R5
%R5 05/25/2018 Store nmov as a class properties (for direct input to the 
%function AverageAcrossMovies in superclass movieData 
%R5 To better deal with file-separator problem in Unix v.s. Windows use
%fullfile function 
%R6 05/31/18 To adapt with different baphy configurations (different baphy 
%trail blocks and movie frames)  
%R7 06/06/18 Compatiable with ROI class R2, movieData R9           
%R8 06/13/18 Modify prePipe function, compatiable with movieData R10 
%(rigid registration and photo bleaching correction)                           
%R9 06/19/18 Add SVD denosing comopatiable with movieData R11           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        flag; %0, for spontaneous activity without stimuli; 1, for activity with stimuli; 2, for wavelength switching 
        nmov; %Number of movies
    end
    
    methods
        function obj = Integration(f,flag,nmov)
        %Integration object constructor
            
            if nargin == 1
                %With external stimulation
                flag = 1;
                disp('Assuming Spike2 and Bahpy files are provided')
            end
            
            %Get corresponding file names 
            Namesobj = Names(f,flag);
            filename = Namesobj.filename;
            Spike2name = Namesobj.Spike2name;
            Baphyname = Namesobj.Baphyname;
            Audioname = Namesobj.Audioname;
            ROIname = Namesobj.ROIname;
          
            %Wavelength switching
            if flag == 2
                movieObj = movieData(filename,1); A_ = movieObj.A; %Here if input flag== 2, will not further process
                ROIObj = ROI(); ROI_ = ROIObj.ROIData; %Here if input flag == 2, will not further process
            else
                A_ = []; ROI_ = [];
            end           
            
            %To get blcokN for spike2 object construction
            baphyObj = baphy(Baphyname);
            
            %Subclass objects construction
            obj@Names(f,flag);
            obj@movieData(filename,flag);
            obj@spike2(Spike2name,baphyObj.blockN);
            obj@baphy(Baphyname);
            obj@ROI();
            obj@wlSwitching(flag,A_,ROI_);
          
            obj.flag = flag;
            obj.nmov = nmov;
            
            clear movieObj ROIObj Namesobj baphyObj
          
        end
        
        function FramesByFreqVolu = prePipe(obj,spacialFactor,reg_flag)
        
        %    Preprocessing pipeline to handle a single auditory project
        %    movie. Process includes: photobleaching correction,
        %    rigid registration, applying ROI masks, downsampling,
        %    top-hat filtering, frequency-frame alignment(optional),
        %    frequency-volume maps plotting(optional),gaussian smoothing, 
        %    black-white thresholding, connected component extraction etc.            
        %    
        %    Inputs:
        %        obj               current object
        %        spacialFactor     for spacial downsampling
        %        reg_flag          rigid registration flag
        %   
        %   Outputs:
        %        FramesByFreqVolu  sorted frames by frequency and volume
        
            currentFolder = pwd;
            filename = obj.filename;        
            typeString = 'filtered';
            checkname = [filename(1:length(filename)-4) '_' typeString '.mat'];
            disp(['Processing: ' filename]);
            disp(' ');
 
            %Check the flag to decide whether do frames-frequency-volume
            %alignment
            if obj.flag
                
                try
                    FramesByFreqVolu = obj.FramesByFreqAlign;
                catch ME
                    %Handle damaged spike2 files problem
                    if (strcmp(ME.identifier,'MATLAB:badsubscript'))
                        disp(' ')
                        disp('Index exceeds matrix dimensions. Check Spike2!!')
                        disp('Skip alignment...')
                        disp(' ')
                    end
                    obj.flag = ~obj.flag;
                    FramesByFreqVolu = 0;
                end
            else
                FramesByFreqVolu = 0;
            end
            
            %Correct photobleaching
            A_corrct = Integration.bleachCorrection(obj.A);
            disp('Photobleaching corrected!');
            
            % Decide whether do rigid registration or not
            if reg_flag 
                if ~exist('Fixed_frame.mat','file')
                    A_fixed = A_corrct(:,:,1);
                    save('Fixed_frame.mat','A_fixed');
                else
                    load('Fixed_frame.mat');
                end
                A_corrct = movieData.movieRigReg(A_fixed,A_corrct);
            end
            
            %Apply ROI mask(s)
            A_ROI = Integration.ApplyMask(A_corrct,obj.ROIData);
            clear A_corrct
            
            %Downsampling
            A_DS = Integration.downSampleMovie(A_ROI,spacialFactor);
            clear A_ROI
            
            %SVD denosing of down-sampled A
            de_A = Integration.roiSVD(A_DS);
            disp('SVD denosing is done')
            disp('')
            clear A_DS
            
            %Top-hat filtering
            TH_A = Integration.TopHatFiltering(de_A);
            disp('Top-head filtering is done');
            disp('')
            clear de_A
            
            %Check flag to decide whether to generate frequency/volume maps
            if obj.flag
                %Integration.FreqColorMap(TH_A,filename,obj);
                Integration.FreqVoluColorMap(TH_A,filename,obj.nmov);
            end
            
            %Impose dFOverF to top-hat filtered matrix
            TH_A = Integration.grossDFoverF(TH_A);
            disp('Gloabal dFoverF is done')
            disp(' ')
            
            %Gaussian smoothing
            Ga_TH_A = Integration.GauSmoo(TH_A,2); %set sigma = 2
            disp('Gaussian smoothing is done');
            disp(' ')
            %clear TH_A
            
            %Save filtered matrix
            outputFolder = fullfile(currentFolder,obj.outputFolder);
            mkdir(outputFolder);
            save(fullfile(outputFolder,checkname),'Ga_TH_A','-v7.3');
            
            %Black-white thresholding of pre-processed A
            [BW_ppA,~] = Integration.bwThresholding(Ga_TH_A); %ppA is short for pre-processed A
            clear GA_TH_A;
            
            %Generate connected component
            Integration.GenerateCC(TH_A,BW_ppA,filename)
            clear TH_A BW_ppA

            disp(['Preprocessing done: ' filename]);
            disp('')
        end
        
       
        
        function FramesByFreqVolu = FramesByFreqAlign(obj)
        
        %    Align frames by frequencies being presented.
        %    
        %    Inputs:
        %        obj     Current object
        %    
        %    Outputs:
        %        .mat files that store ouput data 
        %       
         
                               
            filename = obj.filename;
            Idx_baphstar = obj.Idx_baphstar;
            Idx_framesEy = obj.Idx_framesEy;
            savename = [filename(1:length(filename)-4) '_Spike2BaphyIndex.mat'];
          
            %Handle imperfect recording, preserve good trails
            ActualTrailNumber = min([size(Idx_baphstar,2),size(obj.txt_line,1),floor(size(obj.A,3)/40)]);
            disp(['[Baphstar,txt_line,A_size3]= ' num2str([size(Idx_baphstar,2),size(obj.txt_line,1),floor(size(obj.A,3)/40)])])
            disp(['Usable Trails = ' num2str(ActualTrailNumber)]);
            txt_line = obj.txt_line(1:ActualTrailNumber,:);
            [freq,volu] = Integration.GetVoluFreq(txt_line);
            obj.freq = freq;
            obj.volu = volu;
            
            %Remove extra indices from spike2 if it goes beyond the actual
            %movie frames
            Idx_baphstar = Idx_baphstar(:,1:ActualTrailNumber); 
            Idx_framesEy = Idx_framesEy(:,1:40*ActualTrailNumber);
            
            Idx_baphstar(3,:) = freq';
            Idx_baphstar(4,:) = volu';

            %Whether camera frames are enclosed by squrewaves
            for t = 1:size(Idx_framesEy,2)

                if Idx_baphstar(2,1)<Idx_framesEy(1,t)
                    Idx_baphstar(:,1) = [];
                end
                if Idx_baphstar(1,1)<Idx_framesEy(1,t) && Idx_baphstar(2,1)>Idx_framesEy(1,t)
                    Idx_framesEy(2,t) = Idx_baphstar(3,1);
                    Idx_framesEy(3,t) = Idx_baphstar(4,1);
                end

            end
            
            %For some rare cases, seems that camera frames are present
            %prioir to baphstar square waves. Debug this problem here
            tmp_id = find(~Idx_framesEy);
            for i = 1:length(tmp_id)
                Idx_framesEy(tmp_id(i)) = Idx_framesEy(tmp_id(i)+3);
            end

            save(savename,'Idx_framesEy');
            
            %Get list of frequencies and volumes that have been presented
            FreqList = unique(freq);
            VoluList = unique(volu);

            if exist('FreqVoluList.mat')
                load('FreqVoluList.mat')
                FreqVoluList{size(FreqVoluList,1)+1,1} = FreqList;
                FreqVoluList{size(FreqVoluList,1)+1,1} = VoluList;
                save('FreqVoluList.mat','FreqVoluList')
                %fprintf('exsit \n')
            else
                FreqVoluList = {};
                FreqVoluList{1,1} = FreqList;
                FreqVoluList{2,1} = VoluList;
                save('FreqVoluList.mat','FreqVoluList');
                %fprintf('new \n')
            end
            
            %Construct FramesByFreq that store sorted frames (by freq)
            FramesByFreq = cell(2,size(FreqList,1));
            for k = 1 : size(FreqList,1)
                FramesByFreq{1,k} = FreqList(k,1);
                FramesByFreq{2,k} = find(Idx_framesEy(2,:) == FreqList(k,1))';
            end
            FreqFilename = [filename(1:length(filename)-4) '_FramesByFreq.mat'];
            save(FreqFilename,'FramesByFreq');


            %Construct FramesByFreq that store sorted frames (by freq&volu)
            FramesByFreqVolu = cell(size(VoluList,1),size(FreqList,1));
            for i = 1 : size(FreqList,1)
                temp_Idx1 = find(Idx_framesEy(2,:) == FreqList(i,1));
                for j = 1: size(VoluList,1)
                    temp_Idx_framesEy = [temp_Idx1; Idx_framesEy(2:3,temp_Idx1)];
                    temp_Idx2 =(temp_Idx_framesEy(3,:) == VoluList(j,1));
                    FreqVoluSpecIdx_framesEy = temp_Idx_framesEy(:,temp_Idx2);
                    temp_block = [FreqVoluSpecIdx_framesEy(2,1);FreqVoluSpecIdx_framesEy(3,1);...
                        FreqVoluSpecIdx_framesEy(1,:)'];
                    FramesByFreqVolu{j,i} = temp_block;
                end
            end
            FreqVoluFilename = [filename(1:length(filename)-4) '_FramesByFreqVolu.mat'];
            save(FreqVoluFilename,'FramesByFreqVolu');


        end
                     
    end
    
    
    methods(Static)
        
        function A3 = makeCCMovie(filename,CC,sz)
        %make connected component movie
            A3 = false(sz);           
            for i = 1:CC.NumObjects
                A3(CC.PixelIdxList{i}) = true;
            end
            M(size(A3,3)) = struct('cdata',[],'colormap',[]);
            for fr=1:size(A3,3)
                [I2, map] = gray2ind(A3(:,:,fr), 8); %figure; imshow(I2,map)
                M(fr) = im2frame(I2,map);
            end
            writeMovie(M,filename(1:length(filename)-4));
        end
        
        function CC = ccThresh(CC,minSpotSize)
        
        %    Thresholding the connected componnents based on each
        %    component's size and temporal duration
        %    
        %    Inputs: 
        %       CC     Connected components from bwconncomp
        %       miniSpotSize    Thresholding size factor
        %
        %   Outputs:
        %      CC     Renewed CC
        
            CC = Integration.ccSizeThresh(CC,minSpotSize);
            CC = Integration.ccTimeThresh(CC);
        end
        
        
        function CC = ccSizeThresh(CC,minSpotSize)
        
        %    Thresholding the CC based on length of each PixelIdxList
        %    If length is smaller than minimum spot size, kill it.
        %    
        %    Inputs: 
        %    CC     Connected components from bwconncomp
        %    miniSpotSize    Thresholding size factor
        %
        %    Outputs:
        %    CC     Renewed CC
        
            sz = cellfun(@length,CC.PixelIdxList);
            %minSpotSize = prctile(sort(sz),90);
            kill = sz<minSpotSize;
            CC.PixelIdxList(kill) = [];
            CC.NumObjects = length(CC.PixelIdxList);     
        end
        
        
        function CC = ccTimeThresh(CC)
        
        %    Thresholding the CC based on each component's temporal
        %    duration, which is the length of the BoundingBox side along the
        %    time dimension. 
        %   
        %    Inputs:
        %    CC     Connected components from bwconncomp
        %    
        %    Outputs:
        %    CC     Renewed CC
        
            STATS = regionprops(CC,'Area','BoundingBox');
            roiBoundingBox = zeros(length(STATS),6);
            for i = 1:length(STATS)
                roiBoundingBox(i,:) = STATS(i).BoundingBox;
            end
            durations = roiBoundingBox(:,6);
            newPixelIdxList = CC.PixelIdxList;
            
            %duration should be larger/equal to 2 frames
            newPixelIdxList(durations(:)<2) = [];
            CC.PixelIdxList = newPixelIdxList;
            CC.NumObjects = length(newPixelIdxList);                 
        end
        
        
        function fileDetector()
        
        %    Generate txt files based on file detected. Generate files.txt as a list
        %    of names of .tif movies. Generate Spike2files.txt as a list of names of 
        %    .mat spike2 data. Generat baphyfiles.txt as a list of names of baphy
        %    files. 
        %    
        %    Warning: The filenames of spike2/tif moives/baphys should be in
        %    a descending order accordingly. Otherwise the matching will be
        %    incorrect!!!
        

            tmp_list = dir;
            screen_types = {'.mat','.m','.tif'};
            disp(['Working Folder: ' pwd])
            tifID = fopen('files.txt','w+');
            spike2ID = fopen('Spike2files.txt','w+');
            baphyID = fopen('baphyfiles.txt','w+');

            for i = 1 : length(tmp_list)

                curName = tmp_list(i).name;
                
                %skip '.' and '..'
                if length(curName) < 2
                    continue
                end

                last2chars_name = curName(end-1:end);
                count = 0;
                flag = false;
                
                %Write in different .txt files according to the type of the
                %detected files
                while (count < length(screen_types))&&(~flag)
                    count = count + 1;
                    curr_type = screen_types{count};
                    last2chars_type = curr_type(end-1:end);

                    if last2chars_name == last2chars_type
                        flag = true;
                        switch curr_type
                            case '.mat'
                                fprintf(spike2ID,'%s\r\n',curName);
                            case '.m'
                                fprintf(baphyID,'%s\r\n',curName);
                            case '.tif'
                                if isempty(strfind(curName,'@00'))
                                    fprintf(tifID,'%s\r\n',curName);
                                end
                        end
                    end

                end
            end

            fclose('all');
        end
        
        function minSize= minSpotSize(CC,sz)
        
        %    Compute the minimum size of CC to be preserved. 
        %
        %    Inputs:
        %        CC      connected components
        %        sz      size of the video
        %
        %    Outputs:
        %        minSize     minimum size chosen for thresholding
        
            
            PixelIdxList = CC.PixelIdxList;
            sizeCounts = cellfun(@length,PixelIdxList);
            %Choose the samller one btw 95% percentile of all CC sizes and 1% of image size 
            minSize = min(prctile(sizeCounts,95),sz(1)*sz(2)/100);
        
        end
        
        function [Idx,Idx1,Idx2] = processMovies(f,flag,spacialFactor,nmov,reg_flag)
        
        %    Proccess movies (called by audPipe) seuqentially. This can
        %    handle different scenarios providing different flags.
        %    
        %    Inputs:
        %        f      The f-th movie to process
        %        flag   Represents different experimental setups
        %        spacialFactor     For downsampling
        %    
        %    Outputs:
        %        Idx,Idx1,Idx2    Frames indices sorted by frequencies or volume
        
            Idx = []; Idx1 = []; Idx2 = [];
            if flag == 2
                IntgA = Integration(f,flag,nmov);
                IntgA.A = [];
                IntgA.ROIData = [];
                
                %Do channel 1 
                IntgA.A = IntgA.A1;
                IntgA.ROIData = IntgA.ROI1;
                IntgA.ROI1 = [];
                IntgA.A1 = [];
                currentFolder = pwd;
                newFolder = fullfile(currentFolder,'Channel_1');
                if ~exist(newFolder,'dir')
                    mkdir(newFolder);
                end
                cd(newFolder)
                Idx1 = IntgA.prePipe(spacialFactor,reg_flag);
                cd(currentFolder);      
                fprintf(['No.' num2str(f) ' Channel 1 is finished\n']);

                %Do channel 2
                IntgA.A = IntgA.A2;
                IntgA.ROIData = IntgA.ROI2;
                IntgA.ROI2 = [];
                IntgA.A2 = [];
                currentFolder = pwd;
                newFolder = fullfile(currentFolder,'Channel_2');
                if ~exist(newFolder,'dir')
                    mkdir(newFolder);
                end
                cd(newFolder)
                Idx2 = IntgA.prePipe(spacialFactor,reg_flag);
                cd(currentFolder);   
                fprintf(['No.' num2str(f) ' Channel 2 is finished\n']);

                clear IntgA;

            else
                IntgA = Integration(f,flag,nmov);
                Idx = IntgA.prePipe(spacialFactor,reg_flag);              
                fprintf(['No.' num2str(f) ' is finished\n']);
                disp('');
                clear IntgA;
            end
        end

        function doAveragingAcrossMovies(flag,IdxInAll,IdxInAll_1,IdxInAll_2)
        
        %    Average across different movies with the AveragedMatrix_ mat
        %    files. This can handle different scenarios providing different flags.
        %    
        %    Inputs:
        %        flag   Represents different experimental setups  
        %        IdxInALL,IdxInAll_1,IdxInAll_2    Indices of sorted frames
            
        

            if flag == 2
                currentFolder = pwd;
                Folder1 = fullfile(currentFolder,'Channel_1');
                cd(Folder1)
                Integration.AveragedMapsAcrossMovies();
                savename = 'IdxInAll_1';
                save(savename,'IdxInAll_1');

                cd(currentFolder);   

                currentFolder = pwd;
                Folder2 = fullfile(currentFolder,'Channel_2');
                cd(Folder2)
                Integration.AveragedMapsAcrossMovies();
                savename = 'IdxInAll_2';
                save(savename,'IdxInAll_2');
                cd(currentFolder);
            else
                Integration.AveragedMapsAcrossMovies();
                savename = 'IdxInAll';
                save(savename,'IdxInAll')
            end
        end
        
        
        function GenerateCC(TH_A,BW_ppA,filename)
        
        %    Generate connected components from pre-processed matrix
        %    
        %    Inputs:
        %        TH_A    Matrix that has been filtered
        %        BW_A    Matrix that has been black-white thresholded
        %        filename   current filename
        
            
            %Connected components of thresholded A.
            disp('Processing Connected Components and Regionprops')
            disp('')
            sz = size(BW_ppA);
            CC = bwconncomp(BW_ppA);
            minSize = Integration.minSpotSize(CC,sz);
            CC = Integration.ccThresh(CC,minSize);
            STATS = regionprops(CC,TH_A,'Area','BoundingBox', 'Centroid',...
                'MaxIntensity', 'MinIntensity', 'MeanIntensity','MinorAxisLength');
            %Integration.makeCCMovie(filename,CC,sz);

            %Save connected components and regionprops STATS
            region.domainData.CC = CC;
            region.domainData.STATS = STATS;
            typeString = 'region';
            checkname = ['CC_' filename(16:length(filename)-4) '_' typeString '.mat'];
            save(checkname,'region');
        end
        
    end
end