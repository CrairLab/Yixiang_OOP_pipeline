classdef Integration < spike2 & baphy & movieData & Names & ROI & wlSwitching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration is a subclass that inherits from all superclasses to
% analyze auditory experiment data. Specifically, Integration inherits
% from spike2\baphy\movidData\Names\ROI\wlSwitching classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All previous record is saved on Github
% Visit https://github.com/CrairLab/Yixiang_OOP_pipeline for more info
% Author: yixiang.wang@yale.edu
% Latest update:
% R31 02/10/20 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    properties
        flag; %0, for spontaneous activity without stimuli; 1, for activity with stimuli; 2, for wavelength switching 
        nmov; %Number of movies
        smallMask; %Mask that conform to downsampled ROI-applied movie
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
                movieObj = movieData(filename,1); A_ = movieObj.A; %Here if input flag == 2, will not further process
                ROIObj = ROI(); ROI_ = ROIObj.ROIData; %Here if input flag == 2, will not further process
            else
                A_ = []; ROI_ = [];
            end           
            
            %To get blcokN for spike2 object construction
            baphyObj = baphy(Baphyname);
            
            %Subclass objects construction
            obj@Names(f,flag);
            obj@spike2(Spike2name,baphyObj);
            obj@baphy(Baphyname);
            obj@movieData(filename,flag);
            obj@ROI();
            obj@wlSwitching(flag,A_,ROI_);
          
            obj.flag = flag;
            obj.nmov = nmov;
            
            clear movieObj ROIObj Namesobj baphyObj
          
        end
        
        function FramesByFreqVolu = prePipe(obj,param)
        
        %    Preprocessing pipeline to handle a single auditory project
        %    movie. Process includes: photobleaching correction,
        %    rigid registration, applying ROI masks, downsampling,
        %    top-hat filtering, frequency-frame alignment(optional),
        %    frequency-volume maps plotting(optional),gaussian smoothing, 
        %    black-white thresholding, connected component extraction etc.            
        %    
        %    Inputs:
        %        obj               current object
        %        param             parameter struct to feed in different
        %                          function
        %   
        %   Outputs:
        %        FramesByFreqVolu  sorted frames by frequency and volume
        
            currentFolder = pwd;
            filename = obj.filename;        
            disp(['Processing: ' filename]);
            disp(' ');
 
            outputFolder = fullfile(currentFolder,obj.outputFolder); 
            mkdir(outputFolder);                    
            
            if obj.flag 
            %Try to do frame-frequency alignment if flag == 1;
            %For spontaneous cases this step will be skipped.
            %If smallMask property existed, instance.mat file must had
            %already exisited, so does frequency-volume analysis
                try
                    FramesByFreqVolu = obj.FramesByFreqAlign;
                catch ME
                    %Handle damaged spike2 files problem
                    if (strcmp(ME.identifier,'MATLAB:badsubscript'))
                        disp('')
                        disp('Index exceeds matrix dimensions. Check Spike2!!')
                        disp('Skip alignment...')
                        disp('')
                    elseif (strcmp(ME.identifier,'MATLAB:unassignedOutputs'))
                        disp('')
                        disp('Unassigned Outputs, probably baphy file not provided. Skip alignment...')
                        disp('')
                    end
                    %obj.flag = ~obj.flag;
                    FramesByFreqVolu = 0;
                end
            else
                FramesByFreqVolu = 0;
            end
            
                            
            %Detect if frame-discarding registration happened
            checkname = [filename(1:length(filename)-4) '_ori_DS_registered.mat'];
            checkpath = fullfile(outputFolder,checkname);
            if exist(checkpath,'file')
                movTag = 'dsc';
            else
                movTag = '';
            end

            checkname = [filename(1:length(filename)-4) '_filtered' movTag '.mat'];
            checkpath = fullfile(outputFolder,checkname);
            %If detected _filtered .mat file, skip the whole pre-processing
            if exist(checkpath,'file')
                disp('Detect filtered matrix, skip the whole pre proessing!')
                load(checkpath);
            else


                %If instance matrix exsits, load the object/instance
                if ~isempty(obj.smallMask)
                    %If smallMask property existed, instance.mat file must had
                    %already exisited
                    A_DS = obj.A;
                    sz = size(obj.smallMask);

                    %Recover 3D matrix A_DS
                    if length(size(A_DS)) == 2
                        A_DS = reshape(A_DS, [sz(1) sz(2) size(A_DS,2)]);
                    end                
                else                    
                    %If instance matrix not detected, create one

                    %Correct photobleaching
                    A_corrct = Integration.bleachCorrection(obj.A);
                    obj.A = [];
                    disp('Photobleaching corrected!');


                    %Assess movement, doing translation registration at the same time
                    %Run the movement assessment
                    tic;
                    if param.moveAssessFlag
                        %If true, call the movAssess function
                        [A_registered, A_ori, tform_all, NormTform_all, movIdx_saved] = ...
                        movieData.movAssess(A_corrct);

                        if size(A_ori,3) ~= size(A_registered,3)
                        %Save downsampled and registered original movie
                        %Only save the original movie if there are frames being
                        %identified as moving and discarded
                            A_ori_DS = Integration.downSampleMovie(A_ori,param.spacialFactor);
                            A_ori_DS = reshape(A_ori_DS, [size(A_ori_DS,1)*size(A_ori_DS,2),...
                                size(A_ori_DS,3)]);
                            checkname = [filename(1:length(filename)-4) '_ori_DS_registered.mat'];
                            save(fullfile(outputFolder,checkname),'A_ori_DS','-v7.3');
                            clear A_ori
                            clear A_ori_DS
                            movTag = 'dsc';
                        else
                            movTag = '';
                        end                    
                        disp('Movement assessment finished...Time cost = ')
                    else
                        [A_registered, tform_all] = movieData.dftReg(A_corrct, 'forAll');
                        NormTform_all = sqrt(tform_all(:,3).^2 + tform_all(:,4).^2);
                        movIdx_saved = ones(size(A_corrct,3),1);
                        disp(['Mean tform magnitude (minus I) = ' num2str(mean(NormTform_all))]);
                        disp('Not calling movAsess function. Only do rigid registration. Time cost = ')
                        movTag = '';
                    end


                    toc;

                    %Save movemet assessment results
                    checkname = [filename(1:length(filename)-4) '_moveAssess' movTag '.mat'];
                    save(fullfile(outputFolder,checkname),'tform_all','NormTform_all','movIdx_saved');                


                    %Gaussian smoothing
                    A_registered = Integration.GauSmoo(A_registered,1); %set sigma = 1
                    disp('Gaussian smoothing is done');
                    disp(' ')

                    %Apply ROI mask(s)
                    A_ROI = Integration.ApplyMask(A_registered,obj.ROIData);
                    disp('Successfully apply ROI')
                    clear A_corrct
                    
                    %Get averaged dF/F outside of the ROI
                    A_out = A_registered .* (A_ROI == 0);
                    A_out(A_out == 0) = nan;
                    Avg_out = nanmean(A_out,1);
                    Avg_out = nanmean(Avg_out,2);
                    Avg_out_dFoF = Integration.grossDFoverF(Avg_out);
                    if any(isnan(Avg_out_dFoF))
                        Avg_out_dFoF = [];
                        disp('Background not defined! Check if .roi file is provided!')
                    end
                    checkname = [filename(1:length(filename)-4) '_out_dFoF.mat'];
                    save(fullfile(outputFolder,checkname),'Avg_out_dFoF');
                    clear A_registered Avg_out_dFoF

                    %Focusing on just the ROI part of the movie
                    A_ROI = movieData.focusOnroi(A_ROI);

                    %Downsampling
                    A_DS = Integration.downSampleMovie(A_ROI,param.spacialFactor);
                    disp(['Successfully downsample by factor of ' num2str(param.spacialFactor)])
                    clear A_ROI

                    %Get the downsampled roi mask
                    sz = size(A_DS);
                    ds_Mask = repmat((A_DS(:,:,1) ~= 0),[1,1,sz(3)]);
                    obj.smallMask = ds_Mask(:,:,1);              

                    %"Raw" data stored (reshape to 2D to save space)
                    obj.A = reshape(A_DS, [sz(1)*sz(2), sz(3)]);

                    %Save the instance as an object
                    checkname = [filename(1:length(filename)-4) '_instance.mat'];
                    save(fullfile(outputFolder,checkname),'obj','-v7.3');
                    %delete(filename);
                    disp('Instance.mat file saved, consider delete raw .tif movies!')
                end

                %Top-hat filtering
                if ~obj.flag
                    %TH_A = Integration.TopHatFiltering(A_DS);
                    TH_A = A_DS;
                    %disp('Top-hat filtering is done');
                    disp('Not doing top-hat filtering here')
                else
                    sz = size(A_DS); se = strel('rectangle',[sz(1)*2,sz(2)]*2);
                    TH_A = imtophat(A_DS, se);
                    disp('For stimulation experiments, not doing top-hat filtering in time dimmension')
                    disp('')
                end
                clear A_DS

                %ICA analysis of the matrix
                if isfield(param,'ICAflag')
                    if param.ICAflag && ~param.flag
                    %Only do ICA analysis for spontaneous case
                        try
                            [icasig, M, W, corr_map] = movieData.getICA(TH_A);
                            checkname = [filename(1:length(filename)-4) '_ICA.mat'];
                            save(fullfile(outputFolder,checkname),'icasig', 'M', 'W', 'corr_map')
                        catch
                            disp('ICA analysis failed')
                        end
                    else
                        disp('Skip ICA analysis')
                    end
                end

                %Centered data around origins (gross dFoF)
                A_mean = nanmean(TH_A,3);
                TH_A = TH_A./A_mean - 1;

                %SVD denosing of down-sampled A
                try 
                    iniDim = param.iniDim;
                catch
                    iniDim = 1; param.iniDim = iniDim;
                end
                %iniDim = iniDim + iniDimFlag;
                [de_A,U,S,V,iniDim,PC_exp] = Integration.roiSVD(TH_A, iniDim);
                %Reaply downsampled roi mask
                if ~exist('ds_Mask','var')
                    ds_Mask = obj.smallMask;
                    ds_Mask = repmat(ds_Mask, [1,1,size(de_A,3)]);
                end
                de_A = de_A.*ds_Mask;
                de_A(de_A == 0) = nan;
                checkname = [filename(1:length(filename)-4) '_SVD.mat'];
                save(fullfile(outputFolder,checkname),'U','S','V','PC_exp'...
                    ,'param','iniDim');
                disp('SVD denosing is done')
                disp('')
                clear TH_A

                %Recover the data, this is important for later dFoF
                de_A = de_A.*A_mean + A_mean;

               %Impose dFOverF to downsampled matrix
                if ~obj.flag
                    A_dFoF = Integration.grossDFoverF(de_A);
                    A_dFoF = A_dFoF.*ds_Mask;
                    clear A_registered
                    disp('Gloabal dFoverF is done')
                    disp(' ')
                else 
                    A_dFoF = de_A;
                    disp('Do not do gross dFoF for stimulation experiment!')
                    disp('')
                end
                clear de_A;

                %Z-scoring de_A along the time dimension
                %de_A = zscore(de_A,1,3);
                %disp('Z-scored reconstructed matrix')

                %Save filtered matrix
                checkname = [filename(1:length(filename)-4) '_filtered' movTag '.mat'];
                save(fullfile(outputFolder,checkname),'A_dFoF','-v7.3');
            end
                    
            exptparam.PreStimSilence = obj.PreStimSilence;
            exptparam.PostStimSilence = obj.PostStimSilence;
            exptparam.Duration = obj.Duration;
            exptparam.BlockDura = obj.BlockDura;
                      
            %Check flag to decide whether to generate frequency/volume maps
            if obj.flag
                %Integration.FreqColorMap(TH_A,filename,obj);
                A_dFoF(A_dFoF == 0) = nan; % for better intra-Video dFOverF
                Integration.FreqVoluColorMap(A_dFoF,filename,obj.nmov,exptparam);
                %de_A(isnan(de_A)) = 0; %resume for later filtering
            end

            %if obj.flag == 0 %compute connected components only for spontaneous case
                
            %Generate connected components using renewCC function
            Integration.renewCC(A_dFoF,param.CCthreshold, outputFolder,filename)
            %end
            clear A_dFoF ppA_roi
            
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
            savename = [filename(1:length(filename)-4),'_Spike2BaphyIndex.mat'];
            
            %Reconstruct baphy and spike2 properties if not defined
            if isempty(Idx_baphstar) || isempty(Idx_framesEy)
                
                    %Find corresponding filenames
                    outputFolder = obj.outputFolder;            
                    f = str2double(outputFolder(7:end));
                    Spike2filelist = readtext('Spike2files.txt','\n');
                    Spike2name = Spike2filelist{f};
                    baphyfilelist = readtext('baphyfiles.txt','\n');
                    Baphyname = baphyfilelist{f};
                
                    %To get blcokN for spike2 object construction
                    baphyObj = baphy(Baphyname);          
                    %Construct spike2 class
                    spike2Obj = spike2(Spike2name,baphyObj);
                    obj.Idx_baphstar = spike2Obj.Idx_baphstar;
                    obj.Idx_framesEy = spike2Obj.Idx_framesEy;
                    obj.BlockDura = baphyObj.BlockDura;
                    obj.evt_line = baphyObj.evt_line;
                    Idx_baphstar = obj.Idx_baphstar;
                    Idx_framesEy = obj.Idx_framesEy;
            end
          
            %Handle imperfect recording, preserve good trails
            BlockDura = obj.BlockDura; %Duration of a single block
            
            if length(size(obj.A)) == 2
                frameN = size(obj.A,2);
            else
                frameN = size(obj.A,3);
            end
            
            ActualTrailNumber = min([size(Idx_baphstar,2),size(obj.evt_line,1),floor(frameN/BlockDura)]);
            disp(['[Baphstar,evt_line,A_size3]= ' num2str([size(Idx_baphstar,2),size(obj.evt_line,1),floor(frameN/BlockDura)])])
            disp(['Usable Trails = ' num2str(ActualTrailNumber)]);
            evt_line = obj.evt_line(1:ActualTrailNumber,:);
            [freq,volu] = Integration.GetVoluFreq(evt_line);
            obj.freq = freq;
            obj.volu = volu;
            
            %Remove extra indices from spike2 if it goes beyond the actual
            %movie frames
            Idx_baphstar = Idx_baphstar(:,1:ActualTrailNumber); 
            Idx_framesEy = Idx_framesEy(:,1:BlockDura*ActualTrailNumber);
            
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
            FreqFilename = [filename(1:length(filename)-4),'_FramesByFreq.mat'];
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
            FreqVoluFilename = [filename(1:length(filename)-4),'_FramesByFreqVolu.mat'];
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
        
        function CC = ccThresh_3D(CC,minSpotSize)
        
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
        
        
        function CC = ccTimeThresh(CC, durT)
        
        %    Thresholding the CC based on each component's temporal
        %    duration, which is the length of the BoundingBox side along the
        %    time dimension. A component should last for at least 3 frames
        %   
        %    Inputs:
        %    CC     Connected components from bwconncomp
        %    durT   duration threshold
        %
        %    Outputs:
        %    CC     Renewed CC
        
            if nargin < 2
                durT = 3;
            end
            %disp(['Duration threshold = ' num2str(durT)])
        
            STATS = regionprops(CC,'Area','BoundingBox');
            roiBoundingBox = zeros(length(STATS),6);
            for i = 1:length(STATS)
                roiBoundingBox(i,:) = STATS(i).BoundingBox;
            end
            durations = roiBoundingBox(:,6);
            newPixelIdxList = CC.PixelIdxList;
            
            %duration should be larger/equal to 3 frames
            newPixelIdxList(durations(:) < durT) = [];
            CC.PixelIdxList = newPixelIdxList;
            CC.NumObjects = length(newPixelIdxList);                 
        end
        
        
        function fileDetector()
        
        %    Generate txt files based on file detected. Generate files.txt as a list
        %    of names of .tif movies or existed instance .mat files. 
        %    Generate Spike2files.txt as a list of names of .mat spike2 data. 
        %    Generat baphyfiles.txt as a list of names of baphy files. 
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
                if length(curName) <= 2
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
                                if ~contains(curName,'@00')
                                    fprintf(tifID,'%s\r\n',curName);
                                end
                        end
                    end
                end                               
            end
            
            %Detect existed instace .mat files and write to files.txt
            filelist = readtext('files.txt',' ');
            nmov = size(filelist,1);
            if nmov == 0
                disp('Did not detect .tif files, try to detect instance.mat files')
                for i = 1 : length(tmp_list)
                    curName = tmp_list(i).name;
                    %Detect subfolders containing keyword 'output'
                    if contains(curName, 'output')
                        %Search the sub folders
                        subFileList = dir(fullfile(cd, curName));
                        for j = 1:length(subFileList)
                            subFileName = subFileList(j).name;
                            %Detect instance .mat files containing keyword 'instance'
                            if contains(subFileName, 'instance.mat')
                                %Save the entire path for future loading
                                writeName = [subFileName(1:end-13), '.tif'];
                                fprintf(tifID,'%s\r\n',writeName);
                            end
                        end
                    elseif contains(curName, 'instance.mat')
                        %Search the root folder
                        writeName = [curName(1:end-13), '.tif'];
                        fprintf(tifID,'%s\r\n',writeName);                        
                    end
                end
            end
            
            fclose('all');
            
        end
        
        function filelist = fileDetector_keyword(keyword)
        % Detect files whose names contain the provided keyword
        % Inputs:
        %   keyword     string for regexp function
        %
        % Outputs:
        %   filelist    List of qualified files
        %
            if nargin == 0
                disp('No key string provided!')
                keyword = '.';
            end

            temp_info = dir;
            filelist = {};
            n = 0;
            for i = 1:size(temp_info,1)  
                if regexp(temp_info(i,1).name,keyword)
                    n = n+1;
                    filelist{n,1} = temp_info(i,1).name;
                end
            end
        end
        
        function minSize= minSpotSize(CC)
        
        %    Compute the minimum size of CC to be preserved. 
        %
        %    Inputs:
        %        CC      connected components
        %                %
        %    Outputs:
        %        minSize     minimum size chosen for thresholding
        
            
            PixelIdxList = CC.PixelIdxList;
            sizeCounts = cellfun(@length,PixelIdxList);
            %Choose the bigger one btw 20% percentile of all CC sizes and 0.1% of image size 
            sz = CC.ImageSize;
            minSize = max(prctile(sizeCounts,20),sz(1)*sz(2)/1000);
        
        end
        
        function [Idx,Idx1,Idx2] = processMovies(f,nmov,param)
        
        %    Proccess movies (called by audPipe) seuqentially. This can
        %    handle different scenarios providing different flags.
        %    
        %    Inputs:
        %        f      The f-th movie to process
        %        param  parameters that users feed in
        %        
        %    Outputs:
        %        Idx,Idx1,Idx2    Frames indices sorted by frequencies or volume
        
            Idx = []; Idx1 = []; Idx2 = [];
            flag = param.flag;
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
                Idx1 = IntgA.prePipe(param);
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
                Idx2 = IntgA.prePipe(param);
                cd(currentFolder);   
                fprintf(['No.' num2str(f) ' Channel 2 is finished\n']);

                clear IntgA;

            else
                %Single channel recording
                IntgA = Integration(f,flag,nmov);
                Idx = IntgA.prePipe(param);              
                fprintf(['No.' num2str(f) ' is finished\n']);
                disp('');
                clear IntgA;
            end
        end

        function doAveragingAcrossMovies(flag,IdxInAll,IdxInAll_1,IdxInAll_2,exptparam)
        
        %    Average across different movies with the AveragedMatrix_ mat
        %    files. This can handle different scenarios providing different flags.
        %    
        %    Inputs:
        %        flag   Represents different experimental setups  
        %        IdxInALL,IdxInAll_1,IdxInAll_2    Indices of sorted frames
        %        exptparam   experiment parameters
            
        if isempty(exptparam)
            exptparam.PreStimSilence = 1;
            exptparam.PostStimSilence = 3.5;
            exptparam.Duration = 0.5;
            exptparam.BlockDura = 60;
        end

            if flag == 2
                %Channel 1
                currentFolder = pwd;
                Folder1 = fullfile(currentFolder,'Channel_1');
                cd(Folder1)
                Integration.AveragedMapsAcrossMovies(exptparam);
                savename = 'IdxInAll_1';
                save(savename,'IdxInAll_1');

                cd(currentFolder);   
                
                %Channel 2
                currentFolder = pwd;
                Folder2 = fullfile(currentFolder,'Channel_2');
                cd(Folder2)
                Integration.AveragedMapsAcrossMovies(exptparam);
                savename = 'IdxInAll_2';
                save(savename,'IdxInAll_2');
                cd(currentFolder);
            else
                %Single channel recording
                Integration.AveragedMapsAcrossMovies(exptparam);
                savename = 'IdxInAll';
                save(savename,'IdxInAll')
            end
        end
        
        
        function BW_ppA = filterCC_byFrame(BW_ppA)
            
        % Use 2D propertity of each frame (size, eccentricity, orientation)
        % to filter connected components and reconstruct the binary movie
        %
        % Input/Output:
        %   BW_ppA          Binary 3D matrix
        
            sz = size(BW_ppA);
                for i = 1:sz(3)

                    cur_img = BW_ppA(:,:,i);
                    cur_CC = bwconncomp(cur_img);

                    %Filter by size
                    szList = cellfun(@length,cur_CC.PixelIdxList);
                    kill_sz= szList< sz(1)*sz(2)/1000;
                    cur_CC.PixelIdxList(kill_sz) = [];
                    cur_CC.NumObjects = length(cur_CC.PixelIdxList);

                    %Filter by ellipse properties
                    leftBound = round(size(cur_img,2)/2);
                    rightBound = round(size(cur_img,2)/2);

                    cur_STATS = regionprops(cur_CC,'Eccentricity','Orientation','Centroid'); 
                    kill_eclipse = zeros(1,length(cur_STATS));
                    for j = 1:length(cur_STATS)
                        cur_struct = cur_STATS(j);

                        if  cur_struct.Eccentricity < 0.6 %filter by eccentricity
                            kill_eclipse(1,j) = 1;
                        else
                            if (cur_struct.Centroid(1) < leftBound) %for bands on the left 1/3, specify reasonable orientation range
                                if (cur_struct.Orientation < -75)||(cur_struct.Orientation > -15)
                                    kill_eclipse(1,j) = 1;
                                end
                            elseif (cur_struct.Centroid(1) > rightBound) %for bands on the right 1/3, specify reasonable orientation range
                                if (cur_struct.Orientation > 75)||(cur_struct.Orientation < 15)
                                    kill_eclipse(1,j) = 1;
                                end
                            end
                        end        
                    end

                    cur_CC.PixelIdxList(logical(kill_eclipse)) = [];
                    cur_CC.NumObjects = length(cur_CC.PixelIdxList); 

                    %Reconstruct filtered image
                    x = [];
                    for k = 1:length(cur_CC.PixelIdxList)
                        cur_cell = cur_CC.PixelIdxList(k); 
                        x = [x;cur_cell{1}];
                    end
                    rcs_img = zeros(sz(1),sz(2));
                    rcs_img(x) = 1;
                    BW_ppA(:,:,i) = rcs_img;
                end

        end
        
        
        function [region,BW_ppA] = GenerateCC(ppA,BW_ppA,frameFlag,timeFlag,durT)
        
        %    Generate connected components from pre-processed matrix
        %    
        %    Inputs:
        %        ppA     Matrix that has been pre-processed
        %        BW_A    Matrix that has been black-white thresholded
        %        filename   current filename
        %        frameFlag  whether filtering by frame properties
        %        timeFlag   whether filtering by duration
        %        durT       duration threshold
        %
        %    Outputs:
        %        region  store all region properties
        %        BW_ppA  filtered blakc-white matrix
            
            if nargin<3
                frameFlag = 0;
                timeFlag = 0;
                durT = 3;
            end
            
            %Connected components of thresholded A.
            disp('Processing Connected Components and Regionprops')
            disp('')
            
            %Filter BW_ppA by frame using 2D properties of the CCs
            if frameFlag
                BW_ppA = Integration.filterCC_byFrame(BW_ppA);
            end
            
            %Filter BW_ppA using the 3D properties of the CCs
            CC = bwconncomp(BW_ppA);
            %minSize = Integration.minSpotSize(CC);
            minSpotSize = 9;
            CC = Integration.ccThresh_3D(CC,minSpotSize);
            if timeFlag
                CC = Integration.ccTimeThresh(CC, durT);
            end
            
            STATS = regionprops(CC, ppA,'Area','FilledArea', 'BoundingBox',...
                'Centroid',...
                'MaxIntensity', 'MinIntensity', 'MeanIntensity');
            %Integration.makeCCMovie(filename,CC,sz);
            
            %Reconstruct filtered image
            x = [];
            for k = 1:length(CC.PixelIdxList)
                cur_cell = CC.PixelIdxList(k); 
                x = [x;cur_cell{1}];
            end
            rcs_movie = zeros(size(BW_ppA));
            rcs_movie(x) = 1;
            BW_ppA = logical(rcs_movie);

            %Store connected components and regionprops STATS
            region.CC = CC;
            region.STATS = STATS;
        end
        
        
        function minuteFreqCC(BW_ppA,ppA_roi,outputFolder,filename)            
        %   Compute connected components within each minutes of the movie
        %
        %   Inputs:
        %       BW_ppA     pre-processed matrix black and white movie
        %       outputFolder     folder that users define to save the oupputs
        %       filename    filename from which current matrix is read in
        %
        %   Outputs:
        %       .mat files that store the information average/standard
        %       deviation of number of bands within 1 minutes of the movie
        
            iterN = ceil(size(BW_ppA,3)/600); %Number of iteration
            CC_freq = [];       
            for i = 1 : iterN
                frame_ini = (i-1)*600 + 1;
                frame_end = i*600;           
                %Sample each minutes of the movie
                try
                    sub_BW_ppA = BW_ppA(:,:,frame_ini:frame_end);
                    sub_ppA = ppA_roi(:,:,frame_ini:frame_end);
                    frameN = 600;
                catch
                    sub_BW_ppA = BW_ppA(:,:,frame_ini:end);
                    sub_ppA = ppA_roi(:,:,frame_ini:end);
                    frameN = size(BW_ppA,3) - frame_ini + 1;
                end
                %Calculate CCs within each minutes
                [cur_region,~] = Integration.GenerateCC(sub_ppA,sub_BW_ppA);
                cur_CCN = cur_region.CC.NumObjects;
                cur_CC_freq = cur_CCN *(600./frameN); %In case there aren't 600 frames in this sub-movie
                CC_freq(i,1) = cur_CC_freq;
            end
            avg_minute_freq = mean(CC_freq);
            std_minute_freq = std(CC_freq);
            minute_stat = [avg_minute_freq,std_minute_freq];
            checkname = ['CC_' filename(1:length(filename)-4) '_1minFreq.mat'];
            save(fullfile(outputFolder,checkname),'CC_freq','minute_stat');
        end
        

        function renewCC(ppA_roi,thresh,outputFolder,filename)
        
        %   Renew connected components in case the standards to define what is
        %   a qualified CC is changed.
        %   
        %   Inputs:
        %   ppA_roi     pre-processed matrix after being chopped to roi
        %   outputFolder     folder that users define to save the oupputs
        %   filename    filename from which current matrix is read in
        %
        %   Outputs:
        %   .mat files that store regionprops, CCs, and BW_ppA matrix after
        %   filtering.
                

        if nargin < 3
            outputFolder = cd;
            filename = 'current.mat';
            if nargin < 2
                thresh = 0.05;
            end
        end

                %Black-white thresholding of pre-processed A
                BW_ppA = imbinarize(ppA_roi,thresh);
                
                %Compute connected components within each minutes of the movie
                %Integration.minuteFreqCC(BW_ppA,ppA_roi,outputFolder,filename) 

                %Generate connected component
                [region,BW_ppA] = Integration.GenerateCC(ppA_roi,BW_ppA);
                checkname = ['CC_' filename(1:length(filename)-4) '_region.mat'];
                save(fullfile(outputFolder,checkname),'region');

                %Save binary movie
                checkname = ['Binary_' filename(1:length(filename)-4) '.mat'];
                save(fullfile(outputFolder,checkname),'BW_ppA','-v7.3');
                clear BW_ppA
        end
        
        
        
        function [curLoad,outputFolder,filename] = readInSingleMatrix(tag,specificF)
        %Read in variables from a .mat file and save them as a struct named
        %curLoad. Read in *_filter.mat files in default.
        %
        %Inputs:
        %   tag          tag specified what kind of .mat files to read
        %   specificF    number specify which .mat file to read
        %
        %Outputs:
        %   curLoad      a struct to store variables read in from .mat
        %   outputFolder     detected outputFolder name
        %   filename     detected filename
        
            disp('Make sure the tag variable can identify .mat files to be read in')
            if ~exist('tag','var')
                tag = 'filtered';
            end

            if ~exist('specificF','var')
                %Did not specify which file to read
               warning('Please specify which file to read!')
            end

            f = specificF;
            cur_Names = Names(f,0);
            filename = cur_Names.filename;
            currentFolder = pwd;
            outputFolder = fullfile(currentFolder,cur_Names.outputFolder);
            checkname = [filename(1:length(filename)-4) '_' tag '.mat'];               
            if exist(fullfile(outputFolder,checkname),'file')
                %Check whether pre-processing has been done before in
                %subfolders
                disp([tag ' matrix detected, loading .mat file...'])
                curLoad = load(fullfile(outputFolder,checkname));
                disp('')
            else
                disp('')
                disp([tag ' No matrix detected!!!'])
                disp('Try loading matrix from pwd')
                outputFolder = fullfile(currentFolder);
                if exist(fullfile(outputFolder,checkname),'file')
                    %Check whether pre-processing has been done before in
                    %the root folder
                    disp([tag ' matrix detected, loading .mat file...'])
                    curLoad = load(fullfile(outputFolder,checkname));
                    disp('')
                else
                    disp([tag ' No matrix detected!!!'])
                    curLoad = [];
                    warning('Inquired matrix does not exist...')
                end
                disp('')
            end
        end
        
        
        
        function Avg_out_dFoF = GenerateOutsideAverageFromScratch(f)
        %Detect or generate background trace (averaged trace of pixels
        %outside of given roi) from scratch
        %Inputs:
        %     f    handles(number) of .tif file
        %Outputs:
        %     Avg_out_dFoF    averaged trace of pixels outside of given roi
        
            try
                [out_dFoF,outputFolder,filename]  = Integration.readInSingleMatrix('out_dFoF', f);
                disp(['Current folder: ' outputFolder])
                checkname = [filename(1:length(filename)-4) '_out_dFoF.mat'];
                if isempty(out_dFoF)
                    %Read the movie 
                    curA = movieData.inputMovie(filename);
                    %Photobleaching
                    A_corrct = Integration.bleachCorrection(curA);
                    %Gaussian smoothing
                    A_smoothed = Integration.GauSmoo(A_corrct,1);
                    sz = size(curA);
                    %Get ROI
                    ROIobj = ROI(); ROIData = ROIobj.ROIData;
                    if ~isempty(ROIData)
                        x = ROIData.mnCoordinates(:,1);
                        y = ROIData.mnCoordinates(:,2);
                        %Generate mask for pixels outside the ROI
                        BW = poly2mask(x,y,sz(1),sz(2)); BW_out = ~BW;
                        BW_out = repmat(BW_out, [1,1,sz(3)]);
                        A_out = A_smoothed .* BW_out;
                        A_out(A_out == 0) = nan;
                        %Get the averaged trace outside the ROI
                        Avg_out = nanmean(A_out,1);
                        Avg_out = nanmean(Avg_out,2);
                        Avg_out_dFoF = Integration.grossDFoverF(Avg_out);
                    else 
                        warning('Can not find .roi file!')
                        Avg_out_dFoF = [];
                    end
                    %Save the result
                    save(fullfile(outputFolder,checkname),'Avg_out_dFoF');
                else
                    Avg_out_dFoF = out_dFoF.Avg_out_dFoF;
                    disp('Detected and loaded background trace!')
                end
                
                [moveAssessLoad,~,~]  = Integration.readInSingleMatrix('moveAssess', f);
                if isempty(moveAssessLoad)
                    [moveAssessLoad,~,~]  = Integration.readInSingleMatrix('moveAssessdsc', f);
                    if isempty(moveAssessLoad)
                        disp('Unable to detect movAssess file!')
                    else
                        [Avg_out_dFoF,~,~] = movieData.discardFrames(Avg_out_dFoF, moveAssessLoad.NormTform_all);
                        save(fullfile(outputFolder,checkname),'Avg_out_dFoF');
                    end
                else
                        [Avg_out_dFoF,~,~] = movieData.discardFrames(Avg_out_dFoF, moveAssessLoad.NormTform_all);
                        save(fullfile(outputFolder,checkname),'Avg_out_dFoF');
                end
                
            catch
                warning('Something wrong read in/generate the background trace!')
                Avg_out_dFoF = [];
            end
        end
    end
end