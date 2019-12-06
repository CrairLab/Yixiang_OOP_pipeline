classdef spike2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spike2 stores the information of different channels acquaired from the
% Spike2 software. It also provides some functions to pre-process the data
% and extract indicies information from different channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All previous record is saved on Github
% Visit https://github.com/CrairLab/Yixiang_OOP_pipeline for more info
% Author: yixiang.wang@yale.edu
% Latest update:
% R6 12/06/19 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    properties
        Squarewaves; % baphy square waves raw data
        Cameraframes; % camera frames raw data
        Idx_baphstar; % baphy square waves indices
        Idx_framesEy; % camera frames indices
        Idx_framesEy_trial; % canera franes trail indices
        %Spike2filename;
        %filename;
    end
    
    
    methods
        
        function obj = spike2(Spike2filename,baphyObj)
        %spike2 object constructor
        %Inputs:
        %   Spike2filename    Spike2 filename
        %   baphyObj          corresponding baphy object
            
            %obj.filename = filename;
            %obj.Spike2filename = Spike2filename;
            
            blockN = baphyObj.blockN;
            BlockDura = baphyObj.BlockDura;
            
            try 
                disp(Spike2filename);
                load(Spike2filename);

                try
                    for i = 1:length(data)
                        cur_title = data{i}.title;
                        switch cur_title 
                            case 'baphstar'
                                obj.Squarewaves = data{i}.values;
                            case 'framesEy'
                                obj.Cameraframes = data{i}.values;
                        end
                    end
                        
                    Idx_baphstar = spike2.BaphystarInfo(obj.Squarewaves);
                    [Idx_framesEy, Idx_framesEy_trial] = spike2.FramesEyInfo(obj.Cameraframes,blockN,BlockDura);

                    %SaveName = [filename(1:length(filename)-4) '_Spike2Index.mat'];
                    %save(SaveName,'Idx_baphstar','Idx_framesEy','Idx_framesEy_trial');

                    obj.Idx_baphstar = Idx_baphstar;
                    obj.Idx_framesEy = Idx_framesEy;
                    obj.Idx_framesEy_trial = Idx_framesEy_trial;
                    
                catch
                    disp('')
                    disp('!!!Warning Spike2 Info Extraction Failed!!!')
                    disp('Information missing during spike2 recording!')
                    disp('')
                    disp('Bypassing the issue.Generating Dummy Index Trains...')
                    disp('')
                    [obj.Idx_baphstar,obj.Idx_framesEy] = spike2.DummyTrains(blockN,BlockDura);
                end

            catch ME
                if (strcmp(ME.message,'Argument must contain a character vector.'))
                    disp('')
                    disp('As no spike2 is provided, generating Dummy Index Trains...')
                    disp('')
                    [obj.Idx_baphstar,obj.Idx_framesEy] = spike2.DummyTrains(blockN,BlockDura);
                end
                disp('Spike2 files not provided');
            end
        end

    end
    
    
    
    methods (Static)
        
        function [Idx_baphstar,Idx_framesEy] = DummyTrains(blockN, BlockDura)            
        
        %    This function will generate dummy framesEy and baphstar
        %    trains if either the information is not successfully
        %    extracted from the spike2 files or the spike2 file is not
        %    provided. The dummy trains might work if:
        %    1.Baphy files are provided
        %    2.Assuming frames and baphy square waves were intrinsically aligned
        %    3.Each movie contains 6000 frames
        %    4.Each trail contains BlockDura frames
        %    
        %    Inputs:
        %        blockN           number of trails
        %        BlockDura        Duration of a single block         
        %
        %    Outputs:
        %        Idx_baphstar     Dummy baphy square waves train
        %        Idx_framesEy     Dummy frame indices
           
            
            Idx_baphstar = zeros(2,blockN);
            Idx_framesEy = zeros(1,blockN*BlockDura);
            for i = 1:blockN
                Idx_baphstar(1,i) = (i-1)*(BlockDura+2) + 1;
                Idx_baphstar(2,i) = Idx_baphstar(1,i) + BlockDura + 1;
                Idx_framesEy((i-1)*BlockDura+1:i*BlockDura) = [Idx_baphstar(1,i)+1 : Idx_baphstar(2,i)-1];
            end                         
        end
                
        function Idx_baphstar = BaphystarInfo(Baphystar_values)
        
        %    Extract startponts and endpoints of the square waves
        %    Evaluate the quality of the Spike2 files
        %
        %    Inputs:
        %        Baphystar_values    Channel of the spike2 file recording
        %                            Baphystar values
        %
        %    Outputs:
        %        Idx_baphstar   Indices of the startpoints and endpoints
           
            disp('Baphstar Indices Processing');

            Baphystar_num = length(Baphystar_values);
            Idx_baphstar_ini = [];
            Idx_baphstar_end = [];
            temp_values = Baphystar_values>1;

            for i = 1:Baphystar_num-1
                a1 = temp_values(i);
                a2 = temp_values(i+1);
                if a2>a1
                    Idx_baphstar_ini = [Idx_baphstar_ini i+1];
                end
                if a2<a1
                    Idx_baphstar_end = [Idx_baphstar_end i];
                end
            end

            if temp_values(length(temp_values)) == 1
                Idx_baphstar_end = [Idx_baphstar_end length(temp_values)];
            end

            Idx_baphstar = [Idx_baphstar_ini;Idx_baphstar_end];

            MeanLength = mean(Idx_baphstar(2,:)-Idx_baphstar(1,:));
            StdLength = std(Idx_baphstar(2,:)-Idx_baphstar(1,:));
            disp('Number of trials(last trial included) is:')
            disp(num2str(length(Idx_baphstar)));
            disp('Mean length of trials(last trial included) is:')
            disp(num2str(MeanLength/1000));
            disp('Std of length of trials(last trial included) is:')
            disp(num2str(StdLength/1000));

            disp(' ');

            Idx_baphstar_m = [Idx_baphstar(1,1:length(Idx_baphstar)-1);Idx_baphstar(2,1:length(Idx_baphstar)-1)];
            MeanLength = mean(Idx_baphstar_m(2,:)-Idx_baphstar_m(1,:));
            StdLength = std(Idx_baphstar_m(2,:)-Idx_baphstar_m(1,:));
            disp('Number of trials(last trial excluded) is:')
            disp(num2str(length(Idx_baphstar_m)));
            disp('Mean length of trials(last trial excluded) is:')
            disp(num2str(MeanLength/1000));
            disp('Std of length of trials(last trial excluded) is:')
            disp(num2str(StdLength/1000));
        end
        
        
        function [Idx_framesEy, Idx_framesEy_trial] = FramesEyInfo(FramesEy_values,blockN, BlockDura)
        
        %    Extract startponts and endpoints of the cameraframes
        %    Evaluate the quality of the Spike2 files
        %
        %    Inputs:
        %       FramesEy_values   Channel of the spike2 file recording framesEy
        %
        %    Outputs:
        %       Idx_framesEy   Indices of the startpoints and endpoints
        %       Idx_framesEy_trial   Indices of the camera frames
        
            disp(' ');
            disp('FramesEy Indices Processing');
            
            normalFrameN = blockN * BlockDura;
            FramesEy_num = length(FramesEy_values);
            Idx_framesEy_ini = [];
            Idx_framesEy_end = [];
            Idx_framesEy = [];
            temp_values = FramesEy_values>0.5;

            for i = 1:FramesEy_num-1
                a1 = temp_values(i);
                a2 = temp_values(i+1);
                if a2>a1
                    Idx_framesEy_ini = [Idx_framesEy_ini i+1];
                end
                if a2<a1
                    Idx_framesEy_end = [Idx_framesEy_end i];
                end
            end
            
            %Tackle inperfect frames number 
            try
                Idx_framesEy_ini(length(Idx_framesEy_end)+1:length(Idx_framesEy_ini)) = [];
                              
                disp(' ')
                disp(['Initially identified frames: ' num2str(length(Idx_framesEy_ini))])
                disp(['Default setting: ' num2str(normalFrameN) ' frames per movie!']);                
                disp(' ')
                
                if (length(Idx_framesEy_ini) > normalFrameN) && (length(Idx_framesEy_end) > normalFrameN)
                    disp(['Identified more than ' num2str(normalFrameN) ' frames in spike2, categorized as an error. Omitting extra frames...']);
                    disp('')
                    Idx_framesEy_ini = Idx_framesEy_ini(1:normalFrameN);
                    Idx_framesEy_end = Idx_framesEy_end(1:normalFrameN);
                elseif (length(Idx_framesEy_ini) < normalFrameN) && (length(Idx_framesEy_end) < normalFrameN)
                    disp(['Fewer than ' num2str(normalFrameN) ' frames in spike2 identified!'])
                    disp('')
                else
                    disp([num2str(normalFrameN) ' Frames detected in spike2...'])
                    disp('')
                end
            catch
                disp('Imperfect frames, check spike2 files')
            end
            
            if Idx_framesEy_ini == Idx_framesEy_end
                Idx_framesEy = Idx_framesEy_ini;
                Length_temp1 = Idx_framesEy(2:length(Idx_framesEy));
                Length_temp2 = Idx_framesEy(1:length(Idx_framesEy)-1);
                Length_temp = Length_temp1 - Length_temp2;
                %Usually interval between two consecutive block > 150
                Idx_framesEy_trial_ini = [Idx_framesEy_ini(1) Idx_framesEy_ini(Length_temp>150)+1];
                Idx_framesEy_trial_end = [Idx_framesEy_ini(Length_temp>150) Idx_framesEy_end(length(Idx_framesEy_end))];
                Idx_framesEy_trial = [Idx_framesEy_trial_ini;Idx_framesEy_trial_end];
                Length_trial = Idx_framesEy_trial(2,:) - Idx_framesEy_trial(1,:);

                MeanLength = mean(Length_trial);
                StdLength = std(Length_trial);
                disp('Number of trials(last trial included) is:')
                disp(num2str(length(Length_trial)));
                disp('Mean length of trials(last trial included) is:')
                disp(num2str(MeanLength/1000));
                disp('Std of length of trials(last trial included) is:')
                disp(num2str(StdLength/1000));

                disp(' ');

                Length_trial_m = Length_trial(1,1:length(Length_trial)-1);
                MeanLength = mean(Length_trial_m);
                StdLength = std(Length_trial_m);
                disp('Number of trials(last trial excluded) is:')
                disp(num2str(length(Length_trial_m)));
                disp('Mean length of trials(last trial excluded) is:')
                disp(num2str(MeanLength/1000));
                disp('Std of length of trials(last trial excluded) is:')
                disp(num2str(StdLength/1000));

            else 
                Idx_framesEy = [Idx_framesEy_ini;Idx_framesEy];
                warning('Camera frames_ini is incosistent with Camera frames_end. Pls check the data');
            end
            
        end
        
        function Spike2Matlab(spike2path,cedpath)
        
        %    Convert spike2 files to .mat files
        %    Warning: original spike2 files will be deleted!
        %    
        %    Inputs:
        %        cedpath       CED library path
        %        spike2path    spike2 files path
        %    
        %    Outputs:
        %        .mat files with 1000Hz binning rate
            
        
            if nargin<2
                cedpath='E:\CEDMATLAB\CEDS64ML'; % path to the CED Matlab toolbox
            end
            
            try
                CEDS64LoadLib(cedpath);
                if nargin >= 1
                    cd(spike2path)
                end

                spike2.fileDetector()
                spike2List = readtext('RawSpike2files.txt',' ');

                for i = 1:length(spike2List)
                    spike2.convertSpike2(spike2List{i});
                end

                delete('RawSpike2files.txt');
                CEDS64CloseAll();
            catch
                disp('Can not detect CEDS64 libray or can not find spike2 files!')
            end
         
            
        end
        
        
        
        function fileDetector()
        
        %    Generate txt files based on file detected. Generate RawSpike2files.txt 
        %    as a list of names of .smr/.smrx spike2 data. 
        %
        %    Warning: The filenames should be in descending order accordingly. 
        %    Otherwise the matching will be incorrect!!!
        

            tmp_list = dir;
            screen_types = {'.smrx','.smr'};
            disp(['Working Folder: ' pwd])
            spike2ID = fopen('RawSpike2files.txt','w+');

            for i = 1 : length(tmp_list)

                curName = tmp_list(i).name;

                if length(curName) <= 4
                    continue
                end

                switch curName(end-3:end)
                    case '.smr'
                        fprintf(spike2ID,'%s\r\n',curName);
                    case 'smrx'
                        fprintf(spike2ID,'%s\r\n',curName);
                end

            end
            fclose('all');

        end

        
        
        function convertSpike2(filename)
        
        %    Convert .smr or .smrx file to .mat files, provided binRate (default =
        %    1000Hz). Save the .mat file as -v7.3.
        %
        %    Inputs: 
        %        filename    current .smr or .smrx filename
        %    Outputs:
        %        Corresponding .mat file
        

            display(strcat('loading in smr: ',filename));
            fhand = CEDS64Open(filename);
            ichannum = min(CEDS64MaxChan(fhand),20); %Why 20? Oh it's just a random & resonably large number

            maxTimeTicks = CEDS64ChanMaxTime(fhand,1)+1;
            timeBase = CEDS64TimeBase(fhand); % Actual time for each time unit
            binRate = 1000; % Can be changed if higher resolution is required
            finalSize = ceil(maxTimeTicks*timeBase*binRate); %Final size of the downsampled data
            %data=nan(maxTimeTicks,ichannum);

            % get waveform data from each channel
            counter = 0;

            goodChans = [];
            for ichan=1:ichannum % Identify used channe
                %file name, channel num, max num of points, start and end time in ticks
                iType = CEDS64ChanType(fhand,ichan);
                if iType > 0
                    counter = counter+1;
                    goodChans(counter) = ichan;
                end
            end

            data = cell(counter,1);
            CamFramesTrain = zeros(finalSize,1);

            for ichan = 1:counter %only loop through the used channel

                [fRead,fVals,~] = CEDS64ReadWaveF(fhand,goodChans(ichan),maxTimeTicks,0,maxTimeTicks);

                if isempty(fVals)
                    [ iRead, vi64T ] = CEDS64ReadEvents(fhand,goodChans(ichan),1e6,0);
                    vi64T = ceil(vi64T.*(timeBase*binRate));
                    CamFramesTrain(vi64T) = 1;
                end

                [~,iTitle] =  CEDS64ChanTitle(fhand,goodChans(ichan)); % Get channel title

                if ~isempty(fVals)
                    skipIntvl = round(fRead/finalSize);
                    dataStruct = struct('title',iTitle,'values',fVals(1:skipIntvl:fRead));
                else
                    dataStruct = struct('title',iTitle,'values',CamFramesTrain);
                end

                data{ichan} = dataStruct;

            end

            display(strcat('Finished: ',filename));

            tmpString = [];
            i = 1;
            while true && i<1000   
                if filename(i) == '.'
                    break
                end
                tmpString = [tmpString filename(i)];
                i = i+1;
            end

            save(strcat(tmpString,'.mat'),'data','-v7.3');
            CEDS64Close(fhand);
            delete(filename);


        end
        
    end
    
end