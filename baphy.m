classdef baphy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This class stores inforamtion extracted from the baphy files, which 
% includes arrays of frequencies/volumes, and other experiment parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All previous record is saved on Github
% Visit https://github.com/CrairLab/Yixiang_OOP_pipeline for more info
% Author: yixiang.wang@yale.edu
% Latest update:
% R4 09/11/19 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    properties
        evt_line;%text extracted from .m files
        freq; %frequency
        volu; %volume 
        blockN; % Number of blocks(trails)
        PreStimSilence; %Silent period before stimulation onset
        PostStimSilence %Silent period after stimulation ends
        Duration % Duration of stimulation
        BlockDura % Duration of a single block (in frames: default 10Hz)
    end
    
    methods
        
        function obj = baphy(baphyfn)
        %Baphy object constructor
        
            try
                %Read in baphy (.m) files
                fid = fopen(baphyfn);
                tline = fgetl(fid);
                evt_line = {};
                
                %Detect key words 'PreStimSilence' and 'exptevents'
                while ischar(tline) %&& (length(evt_line) <= 160)
                    if logical(contains(tline,'exptevents')) && logical(contains(tline,'PreStimSilence'))
                        evt_line = [evt_line;tline];                        
                    elseif logical(contains(tline,'exptparams(1).TrialObject(1)')) && logical(contains(tline,'PreStimSilence')) && ~logical(contains(tline,'UserDefinableFields'))
                    %Get duration, pre/post-stimulation silence from txt lines
                        obj.PreStimSilence= baphy.GetExpInfo(tline);
                    elseif logical(contains(tline,'exptparams(1).TrialObject(1)')) && logical(contains(tline,'PostStimSilence')) && ~logical(contains(tline,'UserDefinableFields'))
                    %Get duration, pre/post-stimulation silence from txt lines
                        obj.PostStimSilence= baphy.GetExpInfo(tline);
                    elseif logical(contains(tline,'exptparams(1).TrialObject(1)')) && logical(contains(tline,'Duration')) && ~logical(contains(tline,'UserDefinableFields'))
                    %Get duration, pre/post-stimulation silence from txt lines
                        obj.Duration= baphy.GetExpInfo(tline);       
                    end
                    tline = fgetl(fid);
                end
                
                obj.BlockDura = round(obj.PreStimSilence + obj.PostStimSilence + obj.Duration).* 10;
                obj.blockN = length(evt_line);
                disp([num2str(length(evt_line)) ' blocks detected in baphy']);
                
                %Get frequencies and volumes info from txt lines
                obj.evt_line = char(evt_line);
                [freq,volu] = baphy.GetVoluFreq(obj.evt_line);
                obj.freq = freq;
                obj.volu = volu;
                
                fclose(fid);
            catch
                disp('Baphy files not provided');
            end
        end
        
        function obj = set_freq(obj,freq)
        %Not used for future development
            obj.freq = freq;
        end
        
        function obj = set_volu(obj,volu)
        %Not used for future development
            obj.volu = volu;
        end
        
    end

    methods(Static)
        
        function [freq,volu] = GetVoluFreq(evt_line)
        %Extract frequency and volume values from txt lines
            ntrial = size(evt_line,1);

            for i = 1:ntrial
                %Extract current line
                cur_line = evt_line(i,:);
                mk1 = strfind(cur_line,',');
                mk2 = strfind(cur_line,':');
                mk3 = strfind(cur_line,'dB');
                freq(i,1) = str2double(evt_line(i,mk1(1)+2:mk2-1));
                volu(i,1) = str2double(evt_line(i,mk2+1:mk3-1));
            end
        end
        
        
        function info = GetExpInfo(tline)
        %Extract experiment infomation from txt lines
            mk1 = strfind(tline, '=');
            mk2 = strfind(tline, ';');
            info = str2double(tline(mk1+2:mk2-1));           
        end
        
        
    end
    
end