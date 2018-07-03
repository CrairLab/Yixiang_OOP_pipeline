classdef baphy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Baphy is a class that stores inforamtion extracted from the baphy
%files, which includes arrays of frequencies and volumes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 12/14/17 Modify the code to tackle the problem when more than 150
%blocks of stimuli were presebted during experiment. 
%Add 'length(txt_line) <= 150' to the first while loop in function baphy
%
%R2 05/31/18 Lift the previous 150-block constrain due to new baphy
%configuration, new properties element blockN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        txt_line;%text extracted from .m files
        freq; %frequency
        volu; %volume 
        blockN; % Number of blocks(trails)
    end
    
    methods
        
        function obj = baphy(baphyfn)
        %Baphy object constructor
        
            try
                %Read in baphy (.m) files
                fid = fopen(baphyfn);
                tline = fgetl(fid);
                txt_line = {};
                
                %Detect key words 'PreStimSilence' and 'Note'
                while ischar(tline) %&& (length(txt_line) <= 160)
                    if logical(~isempty(strfind(tline,'PreStimSilence'))) && logical(~isempty(strfind(tline,'Note')))
                        txt_line = [txt_line;tline];
                    end
                    tline = fgetl(fid);
                end
                
                obj.blockN = length(txt_line);
                disp([num2str(length(txt_line)) ' blocks detected in baphy']);
                
                %Get frequencies and volumes info from txt lines
                obj.txt_line = char(txt_line);
                [freq,volu] = baphy.GetVoluFreq(obj.txt_line);
                obj.freq = freq;
                obj.volu = volu;

                %savename = [filename(1:length(filename)-4) '_BaphyIndex.mat'];
                %save(savename,'freq','volu');

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
        
        function [freq,volu] = GetVoluFreq(txt_line)
        %Extract frequency and volume values from txt lines
          
            ntrial = size(txt_line,1);

            for i = 1:ntrial
                mk1 = strfind(txt_line(i,:),',');
                mk2 = strfind(txt_line(i,:),':');
                freq(i,1) = str2double(txt_line(i,mk1(1)+2:mk2-1));
                volu(i,1) = str2double(txt_line(i,mk2+1:mk2+2));
            end
        end
    end
    
end