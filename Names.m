classdef Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Names is a class that restore filenames of movie\spike2\baphy\Audio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R2 11/27/17   Add a new property outputFolder    
%R2 05/09/18 Code is modified to handle the situation in which spike2 files 
%are not provided. 
%R2 05/23/18 Code is modified to be compatible with different recording
%configurations (spontaneous, with stimulation, wavelength switching) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        filename; %filename of the movie
        Spike2name; %filename of the corresponding spike2 
        Baphyname; %filename of the corresponding baphy
        Audioname; %filename of the correspondign audio (not used now)
        ROIname; %filename of the roi
        outputFolder; %output path
    end
    
    methods
        function obj = Names(f,flag)
        %constructor function read in filenames from various txt files
        %   Inputs:
        %        f      file handle
        %        flag   flag for generating txt files
        
            fileList = readtext('files.txt',' ');
            obj.filename = fileList{f,1};
            
            if flag
                spike2List = readtext('Spike2files.txt',' ');
                baphyList = readtext('baphyfiles.txt',' ');
                try
                    obj.Spike2name = spike2List{f,1};
                catch
                    disp('')
                    disp('Spike2 files not provided!');
                    disp('')
                    obj.Spike2name = [];
                end
                obj.Baphyname = baphyList{f,1};
            else
                obj.Spike2name = [];
                obj.Baphyname = [];
            end
            
            obj.ROIname = 'RoiSet.zip';
            obj.outputFolder = ['output' num2str(f)];
        end
    end
end