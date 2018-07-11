classdef ROI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Class ROI stores ROIData from 'RoiSet.zip' (default) using ReadImagJROI
%function. It provides several static functions to get 3D ROI masks
%and apply the masks to movies.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 05/25/18 
%R2 06/06/18 Detect .roi or .zip files automatically compatiable with 
%integration R7   
%R2 07/11/18 New function generateROIArray (for seed-based corr analysis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
   properties
       ROIData; %cell array containing ROI structures
   end
    
   methods
       function obj = ROI(filename)        
       %     Read in ROI data from .zip or .roi files
       %     object constructor
       %
       %     Inputs:
       %         ROIname    filename of the ROI zip file
       %
       %     Outputs:
       %         ROIData    cell array containing ROI structures
           
            if nargin == 0
                if ~isempty(dir('*.roi'))              
                   ROIname = dir('*.roi');
                   obj.ROIData = ReadImageJROI(ROIname.name);
                else
                   disp('No .roi file detected');
                   disp('Try to capture .zip ROI files')
                   filename = 'RoiSet.zip';
                   try
                      obj.ROIData = ReadImageJROI(filename);
                   catch
                      disp('Unable to find/open the roi files.')
                   end
                end
                else
                try
                   obj.ROIData = ReadImageJROI(filename);
                catch
                   disp('Unable to find/open the roi files.')
                end
            end
       end
   end
   
   
   methods(Static)    
       
       function A2 = ApplyMask(A,ROIData)
       
       %   Apply 3D mask to matrix A (movies). First call the ROIMask
       %   function to generate 2D mask, then extend 2D mask to 3D mask.
       %    
       %    Inputs:
       %        A      Input movie (3D)
       %        ROIData     Cell array containing ROI structures
       
        
             sz = size(A);
             if ~isempty(ROIData)
                 [~, Mask] = ROI.ROIMask(ROIData,sz(1:2));
                 Masks3D = repmat(Mask,[1 1 sz(3)]);
                 A2 = Masks3D.* A;
             else
                 A2 = A;
             end
            
       end
       
       
       function [index,Mask] = ROIMask(ROIData,sz)
       
       %    Get 2D ROI mask from ROIData. If there are more than 1 ROI,
       %    merge the ROIs to a single mask. ROI is considered a polycon
       %    area defined in ImageJ
       %    
       %    Inputs:
       %        ROIData     Cell array containing ROI structures
       %        sz          Matrix(frame) size (2D)
       %    
       %    Ouputs:
       %         index      Indices of the vertices of polycons
       %         Mask       2D mask that contains all ROIs
          
       
           if nargin == 1
               sz(1) = 540;
               sz(2) = 640;
           end
           Mask = zeros(sz);
           for i = 1:length(ROIData)
               try
                   ThisROIStruct = ROIData{i};
               catch
                   ThisROIStruct = ROIData;
               end
               Coordinates = ThisROIStruct.mnCoordinates;
               Mask = Mask + poly2mask(Coordinates(:,1),Coordinates(:,2),sz(1),sz(2));                            
           end
           index = find(Mask);
       end
       
       function [roi,roiPolygon] = generateROIArray(ROI_all,sz)
        % This function generate two cell arrays: roi and roiPolygon to store 
        % cooresponding roi mask and vertex coordinates. Should be able to handle
        % polygon or rectangle roi.
        % 
        % Inputs:
        %   ROI_all    cell array contain ROI structs
        %   sz         size of original movie
        % 
        % Outputs:
        %   roi        cell array of corresponding masks
        %   roiPolygon cell array of vertex coordinates
        %  
            for i = 1:length(ROI_all)
                if isfield(ROI_all{i}, 'mnCoordinates')
                    roi{i} = poly2mask(ROI_all{i}.mnCoordinates(:, 1), ROI_all{i}.mnCoordinates(:, 2), sz(1), sz(2));   
                    roiPolygon{i} = ROI_all{i}.mnCoordinates;
                    %roiName{i} = ROI_all{i}.strName;
                else
                    roi{i} = zeros(sz(1), sz(2));
                    roi{i}(ROI_all{i}.vnRectBounds(1):ROI_all{i}.vnRectBounds(3), ROI_all{i}.vnRectBounds(2):ROI_all{i}.vnRectBounds(4)) = 1;
                    %roiName{i} = ROI_all{i}.strName;
                    roiPolygon{i} = ...
                    [ROI_all{i}.vnRectBounds(2),ROI_all{i}.vnRectBounds(1);...
                     ROI_all{i}.vnRectBounds(2),ROI_all{i}.vnRectBounds(3);...   
                     ROI_all{i}.vnRectBounds(4),ROI_all{i}.vnRectBounds(3);...
                     ROI_all{i}.vnRectBounds(4),ROI_all{i}.vnRectBounds(1)];
                end
                %roiName{i} = ROI_all{i}.strName;
            end

        end
             
   end
       
end