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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
   properties
       ROIData; %cell array containing ROI structures
   end
    
   methods
       function obj = ROI()
           
       %     Read in ROI data from .zip or .roi files
       %     object constructor
       %
       %     Inputs:
       %         ROIname    filename of the ROI zip file
       %
       %     Outputs:
       %         ROIData    cell array containing ROI structures
           
      
           if ~(isempty(dir('*.zip')) && isempty(dir('*.roi')))              
               ROIname = dir('*.zip');
               if isempty(ROIname)
                   ROIname = dir('*.roi');
               end
               obj.ROIData = ReadImageJROI(ROIname.name);
           else
               disp('No RoiSet.zip or .roi file detected');
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
       
       
       
       
       
   end
       
end