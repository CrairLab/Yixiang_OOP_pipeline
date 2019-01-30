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
%R3 09/15/18 New function genSeedsMap and genSeedingROIs (for automatically
%generating seeds for seed-based corr analysis) Compatiable with movieData
%class R13 or higher
%R4 09/19/18 Previous algo to calculate coordinates of downsampled seeds is
%wrong. Corrected in related functions in both ROI and movieData Class
%compatible with movieData R14 or higher 
%R4 09/30/18 Modify the genSeedingROIs function to allow seeds sampling in 
%sub region using input roi file named 'SubRegions.zip'
%R4 10/01/18 Tackle the situation where num_rois is already smaller than 
%total_seeds from the very beginning in function genSeedingROIs
%R4 10/15/18 Improve the ApplyMask function
%R4 01/18/19 Modify the ROIMask function
%R5 01/21/19 Modify the genSeedingROIs function. Only compatible with
%movieData R19+
%R5 01/23/19 Improve ROIMask function
%R5 01/24/19 Improve ROIMask function
%R6 01/26/19 Add function makeSizeConform. Modify several functions to allow
%automatic size conformation (mask size and matrix size). Only compatible
%with movieData R21+ and Integration R 17+
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
                 try
                    [~, Mask] = ROI.ROIMask(ROIData,sz(1:2));
                    Masks3D = repmat(Mask,[1 1 sz(3)]);
                    A2 = Masks3D.* A;
                 catch
                    disp('Mask size does not conform with movie size');
                    disp('Try to conform movie and mask');
                    [~, Mask] = ROI.ROIMask(ROIData);
                    [flag, Mask] = ROI.makeSizeConform(Mask,sz(1:2)); 
                    if flag
                        Masks3D = repmat(Mask,[1 1 sz(3)]);                       
                        A2 = Masks3D.* A;
                        disp('Succeeded in complying the mask size extracted from roi file with matrix size!')
                    else
                        warning('Faied to comply the mask size extracted from roi file with matrix size!')
                        A2 = A;
                    end
                 end
             else
                 A2 = A;
             end
            
       end
       
       
       function [subs,Mask] = ROIMask(ROIData,sz)
       
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
           
           if ~isempty(ROIData)
               %Default size of the image
               default_sz = [540,640];
               Mask = zeros(default_sz);
               for i = 1:length(ROIData)
                   try
                       ThisROIStruct = ROIData{i};
                   catch
                       ThisROIStruct = ROIData;
                   end
                   Coordinates = ThisROIStruct.mnCoordinates;
                   Mask = Mask + poly2mask(Coordinates(:,1),Coordinates(:,2),default_sz(1),default_sz(2));                            
               end
               %Try to conform mask size and matrix size
               [flag, Mask] = ROI.makeSizeConform(Mask,sz(1:2));
               if ~flag
                   Mask = ones(sz);
                   disp('Dimensions do not agree. Current size is: '); sz
                   disp('Make Mask = ones(sz);');
               end
           else
               Mask = ones(sz);
           end
           [subs(:,1),subs(:,2)] = ind2sub(sz,find(Mask));
           
           
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
        
       
       
       function Seeds = genSeedingROIs(total_seeds,sz)

        % Generate rois that serve as seeds for seed-based correlation maps
        %
        % Input:
        %   total_seeds     total number of seeds
        %
        % Output:
        %   Seeds           .mat file that store all the seeds 
        %

            if nargin == 0        
                total_seeds = 850;
                sz = [270 320];
            elseif nargin == 1
                sz = [270 320];
            end

            %Detect if there is pre-defined sub-region for seeds generation
            %If nothing is detected, read in any .roi files or .zip files
            %in current folder
            ROIRegime = ROI('SubRegions.zip');
            if isempty(ROIRegime.ROIData)
                ROIRegime = ROI();
                disp('Did not define sub-region for seeds sampling...')
             end
            
            [~, Mask] = ROI.ROIMask(ROIRegime.ROIData,sz);
            %szM = size(Mask);
            
            %It's important here to first do spatial downsampling than
            %focus on non-zero region of Mask (consistent with order in the
            %Integration pre-processing function
       
            %Mask = movieData.spatialDown(Mask,spatialFactor);
            %if ~(szM(1)*szM(2)/spatialFactor^2 == sz(1)*sz(2))
            %    disp('The input matrxi size (first 2d) = '); sz
            %    disp('The identified Mask size = '); szM
            %    Mask = Mask(any(Mask,2),any(Mask,1));
            %end            
            %disp('The new Mask size = '); size(Mask)
            %Mask = (Mask == 1);
            %Mask = imresize(Mask, downSampleRatio, 'bilinear');
            
            %if isempty(ROIRegime.ROIData)
                %default image size 540*640
            %    Mask = zeros(ceil(540*downSampleRatio),ceil(640*downSampleRatio)) + 1;
            %end
            [m,n] = size(Mask);
            %Define a minimum distance between seeds
            d = ceil(sqrt(m*n./100^2));

            num_rois = ROI.genSeedsMap(Mask,d);
            
            %Iteratively increase the distance in order to get fewer seeds
            %until the number of seeds(rois) aquired is smaller than
            %total_seeds defined by the user
            while num_rois > total_seeds
                d = ceil(d*1.05);
                [num_rois, ~] = ROI.genSeedsMap(Mask,d);
            end
            
            %In case the num_rois is already smaller than total_seeds at 
            %the very beginning
            d = d-1; %Make sure the num_rois > total_seeds;
            if d == 0
                d = 1;
            end
            [num_rois, rois_ini] = ROI.genSeedsMap(Mask,d);
            
            Seeds = rois_ini;
            disp(['Generated ' num2str(num_rois) ' seeds in the defined region(s)'])

        end

        
        
        function [num_rois,rois_ini] = genSeedsMap(Mask,d)

        % Generate an array of seeds evenly cover the Mask with distance between 
        % adjacent rois equals to d
        %  
        % Input:
        %   Mask     Input mask generated by the ROI.ROIMask function
        %   d        Distance spanning adjacent rois
        %
        % Output:
        %   num_rois     number of seeds in the mask 
        %   rois_ini     the initial indices of the seeds 

            [m,n] = size(Mask);
            EM = zeros(m,n);

            roi_id = 1:d:m-1;
            k = 0;
            while 1 + k*d <= n 
                try
                    roi_id_cur = roi_id + k*d*m; 
                    EM(roi_id_cur) = 1;
                    k = k+1;
                catch
                    disp(num2str(k));
                end    
            end

            %imshow(mat2gray(EM))
            Mask_rois = Mask.*EM;
            rois_ini = find(Mask_rois == 1);
            num_rois = floor(sum(Mask_rois(:)));

        end
        
        
        function [flag,Mask] = makeSizeConform(Mask,sz)
        %Estimate spatial factor being applied to the movie. The input Mask
        %is generated based on 540*640 image. However size of matrix (sz)
        %could change after downsampling and focusOnroi function. This
        %function try to make the two sizes comply and return flag == 1 if
        %succeed (0 if fail) as well as the new complied Mask
        %
        %Inputs:
        %   Mask        Current mask generated from ROIData
        %   sz          size of movie matrix
        %
        %Outputs:
        %   flag        1 if succeed; 0 if fail
        %   Mask        renewed Mask
            flag = 0;
            szM = size(Mask);
            
            %Check if Mask size and matrix size are already consistent
            if sum(sz == szM) == 2
                flag = 1;
                return;
            end
            
            %Estimate spatialFactor based on matrix size and mask size
            if exist('parameter.mat','file')
                load('parameter.mat')
                spatialFactor = param.spacialFactor;
                flag = 1;
            else
                Mask_f = movieData.focusOnroi(Mask);
                szMf = size(Mask_f);
                lwRatio = round(szMf./sz);
                if lwRatio(1) == lwRatio(2)
                    spatialFactor = lwRatio(1);
                else
                    warning('Length ratio and width ratio between Mask size and Matrix size')
                    disp('Assuming spatialFactor = 2');
                    spatialFactor = 2;
                end
            end
            
            %Try different order to change mask size to make it comply with
            %matrix size
            Mask_b = Mask;
            Mask = movieData.focusOnroi(Mask);
            Mask = movieData.spatialDown(Mask,spatialFactor);
            disp('Try to conform different sizes. Doing focusOnroi first...')
            szM = size(Mask);
            if sum(sz == szM) ~= 2
                Mask = Mask_b;
                Mask = movieData.spatialDown(Mask,spatialFactor);
                Mask = movieData.focusOnroi(Mask);
                disp('Doing focusOnroi first failed. Doing spatialDown first instead...')
                if sum(sz == szM) ~= 2
                    warning('Failed to conform the matrix size with ROI mask size!')
                    warning('Please checked the roi and the matrix size manually.')
                    flag = 0;
                else
                    flag = 1;
                end
            else
                flag = 1;
            end
            
            %In spatialDown function, fractions can be induced
            Mask = (Mask == 1);
        
        end
        
        
        
   end
       
end