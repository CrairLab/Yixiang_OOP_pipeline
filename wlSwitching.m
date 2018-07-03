classdef wlSwitching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This class is specifically designed to handle wavelength switching data, 
%either with or without ROIs. It includes two major functions to separate
%different channels, to separte different ROIs and correspond them to
%the correct channels.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R1 05/24/18 
%R1 05/27/18 New methods to find corresponding ROI 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        A1; %Matrix for channel 1       
        A2; %Matrix for channel 2
        ROI1;
        ROI2;
    end
        
    methods
        function obj = wlSwitching(flag,A,ROI)
        %Class constructor
            if flag == 2
                [obj.A1,obj.A2] = wlSwitching.separateChannels(A);
                if ~isempty(ROI)
                    disp('You need to provide two ROIs for two channels!');
                    [obj.ROI1,obj.ROI2] = wlSwitching.findCorrROI(A(:,:,1),A(:,:,2),ROI);
                else
                    disp('ROI not provided');
                end
            else
                obj.A1 = []; obj.A2 = []; obj.ROI1 = []; obj.ROI2 = [];
            end
        end
    end
    
    methods (Static)
        
        function [A1,A2] = separateChannels(A)
        
        %    Seperate two channels and create new movies 
        %    by copying each frame twice.
        %    
        %    Input:
        %        A        Original movie
        %    
        %    Outpus:
        %        A1,A1    Separated channels
        
            sz = size(A);
            A1 = [];
            A2 = [];
            if mod(sz(3),2) == 0
               for i = 1:(sz(3)/2)
                    A1 = cat(3,A1,A(:,:,2*i-1));
                    A2 = cat(3,A2,A(:,:,2*i));
               end
            else 
                for i = 1:(sz(3)-1)/2
                    A1 = cat(3,A1,A(:,:,2*i+1));
                    A2 = cat(3,A2,A(:,:,2*i));
                end
            end
            A1 = wlSwitching.extroplateMovie(A1);
            A2 = wlSwitching.extroplateMovie(A2);
        end
        
        
        function extp_A = extroplateMovie(A)
        
        %    Called by the separateChannels function. Duplicate each frame.
        %    
        %    Input: 
        %        A     Input matrix
        %    
        %    Output:
        %        extp_A     Output matrix with duplicated frames
        
            
            extp_A = [];
            for i = 1:size(A,3)
                extp_A = cat(3,extp_A, A(:,:,i),A(:,:,i));
            end
        end
        
        function [ROI1,ROI2] = findCorrROI(F1,F2,ROI)
        
        %    Find correct correspondence between two ROIs and two Channels
        %    by socring overlapping between the given ROIs and >95%
        %    percentile strongest signal in the first two switchign frames.
        %    
        %    Inputs:
        %        F1,F2     First and Second frame
        %        ROI       ROIs (struct)
        %    
        %    Outputs:
        %        ROI1,ROI2    Separated and sorted ROIs
        
            F1_95 =  prctile(F1(:),95);
            F2_95 =  prctile(F2(:),95);
            F1 = F1 > F1_95;
            F2 = F2 > F2_95;

            Mask1 = wlSwitching.maskMatrix(ROI{1,1},size(F1));
            Mask2 = wlSwitching.maskMatrix(ROI{1,2},size(F2));
            FM11 = F1.*Mask1; FM12 = F1.*Mask2;
            FM21 = F2.*Mask1; FM22 = F2.*Mask2;
            Score_11 = sum(FM11(:))/sum(Mask1(:));
            Score_12 = sum(FM12(:))/sum(Mask2(:));
            Score_21 = sum(FM21(:))/sum(Mask1(:));
            Score_22 = sum(FM22(:))/sum(Mask2(:));
            
            ROI1 = cell(1);
            ROI2 = cell(1);
            
            if (Score_11 > Score_12) && (Score_21 < Score_22)
                
                ROI1{1} = ROI{1,1};
                ROI2{1} = ROI{1,2};
            elseif (Score_11 < Score_12) && (Score_21 > Score_22)
                
                ROI1{1} = ROI{1,2};
                ROI2{1} = ROI{1,1};
            else
                disp('Can not distinguish 2 channels...Try another method');
                Score = [Score_11 Score_12 Score_21 Score_22];
                switch find(Score == min(Score))
                    case 1
                        ROI1{1} = ROI{1,2};
                        ROI2{1} = ROI{1,1};
                    case 2
                        ROI1{1} = ROI{1,1};
                        ROI2{1} = ROI{1,2};
                    case 3
                        ROI1{1} = ROI{1,1};
                        ROI2{1} = ROI{1,2};
                    case 4
                        ROI1{1} = ROI{1,2};
                        ROI2{1} = ROI{1,1};
                end
                
                if isempty(ROI1)
                    if exist('Separated_ROI.mat','file')
                        disp('Using previous ROI sorting result...')
                        load('Separated_ROI.mat');
                    end
                end
            end
            
            if ~exist('Separated_ROI.mat','file')
                save('Separated_ROI.mat','ROI1','ROI2');
            end
        end
        
        function Mask = maskMatrix(ThisROIStruct,sz)
        
        %    Making mask matrix from the coordinate.
        %    
        %    Inputs:
        %        ThisROIStruct     The ROI input struct
        %        sz                matrix size
        %   
        %   Outpus:
        %        Mask              Mask matrix    
            
        
             Mask = zeros(sz);
             Coordinates = ThisROIStruct.mnCoordinates;
             Mask = poly2mask(Coordinates(:,1),Coordinates(:,2),sz(1),sz(2));          
        end
    end      

end