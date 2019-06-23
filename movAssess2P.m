function [A,output_all,NormTform_all] = movAssess2P(A, flag)
%   Assess movie and get rid of frames with large movements if flag
%   == 1. Use discrete fourier transform algorithm to compare phase
%   correlation. When filtering frames, use a step function as the
%   convolution kernel (width == 5).
%
%   Inputs:
%     A                   3D matrix
%     flag                flag == 1 discard frames
%
%   Outputs:
%     A                   Registered matrix
%     tform_all           The array that stores all transformation matrices
%     NormTform_all       Norm of each matrices (minus I)
%
%   Author: yixian.wang@yale.edu
%   Version: 0.0.3
%   Date: 06/19/19

    if nargin == 1
        flag = 0;
    end

    sz = size(A);
    mask = isnan(A);
    %Convert nan to 0 if there is any
    A(mask) = 0;
    %Use the median frame as reference
    
    %Using the median frame of the first movie as reference for movies!
    if exist('reference_frame.mat', 'file')
        load('reference_frame.mat', 'A_median');
        warning('Using the median frame of the first movie as reference for movies!')
    else
        A_median = nanmedian(A,3);
        save('reference_frame.mat', 'A_median')
    end

    %Do 2d fft to A_median
    A_mf = fft2(A_median);
    
    parfor i = 1:size(A,3)
        %Get current frame 
        A_c = A(:,:,i);
        %Do 2d fft to current frame
        A_cf = fft2(A_c);
        %Do DFT registration
        [output, Greg] = dftregistration(A_mf,A_cf,2);
        %Output statistics
        output_all(i,:) = output;
        %Construct registered movie Greg_all
        Greg_all(:,:,i) = abs(ifft2(Greg));
    end
    %Only save the portion within the ROI
    A = Greg_all.*~mask;
    %Calculate norm of the transformation matrix
    NormTform_all = sqrt(output_all(:,3).^2 + output_all(:,4).^2);
    
    tic;
    [A, NormTform_all, output_all] = correctArtifact(A, NormTform_all, output_all);
    toc;
    
    
    
    %If the norm is larger than 0.7 (~0.5 pixel at either direction)
    %label this frame as moving. Save indices of stable frames as
    %movIdx_saved
    movIdx_saved = NormTform_all < 0.7;
    saveRatio = sum(movIdx_saved)/sz(3);

    %If more than 5% of the movie has substantial movements, warn the user
    if saveRatio < 0.95
        disp('This movie contains more than 5% moving frames!')
        disp('Please check the quality of the movie!')
    end

    %if flag == 1 discard the neighbouring 5 frames
    if flag
        %Indices of frames to be replaced
        movIdx_replace =  ~movIdx_saved;
        
        %Step function kernel for convolution, width == 5
        filter = [1,1,1,1,1]; %Discard the neighbouring 5 frames 
        movIdx_replace = conv(movIdx_replace, filter, 'same');
        movIdx_replace = movIdx_replace > 0;
        %Delete moving frames
        A(:,:,movIdx_replace) = [];
    end            

    disp(['Mean tform magnitude (minus I) = ' num2str(mean(NormTform_all))]);
    disp(['Relatively stable frames ratio = ' num2str(saveRatio)]);
    

end


function A_ref = selectRef(A)
%Calculate the reference frame from input movie A
    
    %Calculate the averaged fluorescent intensity of each frame
    A_tmp = nanmean(A,1);
    A_tmp = nanmean(A_tmp,2);
    A_int = A_tmp(:);
    
    %Get the subset for top 5% averaged intersity 
    topIdx = A_int >= prctile(A_int,95);
    A_top = A(:,:,topIdx);
    
    %Register the subset movie
    A_top_ref = A_top(:,:,1);
    A_rf = fft2(A_top_ref);
    parfor i = 1:size(A_top,3)
        %Get current frame 
        A_c = A_top(:,:,i);
        %Do 2d fft to current frame
        A_cf = fft2(A_c);
        %Do DFT registration
        [~, Greg] = dftregistration(A_rf,A_cf,2);
        %Construct registered movie Greg_all
        Greg_all(:,:,i) = abs(ifft2(Greg));
    end
    
    %Use the mean of the registered movie as reference
    A_ref = nanmean(Greg_all,3);
    
end


function [A, NormTform_all, output_all] = correctArtifact(A, NormTform_all, output_all)
%Correct artifact movement and update registered movie as well as all
%statistics

    %If the norm is larger than 35 (~25 pixels on each direction), label
    %the frame as an artifact
    artIdx = NormTform_all > 35;
    IdxNumber = find(artIdx);
    round = 0;
    
    while any(artIdx) && round<= 10
        round = round + 1;
        %Search the closest non-artifact frame for each artifact frames
        for i = 1:length(IdxNumber)
            curNumber = IdxNumber(i);
            flag = 1;
            checkpoint = 1;

            %Check the neighbourhood of the current frame, if a non-artifact
            %frame is found, change falg to 0 and exit the while loop
            while flag && checkpoint < 6000

                %Try check in forward direction
                try
                   forwardCheck = artIdx(curNumber - checkpoint) == 0;
                   if forwardCheck == 1
                       flag = 0;
                       refIdx = curNumber - checkpoint;
                   end
                catch
                   forwardCheck = 0;
                end

                %Try check in backward direction            
                try
                   backwardCheck = artIdx(curNumber + checkpoint) == 0;
                   if backwardCheck == 1
                       flag = 0;
                       refIdx = curNumber + checkpoint;
                   end
                catch
                   backwardCheck = 0;
                end           
                checkpoint = checkpoint + 1;
            end

            %Use the closest non-artifact frame as reference frame to register
            %this particular artifact frame
            A_curRef = A(:,:,refIdx);
            A_curArt = A(:,:,curNumber);

            [output, Greg] = dftregistration(fft2(A_curRef),fft2(A_curArt));

            %Upate the registered matrix
            A(:,:,curNumber) = abs(ifft2(Greg));

            %Update the transformation matrix
            output_all(curNumber,:) = output;
        end    

        %Update the movement indicator
        NormTform_all = sqrt(output_all(:,3).^2 + output_all(:,4).^2);
        
        %Update artifact indices
        artIdx = NormTform_all > 35;
        IdxNumber = find(artIdx);
    end
        
end