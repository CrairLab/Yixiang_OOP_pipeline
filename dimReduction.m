classdef dimReduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Class dimReuction get pre-processed movie data and perform various types % 
%of dimensionality reduction analysis on the data including: diffusion    % 
%map, t-SNE.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add new properties. Improve tSNE. Compatible with GUI_dimReduction.
%GUI app. 06/18/19 R2
%Author: Yixiang Wang
%Contact:yixiang.wang@yale.edu

    properties
        fd;        %parameters to further downsample input matrix
        sz_fd;     %size of the further-downsampled movie matrix
        A_rd;      %downsampled movie matrix excluding nan elements
        A_ref;     %Reference frame (mean)
        tflag;     %flag to transpose the matrix or not
        locflag;   %add location index flag
        locfactor; %add location weighted factor
        subIdx;    %Indices of non-nan or non-zero elements
        xy_sub;    %xy indicies of non-nan and non-zero pixels
        Dmap;      %Diffusion map
        Y;         %tSNE result
        tParam;    %Parameters used in tSNE
        dParam;    %Parameters used in diffusion map
    end
    
    methods
        function obj = dimReduction(A, tflag, locflag, locfactor, fd)
        %   Class constructor, read in 3D matrix A from pre-processing
        %
        %   Inputs:
        %   A          Input matrix (pre-processed)
        %   tflag      transpose flag
        %   locflag    add location flag
        %   locfactor  add location weighted factor
        %   par    parameters
        %
        %   Outputs:
        %   obj    dimReduction object
            
            if ~exist('fd','var')
                fd = [2,2]; %further downsample by factors of 2s both spatially and temporally
            end
            
            if ~exist('locfactor','var')
                locfactor = 1; %control the contribution of location indices
            end
            
            if nargin==1
                tflag = 0;
                locflag = 0;
            end
            
            obj.fd = fd;
            obj.tflag = tflag;
            obj.locflag = locflag;
            obj.locfactor = locfactor;
            
            %Further downsample input movie
            A_fd = movieData.downSampleMovie(A, fd(1), fd(2));
            sz = size(A_fd);
            
            %Save the mean as the reference frame
            obj.A_ref = nanmean(A_fd, 3);
            
            %Output the further-downsampled movie for quality check
            movieData.makePseudoColorMovie(A_fd, ['A_fd_' num2str(fd(1))...
                '_' num2str(fd(2)) '.avi'])
            clear A;          
                     
            %Save the size/dimensions of the further-downsampled movie
            obj.sz_fd = size(A_fd);
            
            %add location indices
            if locflag && (tflag == 0)
                A_fd = dimReduction.addLocIdx(A_fd, locfactor);
                sz(3) = sz(3) + 2;
            end 
            
            %Delete pixels that are outside roi (are NaN or 0)
            A_re = reshape(A_fd, [sz(1)*sz(2), sz(3)]);
            subIdx = and((~isnan(A_re(:,1))),(~(A_re(:,1)==0)));
            A_re = A_re(subIdx, :);
            obj.subIdx = subIdx;
            obj.xy_sub = dimReduction.genXYind(sz, subIdx);
            
            if tflag %Transpose the matrix
                disp('Transposed the matrix')
                obj.A_rd = A_re';
            else
                obj.A_rd = A_re;
            end
            
            disp(['Actual shape of the matrix to be processed: ' ...
                num2str(size(obj.A_rd))]);
            clear A_re;
            
            [obj.Dmap, obj.dParam, ~] = dimReduction.diffmap(obj.A_rd);
            [obj.Y, obj.tParam, ~] = dimReduction.doTSNE(obj.A_rd);
            
        end       
    end
    
    methods(Static)
              
        function [Dmap, dParam, vals] = diffmap(xs, t, m, sigma)
        %   Perform diffusion map. 
        %
        %   Inputs:
        %       xs       Input matrix, each row represents a data point
        %       t        Time step
        %       m        Number of coordinate functions/non-trivial dimensions
        %       sigma    Sigma for the gaussian kernal where estimate connectivity
        %
        %   Outputs:
        %       Dmap     Diffusion map (embedded coordinates)
        %       dParam   Output the parameters used in diffusion map
        %       vals     Corresponding eigenvalues
        
            if nargin<2
                t = 2;
                m = 3;
                sigma = [];
            end
            
            %Calcultae pairwise distances
            if any(isnan(xs))
                disp('Detected NaN in input matrix, pdist might fail...')
            end
            pdist_xs = pdist(xs);
            std_dist = std(pdist_xs);

            if isempty(sigma)||isnan(sigma)
                sigma = 2*std_dist;
            elseif sigma < 0.25*std_dist
                warning('Input sigma is too small. Rest sigma to 0.25*std_dist!')
                sigma = 0.25*std_dist;
                warning(['New sigma is = ' num2str(sigma)]);
            end
            disp(['Sigma is: ', num2str(sigma)])
            
            %Plot the histogram of pdist_xs
            figure;
            histogram(pdist_xs)

            %Calculate weight matrix
            distMatrix = squareform(pdist_xs);
            W = exp(-(distMatrix.^2)./(2*sigma^2));

            %Calculate weighted M matrix
            d = sum(W,1);
            M = W./d';
            k = m+1;

            %Get the first k eigenvectors/eigenvalues
            [psi, L] = eigs(M, k);

            %Omit the trivial first dimension
            vals = diag(L);
            vals = vals(2:end);
            Dmap = psi(:,2:end).*vals'.^t;
            
            %save parameters
            dParam.t = 2;
            dParam.m = 3;
            dParam.sigma = sigma;
            
            disp('Diffusion map analysis is done!')

        end
                
        function [Y, par, loss] = doTSNE(xs, par)
        %   Do tSNE analysis using input matrix xs and parameters par
        %
        %   Inptus:
        %       xs     input matrix
        %       par    parameter struct
        %
        %   Outputs:
        %       Y      embeddings, default 3D
        %       par    output the parameters used in tSNE
        %       loss   KL divergence loss after optimization
        
            try
                Distance = par.Distance; %Distance method
                nd = par.nd; % Number of dimension
                np = par.np; % NumberPrint
                px = par.px; % Perplexity
                vflag = par.vflag; % Verbose flag
                stdFlag = par.stdFlag; %Standardize flag
                Exaggeration = par.Exaggeration; % Exaggeration factor 
                options = par.options; % Other option
            catch
                disp('Did not provide/ can not read in parameters...')
                disp('Using default parameters')
                Distance = 'euclidean'; par.Distance = Distance;
                nd = 3; par.nd = nd;
                np = 50; par.np = np;
                px = 30; par.px = px;
                stdFlag = false; par.stdFlag = stdFlag;
                vflag = 1; par.vflag = vflag;
                Exaggeration = 8; par.Exaggeration = Exaggeration;
                options = statset('MaxIter', 500); par.options = options;
            end
            disp(['Number 2 dimensions = ' num2str(nd)])
            disp(['Distance = ' Distance])
            disp(['Standardize flag = ' num2str(stdFlag)])
            disp(['Exaggeration = ' num2str(Exaggeration)])
            
            [Y, loss] = tsne(xs, 'Distance', Distance, 'NumDimensions', nd, 'Standardize', stdFlag, ...
            'Perplexity', px,'Exaggeration', Exaggeration, 'NumPrint', np, 'Verbose', vflag, ...
            'Options', options);
                      
        end
        
        
        function A_rcs = rcsFromSVD(U, V, S, iniDim)
        %   Reconstruct matrix from S,V,D matrix starting from initial 
        %   dimension specified by iniDim
           if nargin < 4 
              iniDim = 2;
           end
           disp(['Initial Dimension = ' num2str(iniDim)])
           sz = size(V);
           if length(sz) == 3
               V = reshape(V, [sz(1)*sz(2) sz(3)]);
               A_rcs = V(:,iniDim:end)*S(iniDim:end,iniDim:end)*U(:,iniDim:end)';
               A_rcs = reshape(A_rcs,[sz(1),sz(2),size(A_rcs,2)]);
           else
               disp('Input V matrix should have 3 dimensions')
           end           
        end
        
        
        function plotEmbedding(Embedding, plotFlag, binSize)
        %   Plot 2D/3D embeddings
        %
        %   Inputs:
        %       Embedding   Embedded/Mapped matrix
        %       plotFlag    flag to subsampling in time direcion
        %       binSize     binning size for subsampling 
        
            if nargin<3
                binSize = 300;
            end
        
            cmap = 1:size(Embedding,1);
            %Combinations for 2D embeddings
            if size(Embedding,2) == 3    
                comb = [1 2;1 3;2 3];
                %Plot 3D embedding
                figure();
                ax3D = gca();
                scatter3(Embedding(:,1),Embedding(:,2),Embedding(:,3),[],cmap,'filled')
                title(ax3D, '3D embedding');
                colorbar;
            else
                comb = [1 2];
            end
            
            %Plot 2D embeddings
            for i = 1:size(comb,1)
                figure();
                ax = gca();
                scatter(Embedding(:,comb(i,1)),Embedding(:,comb(i,2)),[],cmap)
                title(ax, ['2D embedding' num2str(comb(i,:))]);
                colorbar;
            end       
            
            %Subsampling every 300 time points to see whether repeated
            %patterns show up in diffusion maps
            if plotFlag                
                for i = 1:ceil(size(Embedding,1)/binSize)
                    first = (i-1)*binSize+1;
                    last = first+binSize;
                    if last>size(Embedding,1)
                        last = size(Embedding,1);
                    end
                    cur_fig = figure('visible','off'); 
                    scatter3(Embedding(first:last,1),Embedding(first:last,2),...
                        Embedding(first:last,3),[],first:last,'filled')
                    xlim(gca, ax3D.XLim);
                    ylim(gca, ax3D.YLim);
                    zlim(gca, ax3D.ZLim);
                    colorbar;
                    hold on
                    saveas(cur_fig, ['Trajectory_' num2str(first) '-' num2str(last) '.jpg']);
                end
            end                   
        end
        
        
        function showFrameEmbedding(A, Embedding, binSize)
        %   Show embedded points and corresponding frame
        %
        %   Inputs:
        %       A           input movie matrix
        %       Embedding   Embedded low-dimension representation
        %       binSize     binning size for subsampling 
        
            if nargin<3
                binSize = 200;
            end
        
            figure;
            max_x = max(Embedding(:,1));
            min_x = min(Embedding(:,1));
            max_y = max(Embedding(:,2));
            min_y = min(Embedding(:,2));

            for i = 1:size(A,3)
                %Plot embedding (using first two coordinates)      
                subplot(2,1,1)             
                scatter(Embedding(i,1), Embedding(i,2),[],'ro','filled')
                xlim(gca, [0.8*min_x 1.1*max_x])
                ylim(gca, [0.8*min_y 1.1*max_y])
                hold on 
                
                %Update embedding plot every 200 iterations
                if mod(i,binSize)== 0
                    hold off
                end
                                
                %Plot corresponding frame on the right
                subplot(2,1,2)
                imshow(mat2gray(A(:,:,i)))
                pause(0.1)
            end
        end
        
                
        function redKmeans(X)
        %   Use Kmenas plus plus algo to cluster reduced data
        %   
        %   Inputs:
        %       X     embeddings (eg. diffusion maps)
        
            mean_X = mean(X,1);
            F_values = [];   %F measure
            sumSQE = [];     %Normalized intra-cluster error
            
            iniK = 2;        %# of clusters starts from 2
            endK = 30;       %# of clusters should not exceed 30
            Kmeans_idx_all = cell(endK,1);
            Kmeans_C_all = cell(endK,1);
            Kmeans_sumd_all = cell(endK,1);
            Kmeans_D_all = cell(endK,1);
            
            %Calculte F values and normalized intra-cluster erros for
            %different K (number of clusters)
            parfor k = iniK:endK
                %scatter3(pj_CC(:,1),pj_CC(:,2),pj_CC(:,3));
                disp(['Try cluster number = ' num2str(k)])
                %Replicate 50 times for randomized optimization purposes
                [idx,C,sumd,D] = kmeans(X ,k, 'Replicates', 50);

                nPoints = [];  %Number of points in each clusters
                intraVar = []; %Intra cluster variance
                for i = 1:k
                    nPoints = [nPoints; sum(idx == i)];
                    intraVar = [intraVar sum(D(idx == i, i))];
                end
                
                %Calculate F values
                F_values = [F_values sum(nPoints.*sum((C - mean_X).^2,2))/sum(intraVar,2)*(size(X,1)-k)/(k-1)];

                %Calculate normalized intra-cluster error
                sumSQE = [sumSQE sum(intraVar,2)/sum(sum((X - mean_X).^2,2))];

                Kmeans_idx_all{k,1} = idx;
                Kmeans_C_all{k,1} = C;
                Kmeans_sumd_all{k,1} = sumd;
                Kmeans_D_all{k,1} = D;
            end
        
            %plot normalized intra-cluster errors from different Ks
            h = figure;
            plot([iniK:endK],sumSQE,'LineWidth',3)
            xlabel('Cluster numbers')
            ylabel('Normalized intra-cluster error')
            title('Clustering evaluation normalized intra-cluster error')
            saveas(h,'IntraInter.png')

            %plot Calinski-Harabasz index from different Ks
            h = figure;
            plot([iniK:endK],F_values,'LineWidth',3)
            xlabel('Cluster numbers')
            ylabel('Calinski-Harabasz index (F values)')
            title('Clustering evaluation Calinski-Harabasz')
            saveas(h,'CH.png')

            %Choose best K based on empirical standards
            %scale_F = 1:length(F_values);
            K_sqe = find(sumSQE <0.05,1) + iniK - 1;

            F_values_ = F_values;
            max_K_f = find(F_values == max(F_values),1) + iniK - 1;
            if max_K_f == iniK
                F_values_(1:5) = -1;
                max_K_f = find(F_values == max(F_values_),1) + iniK - 1;
            end

            K = max([K_sqe max_K_f]);
            disp(['The program think the besk K is' num2str(K)])
            K = input('What do you think is the best K?');
            disp(['Best cluster number = ' num2str(K)])

            C = Kmeans_C_all{K,1};
            idx = Kmeans_idx_all{K,1};
            sumd = Kmeans_sumd_all{K,1};
            D = Kmeans_D_all{K,1};

            save('Kmeans_result.mat','K','idx','C','sumd','D')
        end
        
        
        function [A_addloc, w] = addLocIdx(A, f)
        %   This function read in movie matrix and attach location indices
        %   to the end of the matrix. Input matrix should be 3D.
        %
        %   Inputs:
        %       A     input movie matrix, must be 3D
        %       f     factor to control the contribution of location indices
        %             to final pairwaise distances
        %   
        %   Outputs:
        %       A_addloc    output matrix with location indicies
        %       w           weighted contribution of location indices     
        
            disp('Input matrix should have each samples on rows')
            
            %Compute pairwise distance to calculate contribution weight
            %Only randomly sample 1000 points to speed up this process            
            sz = size(A);
            A_re = reshape(A, [sz(1)*sz(2), sz(3)]);
            A_re(isnan(A_re(:,1)), :) = [];
            A_re(A_re(:,1)==0, :) = [];
            sz_new = size(A_re);
            
            if sz_new(1)*sz_new(2) < 1000
                subID = (1:sz_new(1)*sz_new(2))';
            else
                subID = randi(sz_new(1),[1000,1]);
            end
            A_sub = A_re(subID,:);
            pdist_ori = pdist(A_sub);
            
            xy_sub = dimReduction.genXYind(sz);
            pdist_xy = pdist(xy_sub);
        
            if nargin == 1
                f = 1;                
            end         
            
            %Contribution weight of spatial indices scaled by factor f
            w = (mean(pdist_ori)./mean(pdist_xy)).^2;
            w_f = f.*w;
            disp(['Weight of location indices = ', num2str(w_f)])
            
            A_addloc = reshape(A, [sz(1)*sz(2), sz(3)]);
            A_addloc(:, end+1:end+2) = w_f.*xy_sub; 
            A_addloc = reshape(A_addloc, [sz(1) sz(2) sz(3)+2]);
               
        end
        
        
        function xy_sub = genXYind(sz, subIdx)
        %   Generate XY indices of a matrix given size sz. Can choose only
        %   output subgroup of these indices specified by subIdx
        %
        %   Inputs:
        %       sz        size of a desired matrix, at least 2D
        %       subIdx    Indices of a subgroup of the matrix
        
            if nargin == 1
                subIdx = (1:sz(1)*sz(2))';
            end
            
            %generate x coordinates
            x_idx = ones(sz(1),sz(2)).*(1:sz(1))';
            x_idx = reshape(x_idx, [sz(1)*sz(2), 1]);
            
            %generate y coordinates
            y_idx = ones(sz(1),sz(2)).*(1:sz(2));
            y_idx = reshape(y_idx, [sz(1)*sz(2), 1]);
            
            xy_idx = [x_idx y_idx];
            xy_sub = xy_idx(subIdx,:);
        
        end
        
    end

end