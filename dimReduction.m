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
        PT;        %PHATE result
        tParam;    %Parameters used in tSNE
        dParam;    %Parameters used in diffusion map
        pParam;    %Parameters used in PHATE
        adaptive;  %Whether use adaptive kernels for diffusion map analysis
    end
    
    methods
        function obj = dimReduction(A, varargin)
        %   Class constructor, read in 3D matrix A from pre-processing
        %
        %   Inputs:
        %   A          Input matrix (pre-processed)
        %   tflag      transpose flag
        %   locflag    add location flag
        %   locfactor  add location weighted factor
        %   adaptive   whether use adaptive kernels for diffusion map 
        %
        %   Outputs:
        %   obj    dimReduction object
            
            fd = [2, 2];
            locfactor = 1;
            adaptive = 1;
            tflag = 0;
            locflag = 0;
            
            % get input parameters
            for i=1:length(varargin)
                % fd for further downsampling
                if(strcmp(varargin{i},'fd'))
                   fd = lower(varargin{i+1});
                end
                
                % add location weighted factor
                if(strcmp(varargin{i},'locfactor'))
                   locfactor = lower(varargin{i+1});
                end
                
                % whether to use an adpative kernel
                if(strcmp(varargin{i},'adaptive'))
                   adaptive = lower(varargin{i+1});
                end
                
                % transpose flag
                if(strcmp(varargin{i},'tflag'))
                   tflag = lower(varargin{i+1});
                end
                
                 % whether use locations as constraints
                if(strcmp(varargin{i},'locflag'))
                   locflag = lower(varargin{i+1});
                end
                
            end
            
            if ~exist('fd','var')
                fd = [2,2]; %further downsample by factors of 2s both spatially and temporally
            end
            
            if ~exist('locfactor','var')
                locfactor = 1; %control the contribution of location indices
            end
            
            if ~exist('adaptive', 'var')
                adaptive = 1;
            end
            
            if nargin==1
                tflag = 0;
                locflag = 0;
            end
            
            %Construct properties
            obj.fd = fd;
            obj.tflag = tflag;
            obj.locflag = locflag;
            obj.locfactor = locfactor;
            obj.adaptive = adaptive;
            
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
            
            [obj.Dmap, obj.dParam] = dimReduction.diffmap(obj.A_rd, 2, 30, 5, adaptive);
            [obj.Y, obj.tParam, ~] = dimReduction.doTSNE(obj.A_rd);
            [obj.PT, obj.pParam, ~] = dimReduction.doPHATE(obj.A_rd, 'ndim', 3);
            
        end       
    end
    
    methods(Static)
              
        function [Dmap, dParam] = diffmap(xs, t, m, sigma, adaptive)
        %   Perform diffusion map. 
        %
        %   Inputs:
        %       xs       Input matrix, each row represents a data point
        %       t        Time step
        %       m        Number of coordinate functions/non-trivial dimensions
        %       sigma    Sigma for the gaussian kernal where estimate connectivity
        %                when use adpative kernels, sigma represents the
        %                kth nearest neighbour
        %       adaptive   whether use adaptive kernels for diffusion map 
        %
        %   Outputs:
        %       Dmap     Diffusion map (embedded coordinates)
        %       dParam   Output the parameters used in diffusion map
        %       vals     Corresponding eigenvalues
        
            if nargin<2
                t = 2;
                m = 30;
                sigma = [];
                adaptive = 1;
            end
            
            %Calcultae pairwise distances
            if any(isnan(xs))
                disp('Detected NaN in input matrix, pdist might fail...')
            end
            pdist_xs = pdist(xs);
            std_dist = std(pdist_xs);
            distMatrix = squareform(pdist_xs); %compose distance matrix
            
            if ~adaptive
                if isempty(sigma)||isnan(sigma)
                    sigma = 2*std_dist;
                elseif sigma < 0.25*std_dist
                    warning('Input sigma is too small. Rest sigma to 0.25*std_dist!')
                    sigma = 0.25*std_dist;
                    warning(['New sigma is = ' num2str(sigma)]);
                end
                disp(['Sigma is: ', num2str(sigma)])
            else
                %Use adaptive kernels
                if isempty(sigma)||isnan(sigma)
                    p = 5;
                else
                    p = ceil(sigma); %kernel size = distances to the pth neighbour
                end
                distMatrix_sorted = sort(distMatrix, 1);
                sigma = distMatrix_sorted(p, :);
                disp(['Adaptive kernel, p = ' num2str(p)])
                dParam.p = p;
            end
            
            %Plot the histogram of pdist_xs
            %figure;
            %histogram(pdist_xs)

            %Calculate weight matrix
            W = 0.5*(exp(-(distMatrix.^2)./(sigma.^2)) + exp(-(distMatrix.^2)./(sigma'.^2)));

            %Calculate weighted M matrix
            d = sum(W,1); %Degree matrix
            d_ = diag(sqrt(1./d)); %D^(-1/2)
            M = d_*W*d_; %Normalized markov matrix
            %M = W./d';
            k = m+1;

            %Get the first k eigenvectors/eigenvalues
            [psi, L] = eigs(M, k);
            psi = d_*psi;
            psi = psi./vecnorm(psi);

            %Omit the trivial first dimension
            vals = diag(L);
            Dmap = psi(:,2:end).*(vals(2:end)'.^t);
            
            %plot eigenvalues of graph Laplacian
            f = figure;
            lv = 1-vals;
            plot(lv); 
            
            %save parameters and eigenvectors
            dParam.t = t;
            dParam.m = m;
            dParam.sigma = sigma;
            dParam.adaptive = adaptive;
            dParam.psi = psi;
            dParam.vals = vals;
            
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
        
        
        function [Y, par, others] = doPHATE(X, varargin)
        % phate  Run PHATE for visualizing noisy non-linear data in lower dimensions
        % PLEASE SEE the original code at https://github.com/KrishnaswamyLab/PHATE/tree/master/Matlab
        %   Y = phate(data) runs PHATE on data (rows: samples, columns: features)
        %   with default parameter settings and returns a 2 dimensional embedding.
        %
        %   If data is sparse PCA without mean centering will be done to maintain
        %   low memory footprint. If data is dense then normal PCA (with mean
        %   centering) is done.
        %
        %   Y = phate(data, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to
        %   specify optional parameter name/value pairs that control further details
        %   of PHATE.  Parameters are:
        %
        %   'ndim' - number of (output) embedding dimensions. Common values are 2
        %   or 3. Defaults to 2.
        %
        %   'k' - number of nearest neighbors for bandwidth of adaptive alpha
        %   decaying kernel or, when a=[], number of nearest neighbors of the knn
        %   graph. For the unweighted kernel we recommend k to be a bit larger,
        %   e.g. 10 or 15. Defaults to 5.
        %
        %   'a' - alpha of alpha decaying kernel. when a=[] knn (unweighted) kernel
        %   is used. Defaults to 40.
        %
        %   't' - number of diffusion steps. Defaults to [] wich autmatically picks
        %   the optimal t.
        %
        %   't_max' - maximum t for finding optimal t. if t = [] optimal t will be
        %   computed by computing Von Neumann Entropy for each t <= t_max and then
        %   picking the kneepoint. Defaults to 100.
        %
        %   'npca' - number of pca components for computing distances. Defaults to
        %   100.
        %
        %   'mds_method' - method of multidimensional scaling. Choices are:
        %
        %       'mmds' - metric MDS (default)
        %       'cmds' - classical MDS
        %       'nmmds' - non-metric MDS
        %
        %   'distfun' - distance function. Default is 'euclidean'.
        %
        %   'distfun_mds' - distance function for MDS. Default is 'euclidean'.
        %
        %   'pot_method' - method of computing the PHATE potential dstance. Choices
        %   are:
        %
        %       'log' - -log(P + eps). (default)
        %
        %       'sqrt' - sqrt(P). (not default but often produces superior
        %       embeddings)
        %
        %       'gamma' - 2/(1-\gamma)*P^((1-\gamma)/2)
        %
        %   'gamma' - gamma value for gamma potential method. Value between -1 and
        %   1. -1 is diffusion distance. 1 is log potential. 0 is sqrt. Smaller
        %   gamma is a more locally sensitive embedding whereas larger gamma is a
        %   more globally sensitive embedding. Defaults to 0.5.
        %
        %   'pot_eps' - epsilon value added to diffusion operator prior to
        %   computing potential. Only used for 'pot_method' is 'log', i.e.:
        %   -log(P + pot_eps). Defaults to 1e-7.
        %
        %   'n_landmarks' - number of landmarks for fast and scalable PHATE. [] or
        %   n_landmarks = npoints does no landmarking, which is slower. More
        %   landmarks is more accurate but comes at the cost of speed and memory.
        %   Defaults to 2000.
        %
        %   'nsvd' - number of singular vectors for spectral clustering (for
        %   computing landmarks). Defaults to 100.
        %
        %   'kernel' - user supplied kernel. If not given ([]) kernel is
        %   computed from the supplied data. Supplied kernel should be a square
        %   (samples by samples) symmetric affinity matrix. If kernel is
        %   supplied input data can be empty ([]). Defaults to [].  
            
            %Scrape colums or rows consist of NaN
            if any(~isnan(sum(X,1)))
                X(:,isnan(X(1,:))) = [];
            else
                X(isnan(X(:,1)),:) = [];
            end
            
            %Default parameter
            npca = 100;
            k = 5;
            nsvd = 100;
            n_landmarks = 2000;
            ndim = 2;
            t = [];
            mds_method = 'mmds';
            distfun = 'euclidean';
            distfun_mds = 'euclidean';
            pot_method = 'log';
            K = [];
            a = 40;
            Pnm = [];
            t_max = 100;
            pot_eps = 1e-7;
            gamma = 0.5;

            % get input parameters
            for i=1:length(varargin)
                % k for knn adaptive sigma
                if(strcmp(varargin{i},'k'))
                   k = lower(varargin{i+1});
                end
                % a (alpha) for alpha decaying kernel
                if(strcmp(varargin{i},'a'))
                   a = lower(varargin{i+1});
                end
                % diffusion time
                if(strcmp(varargin{i},'t'))
                   t = lower(varargin{i+1});
                   if isnan(t)
                       t = [];
                   end
                end
                % t_max for VNE
                if(strcmp(varargin{i},'t_max'))
                   t_max = lower(varargin{i+1});
                end
                % Number of pca components
                if(strcmp(varargin{i},'npca'))
                   npca = lower(varargin{i+1});
                end
                % Number of dimensions for the PHATE embedding
                if(strcmp(varargin{i},'ndim'))
                   ndim = lower(varargin{i+1});
                end
                % Method for MDS
                if(strcmp(varargin{i},'mds_method'))
                   mds_method =  varargin{i+1};
                end
                % Distance function for the inputs
                if(strcmp(varargin{i},'distfun'))
                   distfun = lower(varargin{i+1});
                end
                % distfun for MDS
                if(strcmp(varargin{i},'distfun_mds'))
                   distfun_mds =  lower(varargin{i+1});
                end
                % nsvd for spectral clustering
                if(strcmp(varargin{i},'nsvd'))
                   nsvd = lower(varargin{i+1});
                end
                % n_landmarks for spectral clustering
                if(strcmp(varargin{i},'n_landmarks'))
                   n_landmarks = lower(varargin{i+1});
                end
                % potential method: log, sqrt, gamma
                if(strcmp(varargin{i},'pot_method'))
                   pot_method = lower(varargin{i+1});
                end
                % kernel
                if(strcmp(varargin{i},'kernel'))
                   K = lower(varargin{i+1});
                end
                % kernel
                if(strcmp(varargin{i},'gamma'))
                   gamma = lower(varargin{i+1});
                end
                % pot_eps
                if(strcmp(varargin{i},'pot_eps'))
                   pot_eps = lower(varargin{i+1});
                end
            end
            
            [Y, others.P, others. K] = phate(X, 'ndim', ndim, 'k', k, 'a', ...
                a, 't', t, 't_max', t_max, 'npca', npca,...
                'mds_method', mds_method, 'distfun', 'distfun', distfun,...
                'distfun_mds', distfun_mds, 'nsvd', nsvd, 'n_landmarks', ...
                n_landmarks, 'pot_method', pot_method, 'K', K, 'gamma', gamma,...
                'pot_eps', pot_eps);
            
            %Save all parameters
            par.k = k;
            par.ndim = ndim;
            par.a = a;
            par.t = t;
            par.t_max = t_max;
            par.npca = npca;
            par.mds_method = mds_method;
            par.distfun = distfun;
            par.distfun_mds = distfun_mds;
            par.nsvd = nsvd;
            par.n_landmarks = n_landmarks;
            par.pot_method = pot_method;
            par.K = K;
            par.gamma = gamma;
            par.pot_eps = pot_eps;    
            
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