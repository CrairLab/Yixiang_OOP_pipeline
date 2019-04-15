function [Dmap, vals] = diffmap(xs, t, m, sigma)
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
%       vals     Corresponding eigenvalues

    pdist_xs = pdist(xs);
    
    if isempty(sigma)
        sigma = 0.25*std(pdist_xs);
    end
    disp(['Sigma is: ', num2str(sigma)])
    
    distMatrix = squareform(pdist_xs);
    W = exp(-distMatrix./(2*sigma^2));
    
    d = sum(W,1);
    M = W./d';
    k = m+1;
    
    [psi, L] = eigs(M, k);
    
    vals = diag(L);
    vals = vals(2:end);
    Dmap = psi(:,2:end).*vals'.^t;
    
end
