function [u,w,E] = irls(u0,T,M,N,varargin)
%%
% Usage: [u,w,E] = irls(u0,T,M,N,Name,Value)
%
% Input(s)/Output(s):
%
%   u0     : (hypermatrix of double) sequence of low-resolution images
%   T      : (matrix of double) translation vector, must have exactly two
%            columns, i.e., T = [tx,ty], where tx and ty denote the
%            horizontal and vertical components of the translation vectors
%   M      : (scalar >= size(u0,2)) width of the high-resolution domain
%   N      : (scalar >= size(u0,1)) height of the high-resolution domain
%
%   u      : (matrix of double) output high-resolution image (minimizer of
%            the l1-l2 energy)
%   w      : (vector of double) output weights (w(j) = 1/eta(j))
%   E      : (vector of double) l1-l2 energy evolution (E(k) = l1-l2 energy
%            at iteration k) 
%
% Optional Name-Value pair arguments:
%
%   ['verbose',v]   :  (scalar logical, default v = false), set v = true or
%                      false to enable or disable the verbose mode
%
%   ['niter',niter] : (scalar integer, default niter = 100) maximal number
%                     of iterations
% 
%   ['weights',w]   : (vector of double containing size(u0,3) elements,
%                     default weights = ones(1,1,size(u0,3)) initial
%                     weighting coefficients
%
%   ['eps',e]       : (scalar double, default e = 0) impose eta(j) >= eps
%                     (for all entries j) in the IRLS scheme
%   
%   ['stop',r]      : (scalar double, default r = 1e-5) stop iterations
%                     when the relative gap between E(k) and E(k-1) is less
%                     than r, i.e., when |E(k) - E(k-1)| <= r * E(k) 
%
% Description: iteratively reweighted least-squares procedure for the
% minimization of the l1-l2 energy
%
%      El1l2 = u -> sum_{j = 1}^{L} ||Aj(u) - u0(:,:,j)||_2
%

%% Control number of inputs
if(nargin < 4)
    help irls;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('u0');
p.addRequired('T');
p.addRequired('M');
p.addRequired('N');
p.addParameter('verbose',false);
p.addParameter('weights',ones(1,1,size(u0,3)));
p.addParameter('niter',100);
p.addParameter('eps',0);
p.addParameter('stop',1e-5);
parse(p,u0,T,M,N,varargin{:});
verbose = p.Results.verbose;
w = p.Results.weights;
niter = p.Results.niter;
e = p.Results.eps; 
r = p.Results.stop; 

%% consistency checks
% input u0 (hypermatrix of double real numbers)
if(~isreal(u0) || numel(size(u0)) ~= 3)
    help irls;
    error('input ''u0'' must be an hypermatrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help irls;
    error('input ''T'' must have exactly two columns of double real numbers');
end
if(size(T,1) ~= size(u0,3))
    help irls;
    error(['input ''T'' must have the same number of lines as the number of ' ...
           'low-resolution images in the input sequence (i.e. size(T,1) == ' ...
           'size(u0,3))']);
end
% input M (scalar >= size(u0,2), no decimal part)
if(~isreal(M) || ~isscalar(M) || M ~= floor(M) || M <= size(u0,2))
    help irls;
    error('input M must be a real scalar number, without decimal part (M == floor(M)), larger than or equal to the width of the input sequence (M >= size(u0,2))');
end
% input N (scalar >= size(u0,1), no decimal part)
if(~isreal(N) || ~isscalar(N) || N ~= floor(N) || N <= size(u0,1))
    help irls;
    error('input N must be a real scalar number, without decimal part (N == floor(N)), larger than or equal to the height of the input sequence (N >= size(u0,1))');
end
% input 'verbose'
if(~islogical(verbose) || ~isscalar(verbose))
    help irls;
    error('you must set v=true or v=false for the optional Name-Value pair argument [''verbose'',v]');
end
% input 'weights'
if(~isreal(w) || numel(w) ~= size(u0,3))
    help irls;
    error('optional input ''weights'' must contain %d (real-valued) elements',size(u0,3));
end

%% CORE OF THE MODULE
[n,m] = size(u0,[1,2]); 
E = zeros(niter,1); 
for iter = 1:niter 

    % IRLS iteration 
    u = leastsquares_superres(u0,T,M,N,'weights',w); 
    z = sqrt(sum((u0-simulator(u,T,m,n)).^2,[1,2])); 
    eta = max(e,z); 
    w = 1./eta;
    
    % compute l1-l2 energy and deal with verbose mode 
    E(iter) = sum(eta(:)); 
    if(verbose)
        fprintf("iteration %d: l1-l2 energy = %.17e\n",iter,E(iter));
    end
    
    % deal with stopping criterion
    if(iter >= 2) && (abs(1-E(iter-1)/E(iter)) <= r)
        E = E(1:iter); 
        break; 
    end
    
end

end