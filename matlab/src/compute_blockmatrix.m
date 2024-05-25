function M = compute_blockmatrix(T,gam,size_p,size_q,varargin)
%%
% Usage: [M] = compute_blockmatrix(T,gam,size_p,size_q,Name,Value)
%
% Input(s)
%
%   T      : (input 2-columns matrix of double) sequence of translation
%            vectors, T = [tx,ty], where tx denotes the horizontal
%            component of the translation vector, and ty the vertical
%            component of the translation vector
%   gam    : scalar double
%   size_p : scalar double (with integer value)
%   size_q : scalar double (with integer value)
%
% Optional Name-Value pair arguments:
%
%   ['weights',w] : (vector of double containing size(T,1) elements,
%                   default weights = ones(1,1,size(T,1)) computed weighted
%                   bloc matrix (see details below)
% 
% Output:
%
% a square matrix M of complex double, with size Z x Z (denoting
% Z=size_p*size_q), and such as, the entry of M at line k (1 <= k <= Z) and
% column l (1 <= l <= Z) is given by
%
%  M(k,l) = gam * sum_{id=1}^{L} weights(id) * exp(2*1i*pi*(mux*T(id,1) + muy*T(id,2))),
%
% denoting
%
% L = size(T,1)
% mux = floor((l-1)/size_q) - floor((k-1)/size_q)
% muy = mod((l-1),size_q) - mod((k-1),size_q)
%

%% Control number of inputs
if(nargin < 4)
    help compute_blockmatrix;
    error('Incorrect number of input(s)');
end

%% Parser (no error check)
p = inputParser;
p.addRequired('T');
p.addRequired('gam');
p.addRequired('size_p');
p.addRequired('size_q');
p.addParameter('weights',ones(1,1,size(T,1)));
parse(p,T,gam,size_p,size_q,varargin{:});
weights = p.Results.weights;

%% CORE OF THE MODULE: compute block matrix
Z = size_p*size_q;
[l,k] = meshgrid(1:Z,1:Z); % k = row index, l = column index
mux = floor((l-1)/size_q) - floor((k-1)/size_q);
muy = mod(l-1,size_q) - mod(k-1,size_q);
M = gam * sum(reshape(weights,[1,1,size(T,1)]) .* reshape(exp(2*1i*pi*(mux(:)*T(:,1)' + muy(:)*T(:,2)')),[Z,Z,size(T,1)]),3);
end
