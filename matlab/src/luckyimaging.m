function out = luckyimaging(u0,T,M,N,w,n)
%%
% Usage: out = luckyimaging(u0,T,M,N,w,n)
%
% Input(s)/Output(s):
%
%   u0  : (hypermatrix of double) sequence of low-resolution images
%   T   : (matrix of double) translation vector, must have exactly two
%         columns, i.e., T = [tx,ty], where tx and ty denote the horizontal
%         and vertical components of the translation vectors 
%   M   : (scalar >= size(u0,2)) width of the high-resolution domain
%   N   : (scalar >= size(u0,1)) height of the high-resolution domain
%   w   : (vector of double, size(w) == [1,1,size(u0,3)]) IRLS weights (use
%         the irls module to compute those weights) 
%   n   : (scalar, 1 <= n < size(u0,3)) number of images from the input
%         low-resolution sequence u0 to use to compute the output high
%         resolution image (size(u0,3)-n from u0 will be ignored)
%
%   out : (matrix of double) output image (out = u_lucky^{n})
%
% Description: lucky imaging procedure for least-squares super-resolution

%% Control number of inputs
if(nargin < 6)
    help luckyimaging;
    error('Incorrect number of input(s)');
end

%% consistency checks
% input u0 (hypermatrix of double real numbers)
if(~isreal(u0) || numel(size(u0)) ~= 3)
    help luckyimaging;
    error('input ''u0'' must be an hypermatrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help luckyimaging;
    error('input ''T'' must have exactly two columns of double real numbers');
end
if(size(T,1) ~= size(u0,3))
    help luckyimaging;
    error(['input ''T'' must have the same number of lines as the number of ' ...
           'low-resolution images in the input sequence (i.e. size(T,1) == ' ...
           'size(u0,3))']);
end
% input M (scalar >= size(u0,2), no decimal part)
if(~isreal(M) || ~isscalar(M) || M ~= floor(M) || M <= size(u0,2))
    help luckyimaging;
    error('input M must be a real scalar number, without decimal part (M == floor(M)), larger than or equal to the width of the input sequence (M >= size(u0,2))');
end
% input N (scalar >= size(u0,1), no decimal part)
if(~isreal(N) || ~isscalar(N) || N ~= floor(N) || N <= size(u0,1))
    help luckyimaging;
    error('input N must be a real scalar number, without decimal part (N == floor(N)), larger than or equal to the height of the input sequence (N >= size(u0,1))');
end
% input w (vector of double, size(w) == [1,1,size(u0,3)])
if(~isreal(w) || ~all(size(w) == [1,1,size(u0,3)]))
    help luckyimaging;
    error('input w must be made of real number and such as size(w) == [1,1,size(u0,3)]');
end
% input n (scalar in [1,size(u0,3)], no decimal part)
if(~isreal(n) || ~isscalar(n) || n ~= floor(n) || n < 1 || n > size(u0,3))
    help luckyimaging;
    error('input n must be a real scalar number, with no decimal part (n == floor(n)), and such as 1 <= n <= size(u0,3)');
end

%% CORE OF THE MODULE
[~,I] = sort(w);
nsuppr = size(u0,3)-n; % number of images to remove from the input sequence 
out = leastsquares_superres(u0(:,:,I(nsuppr+1:end)),T(I(nsuppr+1:end),:),M,N);

end