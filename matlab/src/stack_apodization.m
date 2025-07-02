function [u0_apod,apod_hr,apod_lr] = stack_apodization(u0,T,M,N,varargin)
%%
% usage [u0_apod,apod_hr,apod_lr] = stack_apodization(u0,T,M,N,Name,Value)
%
% Input(s)/Output(s):
%
%   u0      : (hypermatrix of double) sequence of low-resolution images
%   T       : (matrix of double) translation vector, must have exactly two
%             columns, i.e., T = [tx,ty], where tx and ty denote the
%             horizontal and vertical components of the translation vectors
%   M      : (scalar >= size(u0,2)) width of the high-resolution domain
%   N      : (scalar >= size(u0,1)) height of the high-resolution domain
%
%   u0_apod : (hypermatrix of double) apodized stack of low-resolution images
%   apod_hr : (matrix of double) high-resolution apodization filter
%   apod_lr : (hypermatrix of double) low-resolution apodization filters
%             (such as u0_apod = u0 .* apod_lr)
%
% Optional Name-Value pair arguments:
%
%   ['r',r] : (scalar positive double, default r = .025) smoothness
%             parameter of the Tukey apodization profile
%
% Description: compute low/high resolution multiplicative extended
%              apodization filters.
%

%% Control number of inputs
if(nargin < 4)
    help stack_apodization;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('u0');
p.addRequired('T');
p.addRequired('M');
p.addRequired('N');
p.addParameter('r',.025);
parse(p,u0,T,M,N,varargin{:});
r = p.Results.r;

%% consistency checks
% input u0 (hypermatrix of double real numbers)
if(~isreal(u0) || numel(size(u0)) ~= 3)
    help stack_apodization;
    error('input ''u0'' must be an hypermatrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help stack_apodization;
    error('input ''T'' must have exactly two columns of double real numbers');
end
if(size(T,1) ~= size(u0,3))
    help stack_apodization;
    error(['input ''T'' must have the same number of lines as the number of ' ...
           'low-resolution images in the input sequence (i.e. size(T,1) == ' ...
           'size(u0,3))']);
end
% input M (scalar >= size(u0,2), no decimal part)
if(~isreal(M) || ~isscalar(M) || M ~= floor(M) || M <= size(u0,2))
    help leastsquares_superres;
    error('input M must be a real scalar number, without decimal part (M == floor(M)), larger than or equal to the width of the input sequence (M >= size(u0,2))');
end
% input N (scalar >= size(u0,1), no decimal part)
if(~isreal(N) || ~isscalar(N) || N ~= floor(N) || N <= size(u0,1))
    help leastsquares_superres;
    error('input N must be a real scalar number, without decimal part (N == floor(N)), larger than or equal to the height of the input sequence (N >= size(u0,1))');
end
% input r (scalar positive number)
if(~isreal(r) || ~isscalar(r) || r <= 0)
    help stack_apodization;
    error('input sigma (sharpness parameter of the apodization profile) must be a real number > 0');
end

%% CORE OF THE MODULE
    
% retrieve dimensions of the low-resolution sequence and compute the
% super-resolution factors (zx,zy) 
[n,m,L] = size(u0);
zx = M/m; 
zy = N/n;

% retrieve the maximum displacements along the horizontal & vertical
% directions 
Dx = max(abs(zx*T(:,1)));
Dy = max(abs(zy*T(:,2)));
dx = Dx/(M-1); 
dy = Dy/(N-1); 

% compute the high-resolution apodization filter
x = (0:M-1)/(M-1);
y = (0:N-1)'/(N-1);
apod_hr = tukey(x, 'r', r, 'd', dx) .* tukey(y, 'r', r, 'd', dy);

% compute the low-resolution apodization filters
x = (0:m-1);
y = (0:n-1)';
x = zx * (x + reshape(T(:, 1), [1, 1, L])) / (M-1);
y = zy * (y + reshape(T(:, 2), [1, 1, L])) / (N-1);
apod_lr = reshape(tukey(x, 'r', r, 'd', dx) .* tukey(y, 'r', r, 'd', dy), [n,m,L]);

% compute the apodized sequence
u0_apod = u0.*apod_lr;

end
