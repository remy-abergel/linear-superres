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
%   ['sigma',s] : (scalar positive double, default s = 1.0), sharpness
%                 parameter for the profile such as 10*sigma = length of
%                 the transition intervals where the profile increases from
%                 0 to 1 or decrease from 1 to 0.
%
% Description: compute low/high resolution multiplicative apodization
%              filters.
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
p.addParameter('sigma',1);
parse(p,u0,T,M,N,varargin{:});
sigma = p.Results.sigma;

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
% input sigma (scalar positive number)
if(~isreal(sigma) || ~isscalar(sigma) || sigma <= 0)
    help stack_apodization;
    error('input sigma (sharpness parameter of the apodization profile) must be a real number > 0');
end

%% CORE OF THE MODULE

% retrieve dimensions of the low-resolution sequence and compute the
% super-resolution factors (zx,zy) 
[n,m,L] = size(u0);
zx = M/m; 
zy = N/n;

% macro for the 1D profile
f = @(Q,D,t) 0.5*erfc((abs((Q-1)/2-t) - ((Q-1)/2-D-5*sigma))/(sigma*sqrt(2)));

% retrieve the maximum displacements along the horizontal & vertical
% directions 
Dx = max(abs(zx*T(:,1)));
Dy = max(abs(zy*T(:,2)));

% compute the high-resolution apodization filter
[x,y] = meshgrid(0:M-1,0:N-1);
apod_hr = f(M,Dx,x).*f(N,Dy,y);

% compute the low-resolution apodization filters
[x,y] = meshgrid(0:m-1,0:n-1);
x = zx*(x(:)+T(:,1)');
y = zy*(y(:)+T(:,2)');
apod_lr = reshape(f(M,Dx,x).*f(N,Dy,y),[n,m,L]);

% compute the apodized sequence
u0_apod = u0.*apod_lr;

end
