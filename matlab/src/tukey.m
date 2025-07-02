function y = tukey(x, varargin)
%%
% usage y = tukey(x, r, Name, Value)
%
% Input(s)/Output(s):
%
%   x : (vector, matrix or hypermatrix of double) sampling nodes in [0,1]    
%   y : (vector, matrix or hypermatrix of double) output profile (same size
%       as x)
%
% Optional Name-Value pair arguments:
%
%   r : (scalar positive double, default = 0.05) smoothness parameter,
%       controls the speed of the transitions from 0 to 1 and from 1
%       to 0 of the computed profile (see profile definition below)
%   
%   ['d', d] : (scalar non-negative double, default = 0) rescale the
%              profile along the x-axis from [0,1] to [d, 1-d], set
%              the output profile values to zero over [0,d) and
%              (1-d,1], i.e., set
%              
%                      / 0            if 0 <= x < d
%              y(x) = <  t(a * x + b) if d <= x <= 1-d
%                      \ 0            if 1-d < x <= 1
%
%              where the rescaling parameters a=-1/(2*d-1) and b=-a*d
%              are such that a*[0,1]+b = [d,1-d], and where t denotes
%              the tukey apodization profile defined by
%
%                      / 1/2*(1-cos(2*pi*x/r))     if 0 <= x < r/2
%              t(x) = <  1                         if r/2 <= x <= 1-r/2
%                      \ 1/2*(1-cos(2*pi*(x-1)/r)) if 1-r/2 < x <= 1
%
% Description: compute (rescaled & extended by zero) Tukey apodization
%              profile.

%% Control number of inputs
if(nargin < 1)
    help tukey;
    error('Incorrect number of input(s)');
end

%% Parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('x');
p.addParameter('r', .05);
p.addParameter('d', 0);
parse(p,x,varargin{:});
r = p.Results.r;
d = p.Results.d;

%% Consistency checks

% input x (vector, matrix or hypermatrix of double real numbers)
if(~isreal(x))
    help tukey;
    error('input ''x'' must be a vector, matrix or hypermatrix of double numbers');
end

% input r (scalar positive number)
if(~isreal(r) || ~isscalar(r) || r <= 0)
    help tukey;
    error('input r (sharpness parameter of the apodization profile) must be a real number > 0');
end

% input d (scalar number in [0,1])
if(~isreal(d) || ~isscalar(d) || d < 0 || d > 1)
    help tukey;
    error('input d must be a scalar real number in [0,1]');
end

%% CORE OF THE MODULE

% original tukey function
function y = t(x, r)
    y = ones(size(x));
    idx = find(x < .5*r);
    y(idx) = .5 - .5 * cos(2*pi*x(idx)/r);
    idx = find(1 - .5*r < x);
    y(idx) = .5 - .5 * cos(2*pi*(x(idx)-1)/r);
end

% rescale and extend output
alf = -1. / (2*d - 1);
bet = - alf * d;
y = zeros(size(x));
id = find((d <= x) & (x <= 1.-d));
y(id) = t(alf * x(id) + bet, r);

end
