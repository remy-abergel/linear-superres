function y = modified_tukey(x, varargin)
%%
% usage y = modified_tukey(x, r, Name, Value)
%
% Input(s)/Output(s):
%
%   x : (vector, matrix or hypermatrix of double) sampling nodes
%
%   y : (vector, matrix or hypermatrix of double) output profile (same
%       size as x)
%
% Optional Name-Value pair arguments:
%
%   ['r', r] : (scalar positive double, default = 0.025) smoothness
%              parameter
%
%   ['d', d] : (scalar non-negative double, default = 0)
%              zero-extension parameter
%
% Description: compute the modified Tukey apodization profile defined by
%
%          /
%         |             0                  if          x <= d
%         |
%         |  .5 - .5 * cos(2*pi*(x-d)/r)   if     d < x <= d+r/2
%         |
% T(x) = <              1                  if  d+r/2 < x <= 1-d-r/2
%         |
%         | .5 - .5 * cos(2*pi*(1-d-x)/r)  if   1-d-r/2 < x <= 1-d
%         |
%         |             0                  if         1-d < x
%          \
%

%% Control number of inputs
if(nargin < 1)
    help modified_tukey;
    error('Incorrect number of input(s)');
end

%% Parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('x');
p.addParameter('r', .025);
p.addParameter('d', 0);
parse(p,x,varargin{:});
r = p.Results.r;
d = p.Results.d;

%% Consistency checks

% input x (vector, matrix or hypermatrix of double real numbers)
if(~isreal(x))
    help modified_tukey;
    error('input ''x'' must be a vector, matrix or hypermatrix of double numbers');
end

% input d (scalar number in [0,1])
if(~isreal(d) || ~isscalar(d) || d < 0 || d > 1)
    help modified_tukey;
    error('input d must be a scalar real number in [0,1]');
end

% input r (scalar positive number in (0,1-2*d])
if(~isreal(r) || ~isscalar(r) || r <= 0 || r > 1-2*d)
    help modified_tukey;
    error('input r (sharpness parameter of the apodization profile) must be a real number in (0,1-2*d]');
end

%% CORE OF THE MODULE
id1 = find(x <= d);
id2 = find((x > d) & (x < d+r/2));
id3 = find((x > 1-d-r/2) & (x <= 1-d));
id4 = find(x > 1-d);
y = ones(size(x));
y(id1) = 0.;
y(id2) = .5 - .5*cos(2*pi*(x(id2)-d)/r);
y(id3) = .5 - .5*cos(2*pi*(1-d-x(id3))/r);;
y(id4) = 0.;

end
