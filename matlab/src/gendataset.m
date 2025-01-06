function [u0,ref] = gendataset(in,T,zx,zy,varargin)
%%
% Usage: [u0,ref] = gendataset(in,T,zx,zy,Name,Value)
%
% Input(s)/Output(s):
%
%   in : (input matrix of double) input image
%   T  : (input 2-columns matrix of double) translation vector, must have
%        exactly two components, i.e., T = [tx,ty], where tx denotes the
%        horizontal component of the translation vector, and ty the
%        vertical component of the translation vector
%   zx : (scalar >= 1) subsampling factor along the horizontal direction
%   zy : (scalar >= 1) subsampling factor along the vertical direction
%   u0 : (output hypermatrix or matrix of double) output low-resolution
%        image(s) 
%   ref: (output matrix of double) output high-resolution reference image
%        (corresponds to a cropping of in)
%
% Optional Name-Value pair arguments:
%
%   ['complex',c] : (scalar logical, default c = false), set c = true to
%                   return the complex signal (do not take the real part)
%
%   ['cx',cx]     : (scalar with integer value, default cx = 10), maximal
%                   number of columns to be discarded from the input image
%                   during the first cropping step
%
%   ['cy',cy]     : (scalar with integer value, default cy = 10), maximal
%                   number of row to be discarded from the input image
%                   during the first cropping step
%
%   ['p',p]       : (scalar with integer value, default p = 5) he final
%                   high-resolution ROI size will be equal to (q/p) * the
%                   dimensions of the cropped HR image
%
%   ['q',q]       : (scalar with integer value, default q = 4), see above
%
% Description: Compute a realistic low-resolution stack from a high
%              resolution image (avoiding periodic-like boundaries).

%% Control number of inputs
if(nargin < 4)
    help gendataset;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
P = inputParser;
P.addRequired('in');
P.addRequired('T');
P.addRequired('zx');
P.addRequired('zy');
P.addParameter('cx',10);
P.addParameter('cy',10);
P.addParameter('p',5);
P.addParameter('q',4);
P.addParameter('complex',false);
parse(P,in,T,zx,zy,varargin{:});
cx = P.Results.cx;
cy = P.Results.cy;
p = P.Results.p;
q = P.Results.q;
cmode = P.Results.complex;

%% consistency checks
% input in (matrix of double real numbers)
if(~ismatrix(in) || isscalar(in) || ~isreal(in))
    help gendataset;
    error('input image in must be a (non-scalar) matrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help gendataset;
    error('input translation vector T must have exactly two columns of double real numbers');
end
% input zx (positive scalar >= 1)
if(~isreal(zx) || ~isscalar(zx) || zx < 1)
    help gendataset;
    error('input zx must be a scalar number >= 1');
end
% input zy (positive scalar >= 1)
if(~isreal(zy) || ~isscalar(zy) || zy < 1)
    help gendataset;
    error('input zy must be a scalar number >= 1');
end
% input cx (scalar in [0,size(in,2)] without decimal part)
if(~isreal(cx) || ~isscalar(cx) || cx ~= floor(cx) || cx < 0 || cx > size(in,2))
    help gendataset;
    error('input cx must be a scalar number in [0,size(in,2)] without decimal part');
end
% input cy (scalar in [0,size(in,1)] without decimal part)
if(~isreal(cy) || ~isscalar(cy) || cy ~= floor(cy) || cy < 0 || cy > size(in,1))
    help gendataset;
    error('input cy must be a scalar number in [0,size(in,1)] without decimal part');
end
% input p (scalar >= 2 without decimal part)
if(~isreal(p) || ~isscalar(p) || p ~= floor(p) || p < 2)
    help gendataset;
    error('input p must be a scalar number >= 1 without decimal part');
end
% input q (positive scalar in [1,p-1], without decimal part)
if(~isreal(q) || ~isscalar(q) || q ~= floor(q) || q < 1 || q >= p)
    help gendataset;
    error('input q must be a scalar number in [1, p-1] without decimal part');
end
% input 'complex'
if(~islogical(cmode) || ~isscalar(cmode))
    help gendataset;
    error('you must set c=true or c=false for the optional Name-Value pair argument [''complex'',c]');
end

%% CORE OF THE MODULE: compute a realistic low-resolution stack

% first cropping step along the X & Y axes (we will remove up to cx columns
% and cy rows from the input HR image in order fulfill as best as possible
% the zx and zy requirements)
zx_tmp = zx; zy_tmp = zy; 
[N0, M0] = size(in);
Ax = floor((M0 - (0:cx)) / p);
Ay = floor((N0 - (0:cy)) / p);
ax = round(Ax/zx_tmp);
ay = round(Ay/zy_tmp);
Cx = (zx_tmp - Ax./ax).^2;
Cy = (zy_tmp - Ay./ay).^2;
idx = find(Cx == min(Cx), 1);
idy = find(Cy == min(Cy), 1);
Ax = Ax(idx);
ax = ax(idx);
Ay = Ay(idy);
ay = ay(idy);
zx = Ax/ax;
zy = Ay/ay;
if(zx_tmp ~= zx)
    warning("zx has been changed into %.17g (instead of zx=%g) to allow integer width for the output image",zx,zx_tmp);
end
if(zy_tmp ~= zy) 
    warning("zy has been changed into %.17g (instead of zy=%g) to allow integer height for the output image",zy,zy_tmp);
end
M0prime = Ax * p;
N0prime = Ay * p;
X0 = floor((M0 - M0prime)/2);
Y0 = floor((N0 - N0prime)/2);
ref = in(Y0+(1:N0prime), X0+(1:M0prime));

% perform the periodic + smooth decomposition of the cropped HR image
[P,S] = perdecomp(ref);

% generate (unrealistic) LR stack: use the simulator for the periodic
% component and use bilinear interpolation for the smooth component
m = q*ax;
n = q*ay;
M = q*Ax;
N = q*Ay;
X0 = floor((M0prime - M)/2);
Y0 = floor((N0prime - N)/2);
x0 = X0/zx;
y0 = Y0/zy;
k0 = floor(x0);
l0 = floor(y0);
dx = x0 - k0;
dy = y0 - l0;
T2 = T+ [dx, dy];
m0 = ax*p;
n0 = ay*p;
p0 = simulator(P, T2, m0, n0);
[x, y] = meshgrid(0:m0-1, 0:n0-1);
s0 = zeros(n0,m0,size(T,1));
for id = 1:size(T,1)
    s0(:,:,id) = interp2(S, 1 + zx*(x+T2(id,1)), 1 + zy*(y+T(id,2)), 'linear', 0);
end
u0 = p0 + s0;

% crop the generated stack to get rid of unrealistic boundaries and compute
% the associated HR reference image
ref = ref(Y0+(1:N), X0+(1:M));
u0 = u0(l0+(1:n), k0+(1:m), :);
if(~cmode); u0 = real(u0); end

end