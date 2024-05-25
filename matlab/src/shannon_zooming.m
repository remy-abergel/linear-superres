function out = shannon_zooming(in,W,H,varargin)
%%
% Usage: out = shannon_zooming(in,M,N,Name,Value)
%
% Input(s)/Output(s):
%
%   in     : (matrix of double) input image
%   W      : (scalar) width of the output image (horizontal zoom factor :
%            zx = M/size(in,2)) 
%   H      : (scalar) height of the output image (vertical zoom factor :
%            zy = N/size(in,1)) 
% 
%  out     : (matrix of double) output zoomed image
%
% Optional Name-Value pair arguments:
%
%   ['complex',c] : (scalar logical, default c = false), set c = true to
%                   return the complex signal (do not take the real part) 
%
% Description: image zooming using the complex variant of the Shannon
% interpolation 

%% Control number of inputs
if(nargin < 3)
    help shannon_zooming;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('in');
p.addRequired('W');
p.addRequired('H');
p.addParameter('complex',false);
parse(p,in,W,H,varargin{:});
cmode = p.Results.complex;

%% consistency checks
% input in (matrix double real numbers)
if(~isreal(in))
    help shannon_zooming;
    error('input ''in'' must be matrix of double real numbers');
end
% input W (scalar, no decimal part)
if(~isreal(W) || ~isscalar(W) || W ~= floor(W))
    help shannon_zooming;
    error('input ''W'' must be a real scalar number without decimal part (W == floor(W))');
end
% input height (scalar, no decimal part)
if(~isreal(H) || ~isscalar(H) || H ~= floor(H))
    help shannon_zooming;
    error('input ''H'' must be a real scalar number without decimal part (W == floor(H))');
end
% optional input 'complex'
if(~islogical(cmode) || ~isscalar(cmode))
    help shannon_zooming;
    error('you must set c=true or c=false for the optional Name-Value pair argument [''complex'',c]');
end

%% CORE OF THE MODULE

% retrive sizes
[h,w] = size(in); 
hmin = min(H,h);
wmin = min(W,w); 

% compute frequency coordinates
a = -floor(wmin/2) + (0:wmin-1); 
b = -floor(hmin/2) + (0:hmin-1); 

% compute output image 
dft_out = zeros(H,W); 
dft_in = fft2(in); 
dft_out(1+mod(H+b,H),1+mod(W+a,W)) = ((W*H)/(w*h)) * dft_in(1+mod(h+b,h),1+mod(w+a,w));
out = ifft2(dft_out); 
if(~cmode); out = real(out); end
    
end

