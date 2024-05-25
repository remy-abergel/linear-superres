function out = register_stack(in,T,varargin)
%%
% usage out = register_stack(in,T,Name,Value)
%
% Input(s)/Output(s):
%
%   in      : (hypermatrix of double) image sequence
%   T       : (matrix of double) sequence of registration vectors, must
%             have exactly two columns, i.e., T = [tx,ty], where tx and ty
%             denote the horizontal and vertical components of the
%             translation vectors
% 
%  out      : (hypermatrix of double) output registered image sequence
%
% Optional Name-Value pair arguments:
%
%   ['complex',c] : (scalar logical, default c = false), set c = true to
%                   return the complex signal (do not take the real part)
%
% Description: register a stack of images from a sequence of registration
%              displacements using the complex variant (Uc) of the Shannon
%              interpolate
%

%% Control number of inputs
if(nargin < 2)
    help register_stack;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('in');
p.addRequired('T');
p.addParameter('complex',false);
parse(p,in,T,varargin{:});
cmode = p.Results.complex;

%% consistency checks
% input in (hypermatrix of double real numbers)
if(~isreal(in) || numel(size(in)) ~= 3)
    help register_stack;
    error('input ''in'' must be an hypermatrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help register_stack;
    error('input ''T'' must have exactly two columns of double real numbers');
end
if(size(T,1) ~= size(in,3))
    help register_stack;
    error(['input ''T'' must have the same number of lines as the number of ' ...
           'low-resolution images in the input sequence (i.e. size(T,1) == ' ...
           'size(in,3))']);
end
% input 'complex'
if(~islogical(cmode) || ~isscalar(cmode))
    help register_stack;
    error('you must set c=true or c=false for the optional Name-Value pair argument [''complex'',c]');
end

%% CORE OF THE MODULE

% retrieve dimensions 
[n,m,L] = size(in); 

% compute the frequency coordinate grid 
a = ifftshift(-floor(m/2)+(0:m-1));
b = ifftshift(-floor(n/2)+(0:n-1));
[A,B] = meshgrid(a,b);
A = A(:); B = B(:);

% compute ramp-phasis terms
phi = @(a,b,dx,dy) exp(-2*1i*pi*(a*dx/m + b*dy/n));
F = reshape(phi(A,B,T(:,1)',T(:,2)'),[n,m,L]);

% register stack
out = ifft2(fft2(in) .* F); 
if(~cmode); out = real(out); end

end
