function v = sharpening(u,f)
%%
% Usage: v = sharpening(u,f)
%
% Input(s)/Output(s):
%
%   u   : (matrix of double) input image
%   f   : (function handle or macro) a radial profile y = f(rho) 
%         defined on [0,sqrt(2)], default f = @(rho) 1+5*(1-exp(-rho))
%   v   : (matrix of double) output sharpened image
%
% Description : sharpening by frequency amplification with a given radial
% profile 
%

%% Control number of inputs
if(nargin < 1)
    help sharpening;
    error('Incorrect number of input(s)');
end

%% check consistency

% input u 
if(~ismatrix(u))
    help sharpening;
    error("input image 'u' must be a matrix of double numbers");
end

% input f (function handle)
if(nargin == 1)
    f = @(rho) 1+5*(1-exp(-rho));
else if(~isa(f,'function_handle'))
    help sharpening
    error("input 'f' must be a function handle");    
end

%% CORE OF THE MODULE
[N,M] = size(u); 
alf = ifftshift(-floor(M/2) + (0:M-1)); 
bet = ifftshift(-floor(N/2) + (0:N-1));
[ALF,BET] = meshgrid(alf/(M/2),bet/(N/2)); 
rho = hypot(ALF,BET);
v = real(ifft2(fft2(u).*f(rho))); 

end