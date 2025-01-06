function [p,s] = perdecomp(u)
%%
% Usage: [p,s] = perdecomp(u)
%
% Input(s)/Output(s):
%
%   u : (matrix of double) input image
%   p : (matrix of double) output periodic component of u
%   s : (matrix of double) output smooth component of u
%
% Description: periodic + smooth decomposition of an image (u = p + s)
%
[ny,nx] = size(u); 
u = double(u);
X = 1:nx; Y = 1:ny;
v = zeros(ny,nx);
v(1,X)  = u(1,X)-u(ny,X);
v(ny,X) = -v(1,X);
v(Y,1 ) = v(Y,1 )+u(Y,1)-u(Y,nx);
v(Y,nx) = v(Y,nx)-u(Y,1)+u(Y,nx);
fx = repmat(cos(2.*pi*(X -1)/nx),ny,1);
fy = repmat(cos(2.*pi*(Y'-1)/ny),1,nx);
fx(1,1)=0.;   % avoid division by 0 in the line below
s = real(ifft2(fft2(v)*0.5./(2.-fx-fy)));
p = u-s;
end
