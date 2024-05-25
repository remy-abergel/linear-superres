function [uls,ampl] = leastsquares_superres(u0,T,M,N,varargin)
%%
% Usage: [uls,ampl] = leastsquares_superres(u0,T,M,N,Name,Value)
%
% Input(s)/Output(s):
%
%   u0     : (hypermatrix of double) sequence of low-resolution images
%   T      : (matrix of double) translation vector, must have exactly two
%            columns, i.e., T = [tx,ty], where tx and ty denote the
%            horizontal and vertical components of the translation vectors
%   M      : (scalar >= size(u0,2)) width of the high-resolution domain
%   N      : (scalar >= size(u0,1)) height of the high-resolution domain
%
%   uls    : (matrix of double) output high-resolution image
%   ampl   : (matrix of double) error amplification coefficient (Fourier
%            domain)
%
% Optional Name-Value pair arguments:
%
%   ['complex',c] : (scalar logical, default c = false), set c = true to
%                   return the complex signal (do not take the real part) 
%
%   ['weights',w] : (vector of double containing size(u0,3) elements,
%                   default weights = ones(1,1,size(u0,3)) use this
%                   optionnal input to replace each Aj operator by
%                   Aj * sqrt(weights(j)) and u0(:,:,j) by
%                   u0(:,:,j) * sqrt(weigths(j)) in the least-squares
%                   problem (useful for computing the IRLS algorithm
%                   iterations)
%
%   ['eps',e]     : (scalar double, default e = 1e-9) threshold for the
%                   singular values in the pseudo-inversion routine
%                   (treat as zero all singular values less than e when
%                   computing the pseudo inverse of a matrix)
%
% Description: Super-resolution using the least-squares estimator 
%

%% Control number of inputs
if(nargin < 4)
    help leastsquares_superres;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('u0');
p.addRequired('T');
p.addRequired('M');
p.addRequired('N');
p.addParameter('complex',false);
p.addParameter('weights',ones(1,1,size(u0,3)));
p.addParameter('eps',1e-9);
parse(p,u0,T,M,N,varargin{:});
cmode = p.Results.complex;
weights = p.Results.weights; 
e = p.Results.eps;

%% consistency checks
% input u0 (hypermatrix of double real numbers)
if(~isreal(u0) || numel(size(u0)) ~= 3)
    help leastsquares_superres;
    error('input ''u0'' must be an hypermatrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help leastsquares_superres;
    error('input ''T'' must have exactly two columns of double real numbers');
end
if(size(T,1) ~= size(u0,3))
    help leastsquares_superres;
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
% optional input 'complex'
if(~islogical(cmode) || ~isscalar(cmode))
    help leastsquares_superres;
    error('you must set c=true or c=false for the optional Name-Value pair argument [''complex'',c]');
end
% optional input 'weights'
if(~isreal(weights) || numel(weights) ~= size(u0,3))
    help leastsquares_superres;
    error('optional input ''weights'' must contain %d (real-valued) elements',size(u0,3));
end
% optional input 'eps'
if(~isreal(e) || ~isscalar(e))
    help leastsquares_superres;
    error('optional input ''eps'' must be a scalar real number');
end

%% CORE OF THE MODULE

% retrieve dimensions of the low-resolution sequence and compute the
% super-resolution factors (zx,zy) 
[n,m,L] = size(u0);
zx = M/m; 
zy = N/n;

% compute the frequency coordinate grid associated to the low-resolution
% domain 
a = ifftshift(-floor(m/2)+(0:m-1));
b = ifftshift(-floor(n/2)+(0:n-1));
[A,B] = meshgrid(a,b);
A = A(:); B = B(:);

% compute the frequency coordinate grid associated to the high-resolution
% domain 
alf = ifftshift(-floor(M/2)+(0:M-1));
bet = ifftshift(-floor(N/2)+(0:N-1));
[ALF,BET] = meshgrid(alf,bet);
ALF = ALF(:); BET = BET(:);

% compute phase ramping terms 
phi = @(alf,bet,dx,dy) exp(-2*1i*pi*(alf*dx/M + bet*dy/N));
F = reshape(phi(ALF,BET,zx*T(:,1)',zy*T(:,2)'),[N,M,L]);

% compute bounds for the sets P(A) and Q(B) 
PMIN = ceil(-(M+2*A)/(2*m)); 
PMAX = ceil((M-2*A)/(2*m))-1;
QMIN = ceil(-(N+2*B)/(2*n)); 
QMAX = ceil((N-2*B)/(2*n))-1;

% precompute block matrices
M1_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),floor(zx),floor(zy),'weights',weights),e);
if(zx ~= floor(zx))
    M2_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),ceil(zx),floor(zy),'weights',weights),e);
end
if(zy ~= floor(zy))
    M3_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),floor(zx),ceil(zy),'weights',weights),e);
end
if(zx ~= floor(zx) && zy ~= floor(zy))
    M4_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),ceil(zx),ceil(zy),'weights',weights),e); 
end
     
% compute dft_v = DFT(adjA(u0)), denoting by adjA the adjoint of the A operator)
proj1d = @(gam,size)(mod(gam,size)<(size/2)).*mod(gam,size)+(mod(gam,size)>=(size/2)).*(mod(gam,size)-size);
AA = proj1d(ALF,m); % horizontal aliased position of ALF in the low-frequency domain
BB = proj1d(BET,n); % vertical aliased position of BET in the low-frequency domain
weights = reshape(weights,[1,1,L]);
dft_tmp = fft2(u0.*weights); 
dft_v = sum(reshape(dft_tmp((0:L-1)*m*n + mod(m+AA,m)*n + mod(n+BB,n) + 1),[N,M,L]).*F,3); % DFT of adjA(u0)

% main loop 
dft_uls = zeros(N,M);
ampl = zeros(N,M); 
for idfreq = 1:m*n
    
    % extract frequency (a,b) and its corresponding
    % [pmin(a),pmax(a),qmin(b),qmax(b)] 
    a = A(idfreq);
    b = B(idfreq);
    pmin = PMIN(idfreq);
    pmax = PMAX(idfreq);
    qmin = QMIN(idfreq);
    qmax = QMAX(idfreq);
    size_p = pmax - pmin + 1;
    size_q = qmax - qmin + 1;
    Z = size_p * size_q;
    
    % store in lexicographical order the elements (p,q) of P(a) x Q(b)
    [P,Q] = meshgrid(pmin:pmax,qmin:qmax);
    P = P(:); Q = Q(:);
    
    % retrieve all frequencies (alf,bet) that are aliased at position (a,b)
    % in the low-resolution frequency domain
    alf = a+P*m;
    bet = b+Q*n;
    
    % compute right-hand term of the linear system in (32)
    right = reshape(dft_v(mod(M+alf,M)*N+mod(N+bet,N)+1),[Z,1]);
    
    % retrieve DFT coefficients of the output image at thos locations
    % (alf,bet) by means of a matrix-vector multiplication between the
    % appropriate pseudo-inverse block matrix and right-hand side term
    if(size_p == floor(zx) && size_q == floor(zy))
        dft_uls(mod(M+alf,M)*N+mod(N+bet,N)+1) = M1_PINV*right;
        ampl(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M1_PINV*F(1+ ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy)); 
    elseif(size_p == ceil(zx) && size_q == floor(zy))
        dft_uls(mod(M+alf,M)*N+mod(N+bet,N)+1) = M2_PINV*right;
        ampl(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M2_PINV*F(1+ ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy)); 
    elseif(size_p == floor(zx) && size_q == ceil(zy))
        dft_uls(mod(M+alf,M)*N+mod(N+bet,N)+1) = M3_PINV*right;
        ampl(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M3_PINV*F(1+ ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy)); 
    elseif(size_p == ceil(zx) && size_q == ceil(zy))
        dft_uls(mod(M+alf,M)*N+mod(N+bet,N)+1) = M4_PINV*right;
        ampl(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M4_PINV*F(1+ ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy)); 
    end
    
end
uls = ifft2(dft_uls);
if(~cmode); uls = real(uls); end

end
