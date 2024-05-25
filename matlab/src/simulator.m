function u0 = simulator(u,T,m,n,varargin)
%%
% Usage: u0 = simulator(u,T,m,n,Name,Value)
%
% Input(s)/Output(s):
%
%   u  : (input matrix of double) input image
%   T  : (input 2-columns matrix of double) translation vector, must have
%        exactly two components, i.e., T = [tx,ty], where tx denotes the
%        horizontal component of the translation vector, and ty the
%        vertical component of the translation vector
%   m : (positive scalar <= size(u,2)) width of the low-resolution domain
%   n : (positive scalar <= size(u,1)) height of the low-resolution domain
%   u0 : (output hypermatrix or matrix of double) output image(s)
%
% Optional Name-Value pair arguments:
%
%   ['complex',c] : (scalar logical, default c = false), set c = true to
%                   return the complex signal (do not take the real part)
%
%   ['verbose',v] : (scalar logical, default v = false), set v = true to
%                   enable verbose mode
%
% Description:
%
%  Compute a stack of shifted and subsampled images, such as
%
%  for all (k,l) in {0,...,m-1} x {0,...,n-1},
%  for all j in {1,...,L},
%
%    u0(1+l,1+k,j) = Real part of Uc(zx*(k+T(j,1)),zy*(l+T(j,2))),
%
%  denoting by Uc the "complex variant" of the Shannon interpolate of
%  'u', and denoting zx = M/m, zy = N/n, L = size(T,1).
%

%% Control number of inputs
if(nargin < 4)
    help simulator;
    error('Incorrect number of input(s)');
end

%% parser (consistency checks are done after, to allow precise error messages)
p = inputParser;
p.addRequired('u');
p.addRequired('T');
p.addRequired('m');
p.addRequired('n');
p.addParameter('verbose',false);
p.addParameter('complex',false);
parse(p,u,T,m,n,varargin{:});
verbose = p.Results.verbose;
cmode = p.Results.complex;

%% consistency checks
% input u (matrix of double real numbers)
if(~ismatrix(u) || isscalar(u) || ~isreal(u))
    help simulator;
    error('input image u must be a (non-scalar) matrix of double real numbers');
end
% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help simulator;
    error('input translation vector T must have exactly two columns of double real numbers');
end
% input m (positive scalar <= size(u,2) without decimal part)
if(~isreal(m) || ~isscalar(m) || m ~= floor(m) || m <= 0 || m > size(u,2))
    help simulator;
    error('input m must be a positive scalar number, without decimal part (m == floor(m)), smaller than or equal to width of the input image (m <= size(u,2))');
end
% input n (positive scalar <= size(u,1) without decimal part)
if(~isreal(n) || ~isscalar(n) || n ~= floor(n) || n <= 0 || n > size(u,1))
    help simulator;
    error('input n must be a positive scalar number, without decimal part (n == floor(n)), smaller than or equal to width of the input image (n <= size(u,1))');
end
% input 'verbose'
if(~islogical(verbose) || ~isscalar(verbose))
    help simulator;
    error('you must set v=true or v=false for the optional Name-Value pair argument [''verbose'',v]');
end
% input 'complex'
if(~islogical(cmode) || ~isscalar(cmode))
    help simulator;
    error('you must set c=true or c=false for the optional Name-Value pair argument [''complex'',c]');
end


%% CORE OF THE MODULE: compute the low-resolution stack of images

% retrieve dimensions and compute the super-resolution factors 
[N,M] = size(u);
L = size(T,1);
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
F = reshape(phi(ALF,BET,-zx*T(:,1)',-zy*T(:,2)'),[N,M,L]);

% compute bounds for the sets P(A) and Q(B)
PMIN = ceil(-(M+2*A)/(2*m));
PMAX = ceil((M-2*A)/(2*m))-1;
QMIN = ceil(-(N+2*B)/(2*n));
QMAX = ceil((N-2*B)/(2*n))-1;

% compute the output low-resolution sequence in the Fourier domain
dft_u_shifted = fft2(u) .* F;
dft_u0 = zeros(n,m,L);

% deal with frequencies related to the first block
size_p = floor(zx); size_q = floor(zy);
ID = find((PMAX-PMIN+1 == size_p) & (QMAX-QMIN+1 == size_q));
AA = A(ID);
BB = B(ID);
ADR = (0:L-1)*m*n + mod(m+AA,m)*n+mod(n+BB,n) + 1;
PPMIN = PMIN(ID);
QQMIN = QMIN(ID);
for idk = 0:size_p-1
    ALF = AA + m*(PPMIN+idk);
    for idl = 0:size_q-1
        BET = BB + n*(QQMIN+idl);
        adr = (0:L-1)*M*N + mod(M+ALF,M)*N + mod(N+BET,N) + 1;
        dft_u0(ADR) = dft_u0(ADR) + dft_u_shifted(adr);
    end
end

% deal with frequencies related to the second block
if(zx ~= floor(zx))
    size_p = ceil(zx); size_q = floor(zy);
    ID = find((PMAX-PMIN+1 == size_p) & (QMAX-QMIN+1 == size_q));
    AA = A(ID);
    BB = B(ID);
    ADR = (0:L-1)*m*n + mod(m+AA,m)*n+mod(n+BB,n) + 1;
    PPMIN = PMIN(ID);
    QQMIN = QMIN(ID);
    for idk = 0:size_p-1
        ALF = AA + m*(PPMIN+idk);
        for idl = 0:size_q-1
            BET = BB + n*(QQMIN+idl);
            adr = (0:L-1)*M*N + mod(M+ALF,M)*N + mod(N+BET,N) + 1;
            dft_u0(ADR) = dft_u0(ADR) + dft_u_shifted(adr);
        end
    end
end

% deal with frequencies related to the third block
if(zy ~= floor(zy))
    size_p = floor(zx); size_q = ceil(zy);
    ID = find((PMAX-PMIN+1 == size_p) & (QMAX-QMIN+1 == size_q));
    AA = A(ID);
    BB = B(ID);
    ADR = (0:L-1)*m*n + mod(m+AA,m)*n+mod(n+BB,n) + 1;
    PPMIN = PMIN(ID);
    QQMIN = QMIN(ID);
    for idk = 0:size_p-1
        ALF = AA + m*(PPMIN+idk);
        for idl = 0:size_q-1
            BET = BB + n*(QQMIN+idl);
            adr = (0:L-1)*M*N + mod(M+ALF,M)*N + mod(N+BET,N) + 1;
            dft_u0(ADR) = dft_u0(ADR) + dft_u_shifted(adr);
        end
    end
end

% deal with frequencies related to the fourth block
if(zx ~= floor(zx) && zy ~= floor(zy))
    size_p = ceil(zx); size_q = ceil(zy);
    ID = find((PMAX-PMIN+1 == size_p) & (QMAX-QMIN+1 == size_q));
    AA = A(ID);
    BB = B(ID);
    ADR = (0:L-1)*m*n + mod(m+AA,m)*n+mod(n+BB,n) + 1;
    PPMIN = PMIN(ID);
    QQMIN = QMIN(ID);
    for idk = 0:size_p-1
        ALF = AA + m*(PPMIN+idk);
        for idl = 0:size_q-1
            BET = BB + n*(QQMIN+idl);
            adr = (0:L-1)*M*N + mod(M+ALF,M)*N + mod(N+BET,N) + 1;
            dft_u0(ADR) = dft_u0(ADR) + dft_u_shifted(adr);
        end
    end
end

% perform renormalization and go back to the spatial domain
u0 = ifft2((n*m)/(M*N) * dft_u0);
if(~cmode); u0 = real(u0); end

end
