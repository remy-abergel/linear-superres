function [AMPL,MSE,PSNR] = error_prediction(T,m,n,M,N,varargin)
%%
% Usage: [AMPL,MSE,PSNR] = error_prediction(T,m,n,M,N,Name,Value)
% 
% Input(s)/Output(s):
%
%   T      : (matrix of double) translation vector, must have exactly two
%            columns, i.e., T = [tx,ty], where tx and ty denote the
%            horizontal and vertical components of the translation vectors
%   m      : (scalar positive) width of the low-resolution domain 
%   n      : (scalar positive) height of the low-resolution domain
%   M      : (scalar >= m) width of the high-resolution domain
%   N      : (scalar >= n) height of the high-resolution domain
% 
%   AMPL   : (matrix of double) error amplification coefficient (Fourier
%            domain)
%   MSE    : (scalar, this output is only available when the optional input
%            'sigma' is provided) predicted mean-square-error
%   PSNR   : (scalar, this output is only available when the optional input
%            'sigma' is provided) peak signal-to-noise-ratio computed as
%            PSNR = 10*log10(peakval^2/MSE)
% 
% Optional Name-Value pair arguments:
%
%   ['sigma',sigma]     : (scalar double) set the level of noise (standard
%                         deviation) corrupting the low-resolution sequence
%
%   ['peakval',peakval] : (scalar, default peakval = 255.) peak signal
%                         value used for calculating the predicted peak
%                         signal-to-noise ratio from the predicted
%                         mean-square-error 
%   ['eps',e]           : (scalar double, default e = 1e-9) threshold for
%                         the singular values in the pseudo-inversion
%                         routine (treat as zero all singular values less 
%                         than e when computing the pseudo inverse of a
%                         matrix)
%   ['noshift,nflag]    : (scalar logical, default nflag = false) by
%                         default (nflag = false) the (0,0) frequency
%                         corresponds to the center pixel of AMPL, set 
%                         nflag = true to keep the (0,0) frequency on the
%                         top-left corner of the image (FFTW convention)
% 
% Description: prediction of the super-resolution reconstruction quality
% 

%% Control number of inputs
if(nargin < 5)
    help error_prediction;
    error('Incorrect number of input(s)');
end

%% Parser (manual implementation because MATLAB parser is not case sensitive (inputs (m,M), (n,N) not allowed)

% set default values for optional arguments
sig = []; 
AMPL = []; 
peakval = 255; 
e = 1e-9; 
nflag = false; 

% retrieve optional arguments
for k = 0:floor(numel(varargin)/2)-1
    switch lower(varargin{2*k+1})
        case 'sigma'
            sig = varargin{2*k+2};
        case 'peakval'
            peakval = varargin{2*k+2}; 
        case 'eps'
            e = varargin{2*k+2}; 
        case 'noshift'
            nflag = varargin{2*k+2}; 
        otherwise
            if(isstring(varargin{2*k+1}) || ischar(varargin{2*k+1}))
                error("'%s' is not a recognized parameter. For a list of valid name-value pair arguments, see the documentation for this function.",varargin{2*k+1}); 
            else
                error("Expected a string scalar or character vector for the parameter name."); 
            end
    end
end

%% consistency checks

% check number of elements in varargin (Name-value pair arguments)
if(bitand(numel(varargin),1))
    if(isstring(varargin{end}) || ~ischar(varargin{end}))
        error("No value was given for '%s'. Name-value pair arguments require a name followed by a value.",varargin{end});
    else
        error("Expected a string scalar or character vector for the parameter name."); 
    end
end

% input T (matrix of two double real numbers)
if(~isreal(T) || size(T,2) ~= 2)
    help error_prediction;
    error('input ''T'' must have exactly two columns of double real numbers');
end

% input m (scalar positive, no decimal part)
if(~isreal(m) || ~isscalar(m) || m ~= floor(m))
    help error_prediction;
    error('input m must be a real positive scalar number without decimal part (m == floor(m), m > 0)');
end

% input n (scalar positive, no decimal part)
if(~isreal(n) || ~isscalar(n) || n ~= floor(n))
    help error_prediction;
    error('input n must be a real positive scalar number without decimal part (n == floor(n), n > 0)');
end

% input M (scalar >= m)
if(~isreal(M) || ~isscalar(M) || M ~= floor(M) || m > M)
    help error_prediction;
    error('input M must be a real scalar number, without decimal part (M == floor(M)), grater than or equal to m (M >= m)');
end

% input N (scalar, no decimal part, >= n)
if(~isreal(N) || ~isscalar(N) || N ~= floor(N) || n > N)
    help error_prediction;
    error('input N must be a real scalar number, without decimal part (N == floor(N)), greater than or equal to n (N >= n)');
end

% optional input A (matrix of double, size(A) == [N,M])
if(~isempty(AMPL) && (~isreal(AMPL) || size(AMPL,1) ~= N || size(AMPL,2) ~= M))
    help error_prediction;
    error('optional input ''AMPL'' must be a real matrix of double such as size(AMPL) == [N,M]');
end

% optional input peakval (scalar double)
if(~isscalar(peakval) || ~isreal(peakval))
    help error_prediction;
    error('optional input ''peakval'' must be a scalar double');
end

% optional input sigma (scalar double >= 0)
if(~isempty(sig) && (~isscalar(sig) || ~isreal(sig) || sig < 0))
    help error_prediction;
    error('optional input ''sigma'' must be a nonnegative scalar double');
end

% optional output MSE & PSNR (need to check that optional input 'sigma' is provided)
if((nargout >= 2) && isempty(sig))
    help error_prediction;
    error('optional input ''sigma'' is mandatory to compute the MSE and PSNR outputs');
end

% optional input 'eps'
if(~isreal(e) || ~isscalar(e))
    help error_prediction;
    error('optional input ''eps'' must be a scalar real number');
end

% optional input 'noshift'
if(~islogical(nflag) || ~isscalar(nflag))
    help error_prediction;
    error('optional input ''noshift'' must be a scalar logical');
end

%% CORE OF THE MODULE

% reduce as much as possible domain dimensions (M,m) and (N,n) while
% keeping constant the super-resolution factors zx = M/m and zy = N/n
gx = gcd(M,m); 
gy = gcd(N,n); 
M = M/gx; m = m/gx; % reduce M and m by a factor gx (does not change the ratio zx = M/m)
N = N/gy; n = n/gy; % reduce N and n by a factor gy (does not change the ratio zy = N/n)

% if not provided, compute the error amplification coefficients
if(isempty(AMPL))
    
    AMPL = zeros(N,M);
    
    % compute super-resolution factors 
    zx = M/m; 
    zy = N/n;
    
    % compute the frequency coordinate grid associated to the low-resolution domain
    a = ifftshift(-floor(m/2)+(0:m-1));
    b = ifftshift(-floor(n/2)+(0:n-1));
    [A,B] = meshgrid(a,b);
    A = A(:); B = B(:);
    
    % compute the frequency coordinate grid associated to the high-resolution domain
    alf = ifftshift(-floor(M/2)+(0:M-1));
    bet = ifftshift(-floor(N/2)+(0:N-1));
    [ALF,BET] = meshgrid(alf,bet);
    ALF = ALF(:); BET = BET(:);
    
    % compute ramp-phasis terms
    L = size(T,1); 
    phi = @(alf,bet,dx,dy) exp(-2*1i*pi*(alf*dx/M + bet*dy/N));
    F = reshape(phi(ALF,BET,zx*T(:,1)',zy*T(:,2)'),[N,M,L]);
    
    % compute bounds for the sets P(A) and Q(B)
    PMIN = ceil(-(M+2*A)/(2*m));
    PMAX = ceil((M-2*A)/(2*m))-1;
    QMIN = ceil(-(N+2*B)/(2*n));
    QMAX = ceil((N-2*B)/(2*n))-1;
    
    % precompute block matrices
    M1_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),floor(zx),floor(zy)),e);
    if(zx ~= floor(zx))
        M2_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),ceil(zx),floor(zy)),e);
    end
    if(zy ~= floor(zy))
        M3_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),floor(zx),ceil(zy)),e);
    end
    if(zx ~= floor(zx) && zy ~= floor(zy))
        M4_PINV = pinv(compute_blockmatrix(T,1/(zx*zy),ceil(zx),ceil(zy)),e);
    end
    
    % main loop
    for idfreq = 1:m*n
        
        % extract frequency (a,b) and its corresponding [pmin(a),pmax(a),qmin(b),qmax(b)]
        a = A(idfreq);
        b = B(idfreq);
        pmin = PMIN(idfreq);
        pmax = PMAX(idfreq);
        qmin = QMIN(idfreq);
        qmax = QMAX(idfreq);
        size_p = pmax - pmin + 1;
        size_q = qmax - qmin + 1;
        
        % store in lexicographical order the elements (p,q) of P(a) x Q(b)
        [P,Q] = meshgrid(pmin:pmax,qmin:qmax);
        P = P(:); Q = Q(:);
        
        % retrieve all frequencies (alf,bet) that are aliased at position (a,b)
        % in the low-resolution frequency domain
        alf = a+P*m;
        bet = b+Q*n;
        
        % compute amplification coefficients
        if(floor(zx) == size_p && floor(zy) == size_q)
            AMPL(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M1_PINV*F(1 + ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy));
        elseif(floor(zx) + 1 == size_p && floor(zy) == size_q)
            AMPL(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M2_PINV*F(1 + ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy));
        elseif(floor(zx) == size_p && floor(zy) + 1 == size_q)
            AMPL(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M3_PINV*F(1 + ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy));
        elseif(floor(zx) + 1 == size_p && floor(zy) + 1 == size_q)
            AMPL(mod(M+alf,M)*N+mod(N+bet,N)+1) = sqrt(sum(abs(M4_PINV*F(1 + ((0:L-1)*M*N + mod(M+alf,M)*N + mod(N+bet,N)))).^2,2)/(zx*zy));
        end
        
    end
end

% compute the predicted MSE and PSNR
if(nargout >= 2)
    MSE = sig^2 * sum(AMPL(:).^2) / (M*N);
end
if(nargout >= 3)
    PSNR = 10*log10(peakval^2/MSE);
end

% perform nearest neighbor zooming with factor (gx,gy) to get the
% amplification map at the initial size
AMPL = kron(fftshift(AMPL),ones(gy,gx)); 
if(nflag) 
    AMPL = ifftshift(AMPL); 
end

end

