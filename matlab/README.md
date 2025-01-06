# Linear Super-Resolution through Translational Motion (Matlab package)
  
## Installation of the MATLAB modules

Simply add the [src](src) directory of this MATLAB pachage to your
MATLAB path. All provided modules will be immediately available from
your MATLAB console (**no compilation/mex-interfacing needed**).

From the MATLAB console, use the `help` command to display the
documentation related to each provided module, for instance typing

```matlab
help leastsquares_superres
```

should result in the following message in your MATLAB console:
  
```
Usage: [uls,ampl] = leastsquares_superres(u0,T,M,N,Name,Value)

Input(s)/Output(s):

  u0     : (hypermatrix of double) sequence of low-resolution images
  T      : (matrix of double) translation vector, must have exactly two
           columns, i.e., T = [tx,ty], where tx and ty denote the
           horizontal and vertical components of the translation vectors
  M      : (scalar >= size(u0,2)) width of the high-resolution domain
  N      : (scalar >= size(u0,1)) height of the high-resolution domain

  uls    : (matrix of double) output high-resolution image
  ampl   : (matrix of double) error amplification coefficient (Fourier
           domain)
 
Optional Name-Value pair arguments:
 
  ['complex',c] : (scalar logical, default c = false), set c = true to
                  return the complex signal (do not take the real part) 
 
  ['weights',w] : (vector of double containing size(u0,3) elements,
                  default weights = ones(1,1,size(u0,3)) use this
                  optionnal input to replace each Aj operator by
                  Aj * sqrt(weights(j)) and u0(:,:,j) by
                  u0(:,:,j) * sqrt(weigths(j)) in the least-squares
                  problem (useful for computing the IRLS algorithm
                  iterations)
 
  ['eps',e]     : (scalar double, default e = 1e-9) threshold for the
                  singular values in the pseudo-inversion routine
                  (treat as zero all singular values less than e when
                  computing the pseudo inverse of a matrix)

Description: Super-resolution using the least-squares estimator
```
Now, you may want to jump to [practical
examples](#examples-reproduce-several-experiments-of-the-companion-article)
(and reproduce some experiments presented the companion article) or to
have a closer look to the [modules contained into this package](#modules-description).

## Modules description

  | SOURCE FILE                                            | CONTENT / PURPOSE                                                                                                                                      | RELATIONS WITH THE COMPANION ARTICLE                             |
  |--------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------|
  | [simulator.m](src/simulator.m)                         | Compute a stack of (shifted and subsampled) low-resolution images from an input high resolution image                                                  | implements the operator A defined in Equation (12)               |
  | [stack_apodization.m](src/stack_apodization.m)         | Apodization of a low-resolution sequence (used to avoid reconstruction artifacts when processing real-life sequences)                                  | implements the apodization procedure described in Section 2.3    |
  | [compute_blockmatrix.m](src/compute_blockmatrix.m)     | Computation of a block matrix                                                                                                                          | implements the pseudocode Algorithm 1 described in Section 3.4   |
  | [leastsquares_superres.m](src/leastsquares_superres.m) | Super-resolution using the least-squares estimator                                                                                                     | implements the pseudocode Algorithm 2 described in Section 3.4   |
  | [irls.m](src/irls.m)                                   | Iteratively Reweighted Least-Squares                                                                                                                   | implements the IRLS procedure described in Section 6             |
|                                                        |                                                                                                                                                        |                                                                  |
  | [luckyimaging.m](src/luckyimaging.m)                   | Lucky-imaging procedure for least-squares image super-resolution                                                                                       | implements the lucky-imaging procedure described in Section 6    |
  | [sharpening.m](src/sharpening.m)                       | Image sharpening using a frequency amplification filter                                                                                                | implements the image sharpening procedure described in Section 7 |
  | [imview.m](src/imview.m)                               | Image displayer (display an image into a MATLAB Figure with tight borders AND without interpolation: one pixel of the screen = one pixel of the image) | None                                                             |
  | [mview.m](src/mview.m)                                 | Frame-by-frame movie displayer                                                                                                                         | None                                                             |
  | [gendataset.m](src/gendataset.m)                       | synthesizing realistic low-resolution sequences (without periodic-like boundaries) from an input high-resolution image                                 | None                                                             |

  **Additional note**

  The modules [`compute_blockmatrix`](src/compute_blockmatrix.m) and
  [`leastsquares_superres`](src/leastsquares_superres.m) come with an
  optional `weight` input that slighlty generalizes the corresponding
  pseudocodes algorithms described in the companion article (Algorithm
  1 and Algorithm 2). This optional input is used inside the
  [`irls`](src/irls.m) module (see the explanations in Section 6 of
  the companion article).

## Examples (reproduce several experiments of the companion article)
### Simulate realistic low-resolution sequences

We illustrate a procedure for synthesizing realistic sequences of
low-resolution images from a given high-resolution image and a
sequence of displacements.

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:
   
```matlab
% load and display a high-resolution image (Brooklyn bridge)
ref_large = double(imread('bridge.tif')); 
figure('Name','large reference high-resolution image'); 
imview(ref_large);

% generate a sequence containing 10 random shifts in [-5,5] x [-5,5]
T = -5 + 10*rand(10,2);

% first step: set width and height of the intermediate unrealistic
% low-resolution stack (this sequence will exhibit prediodic-like
% unrealistic boundaries)
m = size(ref_large,2) / 3; % width of the low-resolution domain
n = size(ref_large,1) / 3; % height of the low-resolution domain
zx = size(ref_large,2) / m; % (= 3) subsampling factor along the X-axis
zy = size(ref_large,1) / n; % (= 3) subsampling factor along the Y-axis

% second step: use the simulator to compute a low-resolution stack
% (the unrealistic periodical behavior of the simulator can be
% observed near the boundaries of the low-resolution images)
u0_nonrealistic = simulator(ref_large,T,m,n);
figure(); mview(u0_nonrealistic);

% third step: compute a realistic sequence by cropping the
% non-realistic sequence generated by the 'simulator' module (the crop
% is used to get rid of the unrealistic periodic-like boundaries)
x0 = 81; y0 = 36; % coordinates of the top-left corner of the low-resolution cropping domain
m = 150; n = 108; % width and height of the cropped low-resolution domain
u0_realistic = u0_nonrealistic(y0+(1:n),x0+(1:m),:);

% crop to the high-resolution image in the same area (generate the
% high-resolution reference (ground-truth) image associated to the
% simulated realistic low-resolution sequence)
X0 = zx*x0; Y0 = zy*y0; % coordinates of the top-left corner of the high-resolution cropping domain
M = zx*m; % width of the (cropped) high-resolution domain
N = zy*n; % height of the (cropped) high-resolution domain
ref_crop = ref_large(Y0+(1:N),X0+(1:M));

% display the realistic low-resolution sequence and the corresponding
% high-resolution image
figure('Name','Reference high-resolution image'); imview(ref_crop);
figure(); mview(u0_realistic);
```

Notice that, in the above example, the cropping area of the
low-resolution image yields corners with integer coordinates in the
corresponding high-resolution domain, allowing the extraction of the
reference high-resolution image using a simple crop.

Synthesizing some realistic datasets from a reference image with
different dimensions and using arbitrary subsampling factors
(especially noninteger) can be trickier using the methodology
described above. More generic synthesis of realistic datasets can be
carried out using the gendataset module that added in v1.0.2 (note
that this module was not used in the experiments presented in the
companion research article.

To use this module, you can run the following MATLAB commands:

```matlab
% load a high-resolution image (Brooklyn bridge)
ref_large = double(imread('bridge.tif')); 

% generate a sequence containing 10 random shifts in [-5,5] x [-5,5]
T = -5 + 10*rand(10,2);

% compute a realistic dataset made of a sequence of low-resolution images 
% and the corresponding high-resolution ground-truth from an input 
% high-resolution image (not that the actual subsampling factors may be 
% slightly different from the requested ones, actual values will be printed 
% with 17 digits of precision on the standard output)
zx = 2.3; % requested subsampling factor along the horizontal axis
zy = 2.8; % requested subsampling factor along the vertical axis
[u0_realistic, ref] = gendataset(ref_large,T,zx,zy);

% retrieve actual subsampling factors
[N,M] = size(ref);
[n,m] = size(u0_realistic,[1,2]);
zx = M/m;
zy = N/n;

% display the generated realistic low-resolution sequence and the corresponding
% high-resolution image
figure('Name', 'Reference high-resolution image'); imview(ref); 
figure(); mview(u0_realistic); 
```

### Using apodization to avoid boundary artifacts (reproduce Figure 1 of the companion article)

This experiments illustrates the importance of properly dealing with
periodization artifacts induced by the periodicity of the Shannon
interpolate when computing a super-resolved image using the
leastsquare-superres module.

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:

```matlab
% load a the realistic low-resolution sequence and the corresponding
% sequence of displacements used in Figure 1 (you can also generate
% your own realistic sequence as shown in the previous section)
load('bridge_sequence_fig1.mat','u0_realistic','T');
figure(); mview(u0_realistic);

% load the corresponding reference (ground-truth) image
load('bridge_sequence_fig1.mat','ref_crop');
figure(); imview(ref_crop);
[N,M] = size(ref_crop);

% perform super-resolution over the realistic sequence of
% low-resolution images: we observe important artefacts in the
% super-resolved image
uls1 = leastsquares_superres(u0_realistic,T,M,N);
figure('Name','reconstruction without apodization'); imview(uls1,'black',0,'white',255); % display (clipped) output image

% apodize the sequence using the 'stack_apodization' module
[u0_apod,apod_hr] = stack_apodization(u0_realistic,T,M,N);
figure(); mview(u0_apod);

% perform super-resolution over the apodized sequence (the
% super-resolved image does not exhibits artifact anymore and is very
% close to the apodized reference image
uls2 = leastsquares_superres(u0_apod,T,M,N);
figure('Name','reconstruction with apodization'); imview(uls2,'black',0,'white',255);
figure('Name','apodized reference image'); imview(ref_crop .* apod_hr,'black',0,'white',255);
```

### Super-resolution using the least-squares (reproduce Figure 7 of the companion article)
   
This experiments intents to perform super-resolution with factor 2
along both dimensions (zx=zy=2) from a sequence containing L = 20
noisy low-resolution images, with no perturbation on the
displacements.
   
Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:
   
```matlab
% compute a realistic sequence containing L=20 low-resolution images
% corrupted by an additive Gaussian noise (stdev = 2)
ref_large = double(imread('bridge.tif'));
T = -5 + 10*rand(20,2);
m = size(ref_large,2) / 2;
n = size(ref_large,1) / 2;
u0_nonrealistic = simulator(ref_large,T,m,n);
u0_realistic = u0_nonrealistic(54+(1:162),121+(1:225),:);
u0_noisy = u0_realistic + 2*randn(size(u0_realistic));
figure('Name','low-resolution sequence'); mview(u0_noisy);

% set the width and height of the high-resolution domain (we set here
% the super-resolution factors equal to 2 in both directions: zx=zy=2)
[n,m] = size(u0_noisy,[1,2]);
M = 2*m; % width of the high-resolution domain
N = 2*n; % height of the high-resolution domain

% apodize this sequence to avoid boundary artifacts in further
% processings (see previous section for more details)
[u0_apod,apod_hr] = stack_apodization(u0_noisy,T,M,N);

% display the spectrum (modulus of the DFT) in logarithmic scale of
% the first low-resolution image (remark the aliasing in the
% frequency domain)
figure('Name','spectrum of a low-resolution image'); 
imview(log(1+abs(fftshift(fft2(u0_apod(:,:,1))))),'black',15,'white',2); % black = large values, white = low values

% perform super-resolution, display the result and its spectrum
uls = leastsquares_superres(u0_apod,T,M,N);
figure('Name','least-squares reconstruction'); imview(uls);
figure('Name','spectrum of the least-squares reconstruction'); imview(log(1+abs(fftshift(fft2(uls)))),'black',15,'white',2);

% compare the reconstruction to the apodized reference image (compute
% PSNR & MSE metrics)
ref = ref_large(108+(1:324),242+(1:450)) .* apod_hr; % cropped and apodized reference image
fprintf("MSE = %.2e, PSNR = %2.3g dB\n",mean((ref(:)-uls(:)).^2),10*log10(255^2/mean((ref(:)-uls(:)).^2)));
```
   
### Prediction of the reconstruction quality (reproduce Figure 4, 5, 6 of the companion article)

This experiment illustrates how the reconstruction quality provided by
the least-square estimator can be blindly (i.e. without reference
image) and efficiently predicted. It also illustrates the
non-stationarity of the reconstruction error in the Fourier domain.

**Remark**: in this experiment, to simplify the image synthesis
process, we generate a (non-realistic) sequence using the `simulator`
module and then we apodize it using the `stack-apodization`
module. This provides a realistic apodized sequence, in the sense that
the apodized sequence obtained in this way is visually
indistinguishable from the apodized sequence that we would have
obtained if we had generated a realistic sequence and apodized it
afterwards (because, by construction, the apodization filter vanishes
where the periodic-like artifact occur in the non-realistic simulated
sequence).

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:
   
```matlab
% compute an apodized low-resolution sequence (L = 10, zx = 2.5, zy = 2.4)
ref_large = double(imread('bridge.tif'));
ref_crop = ref_large(108+(1:324),242+(1:450));
[N,M] = size(ref_crop); % dimensions of the high-resolution domain
m = round(M/2.5); % width of the low-resolution domain 
n = round(N/2.4); % height of the low-resolution domain
zx = M/m; % super-resolution factor along the X-axis (horizontal)
zy = N/n; % super-resolution factor along the Y-axis (vertical)     
L = 10; % number of low-resolution images (you can increase L to see the effect on the reconstruction quality)
T = -5 + 10*rand(L,2); % random displacement sequence
u0_nonrealistic = simulator(ref_crop,T,m,n);
sig = 2; % noise level (standard deviation)
u0_noisy = u0_nonrealistic + sig*randn(size(u0_nonrealistic));
sig_apod = 1; % sigma parameter used to design the apodization filters
[u0_apod,apod_hr] = stack_apodization(u0_noisy,T,M,N,'sigma',sig_apod); % generates a realistic apodized low-resolution sequence (from a non-realistic low-resolution sequence)
figure('Name','apodized low-resolution sequence'); mview(u0_apod);

% perform prediction of the least-square super-resolution
% reconstruction quality (MSE & PSNR)
[ampl,mse_pred,psnr_pred] = error_prediction(T,m,n,M,N,'sigma',sig,'peakval',255);
fprintf("predicted reconstruction quality: MSE = %.3e, PSNR = %.3g dB.\n",mse_pred,psnr_pred);

% perform super-resolution reconstruction and compute observed MSE &
% PSNR
uls = leastsquares_superres(u0_apod,T,M,N);
ref_apod = ref_crop .* apod_hr; % apodized ground-truth
mse_observed = mean((uls(:)-ref_apod(:)).^2);
psnr_observed = 10*log10(255^2/mse_observed);
fprintf("observed reconstruction quality: MSE = %.3e, PSNR = %.3g dB.\n",mse_observed,psnr_observed);

% we remark that the observed MSE & PSNR are a bit better than the
% predictions, this is due to apodization which is not taken into
% account in our theoretical study of the error. Indeed, the
% reconstruction error is very small in the areas affected by
% apodization (black borders).
figure('Name','reconstruction absolute error (spatial domain)'); imview(abs(ref_apod - uls));

% removing the black borders (that do not correspond to useful signal)
% from the apodized high-resolution images, the observed MSE and PSNR
% are even closer to their predicted values
delta_x = ceil((M/m)*max(abs(T(:,1))) + 10*sig_apod);
delta_y = ceil((N/n)*max(abs(T(:,2))) + 10*sig_apod);
ref_apod_noborder = ref_apod((1+delta_y):(N-delta_y),(1+delta_x):(M-delta_x));
uls_noborder = uls((1+delta_y):(N-delta_y),(1+delta_x):(M-delta_x));
mse_observed2 = mean((uls_noborder(:)-ref_apod_noborder(:)).^2);
psnr_observed2 = 10*log10(255^2/mse_observed2);
fprintf("observed reconstruction quality (remove black borders): MSE = %.3e, PSNR = %.3g dB.\n",mse_observed2,psnr_observed2);

% display the observed (normalized) reconstruction error in the
% Fourier domain and compare to its expectation
figure('Name','expected reconstruction error (Fourier domain)'); imview(sig*ampl);
figure('Name','observed reconstruction error (Fourier domain)'); imview(abs(fftshift(fft2(uls-ref_apod)))/sqrt(M*N));

% use mview to flip those two error maps (press 'f' to flip images)
mov = sig*ampl; % expected reconstruction error (Fourier domain, center = (0,0) frequency)
mov(:,:,2) = abs(fftshift(fft2(uls-ref_apod)))/sqrt(M*N); % observed reconstruction error (Fourier domain)
figure(); mview(mov);
```

In the following script, we focus on the reconstruction error in both
Fourier & spatial domains, we reproduce (with larger image domains)
experiments similar to that described in Figure 4 and Figure 5.
   
```matlab
% set size of the low and high resolution domains
m = 200; % width of the low-resolution domain
n = 100; % height of the low-resolution domain
M = round(m*2.3); % width of the high-resolution domain 
N = round(n*2.2); % height of the high-resolution domain

% compute reconstruction error (= least-squares operator applied to a
% pure noise sequence) for several values of L
rng(31415); % for repeatability (ensures nice display), you can comment this line to generate different realizations
T9 = -5 + 10*rand(9,2); 
[eps_prime_L9,ampl_L9] = leastsquares_superres(randn(n,m,9),T9,M,N); % reconstruction error for L = 9

T14 = -5 + 10*rand(14,2); 
[eps_prime_L14,ampl_L14] = leastsquares_superres(randn(n,m,14),T14,M,N); % reconstruction error for L = 14

T20 = -5 + 10*rand(20,2); 
[eps_prime_L20,ampl_L20] = leastsquares_superres(randn(n,m,20),T20,M,N); % reconstruction error for L = 20

% display the (normalized) reconstruction error in the Fourier domain
% and the associated map of error amplification coefficients

% reproduce a similar result as that displayed in the first column of Figure 4 (L = 9) 
black_graylevel = min(ampl_L20(:)); 
white_graylevel = max(ampl_L9(:)); 
figure('name','(L = 9) reconstruction error (Fourier domain)'); imview(fftshift(abs(fft2(eps_prime_L9))/sqrt(M*N)),'black',black_graylevel,'white',white_graylevel,'colormap',jet(256)); 
figure('name','(L = 9) error amplification coefficients'); imview(fftshift(ampl_L9),'black',black_graylevel,'white',white_graylevel,'colormap',jet(256));

% reproduce a similar result as that displayed in the second column of Figure 4 (L = 14)
figure('name','(L = 14) reconstruction error (Fourier domain)'); imview(fftshift(abs(fft2(eps_prime_L14))/sqrt(M*N)),'black',black_graylevel,'white',white_graylevel,'colormap',jet(256)); 
figure('name','(L = 14) error amplification coefficients'); imview(fftshift(ampl_L14),'black',black_graylevel,'white',white_graylevel,'colormap',jet(256));

% reproduce a similar result as that displayed in the last column of Figure 4 (L = 20)
figure('name','(L = 20) realization of sqrt(|eps''|^2/(MN))'); imview(fftshift(abs(fft2(eps_prime_L20))/sqrt(M*N)),'black',black_graylevel,'white',white_graylevel,'colormap',jet(256)); 
figure('name','(L = 20) error amplification coefficients'); imview(fftshift(ampl_L20),'black',black_graylevel,'white',white_graylevel,'colormap',jet(256));

% display the reconstruction error in the spatial domain (reproduce 
% a similar result as that displayed in the first row of Figure 5)
figure('name','(L=9) reconstruction error (spatial domain)'); imview(eps_prime_L9,'black',-10,'white',10); 
figure('name','(L=14) reconstruction error (spatial domain)'); imview(eps_prime_L14,'black',-10,'white',10); 
figure('name','(L=20) reconstruction error (spatial domain)'); imview(eps_prime_L20,'black',-10,'white',10); 

% compute and display the empirical standard deviation of the
% reconstruction error in the spatial domain (this simulation takes
% several minutes) to reproduce a similar result as that displayed 
% in the last row of Figure 5
Nsimu = 100; 
eps_prime_L9 = zeros(N,M,Nsimu); 
eps_prime_L14 = zeros(N,M,Nsimu); 
eps_prime_L20 = zeros(N,M,Nsimu); 
for k = 1:Nsimu 
	
	% compute a realization of the reconstruction error for L = 9
	out = leastsquares_superres(randn(n,m,9),T9,M,N); 
	eps_prime_L9(:,:,k) = out; 
	
	% compute a realization of the reconstruction error for L = 14
	out = leastsquares_superres(randn(n,m,14),T14,M,N); 
	eps_prime_L14(:,:,k) = out; 
	
	% compute a realization of the reconstruction error for L = 20
	out = leastsquares_superres(randn(n,m,20),T20,M,N); 
	eps_prime_L20(:,:,k) = out; 
	
	% display progression
	fprintf("simulation %d/%d: done\n",k,Nsimu);
	
end
emp_std_L9 = std(eps_prime_L9,[],3); 
emp_std_L14 = std(eps_prime_L14,[],3); 
emp_std_L20 = std(eps_prime_L20,[],3); 

figure('Name','(L = 9) empirical standard deviation of the reconstruction error (spatial domain)'); imview(emp_std_L9,'black',0,'white',4,'colormap',jet(256)); 
figure('Name','(L = 14) empirical standard deviation of the reconstruction error (spatial domain)'); imview(emp_std_L14,'black',0,'white',4,'colormap',jet(256)); 
figure('Name','(L = 20) empirical standard deviation of the reconstruction error (spatial domain)'); imview(emp_std_L20,'black',0,'white',4,'colormap',jet(256)); 
```

Now, let us compare the accuracy of the PSNR prediction over several
random experiments (we reproduce an experiment similar to that
described in Figure 6, this experiment takes several minutes).

```matlab
% load reference image
ref_large = double(imread('bridge.tif')); 
ref_crop = ref_large(108+(1:324),242+(1:450));
[N,M] = size(ref_crop); % dimensions of the high-resolution domain

% reproduce Figure 6 (a)
m = M/2; % width of the low-resolution domain (zx = M/m = 2) 
n = N/2; % height of the low-resolution domain (zy = N/n = 2)
Nsimu = 50; % number of simulations per tested value of L (in Figure 6 (a) we used Nsimu = 100)
L_list = [4,5,6,8,10,12,14]; % values of L to be tested
sig_noise = 2; % noise level in the low-resolution sequence
sig_apod = 1; % apodization parameter
PSNR_PRED = zeros(Nsimu,numel(L_list));
PSNR_OBSERVED = zeros(Nsimu,numel(L_list));
for idL = 1:numel(L_list)
	L = L_list(idL); 
	for k = 1:Nsimu 
	
	% generate a random low-resolution stack 
	T = -5 + 10*randn(L,2); 
	[u0_apod,apod_hr] = stack_apodization(simulator(ref_crop,T,m,n) + sig_noise*randn(n,m,L),T,M,N,'sigma',sig_apod); 
	
	% compute reference image associated to the apodized sequence (and
	% remove black borders due to apodization)
	ref_apod = ref_crop .* apod_hr; 
	delta_x = ceil((M/m)*max(abs(T(:,1))) + 10*sig_apod); 
	delta_y = ceil((N/n)*max(abs(T(:,2))) + 10*sig_apod); 
	ref_apod_noborder = ref_apod((1+delta_y):(N-delta_y),(1+delta_x):(M-delta_x));
	
	% predict the PSNR of the least-squares super-resolution reconstruction 
	[~,~,psnr_pred] = error_prediction(T,m,n,M,N,'sigma',sig_noise); 
	
	% perform super-resolution reconstruction (and remove black borders due to apodization)
	uls = leastsquares_superres(u0_apod,T,M,N); 
	uls_noborder = uls((1+delta_y):(N-delta_y),(1+delta_x):(M-delta_x));
	
	% compute the observed PSNR between the reconstruction and the reference
	% image (without black borders)
	psnr_observed = 10*log10(255^2/mean((uls_noborder(:)-ref_apod_noborder(:)).^2));
	
	% store computed values 
	PSNR_PRED(k,idL) = psnr_pred; 
	PSNR_OBSERVED(k,idL) = psnr_observed; 
	fprintf("+ L = %d, simu %d/%d: predicted PSNR = %.3g dB, observed PSNR = %.3g dB\n",L,k,Nsimu,psnr_pred,psnr_observed);
	
	end 
end
     
fg = figure(); hold on; 
leg = {}; 
for idL = 1:numel(L_list)
	scatter(PSNR_PRED(:,idL),PSNR_OBSERVED(:,idL),'LineWidth',1,'SizeData',80); 
	leg{end+1} = sprintf("L = %d",L_list(idL));
end
a = gca(); 
plot(a.XLim,a.XLim,'g--','LineWidth',1); 
grid on; a.GridLineStyle = '--';
leg{end+1} = "y = x";
legend(leg,'Location','southeast'); 
ylabel("observed PSNR (dB)"); 
xlabel("predicted PSNR (dB)");
title(sprintf('zx = %g, zy = %g',M/m,N/n)); 

% reproduce Figure 6 (b)
m = round(M/2.1); % width of the low-resolution domain (zx = M/m close to 2.1) 
n = round(N/2.3); % height of the low-resolution domain (zy = N/n close to 2.3)
Nsimu = 50; % number of simulations per tested value of L (in Figure 6 (b) we used Nsimu = 100)
L_list = [9,10,11,13,15,17,19]; % values of L to be tested
sig_noise = 2; % noise level in the low-resolution sequence
sig_apod = 1; % apodization parameter
PSNR_PRED = zeros(Nsimu,numel(L_list));
PSNR_OBSERVED = zeros(Nsimu,numel(L_list));
for idL = 1:numel(L_list)
	L = L_list(idL); 
	for k = 1:Nsimu 
	
	% generate a random low-resolution stack 
	T = -5 + 10*randn(L,2); 
	[u0_apod,apod_hr] = stack_apodization(simulator(ref_crop,T,m,n) + sig_noise*randn(n,m,L),T,M,N,'sigma',sig_apod); 
	
	% compute reference image associated to the apodized sequence (and
	% remove black borders due to apodization)
	ref_apod = ref_crop .* apod_hr; 
	delta_x = ceil((M/m)*max(abs(T(:,1))) + 10*sig_apod); 
	delta_y = ceil((N/n)*max(abs(T(:,2))) + 10*sig_apod); 
	ref_apod_noborder = ref_apod((1+delta_y):(N-delta_y),(1+delta_x):(M-delta_x));
	
	% predict the PSNR of the least-squares super-resolution reconstruction 
	[~,~,psnr_pred] = error_prediction(T,m,n,M,N,'sigma',sig_noise); 
	
	% perform super-resolution reconstruction (and remove black borders due to apodization)
	uls = leastsquares_superres(u0_apod,T,M,N); 
	uls_noborder = uls((1+delta_y):(N-delta_y),(1+delta_x):(M-delta_x));
	
	% compute the observed PSNR between the reconstruction and the reference
	% image (without black borders)
	psnr_observed = 10*log10(255^2/mean((uls_noborder(:)-ref_apod_noborder(:)).^2));
	
	% store computed values 
	PSNR_PRED(k,idL) = psnr_pred; 
	PSNR_OBSERVED(k,idL) = psnr_observed; 
	fprintf("+ L = %d, simu %d/%d: predicted PSNR = %.3g dB, observed PSNR = %.3g dB\n",L,k,Nsimu,psnr_pred,psnr_observed);
	
	end 
end

fg = figure(); hold on; 
leg = {}; 
for idL = 1:numel(L_list)
	scatter(PSNR_PRED(:,idL),PSNR_OBSERVED(:,idL),'LineWidth',1,'SizeData',80); 
	leg{end+1} = sprintf("L = %d",L_list(idL));
end
a = gca(); 
plot(a.XLim,a.XLim,'g--','LineWidth',1); 
grid on; a.GridLineStyle = '--';
leg{end+1} = "y = x";
legend(leg,'Location','southeast'); 
ylabel("observed PSNR (dB)"); 
xlabel("predicted PSNR (dB)");
title(sprintf('zx = %g, zy = %g',M/m,N/n)); 
```

### Influence of the displacement configuration over the quality of the reconstruction (reproduce Figure 9 of the companion article)
   
In this experiment we reproduce some results similar to that displayed
in the first row of Figure 9, showing how the displacement
configuration may affect the quality reconstruction (especially when L
is close to zx*zy).

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:
   
```matlab
% compute three sequences containing 4 low-resolution images using the
% same displacement sequences as that used to compute the images
% displayed in the first row of Figure 9
ref_large = double(imread('bridge.tif')); 
m = size(ref_large,2) / 2; 
n = size(ref_large,1) / 2; 
T1 = dlmread('shifts_fig9a.txt');
T2 = dlmread('shifts_fig9b.txt');
T3 = dlmread('shifts_fig9c.txt');
u0_nonrealistic1 = simulator(ref_large,T1,m,n);
u0_nonrealistic2 = simulator(ref_large,T2,m,n);
u0_nonrealistic3 = simulator(ref_large,T3,m,n);
u0_realistic1 = u0_nonrealistic1(54+(1:162),121+(1:225),:);
u0_realistic2 = u0_nonrealistic2(54+(1:162),121+(1:225),:);
u0_realistic3 = u0_nonrealistic3(54+(1:162),121+(1:225),:);
sig = 2; % noise level (standard deviation)
u0_noisy1 = u0_realistic1 + sig*randn(size(u0_realistic1)); 
u0_noisy2 = u0_realistic2 + sig*randn(size(u0_realistic2)); 
u0_noisy3 = u0_realistic3 + sig*randn(size(u0_realistic3)); 

% set the width and height of the high-resolution domain (we set here
% the super-resolution factors equal to 2 in both directions: zx=zy=2)
[n,m] = size(u0_noisy1,[1,2]); % update the size of the low-resolution domain (after cropping)
M = 2*m; % width of the high-resolution domain
N = 2*n; % height of the high-resolution domain
[u0_apod1,apod_hr1] = stack_apodization(u0_noisy1,T1,M,N); 
[u0_apod2,apod_hr2] = stack_apodization(u0_noisy2,T2,M,N); 
[u0_apod3,apod_hr3] = stack_apodization(u0_noisy3,T3,M,N); 

% perform super-resolution reconstruction from each sequence of
% low-resolution images and compute the predicted PSNR for each
% reconstruction
uls1 = leastsquares_superres(u0_apod1,T1,M,N);
uls2 = leastsquares_superres(u0_apod2,T2,M,N);
uls3 = leastsquares_superres(u0_apod3,T3,M,N);
[~,~,psnr1] = error_prediction(T1,m,n,M,N,'sigma',sig,'peakval',255);
[~,~,psnr2] = error_prediction(T2,m,n,M,N,'sigma',sig,'peakval',255);
[~,~,psnr3] = error_prediction(T3,m,n,M,N,'sigma',sig,'peakval',255);
figure('Name',sprintf('PSNR = %.3g dB (first decile)',psnr1)); imview(uls1,'black',0,'white',255); 
figure('Name',sprintf('PSNR = %.3g dB (median)',psnr2)); imview(uls2,'black',0,'white',255); 
figure('Name',sprintf('PSNR = %.3g dB (last decile)',psnr3)); imview(uls3,'black',0,'white',255); 
```
### Least-squares reconstruction using erroneous displacements (reproduce Figure 10 of the companion article)
   
In this experiment, we perform least-squares reconstruction from
inexact sequences of displacements, which corresponds to a similar
experiment as that displayed in Figure 10.
   
Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:

```matlab
% compute a realistic sequence containing L=20 low-resolution images
% corrupted by an additive Gaussian noise (stdev = 2)
ref_large = double(imread('bridge.tif'));
T_ref = -5 + 10*rand(20,2);
m = size(ref_large,2) / 2; 
n = size(ref_large,1) / 2; 
u0_nonrealistic = simulator(ref_large,T_ref,m,n);
u0_realistic = u0_nonrealistic(54+(1:162),121+(1:225),:);
sig_noise = 2; % noise level (standard deviation)
u0_noisy = u0_realistic + sig_noise*randn(size(u0_realistic)); 
figure('Name','low-resolution sequence'); mview(u0_noisy);

% generate noisy displacement sequences (to model estimation error)
T1 = T_ref + 0.01*randn(size(T_ref)); % small amount of noise
T2 = T_ref + 0.05*randn(size(T_ref)); 
T3 = T_ref + 0.10*randn(size(T_ref)); 
T4 = T_ref + 0.20*randn(size(T_ref)); % large amount of noise

% set the width and height of the high-resolution domain (we set here
% the super-resolution factors equal to 2 in both directions: zx=zy=2)
[n,m] = size(u0_noisy,[1,2]); % dimensions of the low-resolution domain
M = 2*m; % width of the high-resolution domain
N = 2*n; % height of the high-resolution domain

% perform reconstructions from those noisy displacement sequences
% (apodization & least-squares super-resolution)
u0_apod1 = stack_apodization(u0_noisy,T1,M,N); 
u0_apod2 = stack_apodization(u0_noisy,T2,M,N); 
u0_apod3 = stack_apodization(u0_noisy,T3,M,N); 
u0_apod4 = stack_apodization(u0_noisy,T4,M,N); 

uls1 = leastsquares_superres(u0_apod1,T1,M,N); 
uls2 = leastsquares_superres(u0_apod2,T2,M,N); 
uls3 = leastsquares_superres(u0_apod3,T3,M,N); 
uls4 = leastsquares_superres(u0_apod4,T4,M,N); 

% display images (spatial domain)
figure('Name','reconstruction (sigma_delta = 0.01)'); imview(uls1,'black',0,'white',255); 
figure('Name','reconstruction (sigma_delta = 0.05)'); imview(uls2,'black',0,'white',255); 
figure('Name','reconstruction (sigma_delta = 0.1)'); imview(uls3,'black',0,'white',255); 
figure('Name','reconstruction (sigma_delta = 0.2)'); imview(uls4,'black',0,'white',255); 

% display details (close-up views x5)
X0=324; Y0=84; X1=408; Y1=144;
figure('Name','close-up view x5 (sigma_delta = 0.01)'); imview(uls1(1+(Y0:Y1),1+(X0:X1)),'scale',5); 
figure('Name','close-up view x5 (sigma_delta = 0.05)'); imview(uls2(1+(Y0:Y1),1+(X0:X1)),'scale',5); 
figure('Name','close-up view x5 (sigma_delta = 0.1)'); imview(uls3(1+(Y0:Y1),1+(X0:X1)),'scale',5); 
figure('Name','close-up view x5 (sigma_delta = 0.2)'); imview(uls4(1+(Y0:Y1),1+(X0:X1)),'scale',5); 

% perform shannon resamplings (zoom x5)
uls1_zoom = shannon_zooming(uls1,5*M,5*N); 
uls2_zoom = shannon_zooming(uls2,5*M,5*N); 
uls3_zoom = shannon_zooming(uls3,5*M,5*N); 
uls4_zoom = shannon_zooming(uls4,5*M,5*N); 
figure('Name','Shannon zooming x5 (sigma_delta = 0.01)'); imview(uls1_zoom(5*Y0+(1:5*(Y1-Y0+1)),5*X0+(1:5*(X1-X0+1)))); 
figure('Name','Shannon zooming x5 (sigma_delta = 0.05)'); imview(uls2_zoom(5*Y0+(1:5*(Y1-Y0+1)),5*X0+(1:5*(X1-X0+1)))); 
figure('Name','Shannon zooming x5 (sigma_delta = 0.1)'); imview(uls3_zoom(5*Y0+(1:5*(Y1-Y0+1)),5*X0+(1:5*(X1-X0+1)))); 
figure('Name','Shannon zooming x5 (sigma_delta = 0.2)'); imview(uls4_zoom(5*Y0+(1:5*(Y1-Y0+1)),5*X0+(1:5*(X1-X0+1)))); 

% display Fourier spectra
figure('Name','Fourier spectrum (sigma_delta = 0.01)'); imview(log(1+abs(fftshift(fft2(uls1)))),'black',16,'white',1); 
figure('Name','Fourier spectrum (sigma_delta = 0.05)'); imview(log(1+abs(fftshift(fft2(uls2)))),'black',16,'white',1); 
figure('Name','Fourier spectrum (sigma_delta = 0.1)'); imview(log(1+abs(fftshift(fft2(uls3)))),'black',16,'white',1); 
figure('Name','Fourier spectrum (sigma_delta = 0.2)'); imview(log(1+abs(fftshift(fft2(uls4)))),'black',16,'white',1); 
```

### Least-squares reconstruction over real thermal infrared data (reproduce Figure 12 of the companion article)

This experiments intents to perform super-resolution with factor 1.8
along both dimensions (zx=zy=1.8) from a real-life sequence containing
L = 500 thermal infrared images that we acquired by ourselves (the
associated sequence of displacements was estimated from the
low-resolution sequence using algorithm of [Keren et
al.](https://doi.org/10.1109/CVPR.1988.196317).
   
Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:

```matlab
% load the thermal infrared sequence and the estimated sequence of
% displacements
load('flirT640.mat','flirT640');
T = dlmread('shifts_flirT640.txt');

% set width and height of the reconstruction (high-resolution) domain
[n,m,L] = size(flirT640);
M = 1.8*m; % width of the high-resolution domain (super-resolution factor zx = 1.8)
N = 1.8*n; % height of the high-resolution domain (super-resolution factor zy = 1.8)

% normalize & apodize the FLIR sequence to avoid edge effects in the
% reconstruction
u0 = double(flirT640);
u0 = 255. * (u0 - min(u0(:))) / (max(u0(:))-min(u0(:)));
u0_apod = stack_apodization(u0,T,M,N);

% compute the shift-and-median (i.e., the median of the registered
% low-resolution sequence)
u0_registered = register_stack(u0_apod,T);
out_shiftandmedian = median(u0_registered,3);

% perform least-squares super-resolution
uls = leastsquares_superres(u0_apod,T,M,N);

% display and compare those three images (reproduce the first row of Figure 12)
figure('Name','low-resolution image u0_apod(:,:,1)'); imview(u0_apod(:,:,1),'black',20,'white',140,'scale',2);
figure('Name','shift-and-median'); imview(out_shiftandmedian,'black',20,'white',140);
figure('Name','least-squares reconstruction uls'); imview(uls,'black',20,'white',140);

% perform bicubic resamplings (zoom x3, reproduce the second row of Figure 12)
zoom_factor = 3; 
[Xq1,Yq1] = meshgrid(1+(0:round(zoom_factor*m-1))/zoom_factor,1+(0:round(zoom_factor*n-1))/zoom_factor);
[Xq2,Yq2] = meshgrid(1+(0:round(zoom_factor*M-1))/zoom_factor,1+(0:round(zoom_factor*N-1))/zoom_factor);
figure('Name',sprintf('first image of u0_apod (bicubic zoom %g%%)',100*zoom_factor)); imview(interp2(u0_apod(:,:,1),Xq1,Yq1,'bicubic'),'black',20,'white',140);
figure('Name',sprintf('shift-and-median (bicubic zoom %g%%)',100*zoom_factor)); imview(interp2(out_shiftandmedian(:,:,1),Xq1,Yq1,'bicubic'),'black',20,'white',140);
figure('Name',sprintf('least-squares reconstruction uls (bicubic zoom %g%%)',100*zoom_factor)); imview(interp2(uls,Xq2,Yq2,'bicubic'),'black',20,'white',140);

% display and compare the DFT spectra (reproduce the last row of Figure 12)
figure('Name','spectrum of u0_apod(:,:,1)'); imview(log(1+abs(fftshift(fft2(u0_apod(:,:,1))))),'black',16,'white',0);
figure('Name','spectrum of out_shiftandmedian'); imview(log(1+abs(fftshift(fft2(out_shiftandmedian)))),'black',16,'white',0);
figure('Name','spectrum of uls'); imview(log(1+abs(fftshift(fft2(uls)))),'black',16,'white',0);
```

### IRLS & Lucky imaging (reproduce Figure 14 of the companion article)
   
We reproduce here a similar experiment to that proposed in Figure 14
of the companion article (again, the result will formally slightly
differ due to the random generation of the experiment
parameters). This experiment illustrates how the lucky-imaging
procedure can be used to remove outliers from the initial sequence in
order to improve the reconstruction quality.

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:

```matlab
% compute a synthetic sequence containing L=20 low-resolution images
% corrupted by an additive Gaussian noise (stdev = 2), the sequence
% is not apodized yet
ref_large = double(imread('bridge.tif'));
m = size(ref_large,2) / 2;
n = size(ref_large,1) / 2;
L = 20; % number of low-resolution images 
T_ref = -5 + 10*rand(L,2);
u0_nonrealistic = simulator(ref_large,T_ref,m,n);
u0_realistic = u0_nonrealistic(54+(1:162),121+(1:225),:);
sig = 2; % noise level (standard deviation)
u0_noisy = u0_realistic + sig*randn(size(u0_realistic));

% compute perturbated displacement sequence by adding a small amounts
% of Gaussian noise (stdev=0.01) to 14 of those displacements and a
% large amount of Gaussian noise (stdev=3) to the remaining 6
% displacements
id = randperm(L);
T = T_ref;
T(id(1:L-6),:) = T(id(1:L-6),:) + 0.01*randn(L-6,2); % add small perturbations
T(id(L-5:L),:) = T(id(L-5:L),:) + 3*randn(6,2); % add large perturbations

% apodize & process the low-resolution sequence using the perturbated
% sequence
[n,m] = size(u0_noisy,[1,2]);
M = 2*m; N = 2*n;
u0_apod = stack_apodization(u0_noisy,T,M,N);
uls = leastsquares_superres(u0_apod,T,M,N);

% display reconstruction (reproduce first row of Figure 14, display the whole 
% images instead of cropping areas)
figure('Name','least-squares reconstruction'); imview(uls,'black',0,'white',255);
figure('Name','least-squares reconstruction (Shannon zoom x3)'); imview(shannon_zooming(uls,3*M,3*N),'black',0,'white',255);
figure('Name','least-squares reconstruction (Fourier spectrum)'); imview(log(1+abs(fftshift(fft2(uls)))),'black',16,'white',1);

% use the IRLS Algorithm to compute the minimizer of the l1-l2 residual
[ul1l2,weights,El1l2] = irls(u0_apod,T,M,N,'verbose',true);
figure(); plot(El1l2); xlabel('IRLS iterations'); ylabel('l1-l2 energy'); title('Energy decrease');

% display l1-l2 energy minimizer (reproduce second row of Figure 14, display 
% the whole images instead of cropping areas)
figure('Name','l1-l2 reconstruction'); imview(ul1l2,'black',0,'white',255);
figure('Name','l1-l2 reconstruction (Shannon zoom x3)'); imview(shannon_zooming(ul1l2,3*M,3*N),'black',0,'white',255);
figure('Name','l1-l2 reconstruction (Fourier spectrum)'); imview(log(1+abs(fftshift(fft2(ul1l2)))),'black',16,'white',1);

% display IRLS output weights, outliers should be associated to
% significantly low weights
figure(); plot(weights(:),'*-'); hold on; a = gca(); a.XLim = [0,L+1]; a.XTick = 1:L;
for k = id(L-5:L)
	plot([k,k],a.YLim,'r--');
end
legend('weights 1/\eta_j','outliers');
xlabel('index j'); ylabel('weights');

% Lucky-imaging: eliminate N=6 images (associated with the N-smallest
% weights) and perform least-squares reconstruction
[~,I] = sort(weights);
ulucky = leastsquares_superres(u0_apod(:,:,I(7:end)),T(I(7:end),:),M,N); % u_lucky^{14} : image obtained using L-6 = 14 images from the input low-resolution sequence

% Equivalently, you can use the luckyimaging module (generates exactly
% the same output)
ulucky = luckyimaging(u0_apod,T,M,N,weights,14); % u_lucky^{14} : image obtained using L-6 = 14 images from the input low-resolution sequence

% display lucky-imaging result (reproduce last row of Figure 14, display 
% the whole images instead of cropping areas)
figure('Name','lucky reconstruction'); imview(ulucky,'black',0,'white',255);
figure('Name','lucky reconstruction (Shannon zoom x3)'); imview(shannon_zooming(ulucky,3*M,3*N),'black',0,'white',255);
figure('Name','lucky reconstruction (Fourier spectrum)'); imview(log(1+abs(fftshift(fft2(ulucky)))),'black',16,'white',1);

% when the number of outliers is unknown, we can compute a movie, by
% removing 0, 1, 2, ... images from the sequence in the same order
% than the IRLS weights, play the movie, and select the most
% satisfying result (notice that reconstruction artifacts are
% particularly visible in the Fourier domain)
ulucky_mov = zeros(N,M,L-3);
for Nsuppr = 0:L-4
	Nkeep = L-Nsuppr;
	ulucky_mov(:,:,Nsuppr+1) = luckyimaging(u0_apod,T,M,N,weights,Nkeep);
	fprintf("lucky-imaging movie (frame #%d): used %d images from the initial sequence\n",Nsuppr+1,Nkeep);
end
figure(); mview(ulucky_mov,'black',0,'white',255); % display movie in the spatial domain
figure(); mview(log(1+abs(fftshift(fftshift(fft2(ulucky_mov),1),2))),'black',16,'white',1); % display movie in the Fourier domain
```

### Least-squares reconstruction over real microscopy data (reproduce Figure 11 and 15 of the companion article)

The dataset used in the following experiments corresponds to an
in-vivo recording of Purkinje cells of a living rat, which was
acquired with a 2-photons microscope. This sequence was kindly
provided to us by [Jorge Enrique
Ramírez-Buriticá](https://scholar.google.com/citations?user=cQFpBPQAAAAJ&hl=es)
and [Brandon
Stell](https://scholar.google.com/citations?user=WLggr90AAAAJ&hl=es),
who performed this challenging video acquisition in order to study the
temporal spiking activity of in-vivo brain cells. Achieving a
satisfactory temporal sampling (FPS) rate for this sequence
constrained the spatial sampling rate to be set large (4 micrometer
per pixel), leading to strongly aliased images.

The original sequence contains 3921 images with size 45 x 75 each. We
register this sequence using the [Keren et
al.](https://doi.org/10.1109/CVPR.1988.196317) algorithm. We removed
35 images from the original sequence which were associated with to
large displacement (>20 micrometer) in order to keep a large-enough
field of view for the super-resolved image to reconstruct. This leads
to the sequence of 3186 images shared [here](../data/2photons.tif).

In the following, we partially reproduce the experiment presented in
Figure 11 and Figure 15 of the companion article.

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:

```matlab
% load the 2-photons sequence and the estimated sequence of
% displacements
u0 = zeros(75, 45, 3186); 
for id = 1:3186
	u0(:,:,id) = double(imread('2photons.tif', id));
end
T = dlmread('shifts_2photons.txt');

% display the low-resolution sequence: we use the 'scale' 
% parameter of the mview module to upscale the sequence 
% (using nearest-neighbor zoom) with a factor 10 to improve 
% the display
z = 10;
mview(u0, 'scale', z); % set z=1 to disable upscaling

% set width and height of the reconstruction (high-resolution) domain
[n,m,L] = size(u0);
M = 2*m; % width of the high-resolution domain (super-resolution factor zx = 2)
N = 2*n; % height of the high-resolution domain (super-resolution factor zy = 2)

% apodize the FLIR sequence to avoid edge effects in the 
% reconstruction
u0_apod = stack_apodization(u0,T,M,N,'sigma',0.5);

% compute the shift-and-add (i.e., the temporal of the registered
% low-resolution sequence)
u0_registered = register_stack(u0_apod,T);
out_shiftandadd = mean(u0_registered,3);

% perform least-squares super-resolution using the whole sequence
uls = leastsquares_superres(u0_apod,T,M,N);

% perform least-squares super-resolution using the lucky-imaging 
% strategy
[ul1l2,weights,El1l2] = irls(u0_apod,T,M,N,'verbose',true);
ulucky = luckyimaging(u0_apod,T,M,N,weights,L-300); % u_lucky^{2886} : image obtained using L-300 = 2886 images from the input low-resolution sequence

% compare one low-resolution image (the closest to out_shiftandadd 
% in l2 distance) of the apodized sequence, the shift-and-add image,
% and the super-resolved image obtained using the lucky-imaging 
% strategy (reproduce first row of Figure 11 of the companion 
% article)
zz = reshape(sum((u0_apod - repmat(out_shiftandadd,[1,1,L])).^2,[1,2]),[L,1]);
id = find(zz == min(zz));
figure('Name',sprintf('low-resolution image u0_apod(:,:,%d)',id)); imview(u0_apod(:,:,id),'black',0,'white',140);
figure('Name',sprintf('shift-and-add image',id)); imview(out_shiftandadd,'black',0,'white',140);
figure('Name',sprintf('super-resolved image (lucky)',id)); imview(ulucky,'black',0,'white',140);

% same using upscaling to improve the display (display the low-resolution 
% images using factor 10 upscaling and the super-resolved using a 
% factor 5 upscaling so that the displayed image have the same size (but 
% not the same size))
figure('Name',sprintf('low-resolution image u0_apod(:,:,%d) (x10 scale)',id)); imview(u0_apod(:,:,id),'black',0,'white',140,'scale',10);
figure('Name',sprintf('shift-and-add image (x10 scale)',id)); imview(out_shiftandadd,'black',0,'white',140,'scale',10);
figure('Name',sprintf('super-resolved image (lucky) (x5 scale)',id)); imview(ulucky,'black',0,'white',140,'scale',5);

% compare the super-resolved images reconstructed with of without 
% lucky imaging strategy (reproduce the two first columns of 
% Figure 15 of the companion article)
figure('Name',sprintf('super-resolved image (no lucky)',id)); imview(uls,'black',0,'white',140);
figure('Name',sprintf('super-resolved image (lucky)',id)); imview(ulucky,'black',0,'white',140);

% same using x5 upscaling for both high resolution images 
figure('Name',sprintf('super-resolved image (no lucky)',id)); imview(uls,'black',0,'white',140,'scale',5);
figure('Name',sprintf('super-resolved image (lucky)',id)); imview(ulucky,'black',0,'white',140,'scale',5);

```

### Super-resolution and deconvolution (reproduce Figure 17 of the companion article)

This experiments illustrates the benefit of applying frequency
enhancement procedure after the least-squares super-resolution
reconstruction process over a real data sequence (the FLIR T640
thermal infrared image sequence). This experiments reproduces Figure
17 of the companion article.

Make sure you added the [`src`](src) and [`data`](../data) directories
of this package to your MATLAB path and run the following MATLAB
commands:

```matlab
% load the thermal infrared sequence and the estimated sequence of
% displacements
load('flirT640.mat','flirT640'); 
T = dlmread('shifts_flirT640.txt'); 

% set width and height of the reconstruction (high-resolution) domain
[n,m,L] = size(flirT640); 
M = 1.8*m; % width of the high-resolution domain (super-resolution factor zx = 1.8)
N = 1.8*n; % height of the high-resolution domain (super-resolution factor zy = 1.8)

% normalize & apodize the FLIR sequence to avoid edge effects in the
% reconstruction
u0 = double(flirT640); 
u0 = 255. * (u0 - min(u0(:))) / (max(u0(:))-min(u0(:)));
u0_apod = stack_apodization(u0,T,M,N);

% compute the shift-and-median (i.e., the median of the registered
% low-resolution sequence), apply bicubic zooming x1.8 and perform
% sharpening
u0_registered = register_stack(u0_apod,T);
out_sam = median(u0_registered,3); % shift-and-median image
[Xq,Yq] = meshgrid(1+(0:(1.8*m-1))/1.8,1+(0:(1.8*n-1))/1.8); 
out_sam_bicubic = interp2(out_sam,Xq,Yq,'bicubic',0); % bicubic zoom x1.8 of out_sam 
out_sam_bicubic_sharpen = sharpening(out_sam_bicubic,@(rho)1+5*(1-exp(-rho)));

% perform least-squares super-resolution & Sharpening
uls = leastsquares_superres(u0_apod,T,M,N);
uls_sharpen = sharpening(uls,@(rho)1+5*(1-exp(-rho))); 

% display images (reproduce Figure 17, without cropping)
figure('Name','shift-and-median (out_sam)'); imview(out_sam,'black',20,'white',140); 
figure('Name','super-resolution (uls)'); imview(uls,'black',20,'white',140); 
figure('Name','(out_sam) bicubic zoom x 1.8 + sharpening'); imview(out_sam_bicubic_sharpen,'black',20,'white',140); 
figure('Name','(uls) sharpening'); imview(uls_sharpen,'black',20,'white',140); 

% compare close-up views of high-resolution images (out_sam_bicubic_sharpen, uls and uls_sharpen)
display_factor = 2; 
figure('Name','(out_sam) bicubic zoom x 1.8 + sharpening'); imview(out_sam_bicubic_sharpen,'black',20,'white',140,'scale',display_factor); 
figure('Name','uls'); imview(uls,'black',20,'white',140,'scale',display_factor); 
figure('Name','uls + sharpening'); imview(uls_sharpen,'black',20,'white',140,'scale',display_factor); 
```
