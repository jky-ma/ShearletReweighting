%% Exmple code to run the proposed method without a T(G)V regularizer for
%% Fourier measurements and Shearlet/Wavelet Frame

% initialize data
m       = [256, 256];
% DefineBrain;
f       = RasterizePhantom(Brain,m(1));

m       = [256, 256];
f       = double(imresize(rgb2gray(f),m));
f       = f./max(f(:));

% mask in k-space
mask    =   RadialMask(30,m(1),m(2),1)==1;

% fetch Operators
A = getFourierOperator(m,mask);

%% Wavelet
% D = getWaveletOperator(m,2,3);
%% Shearlet
D = getShearletOperator(m,[1 1 2]);

% fetch measurement operator
y       = A.times(f(:));

% set parameters
beta        = 1e5;
alpha       = [1 1];
mu          = [5e3, 1e1, 2e1];
epsilon     = 1e-5;

maxIter     = 5e1;
adaptive    = 'NewIRL1';
correct     = @(x) real(x);
doTrack     = true;

%% solve
out = TGVsolver(y,m,A,D,alpha,beta,mu,'maxIter',maxIter,'adaptive',adaptive,'f',f,'correct',correct,'epsilon',epsilon,'doTrack',doTrack);

% visualize convergence
figure;
plot(out.error);
title('Error')
