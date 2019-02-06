%% Exmple code to run the proposed method without a T(G)V regularizer for 
%% Fourier measurements and Wavelet Frame

% initialize data
m       = [256, 256];
DefineBrain;
f       = RasterizePhantom(Brain,m(1));

% mask in k-space
mask    =   RadialMask(24,m(1),m(2),1)==1;

% fetch Operators
A = getFourierOperator(m,mask);
D = getWaveletOperator(m,2,3);
    
% fetch measurement operator
y       = A.times(f(:));

% set parameters
alpha       = 0; % switch of TV
beta        = 1e4;
mu1         = 0; % switch of TV
mu2         = 6e2;
maxIter     = 1.25e2;
adaptive    = 'NewIRL1';
correct     = @(x) real(x);
doTrack     = true;

%% solve
out = TVsolver(y,m,A,D,alpha,beta,mu1,mu2,'maxIter',maxIter,'adaptive',adaptive,'f',f,'correct',correct,'doTrack',doTrack);

% visualize convergence
figure;
plot(out.error);
title('Error')