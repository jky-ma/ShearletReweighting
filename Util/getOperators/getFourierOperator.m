function A = getFourierOperator(m,mask)
%% function to construct Fourier measurement operator, i.e.
%   M*F, for F unitary 2D/3D Fourier-matrix, and M mask-operator
% 
%% Input: 
%         m     -   dimensions
%         mask  -   mask in k-space
%
%% Output:
%        A      - contains:
%                   A.adj()         -   adjoint
%                   A.times()       -   measurement operator
%                   A.mask          -   mask
% 
% written by Maximilian März, April 2016

    vec = @(x) x(:);
    if length(m) == 2
        maskOp  =   opMask(prod(m),mask);
        F       =   opDFT2(m(1),m(2),true);
    
        A.times = @(x) maskOp*F*x;
        A.adj   = @(x) F'*maskOp*x;
        A.mask  = mask;
        
    elseif length(m) == 3
        A.mask  = mask;
        % create Mask as operator, mask is centered in k-space
        % mask = opMask(prod(m),mask);

        % sampling in k-space
        A.times = @(x) vec((A.mask).*(fftshift(fftn((reshape(x,m))./sqrt(prod(m))))));
        A.adj   = @(x) vec((ifftn(ifftshift(reshape((A.mask).*(reshape(x,m)),m))).*sqrt(prod(m))));
end
        