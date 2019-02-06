function [gradX,gradY,Dx,Dy] = getPeriodicGradBack(m)

% Function that constructs a 2D Gradient with periodic boundary conditions,
% based of backward differences.
%
% Input:    m - vector of dimensions
%
% Output:   gradX - gradient in x-direction
%           gradY - gradient in y-direction
%           Dx    - diagonal gradient in k-space, F*gradX*F'
%           Dy    - diagonal gradient in k-space, F*gradY*F'
%
% written by Maximilian März, April 2016


% periodic backward gradient in 1D 
per1Dx = spdiags([-1*ones(m(1),1),ones(m(1),1)],[-1,0],m(1),m(1));
per1Dy = spdiags([-1*ones(m(2),1),ones(m(2),1)],[-1,0],m(2),m(2));

per1Dx(1,end) = -1;
per1Dy(1,end) = -1;

% Identity matrices
I1 = speye(m(1));
I2 = speye(m(2));

% build gradient with kronecker product
gradX = kron(per1Dx,I2);
gradY = kron(I1,per1Dy);

% get diagonal in k-space with formula for bccb matrices
cX = reshape(gradX(:,1),m);
cY = reshape(gradY(:,1),m);
Dx = fft2(full(cX));
Dy = fft2(full(cY));