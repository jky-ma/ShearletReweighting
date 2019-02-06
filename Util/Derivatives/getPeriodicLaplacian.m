function [L,Lx,Ly,D,Dx,Dy] = getPeriodicLaplacian(m)

% Function that constructs a 2D Laplacian with periodic boundary conditions.
%
% Input:    m - vector of dimensions
%
% Output:   L        - Laplacian 
%           Lx       - second derivative in x-direction
%           Ly       - second derivative in y-direction
%           D        - diagonal Laplacian in k-space, i.e. F*L*F'
%           Dx       - diagonal 2nd derivative in k-space, i.e. F*Dx*F'
%           Dy       - diagonal 2nd derivative in k-space, i.e. F*Dy*F'
%
% 
% written by Maximilian März, April 2016

    % 1D Laplacian in x 
    ex = ones(m(1),1);
    Dxx = spdiags([ex -2*ex ex], [-1 0 1], m(1), m(1));
    Dxx(1,end) = 1; Dxx(end,1) = 1;

    % 1D Laplacian in y
    ey = ones(m(2),1);
    Dyy = spdiags([ey -2*ey ey], [-1 0 1], m(2), m(2));
    Dyy(1,end) = 1; Dyy(end,1) = 1;

    %% 2D Laplacian
    Lx = kron(Dxx,speye(m(2)))*(-1);
    Ly = kron(speye(m(1)),Dyy)*(-1);
    L = Lx + Ly;

    % Construct BCCB matrix according to formula
    c = L(:,1);
    D = (fft2(full(reshape(c,m(1),m(2)))));   % (gradX'*gradX + gradY'*gradY) in k-space

    cx = Lx(:,1);
    Dx = (fft2(full(reshape(cx,m(1),m(2))))); % (gradX'*gradX) in k-space

    cy = Ly(:,1);
    Dy = (fft2(full(reshape(cy,m(1),m(2))))); % (gradY'*gradY) in k-space



    