function [L,Lx,Ly,Lz,D,Dx,Dy,Dz] = getPeriodicLaplacian3D(m)

% Function that constructs a 3D Laplacian with periodic boundary conditions.
%
% Input:    m - vector of dimensions
%
% Output:   L        - Laplacian 
%           Lx       - second derivative in x-direction
%           Ly       - second derivative in y-direction
%           D        - diagonal Laplacian in k-space, i.e. F*L*F'
%           Dx       - diagonal 2nd derivative in k-space, i.e. F*Dx*F'
%           Dy       - diagonal 2nd derivative in k-space, i.e. F*Dy*F'
% Written by Maximilian März, April 2016   

    % 1D Laplacian in x 
    ex = ones(m(1),1);
    Dxx = spdiags([ex -2*ex ex], [-1 0 1], m(1), m(1));
    Dxx(1,end) = 1; Dxx(end,1) = 1;

    % 1D Laplacian in y
    ey = ones(m(2),1);
    Dyy = spdiags([ey -2*ey ey], [-1 0 1], m(2), m(2));
    Dyy(1,end) = 1; Dyy(end,1) = 1;
    
    % 1D Laplacian in z
    ez = ones(m(3),1);
    Dzz = spdiags([ez -2*ez ez], [-1 0 1], m(3), m(3));
    Dzz(1,end) = 1; Dzz(end,1) = 1;
    
    %% 3D Laplacian
    Lx = kron(Dxx,kron(speye(m(2)),speye(m(3))))*(-1);
    Ly = kron(speye(m(1)),kron(Dyy,speye(m(3))))*(-1);
    Lz = kron(kron(speye(m(1)),speye(m(2))),Dzz)*(-1);
    L = Lx + Ly + Lz;
    
    % Construct BCCB matrix according to formula
    c = L(:,1);
    D = (fftn(full(ndSparse(c,[m(1),m(2),m(3)]))));   % -1*(gradX'*gradX + gradY'*gradY) in k-space

    cx = Lx(:,1);
    Dx = (fftn(full(ndSparse(cx,[m(1),m(2),m(3)])))); % -1*(gradX'*gradX) in k-space

    cy = Ly(:,1);
    Dy = (fftn(full(ndSparse(cy,[m(1),m(2),m(3)])))); % -1*(gradY'*gradY) in k-space
    
    cz = Lz(:,1);
    Dz = (fftn(full(ndSparse(cz,[m(1),m(2),m(3)])))); % -1*(gradZ'*gradZ) in k-space