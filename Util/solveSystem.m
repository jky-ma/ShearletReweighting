function [x,y,z] = solveSystem(d,X,FB1,FB2,FB3,iterative,uOld)
% Function that calculates with Gauss elimination the solution of the linear system A*x =
% [FB1,FB2,FB3], where A is defined as below.
% 
% Input:    d           - struct containing blocks of A, d{1} might be a
%                               funtionhandle
%           X           - contains pre-calculated parts of the linear system
%           FB1-3       - RHS
%           iterative 	- flag 
%           uOld        - initial guess in case of iterative solver
% Output:   [x,y,z] - solution of the system.
%
% Written by Maximilian März, April 2016

% vectorization
vec = @(x) x(:);


% Calculate RHS after transform the linear system to an lower-triangular
% system
q = FB3;
p = FB2 - FB3.*conj(d{6})./(d{2}-(d{6}.*conj(d{6}))./d{3});
u = FB1 - FB3.*conj(d{5})./d{3} - (FB2 - (FB3.*conj(d{6}))./d{3}).*(conj(d{4}) - (d{6}.*conj(d{5}))./d{3})...
            ./(d{2} - (d{6}.*conj(d{6})./d{3}));

% solve first line of lower-triangular system
if iterative
    [temp,~] = pcg(X{1,1},vec(u),1e-10,75,[],[],vec(uOld));
    x        = reshape(temp,size(X{2,1}));
else
    x = u./(X{1,1});
end
% solve last two lines of linear system
y = (p-X{2,1}.*x)./X{2,2};
z = (q-d{5}.*x-d{6}.*y)./d{3};
    
    