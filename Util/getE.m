function A = getE(p,m,gradX,gradY)

% Function that calculates E^b(p) for TGV
%
% Input:     p - input argument
%            m - vector of dimensions
%            gradX - derivative in x-direction
%            gradY - derivative in y-direction
%
% Output:    A  -  3D matrix containing the entries of E(p) in
%                       lexicographical order 
%
% Written by Maximilian März, April 2016

A = zeros(m(1),m(2),4);
p1 = p(:,:,1);
p2 = p(:,:,2);

A(:,:,1) = reshape(gradX*p1(:),m);
A(:,:,2) = 0.5*reshape((gradY*p1(:)+gradX*p2(:)),m);
A(:,:,3) = 0.5*reshape((gradY*p1(:)+gradX*p2(:)),m);
A(:,:,4) = reshape(gradY*p2(:),m);

end