function y = IsotropicShrink(x,thresh)
%% Isotropic shrinkage/proximal operator in R^n spaces for l1-norm
% Calculates the solution to the following problem
%    (1)          min_z ||z||_1 +(1/2*thresh) ||z - x||_2^2.     
%
% For an m x n x k-array the isotropic l_1 norm is defined as 
%                  \sum_{i,j}  ||x(i,j,:)||_2.
% Therefore (1) reduces to the subproblems
%    (2)          min_z ||z||_2 +(1/2*thresh) ||z - x||_2^2,   (z \in R^k),
% which are solved by max{||x||_2 - thresh,0} (x/||x||_2).
%
% Input:    x       - m x n x k-array containing the data
%           thresh  - thesholding parameter
% Output    y       - the shrinked data.
%
% Written by Maximilian März, April 2016

% calculate pixelwise 2-norm
k = size(x,3);
r = sqrt(sum(abs(x).^2,3));

y = zeros(size(x));

% solution as given in (2)
for i = 1:k
    y(:,:,i) = max(r-thresh,0).*x(:,:,i)./r;
end
y(isnan(y)) = 0;
