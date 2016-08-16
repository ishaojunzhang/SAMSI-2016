% 2D STRAIGHT INJECTION RESTRICTION OPERATOR
%
% coarse = restrict_strinj_2D(fine)
%
% input parameters:
%    fine........matrix of dimension (2^k-1) x (2^k-1)
% output parameter:
%    coarse......matrix of dimension (2^(k-1)-1) x (2^(k-1)-1)
function coarse = restrict_strinj_2D(fine)

[n,m] = size(fine); coarse = fine(2:2:n-1,2:2:m-1);
