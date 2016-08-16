% 2D LINEAR INTERPOLATION OPERATOR
%
% fine = interpolate02(coarse)
%
% input parameters:
%    coarse......matrix of dimension (2^(k-1)-1)x(2^(k-1)-1)
% output parameter:
%    fine........matrix of dimension (2^k-1)x(2^k-1)
function fine = interpolate02(coarse)

% initialize vector:
[m,n] = size(coarse); fine  = zeros(2*m + 1,2*m + 1);

% set on coarse grid:
fine(2:2:end-1,2:2:end-1) = coarse;
% interpolate in x-direction:
fine(2:2:end-1,1:2:end-2) = 0.5.*coarse; 
fine(2:2:end-1,3:2:end) = fine(2:2:end-1,3:2:end) + 0.5.*coarse;

% interpolate in y-direction:
fine(1:2:end-2,:) = 0.5.*fine(2:2:end-1,:); 
fine(3:2:end,:) = fine(3:2:end,:) + 0.5.*fine(2:2:end-1,:);
