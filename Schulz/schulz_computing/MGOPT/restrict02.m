% 2D FULL WEIGHTING RESTRICTION OPERATOR
%
% coarse = restrict02(fine)
%
% input parameters:
%    fine........matrix of dimension (2^k-1)x(2^k-1)
% output parameter:
%    coarse......matrix of dimension (2^(k-1)-1)x(2^(k-1)-1)
function coarse = restrict02(fine)

coarse = fine(2:2:end-1,2:2:end-1)./4 + ...
    (fine(1:2:end-2,2:2:end-1) + fine(3:2:end,2:2:end-1) + ...
    fine(2:2:end-1,1:2:end-2) + fine(2:2:end-1,3:2:end))./8 + ...
    (fine(1:2:end-2,1:2:end-2) + fine(3:2:end,3:2:end) + ...
    fine(1:2:end-2,3:2:end) + fine(3:2:end,1:2:end-2))./16;
