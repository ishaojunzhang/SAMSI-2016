% 2D FIVE POINT NEGATIVE LAPLACE OPERATOR
%
% v = neg_lap02(u,h_inv2)
%
% input parameters:
%    u........matrix
%    h_inv2...1/h^2, where h is the mesh size in both directions
% output parameter:
%    v........matrix -Delta u, same size as u
function v = neg_lap02(u,h_inv2)

% compute negative laplacean:
vx = [2*u(:,1) - u(:,2), -u(:,1:end-2) + 2*u(:,2:end-1) - ...
    u(:,3:end), -u(:,end-1) + 2*u(:,end)];
vy = [2*u(1,:) - u(2,:); -u(1:end-2,:) + 2*u(2:end-1,:) - ...
    u(3:end,:); -u(end-1,:) + 2*u(end,:)];

v = (vx + vy)*h_inv2;
