% 
% function lap for 5-point-stencil Dirichlet
%
% cf. page 23 of
% Juan Carlos de los Reyes. Numerical PDE-Constrained Optimization. 
% Springer Briefs in Optimization. Springer, 2015.

function [lap]=matrices(n,h)
d(n:n:n^2)=1;
d=d';
e=sparse(n^2,1);
e(1:n:(n-1)*n+1)=1;
b=ones (n^2,1);
a=[b,b-d, -4*b,b-e,b];
lap=-1/(h^2)*spdiags(a, [-n,-1,0,1,n],n^2,n^2);
