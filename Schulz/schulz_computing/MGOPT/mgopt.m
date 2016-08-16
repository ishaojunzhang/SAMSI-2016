% MAIN PROGRAM
% multi-level scheme solving the optimal control problem
%
%           min  1/2*|y-z|^2 + nu/2*|u|^2
%           s.t. laplace y = u + f
%
% on the 2D domain (0,1)x(0,1).

clear fprintf('\n\nMULTILEVEL ALGORITHM:\n') fprintf('build data
...')

% grid:
% -----
% the domain (0,1)x(0,1) is discretized using an equidistant grid
% with (2^k-1) inner grid points per dimension on each level of
% the cycle. the mesh size is given by h = 1/(2^k). the finest
% grid is defined by k_max, the coarsest by k_min.
% the mesh sizes of the grids are therefore defined by:
%           h = [1/(2^k_max),1/(2^(k_max-1)),...,1/(2^k_min)]
k_max = 10;     % finest mesh size: 1/(2^k_max)
k_min = 2;      % coarsest mesh size: 1/(2^k_min)

% parameters:
% -----------
gamma = 1;         % number of recursive multilevel calls
pre = 2;           % number of pre-smoothing steps
post = 2;          % number of post-smoothing steps
tolerance = 1e-9;  % tolerance for stopping criterion
counter_max = 20;  % maximal number of cycles

% initialize values:
% ------------------
% build finest grid:
h = 2^(-k_max);    % mesh size
h_inv = 1/h; grid = linspace(h,1-h,2^k_max -1)'; [xgrid,ygrid] = meshgrid(grid);
% weighting parameter in cost functional:
nu = 1e-3;
% right-hand sides:
f = -4*pi*pi*sin(2*pi*xgrid).*(2*cos(2*pi*ygrid) - 1) - ...
    sin(pi*xgrid).*(ygrid.*ygrid - ygrid);
z = nu*sin(pi*xgrid).*(2 - pi*pi*(ygrid.*ygrid - ygrid)) +...
    sin(2*pi*xgrid).*(cos(2*pi*ygrid) - 1);
% norms of right hand sides (for relative residual):
norm_f_inv = 1/norm(f,'fro'); norm_z_inv = 1/norm(z,'fro');

% exact solutions:
y_exact = sin(2*pi*xgrid).*(cos(2*pi*ygrid)-1); u_exact = sin(pi*xgrid).*(ygrid.*ygrid - ygrid); p_exact = nu*u_exact;
% initial guess:
y = sin(20*pi*xgrid).*(cos(20*pi*ygrid)-1); p = sin(20*pi*xgrid).*(cos(20*pi*ygrid)-1); 
u = sin(20*pi*xgrid).*(cos(20*pi*ygrid)-1);
% initialize residual:
rel_resid_state = norm(neg_lap02(y,h_inv*h_inv) + f + u,...
    'fro')*norm_f_inv;
rel_resid_adj   = norm(neg_lap02(p,h_inv*h_inv) + z - y,...
    'fro')*norm_z_inv;
% initialize counter for cycles:
counter = 0; fprintf('done!\n')

% start multilevel scheme:
% ------------------------
fprintf('start mutilevel scheme ...\n');
start_time = cputime; 
while ((rel_resid_state > tolerance) | (rel_resid_adj > tolerance) &...
        (counter <= counter_max))
   % call recursive multilevel scheme:
   [y,p,u] = multilevel_recursive_oc01(y,p,u,f,z,nu,gamma,...
       pre,post,k_max,k_min);
   % increase counter:
   counter = counter + 1;
   % compute relative residuals:
   rel_resid_state = norm(neg_lap02(y,h_inv*h_inv) + f + u,...
       'fro')*norm_f_inv;
   rel_resid_adj   = norm(neg_lap02(p,h_inv*h_inv) + z - y,...
       'fro')*norm_z_inv;
   fprintf('   iter %4i, relres_state = %6.4e, relres_adj = %6.4e\n',...
       counter,rel_resid_state,rel_resid_adj);
   fprintf('   |y - y_ex| = %6.4e\n',norm(y-y_exact,'fro')*h*h);
   fprintf('   |u - u_ex| = %6.4e\n',norm(u-u_exact,'fro')*h*h);
end
fprintf('done!\n');
elapsed_time = cputime - start_time;

% output:
fprintf('results:\n')
fprintf('  no. of cycles: %4i\n',counter)
fprintf('  relative residual of state eq.: %6.4e\n',rel_resid_state)
fprintf('  relative residual of adj. eq. : %6.4e\n',rel_resid_adj)
fprintf('  elapsed time (s): %8.4f\n',elapsed_time)

fprintf('data:\n')
fprintf('  no. of inner grid points on finest grid:   %6i\n',...
    (2^k_max-1)^2)
fprintf('  no. of inner grid points on coarsest grid: %6i\n',...
    (2^k_min-1)^2)
fprintf('  k_min, k_max : %2i, %2i\n', k_min,k_max)
fprintf('  gamma = %2i, nu1 = %2i, nu2 = %2i\n',gamma,pre,post)

