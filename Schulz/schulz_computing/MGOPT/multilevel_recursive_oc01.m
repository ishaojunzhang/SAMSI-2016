%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% RECURSIVE MULTILEVEL SOLVER
%
%  [y,p,u] = multilevel_recursive_oc01(y0,p0,u0,f,z,nu,gamma,pre,post,k,k_min)
%
% input parameters:
%     y0......initial guess for state
%     p0......initial guess for lagrangean multiplier
%     u0......initial guess for control
%     f.......right hand side of state equation
%     z.......desired state
%     nu......weighting parameter in cost functional
%     gamma...number of recursive calls
%     pre.....number of pre-smoothing steps
%     post....number of post-smoothing steps
%     k.......number corresponding to grid (h = 1/2^k)
%     k_min...coarsest grid
% output parameter:
%     y.......state
%     p.......lagrangean multiplier
%     u.......control
function [y,p,u] = ...
    multilevel_recursive_oc01(y0,p0,u0,f,z,nu,gamma,pre,post,k,k_min)

% mesh size:
h_inv = 2^k; h = 1/h_inv;

% if on coarsest grid, then solve exactly:
if k <= k_min
    % build matrix:
    n_temp = 2^k-1;
    e = ones(n_temp,1)*h_inv;
    A = spdiags([e -2*e e],-1:1,n_temp,n_temp);
    % 2D (positive) Laplace-matrix A:
    I_h = speye(n_temp)*h_inv;
    A = kron(I_h,A) + kron(A,I_h);    % n2 x n2
    clear e I_h

    % solve exactly:
    yp = [A,-1/nu.*speye(n_temp*n_temp);speye(n_temp*n_temp),A]\...
        [f(:);z(:)];
    y = reshape(yp(1:n_temp*n_temp),n_temp,n_temp);
    p = reshape(yp(n_temp*n_temp+1:end),n_temp,n_temp);
    u = p./nu;
    return

    % otherwise smooth, compute residual, restrict, call recursive
    % multilevel scheme, compute coarse-grid-correction and smooth
    % again:
else
    % initialize solution:
    y = y0;
    p = p0;
    u = u0;

    % perform multigrid scheme gamma times:
    for rv = 1:gamma
        % coarse mesh size:
        H_inv = 2^(k-1);
        H = 1/H_inv;
        % pre-smoothing:
        [y,p,u] = smooth_oc01(y,p,u,f,z,nu,pre);

        % compute residual:
        resid_state = f + neg_lap02(y,h_inv*h_inv) + u;
        resid_adj   = z + neg_lap02(p,h_inv*h_inv) - y;

        % restrict to coarser grid:
        resid_state_coarse = restrict02(resid_state);
        resid_adj_coarse = restrict02(resid_adj);

        % compute straight injection:
        y_strinj_coarse = restrict_strinj_2D(y);
        p_strinj_coarse = restrict_strinj_2D(p);
        u_strinj_coarse = restrict_strinj_2D(u);


        % right hand side of coarse problem:
        f_coarse = resid_state_coarse - ...
            neg_lap02(y_strinj_coarse,H_inv*H_inv) - u_strinj_coarse;
        z_coarse = resid_adj_coarse - neg_lap02(p_strinj_coarse,...
            H_inv*H_inv) + y_strinj_coarse;

        % apply recursive multilevel scheme:
        [y_coarse,p_coarse,u_coarse] = ...
        multilevel_recursive_oc01(y_strinj_coarse,...
        p_strinj_coarse,u_strinj_coarse,f_coarse,z_coarse,...
        nu,gamma,pre,post,k-1,k_min);


        % coarse-grid-correction:
        y = y + interpolate02(y_coarse - y_strinj_coarse);
        p = p + interpolate02(p_coarse - p_strinj_coarse);
        u = 1/nu.*p;

        % post-smoothing:
        [y,p,u] = smooth_oc01(y,p,u,f,z,nu,post);
    end
end

