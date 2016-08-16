% COLLECTIVE GAUSS-SEIDEL SMOOTHING
%
% [y,p,u] = smooth_oc01(y,p,u,f,z,nu,cycles)
%
% input parameters:
%    y........state
%    p........lagrange multiplier
%    u........control
%    f........right hand side of state equation
%    z........desired state
%    nu.......weighting parameter in cost functional
%    cycles...number of smoothing cycles
% output parameter:
%    y........state after smoothing
%    p........multiplier after smoothing
%    u........control after smoothing
function [y,p,u] = smooth_oc01(y,p,u,f,z,nu,cycles)

% mesh size h:
[n,m] = size(f); h = 1/(n + 1); h2 = h*h;

% embed variables in zeros:
y = [zeros(1,n+2); zeros(n,1),y,zeros(n,1); zeros(1,n+2)]; 
p = [zeros(1,n+2); zeros(n,1),p,zeros(n,1); zeros(1,n+2)];

% loop over number of cycles:
for rv_cycles = 1:cycles
    % running variable over columns:
    for rvc = 2:n+1
        % running variable over rows:
        for rvr = 2:n+1
            Aij = y(rvr+1,rvc) + y(rvr-1,rvc) + y(rvr,rvc+1) + ...
                y(rvr,rvc-1) - h2*f(rvr-1,rvc-1);
            Bij = p(rvr+1,rvc) + p(rvr-1,rvc) + p(rvr,rvc+1) + ...
                p(rvr,rvc-1) - h2*z(rvr-1,rvc-1);
            u(rvr-1,rvc-1) = (4*Bij + h2*Aij)/(16*nu + h2*h2);
            y(rvr,rvc) = (Aij - h2*u(rvr-1,rvc-1))/4;
            p(rvr,rvc) = nu*u(rvr-1,rvc-1);
        end
    end
end

% eliminate the zero boundary values:
y = y(2:n+1,2:n+1); p = p(2:n+1,2:n+1);
