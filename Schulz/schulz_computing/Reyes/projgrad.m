% 
% Projected Gradient Method for the Optimization of the Poisson Equation
%
% cf. page 79f of
% Juan Carlos de los Reyes. Numerical PDE-Constrained Optimization. 
% Springer Briefs in Optimization. Springer, 2015.

clear all;
n=input('Mesh  points:   ');  h=1/(n+1); 
alpha=input(' Tikhonov regularization parameter: '); 
ua=input(' Lower bound: '); ub=input('Upper bound: '); 
Armijo=input('Line  search (1=yes): ');
[x1,y1]=meshgrid(h:h:1-h,h:h:1-h);	%%%%% Coordinates %%%%%

%%%%% Desired state %%%%% 
desiredstate=inline('x.*y' , 'x' , 'y');
z=feval(desiredstate,x1,y1); z=reshape(z,n^2,1); 

lap=matrices(n,h);	%%%%% Laplacian %%%%%

%%%% Initialization %%%%% 
u=sparse(n^2,1); res=1; iter=0;

while res >= 1e-3 
    iter=iter+1 
    y=lap\u;       %%%%% State equation %%%%%
    p=lap\(y-z);   %%%%% Adjoint solver %%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Armijo line search
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    countarm=0; beta=1;
    gradcost=norm(max(ua,min(u-beta*(p+alpha*u),ub))-u); 
    cost1=1/2*norm(y-z)^2+alpha/2*norm(u)^2;

    if Armijo==1
        armijo=1e5;
        while armijo > -1e-4/beta*gradcost^2 
            beta=1/2^(countarm);
            
            uinc=max(ua,min(u-beta*(p+alpha*u),ub)); 
            yinc=lap\uinc;

            cost2=1/2*norm(yinc-z)^2+alpha/2*norm(uinc)^2; 
            armijo=cost2-cost1;
            gradcost=norm(max(ua,min(u-beta*(p+alpha*u),ub))-u); 
            countarm=countarm+1;
        end
    end
    uprev=u;
    %%%% Projected gradient step %%%% 
    u=max(ua,min(u-beta*(p+alpha*u),ub)); 
    res=norm(u-uprev)
end
