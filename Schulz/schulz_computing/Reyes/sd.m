%
% steepest descent method for optimization of the Poisson equation
%
% cf. page 52f of
% Juan Carlos de los Reyes. Numerical PDE-Constrained Optimization. 
% Springer Briefs in Optimization. Springer, 2015.

clear all;
n=input('Mesh points: ' ); h=1/(n+1); 
alpha=input('Regularization parameter: ' );

[x1,y1]=meshgrid(h:h:1-h,h:h:1-h);	%%%%% Coordinates %%%%%

%%%%% Desired state %%%%% 
desiredstate=inline('x./x' , 'x' , 'y' );
%desiredstate=inline('x.*y' , 'x' , 'y' );
z=feval(desiredstate,x1,y1); z=reshape(z,n^2,1); 

lap=matrices(n,h);	%%%%% Laplacian %%%%%

%%%%% Initialization %%%%% 
u=sparse(n^2,1);
res=1; iter=0; tol=1e-4;

while res >= tol 
    iter=iter+1; 
    y=lap\u;      %%%%% State equation %%%%%
    p=lap\(y-z);  %%%%% Adjoint solver %%%%%
    beta=armijolq(u,y,p,z,alpha,lap);	%% Armijo line search %%
    uprev=u;
    u=u-beta*(p+alpha*u);	%%%%% Gradient step %%%%% 
    res=norm(u-uprev)
end
uu=reshape(u,[n,n]);
yy=reshape(y,[n,n]);