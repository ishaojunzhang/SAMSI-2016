% 
% SQP (full Newton) Method for the Optimization of a Semilinear Equation
%
% cf. page 68 of
% Juan Carlos de los Reyes. Numerical PDE-Constrained Optimization. 
% Springer Briefs in Optimization. Springer, 2015.

clear all;
n=input('Mesh   points:    '); h=1/(n+1); 
alpha=input(' Regularization parameter: '); 
[x1,y1]=meshgrid(h:h:1-h,h:h:1-h);	%%%%% Coordinates  %%%%%

%%%%% Desired state %%%%% 
desiredstate=inline('x.*y' , 'x', 'y');
z=feval(desiredstate,x1,y1); z=reshape(z,n^2,1) ; 
lap=matrices(n,h);	%%%%%  Laplacian %%%%%

%%%%% Initialization %%%%%
u=sparse(n^2,1) ; y=sparse(n^2,1); p=sparse(n^2,1) ; 
res=1; iter=0;

while res >= 1e-3 
    iter=iter+1

    %%%%% SQP step %%%%%
    Y=spdiags(y,0,n^2,n^2); P=spdiags(p,0,n^2,n^2);

    A=[speye(n^2)-6*Y.*P sparse(n^2,n^2) -lap+3*Y.^2 
        sparse(n^2,n^2) alpha*speye(n^2) speye(n^2)
        -lap+3*Y.^2 speye(n^2) sparse(n^2,n^2)];
    
    F=[lap*p+3*Y.^2*p-y+z;-p-alpha*u;lap*y+y.^3-u];
   
    delta=A\F;
    uprev=u; yprev=y; pprev=p; 
    y=y+delta(1:n^2); 
    u=u+delta(n^2+1:2*n^2); 
    p=p+delta(2*n^2+1:3*n^2);
    res=norm(u-uprev) +norm(y-yprev) +norm(p-pprev)
end

