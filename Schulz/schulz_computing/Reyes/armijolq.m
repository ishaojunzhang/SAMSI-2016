% 
% Armijo Line Search for the Linear-Quadratic Problem
%
% cf. page 53 of
% Juan Carlos de los Reyes. Numerical PDE-Constrained Optimization. 
% Springer Briefs in Optimization. Springer, 2015.


function arm=armijolq(u,y,p,z,alpha,lap)
     countarm=0; 
     gradcost=norm(p+alpha*u)^2;
     cost1=1/2*norm(y-z)^2+alpha/2*norm(u)^2; 
     beta=1; armijo=1e5;

     while armijo > -1e-4*beta*gradcost 
         beta=1/2^(countarm);
         uinc=u-beta*(p+alpha*u); 
         yinc=lap\uinc;
         cost2=1/2*norm(yinc-z)^2+alpha/2*norm(uinc)^2; 
         armijo=cost2-cost1;
         countarm=countarm+1; 
     end
             
     arm=beta;
