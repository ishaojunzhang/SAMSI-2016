% Code for the Gibbs Samling Example from SAMSI Optimization Summer School, 2016
% Code written by David Bortz
 
clear all
close all
stream = RandStream.getGlobalStream;
seed = 1244;
reset(stream,seed);

n = 100;
alpha = 0;
beta = 2;
sig2 = 0.5;
true = [alpha,beta,sig2];
x = normrnd(zeros(n,1),1);%randn(n,1);
y = normrnd(alpha+beta.*x,sqrt(sig2));%alpha+beta.*x+sqrt(sig2)*randn(n,1);

% Prior hyperparameters
alpha0=0;
tau2a=10;
beta0=0;
tau2b=10;
nu0=3;
s02=1;
nu0s02=nu0*s02;

% Setting up starting values
alpha=0;
beta=0;
sig2=1;

% Gibbs sampler
M = 1000;
draws = [zeros(M,3)];
for i=1:M
    Var = 1/(1/tau2a+n/sig2);
    Mean = Var*(sum(y-beta*x)/sig2+alpha0/tau2a);
    alpha = normrnd(Mean,sqrt(Var));%Mean+sqrt(var)*randn(1);
    Var = 1/(1/tau2b+sum(x.^2)/sig2);
    Mean = Var*(sum((y-alpha).*x)/sig2+beta0/tau2b);
    beta = normrnd(Mean,sqrt(Var));%Mean+sqrt(var)*randn(1);
    sig2 = 1/gamrnd((nu0+n)/2,1./((nu0s02+sum((y-alpha-beta*x).^2)/2))); %v1
%     sig2 = 1/gamrnd((nu0s02+sum((y-alpha-beta*x).^2)/2),(nu0+n)/2); %v2
    draws(i,:) = [alpha,beta,sig2];
end

% Markov chains + marginal posterior
names = {'alpha','beta','sig2'};
ind = 100:M;
mfrow=[3,3];
figure;clf;


for i=1:3
        subplot(3,3,(i-1)*3+1);
        plot(draws(ind,i),'k');hold on
        plot(draws(:,i),'k');hold on
        plot([1 M],mean(draws(ind,i))*[1 1],'r')
        plot([100 100],[min(draws(ind,i)) max(draws(ind,i))],'b')
        axis([-25 1025 min(draws(ind,i)) max(draws(ind,i))]);
        title(names(i));
        xlabel('iterations')
        
        subplot(3,3,(i-1)*3+2);
        autocorr(draws(ind,i),30);
        axis([-1 31 -.1 1.1])
        
        subplot(3,3,(i-1)*3+3);
        hist(draws(ind,i));
        hold on
        ax = axis;
        plot(true(i)*[1 1], [-10 1000],'r');
        ylabel('Density')
        axis([min(draws(ind,i)) max(draws(ind,i)) 0 300])
end
