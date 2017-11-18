
% This program can be used to reproduce the simulation results of the
% performance of the best subset maximum score approach of Chen and Lee (2017).

clear;

rng(1,'twister');

% initialize simulation configurations

warm_start = 0; % set this variable to 1 for using the warm start strategy
tau = 1.5; % the tuning parameter to construct the refined bound used in the warm start approach
mio = 1; % set mio to 1 for using Method 1 for the MIO formulation
         % set mio to 2 for using Method 2 for the MIO formulation

q=1; % change here for the value of variable selection bound  

N = 100; % size of the training sample
N_val = 5000; % size of the validation sample
R = 100; % simulation repetitions

p=10;  % change here for the auxiliary covariate dimension

type=1;    % type = 1 ==> heteroskedastic error design
           % type = 0 ==> homoskedastic error design

% [X intercept] : focus covariates
beta0=1;
if type==1
beta_s = -1.5;
else
beta_s = -0.35;   
end

beta = [beta0; 0; beta_s; zeros(p-1,1)];

K=length(beta);

bhat=zeros(K-1,R);


rho=0.25;
sigma=ones(K-1,1);
for i=1:K-2
sigma(i+1)=rho^i;
end
sigma=toeplitz(sigma);

gap=zeros(R,1); % MIO gap
rtime=zeros(R,1); % MIO running time
ncount=zeros(R,1); % MIO node count
score=zeros(R,1); % MIO score

DGP_score=zeros(R,1); % in-sample score at the DGP parameter vector
val_score=zeros(R,1); % in-sample score at the estimated parameter vector
DGP_score_test=zeros(R,1); % out-of-sample score at the DGP parameter vector
val_score_test=zeros(R,1); % out-of-sample score at the estimated parameter vector

bnd=[-10*ones(size(bhat,1),1) 10*ones(size(bhat,1),1)];
bnd_h = zeros(size(bhat,1),2,R);

maxT=0;

if q>=1 && p>N
tol=min([0.5*sqrt((1+q)*log(p)*N);0.05*N]);  % early stopping rule
else
tol=0;
end

for i=1:R
disp(i);    
[y,datax] = simulation_data(N,beta,sigma,type);

try

if warm_start == 1 % warm start MIO
[bhat(:,i),score(i),gap(i),rtime(i),ncount(i)]  = warm_start_max_score(y,datax(:,1:2),datax(:,3:end),beta0,q,maxT,tol,bnd,mio,tau);
else % cold start MIO 
[bhat(:,i),score(i),gap(i),rtime(i),ncount(i)]  = max_score_constr_fn(y,datax(:,1:2),datax(:,3:end),beta0,q,maxT,tol,bnd,mio);
end    

catch gurobiError
    fprintf('Error reported\n');
end
 
DGP_score(i) = mean(y == ((datax*beta)>=0)); 
val_score(i) = mean(y == ((datax*[beta0;bhat(:,i)])>=0)); 

if N_val>0
  [y_val,datax_val] = simulation_data(N_val,beta,sigma,type);
  DGP_score_test(i) = mean(y_val == ((datax_val*beta)>=0)); 
  val_score_test(i) = mean(y_val == ((datax_val*[beta0;bhat(:,i)])>=0)); 
 end

end

disp(mean([val_score DGP_score val_score_test DGP_score_test]));
disp(mean([val_score./DGP_score val_score_test./DGP_score_test]));

