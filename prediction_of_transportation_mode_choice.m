% This program can be used to reproduce the empirical results concerning
% an empirical application in the prediction of transportation mode choice
% using the best subset maximum score approach of Chen and Lee (2017).

% The application requires the work-trip mode choice dataset of Horowitz (1993).
% The dataset can be also found in this code repository.

clear;
load('horowitz_data.mat'); % load the work-trip mode choice dataset

% create variables from training and validation samples
% data columins : [Y DCOST CARS DOVTT DIVTT]

q=1;  % variable selection constraint

warm_start = 1; % set this variable to 1 for using the warm start strategy
tau = 1.5; % the tuning parameter to construct the refined bound used in the warm start approach
mio = 1; % set mio to 1 for using Method 1 for the MIO formulation
         % set mio to 2 for using Method 2 for the MIO formulation

series_exp = 0;   % set this variable to 1 for using quadratic expansion terms as covariates

beta0=1;

b=10; % bound value

disp('estimation based on full sample');
Y_tr=data(:,1); temp=data(:,2:end);

n_tr=length(Y_tr);
% [DCOST CARS DOVTT DIVTT]
x_std=(temp-repmat(mean(temp),n_tr,1))./repmat(std(temp),n_tr,1);
x_foc = [x_std(:,1) ones(n_tr,1)];  % [DCOST Intercept]

if series_exp == 1
z2 = x_std(:,2); z3 = x_std(:,3); z4 = x_std(:,4);
x_aux1 = [z2 z3 z4]; % linear terms
x_aux2 = [z2.*z3 z3.*z4 z2.*z4];
x_aux3 = [z2.*z2 z3.*z3 z4.*z4];
x_aux = [x_aux1 x_aux2 x_aux3];
clear temp x_std x_aux1 x_aux2 x_aux3 z2 z3 z4;  
else
x_aux = x_std(:,2:4); % [CARS DOVTT DIVTT]    
clear temp x_std;  
end

k=size(x_foc,2);
d=size(x_aux,2);

bnd=[-b*ones(k-1+d,1) b*ones(k-1+d,1)]; % set the initial parameter bounds

tol = floor(sqrt(log(n_tr)*n_tr)/2); % set the tolerance level value
disp(['tolerance level: ', num2str(tol/n_tr)]);

time_limit = 86400; % set the MIO solver time limit

if warm_start == 1 % warm start MIO
[bhat,score,gap,rtime,ncount]  = warm_start_max_score(Y_tr,x_foc,x_aux,beta0,q,time_limit,tol,bnd,mio,tau);
else % cold start MIO 
[bhat,score,gap,rtime,ncount]  = max_score_constr_fn(Y_tr,x_foc,x_aux,beta0,q,time_limit,tol,bnd,mio);
end

disp('parameter estimates:');
disp(bhat);
disp('avg_score gap time node_count');
disp([score gap rtime ncount]);

