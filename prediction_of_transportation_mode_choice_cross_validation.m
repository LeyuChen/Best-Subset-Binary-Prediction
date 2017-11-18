% This program can be used to reproduce the empirical results concerning
% an empirical application in the prediction of transportation mode choice
% using the best subset maximum score approach of Chen and Lee (2017).
% The program is concerned with the covariate specification using quadratic
% expansion terms as the auxiliary covariates.

% For a given value of q, the program outputs the performance summary statistics
% using training and testing datasets for the best subset prediction approach with
% variable selection bound q using the datasets in the 5-fold cross
% validation procedure.

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

series_exp = 1;

fold=size(tr_ind,2);

beta0=1;

b=10; % bound value
time_limit = 86400; % set the MIO solver time limit

bhat=zeros(10,fold);
    
score=zeros(fold,1);
gap=zeros(fold,1);
in_score=zeros(fold,1);
rtime=zeros(fold,1);
ncount=zeros(fold,1);
p_ratio=zeros(fold,1);

for i=1:fold

data_tr=data(tr_ind(:,i),:);
data_v=data(test_ind(:,i),:);

disp(['estimation based on training sample at fold: ',num2str(i)]);
Y_tr=data_tr(:,1); 
temp=data_tr(:,2:end);    

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

tol = floor(sqrt(log(n_tr)*n_tr)/2);
disp(['tolerance level: ', num2str(tol)]);

if warm_start == 1 % warm start MIO
[bhat(:,i),score(i),gap(i),rtime(i),ncount(i)]  = warm_start_max_score(Y_tr,x_foc,x_aux,beta0,q,time_limit,tol,bnd,mio,tau);
else % cold start MIO 
[bhat(:,i),score(i),gap(i),rtime(i),ncount(i)]  = max_score_constr_fn(Y_tr,x_foc,x_aux,beta0,q,time_limit,tol,bnd,mio);
end

disp('coefficient values: ');
disp(bhat(:,i)');
disp(['gurobi score: ',num2str(score(i))]);
disp(['gurobi absolute gap: ',num2str(gap(i))]);
disp(['gurobi running time: ',num2str(rtime(i))]);
disp(['gurobi node count: ',num2str(ncount(i))]);
in_score(i)=sum(Y_tr==([x_foc x_aux]*[1;bhat(:,i)]>0));
disp(['in-sample score: ',num2str(in_score(i))]);

% validation sample

Y_val=data_v(:,1);
n_val=length(Y_val);
temp=data_v(:,2:end);
x_std=(temp-repmat(mean(temp),n_val,1))./repmat(std(temp),n_val,1);

if series_exp == 1
z2 = x_std(:,2); z3 = x_std(:,3); z4 = x_std(:,4);
x_aux1 = [z2 z3 z4]; % linear terms
x_aux2 = [z2.*z3 z3.*z4 z2.*z4];
x_aux3 = [z2.*z2 z3.*z3 z4.*z4];

x_aux = [x_aux1 x_aux2 x_aux3];
x_v = [x_std(:,1) ones(n_val,1) x_aux];
clear temp x_std x_aux x_aux1 x_aux2 x_aux3 z2 z3 z4;  
else
x_v = [x_std(:,1) ones(n_val,1) x_std(:,2:4)];
clear temp x_std;  
end

y_hat=((x_v*[1;bhat(:,i)])>0);
p_ratio(i)=mean(Y_val==y_hat);
disp(['out-of-sample performance: ',num2str(p_ratio(i))]);

end
disp('average coefficient vector: ');
disp(mean(bhat,2)');
disp(['average score: ',num2str(mean(score))]);
disp(['average gap: ',num2str(mean(gap))]);
disp(['average running time: ',num2str(mean(rtime))]);
disp(['average node count: ',num2str(mean(ncount))]);
disp(['average in-sample score: ',num2str(mean(in_score))]);
disp(['average out-of-sample performance: ',num2str(mean(p_ratio))]);



