
% function input :
% y : vector of binary outcomes
% x_foc : (n by (k+1)) matrix of data for focused covariates, 
%         of which the first column should contain data of 
%         the contiuous regressor with respect to which 
%         scale normalization is imposed
% x_aux : (n by p) matrix of data for auxiliary covariates which will be
%         selected based on the best subset covariate selection procedure
% beta0 : the coefficient taking value either 1 or -1 to normalize the 
%         scale for the first covariate in x_foc       
% q     : the cardinality constraint for the covariate selection
% T     : the time limit specified for the MIO solver
% bnd   : (((k-1)+d) by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients
%         the first (k-1) rows correspond to the bounds of the focused 
%         covariates excluding the first one. 
%         the remaining d rows correspond to the bounds of the auxiliary
%         covariates
% mio   : MIO formulation of the best subset maximum score problem
%         set mio = 1 for Method 1 of the paper
%         set mio = 2 for Method 2 of the paper
% tau   : the tuning parameter for enlarging the estimated bounds

% function output :
% bhat  : the maximum score estimates for the unknown coefficients
% score : the value of maximum score objective function
% gap   : the MIO optimization gap value in case of early termination
%         gap = 0 ==> optimal solution is found within the time limit
% rtime : the time used by the MIO solver in the estimation procedure

function [bhat,score,gap,rtime,ncount] = warm_start_max_score(y,x_foc,x_aux,beta0,q,T,tol,bnd,mio,tau)
bnd_h = get_bnd(y,[x_foc x_aux],beta0,bnd);
bnd_abs = tau*max(abs(bnd_h),[],2);
bnd0 = [max([-bnd_abs bnd(:,1)],[],2) min([bnd_abs bnd(:,2)],[],2)]; % this is the refind bound used for warm start MIO
[bhat,score,gap,rtime,ncount]  = max_score_constr_fn(y,x_foc,x_aux,beta0,q,T,tol,bnd0,mio);
end

