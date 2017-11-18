
% function input :
% x       : (n by k) matrix of covariate data 
% beta0 : the coefficient for the first covariate in x       
% bnd   : ((k-1) by 2) matrix where the first and second columns  
%            respectively store the lower and upper bounds 
%            of the unknown coefficients

% function output :
% bound : ((k-1) by 2) matrix for the refined lower and upper bounds 
%         of the unknown coefficients for the parameter space used
%         in the warm-start approach  

function bound  = get_bnd(y,x,beta0,bnd)

k=size(x,2)-1;

p_hat = 1./(1+exp(-x*logit(y,x)));

constr=repmat(p_hat-0.5,1,size(x,2)).*x;
bound=bnd;

model.sense = '>';
model.A = sparse(constr(:,2:size(x,2)));
model.rhs = -constr(:,1)*beta0;

tol=1e-6;
params.outputflag = 0; 
params.OptimalityTol=tol;
params.FeasibilityTol=tol;
params.IntFeasTol=tol;

model.lb = bound(:,1);
model.ub = bound(:,2); 

for i=1:k
       
objcoef=zeros(1,k);
objcoef(i)=1;
model.obj = objcoef;

model.modelsense = 'min';

try
    result= gurobi(model, params);
    bound(i,1)=result.objval;
    catch gurobiError
    fprintf('Error reported\n');
end

model.modelsense = 'max';
model.lb = bound(:,1);

try
    result = gurobi(model, params);
    bound(i,2)=result.objval;
    catch gurobiError
    fprintf('Error reported\n');
end

model.ub = bound(:,2);

end

end

