
% Given x and beta0, this function solves the following maximization problem
% for each i, max |beta0*x(i,1)+x(i,:)*t| over t confined to the space
% described by bnd

% function input :
% x     : (n by k) matrix of covariate data 
% beta0 : the coefficient for the first covariate in x       
% bnd   : ((k-1) by 2) matrix where the first and second columns  
%         respectively store the lower and upper bounds 
%         of the unknown coefficients

% function output :
% value : the value of the maximized objective function

function value = miobnd_fn(x,beta0,bnd)

n=size(x,1);
k=size(x,2)-1;

model.modelsense = 'max';
model.sense = '>';
model.lb = bnd(:,1);
model.ub = bnd(:,2);

tol=1e-6;

params.outputflag = 0; 
params.OptimalityTol=tol;
params.FeasibilityTol=tol;
params.IntFeasTol=tol;

v=zeros(2,1);

value=zeros(n,1);

for i=1:n
    
alpha =  beta0*x(i,1);

model.obj = x(i,2:k+1);
model.objcon = alpha;

try
    model.A = sparse(x(i,2:k+1));
    model.rhs = -alpha;
    result= gurobi(model, params);
    v(1)=result.objval;
    catch gurobiError
    fprintf('Error reported\n');
end

model.obj = -x(i,2:k+1);
model.objcon = -alpha;

try
    model.A = sparse(-x(i,2:k+1));
    model.rhs = alpha;
    result = gurobi(model, params);
    v(2)=result.objval;
    
  catch gurobiError
    fprintf('Error reported\n');
end
value(i)=max(v);

end

end

