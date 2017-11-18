
function [y,datax] = simulation_data(n,beta,sigma,type)
k=length(beta);
reg = mvnrnd(zeros(n,k-1),sigma,n);
datax=[reg(:,1) ones(n,1) reg(:,2:end)];
if type == 1 % heteroskedasticity
   z=reg(:,1)+reg(:,2);
   e=0.25*(ones(n,1) + 2*(z.^2) + z.^4).*randn(n,1);  
else % homoskedasticity
   e=0.25*randn(n,1);
end
y=((datax*beta)>=e);
end