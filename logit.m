
% This function computes the estimates of the logit regression of y on x. 

function b = logit(y,x)

cnv = 0;

b0 = zeros(size(x,2),1); 


while cnv == 0

ind = x*b0;
P = exp(ind)./(1+exp(ind));

grd = sum(repmat(y-P,1,size(x,2)).*x);
hes = -x'*(repmat(P.*(1-P),1,size(x,2)).*x);

b1 = b0 - inv(hes)*grd';

dev = max(abs(b1-b0));

if dev < 1e-8 
cnv = 1;
end

b0 = b1;

end

b=b1;
end

