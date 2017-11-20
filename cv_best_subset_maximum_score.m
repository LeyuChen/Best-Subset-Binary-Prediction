
% This function computes the optimal q value among a range of q values
% specified by the vector q_range via the cross-validation procedure.
% The indices for the observations to be included in the training and testing 
% sample in the CV folds are specified by tr_ind and test_ind, respectively.


function [best_q, bhat,score,gap,rtime,ncount]=cv_best_subset_maximum_score(tr_ind,test_ind,data,focus_ind,aux_ind,beta0,q_range,T,tol,bnd,mio)

q_num = length(q_range); 
fold=size(tr_ind,2);
score=zeros(fold,q_num);
gap=zeros(fold,q_num);
rtime=zeros(fold,q_num);
ncount=zeros(fold,q_num);
bhat=zeros(length(focus_ind)+length(aux_ind)-1,fold,q_num);
val_score=zeros(q_num,1);

for q=1:q_num
 for i=1:fold
     disp(['(q,fold) : ' num2str(q_range(q)) num2str(i)]);
y=data(tr_ind(:,i),1);
datax=data(tr_ind(:,i),2:end);

[bhat(:,i,q),score(i,q),gap(i,q),rtime(i,q),ncount(i,q)]  = max_score_constr_fn(y,datax(:,focus_ind),datax(:,aux_ind),beta0,q_range(q),T,tol(q),bnd,mio);

y_v=data(test_ind(:,i),1);
datax_v=data(test_ind(:,i),2:end);
val_score(q) = val_score(q)+ mean(y_v == ((datax_v*[beta0;bhat(:,i,q)])>=0)); 
 end
 val_score(q)=val_score(q)/fold;
end
[~,best_q] = max(val_score);
end