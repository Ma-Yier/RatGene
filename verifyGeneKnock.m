function [XX,input_xg] = verifyGeneKnock(input_model,input_x,indGPR,id_biomass,id_target,min_bound)
%UNTITLED6 此处显示有关此函数的摘要
% 

model=input_model;
input_xg=zeros(size(input_x));
for i=1:size(input_x,1)
   if input_x(i,1)>0.1
      input_xg(i,1)=1; 
   end
end
var_x=verifyRatioGene(model,input_xg);

%{
for i=1:size(var_x,1)
   if var_x(i)==1
      index=indGPR(i);
      model.lb(index)=0;
      model.ub(index)=0;
   end
end
%}
index=find(var_x==1);
model.lb(index)=0;
model.ub(index)=0;
model.lb(id_biomass)=0;

% Cplex
[X,FVAL,EXITFLAG]=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);

% Gurobi
%OPTIONS.Display='off';
%[X,FVAL,EXITFLAG]=LINPROG(-model.c,[],[],model.S,model.b,model.lb,model.ub,OPTIONS);

if EXITFLAG==1 
   if (-1)*FVAL>=min_bound
       XX=X(id_target);
       return;
   else
       XX=0;
   end
else
    X=zeros(size(model.rxns,1),1);
    XX=0;
    return;
end

% end funtion
end

