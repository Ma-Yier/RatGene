function [XX,knockout] = verifyKnock(input_model,input_x,id_biomass,id_target,min_bound)
%UNTITLED6 此处显示有关此函数的摘要
% 

model=input_model;
var_x=input_x;
var_x(abs(var_x)<0.0000001)=0;
model.lb(var_x==0)=0;
model.ub(var_x==0)=0;
model.lb(id_biomass)=0;

% knockout
knockout=ones(size(input_x));
knockout(var_x==0)=0;

% Cplex
[X,FVAL,EXITFLAG]=cplexlp(-model.c,[],[],model.S,model.b,model.lb,model.ub);

% Gurobi
%OPTIONS.Display='off';
%[X,FVAL,EXITFLAG]=LINPROG(-model.c,[],[],model.S,model.b,model.lb,model.ub,OPTIONS);

if EXITFLAG==1 
   if (-1)*FVAL>min_bound
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

