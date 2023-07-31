function [XX,knockout] = verifyKnock(input_model,input_x,id_biomass,id_target,min_bound)
%Verify whether the reaction deletion or deletion/addition strategy is valid.
%
%function [XX,knockout] = verifyKnock(input_model,input_x,id_biomass,id_target,min_bound)
%
%INPUTS
%   input_model    The same struct type as the .mat file downloaded from BiGG
%   input_x        The reaction modification strategy
%   id_biomass     The id of biomass reaction
%   id_target      The id of the target met exchange reaction
%   min_bound      Minimum threshold for the biomass growth reaction
%
%OUTPUTS
%   XX        The target reaction rate after the modification
%   knockout  The modification strategy
%
%
% July 31, 2023    Ma Yier
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

