function [model,id_target,TMPR] = introExchange(input_model,id_biomass,id_input,id_met)
%Detect the transport reaction of the target metabolite. If it exists,
%return the reacion id and TMPR. If it not, add a virtual transport 
%reaction for the target metabolite.
%   
%function [model,id_target,TMPR] = introExchange(input_model,id_biomass,id_input,id_met)
%
%INPUTS
%   input_model    The same struct type as the .mat file downloaded from BiGG
%   id_biomass     The id of biomass reaction
%   id_input       The matrix indicate carbon source id and oxygen source id
%   id_met         The id of target metabolite
%
%OUTPUTS
%   model      The input_model with an exchange reaction for the target met
%   id_target  The id of the target met exchanget reaction
%   TMPR       The theoretically maximum production rate for the target met
%
%
% July 31, 2023    Ma Yier
%

model=input_model;

% detect input reaction id(carbon source and oxygen)
num_input=size(id_input,2);
if num_input==1
    id_carbon=id_input(1);
    id_oxygen=0;
elseif num_input==2
    id_carbon=id_input(1);
    id_oxygen=id_input(2);
else
    disp('Input id error. At most 2 allowed.');
    id_target=0;
    TMPR=0;
    return;
end

% find if transport reaction exists or not
name_met=model.mets{id_met};
ex_rxn=char("EX_"+string(name_met));
dx_rxn=char("DM_"+string(name_met));

% exchange rxn is EX_ or DM_
if ismember(ex_rxn,model.rxns) || ismember(dx_rxn,model.rxns)
    
    % target id
    if ismember(ex_rxn,model.rxns)
        id_target=find(strcmp(ex_rxn,model.rxns));
    end
    if ismember(dx_rxn,model.rxns)
        id_target=find(strcmp(dx_rxn,model.rxns));
    end
    
    % target id cannot be input reactions
    if id_target==id_carbon || id_target==id_oxygen
       id_target=0;
       TMPR=0;
       return; 
    end
else
    
    % exchange rxn no exist
    id_target=size(model.rxns,1)+1;
    model.S=[model.S,zeros(size(model.mets,1),1)];
    model.S(id_met,id_target)=-1;
    model.rxns{id_target,1}='EX_target';
    model.grRules{id_target,1}='';
    %model.genes=[model.genes;'target_exchange'];
    %model.grRules=[model.grRules;'target_exchange'];
    model.lb(id_target)=-1000;
    model.ub(id_target)=1000;
    model.c=[model.c;0]; 
end

% calculate TMPR
model.c(id_biomass)=0;
model.c(id_target)=1;

% Cplex
x=cplexlp(-model.c,eye(size(model.rxns,1)),model.ub,model.S,model.b,model.lb,model.ub);

% Gurobi
%OPTIONS.Display='off';
%x=LINPROG(-model.c,eye(size(model.rxns,1)),model.ub,model.S,model.b,model.lb,model.ub,OPTIONS);

TMPR=x(id_target);
model.c(id_biomass)=1;
model.c(id_target)=0;
return;

% end function
end

