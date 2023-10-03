function [minTarget,maxTarget] = boundAnalysis(model,kogene,id_biomass,id_target)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

% init
minTarget=0;
maxTarget=0;

% get rxn knock
var_x=verifyRatioGene(model,kogene);
model1=model;

% delete rxn
index=find(var_x==1);
model1.lb(index)=0;
model1.ub(index)=0;
model1.lb(id_biomass)=0;

% get max biomass and fix
[~,FVAL]=cplexlp(-model1.c,[],[],model1.S,model1.b,model1.lb,model1.ub);
model2=model1;
model2.lb(id_biomass)=(-1)*FVAL;
model2.ub(id_biomass)=(-1)*FVAL;
model2.c(id_biomass)=0;
model2.c(id_target)=1;

% bet case analysis
[~,FVAL1]=cplexlp(-model1.c,[],[],model1.S,model1.b,model1.lb,model1.ub);
maxTarget=(-1)*FVAL1;
if maxTarget<=0.001
    maxTarget=0;
end

% worst case analysis
[~,minTarget]=cplexlp(model1.c,[],[],model1.S,model1.b,model1.lb,model1.ub);
if minTarget<=0.001
    minTarget=0;
end

% end function
end

