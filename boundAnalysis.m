function [minTarget,maxTarget] = boundAnalysis(model,kogene,id_biomass,id_target)
%Analyze the minimum bound and maximum bound of target reaction under a
%condistion that a certain modification strategy is applied to the model
%and the biomass growth reaction is maximized.
%
%function [minTarget,maxTarget] ...,
%   = boundAnalysis(model,kogene,id_biomass,id_target)
%
%INPUTS
%   model       The same struct type as the .mat file downloaded from BiGG
%   kogene      A modification strategy indicates which genes to be added
%               or deleted.
%   id_biomass  Biomass growth reaction ID
%   id_target   Target metabolite production reaction ID
%
%OUTPUTS
%   minTarget   Minimum reaction rate when biomass reaction is maximized
%               in another words, the worst case.
%   maxTarget   Maximum reaction rate when biomass reaction is maximized
%               in another words, the best case.
%
%
% October 5, 2023   Ma Yier
%


% init return variables
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
[~,FVAL1]=cplexlp(-model2.c,[],[],model2.S,model2.b,model2.lb,model2.ub);
maxTarget=(-1)*FVAL1;
if maxTarget<=0.001
    maxTarget=0;
end

% worst case analysis
[~,minTarget]=cplexlp(model2.c,[],[],model2.S,model2.b,model2.lb,model2.ub);
if minTarget<=0.001
    minTarget=0;
end

% end function
end

