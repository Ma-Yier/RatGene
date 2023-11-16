% initilize cobratollbox
%initCobraToolbox;

% set cplex as solver
changeCobraSolver('ibm_cplex');

% download data and load
websave('iMM904.mat','http://bigg.ucsd.edu/static/models/iMM904.mat');
load('iMM904.mat');
model=iMM904;

% prepare important info
id_biomass=find(model.c);
metID=533;
id_carbon=find(strcmp('EX_glc__D_e',model.rxns));
id_oxygen=find(strcmp('EX_o2_e',model.rxns));

% test TMGR>0
model1=model;
model1.lb(id_oxygen)=-20;
model1.lb(id_carbon)=-15;
opt=optimizeCbModel(model1);
if opt.stat~=1
    error('TMGR does not exist');
end
if opt.f<=0
    error('TMGR does not satisfy');
end
fprintf('The TMGR is: %.4f mmol/gDW/h \n', opt.f);
fprintf('continue to RatGene...\n');

% execute RatGene
[x,knockouts]=RatGene(model,model.mets{metID},'LBoxygen',-20, ...,
    'LBcarbon',-15,'timeLimit',500);

% worst analysis
[new_model,id_target,TMPR] = introExchange(model1,id_biomass,[id_carbon,id_oxygen],metID);
kogene=ones(size(new_model.genes));
for i=1:numel(knockouts)
    kogene(strcmp(knockouts{i,1},new_model.genes),1)=0;
end
[minTarget,maxTarget] = boundAnalysis(new_model,kogene,id_biomass,id_target);

% output the results and save
fprintf('The number of gene knockouts is: %d \n', numel(knockouts));
fprintf('The target reaction rate with the strategy applied: %.4f mmol/gDW/h \n', x);
disp("--------------------");
fprintf('The target reaction rate in worst case: %.4f mmol/gDW/h \n', minTarget);
fprintf('The target reaction rate in best case: %.4f mmol/gDW/h \n', maxTarget);
%save('results_test1.mat');
