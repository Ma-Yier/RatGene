% initilize cobratollbox
initCobraToolbox;

% set cplex as solver
changeCobraSolver('ibm_cplex');

% download data and load
websave('iMM904.mat','http://bigg.ucsd.edu/static/models/iMM904.mat');
websave('iJR904.mat','http://bigg.ucsd.edu/static/models/iJR904.mat');
load('iMM904.mat');
load('iJR904.mat');
model=iMM904;

% tset TMGR>0
model1=model;
model1.lb(find(strcmp('EX_o2_e',model.rxns)))=-20;
model1.lb(find(strcmp('EX_glc__D_e',model.rxns)))=-15;
opt=optimizeCbModel(model1);
if opt.stat~=1
    error('TMGR does not exist');
end
if opt.f<0
    error('TMGR does not satisfy');
end
fprintf('The TMGR is: %.4f mmol/gDW/h \n', opt.f);
fprintf('continue to RatGene...\n');

% execute RatGene
[x,knockouts,additions]=RatGene(model,model.mets{608},'LBoxygen',-20, ...,
    'LBcarbon',-15,'timeLimit',500,'addition',iJR904);

% output the results and save
fprintf('The number of gene knockouts is: %d \n', numel(knockouts));
fprintf('The number of gene additions is: %d \n', numel(additions));
fprintf('The target reaction rate with the strategy applied: %.4f mmol/gDW/h \n', x);
save('results_test2.mat');
