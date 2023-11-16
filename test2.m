% initilize cobratollbox
initCobraToolbox;

% set cplex as solver
changeCobraSolver('ibm_cplex');

% download data and load
websave('iMM904.mat','http://bigg.ucsd.edu/static/models/iMM904.mat');
websave('iJR904.mat','http://bigg.ucsd.edu/static/models/iJR904.mat');
load('iMM904.mat');
load('iJR904.mat');
model = uniteGeneModel(iMM904,iJR904);

% prepare important info
id_biomass=find(model.c);
metID=483;
id_carbon=find(strcmp('EX_glc__D_e',model.rxns));
id_oxygen=find(strcmp('EX_o2_e',model.rxns));

% test TMGR>0
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
fprintf('The metabolite is %s \n', model.metNames{metID});
disp("continue to RatGene...");

% execute RatGene
[x,knockouts,additions]=RatGene(iMM904,iMM904.mets{metID},'LBoxygen',-20, ...,
    'LBcarbon',-15,'timeLimit',500,'addition',iJR904);

% worst analysis
[new_model,id_target,TMPR] = introExchange(model1,id_biomass,[id_carbon,id_oxygen],metID);
kogene=[ones(size(iMM904.genes));zeros(numel(new_model.genes)-numel(iMM904.genes),1)];
lendel=numel(knockouts);
lenadd=numel(additions);
for i=1:lendel+lenadd
    if i<=lendel
        kogene(strcmp(knockouts{i,1},new_model.genes),1)=0;
    else
        kogene(strcmp(additions{i-lendel,1},new_model.genes),1)=1;
    end 
end
[minTarget,maxTarget] = boundAnalysis(new_model,kogene,id_biomass,id_target);

% output the results and save
fprintf('The number of gene knockouts is: %d \n', numel(knockouts));
fprintf('The number of gene additions is: %d \n', numel(additions));
fprintf('The target reaction rate with the strategy applied: %.4f mmol/gDW/h \n', x);
disp("--------------------");
fprintf('The target reaction rate in worst case: %.4f mmol/gDW/h \n', minTarget);
fprintf('The target reaction rate in best case: %.4f mmol/gDW/h \n', maxTarget);
save('results_test2.mat');
