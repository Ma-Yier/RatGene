load('toyData/toyCore.mat');
load('toyData/toyEdge.mat');
%load('toyData/toyIntegrate.mat');

% integrate core and edge model
disp('------------------------------------------');
disp('Integrate core and edge model...');
model=uniteGeneModel(toyCore,toyEdge);
id_target=find(strcmp(model.rxns,'R11'));
id_biomass=find(model.c);
model.lb(id_biomass)=0.05;
fprintf('The id of biomass is %i \n',id_biomass);
fprintf('The id of target is %i \n',id_target);

% find strategy
disp('Solve the problem to find the strategy...');
TMGR=2; 
maxLoop=1000;
gap=2/model.lb(id_biomass)/maxLoop;
timeLimit=1000;
loop=10;
[x_target,alpha,knockout]=ratioGene(model,id_biomass,id_target,TMGR,maxLoop,gap,loop,timeLimit);
fprintf('Target reaction rate: %.2f \n',x_target);
disp('------------------------------------------');

% print strategy
rxnSplit=numel(toyCore.rxns);
geneSplit=numel(toyCore.genes);
delList=find(knockout(1:geneSplit)==0);
addList=find([zeros(geneSplit,1);knockout(geneSplit+1:numel(knockout))]);
fprintf('Deletion gene/(s) is/are: \n');
for i=1:numel(delList)
    fprintf(' ');
    disp(model.genes{delList(i),1});
end
fprintf('Addition gene/(s) is/are: \n');
for i=1:numel(addList)
    fprintf(' ');
    disp(model.genes{addList(i),1});
end
disp('------------------------------------------');

% reduce size
disp('Reduce size processing...')
[preNumKnockouts,afterNumKnockouts,newKnockouts]=rmRedundancy(model,knockout,id_target,id_biomass,[rxnSplit;geneSplit]);
fprintf('Reduce deletion size: %i-->%i \n',preNumKnockouts(1),afterNumKnockouts(1));
fprintf('Reduce addition size: %i-->%i \n',preNumKnockouts(2),afterNumKnockouts(2));
disp('------------------------------');
newDelList=find(newKnockouts(1:geneSplit)==0);
newAddList=find([zeros(geneSplit,1);newKnockouts(geneSplit+1:numel(newKnockouts))]);
fprintf('Reduced deletion gene/(s) is/are: \n');
for i=1:numel(newDelList)
    fprintf(' ');
    disp(model.genes{newDelList(i),1});
end
fprintf('Reduced addition gene/(s) is/are: \n');
for i=1:numel(newAddList)
    fprintf(' ');
    disp(model.genes{newAddList(i),1});
end
save('toyData/toyUnion.mat','model');
save('results_test3.mat');
