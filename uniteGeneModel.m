function [universal_model] = uniteGeneModel(model1,model2)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
% model1 core
% model2 out space


% get scale info of models
[m,n]=size(model2.S);
a=size(model2.genes,1);
[x,y]=size(model1.S);
b=size(model1.genes,1);
metPair=zeros(m,2);
rxnPair=zeros(n,2);
genePair=zeros(a,2);

% cmp mets diff
for i=1:m
    for j=1:x
        if string(model2.metNames{i,1})==string(model1.metNames{j,1})
            metPair(i,:)=[i,j];
            break;
        end
    end
end


% cmp rxns diff
for i=1:n
    for j=1:y
        if string(model2.rxnNames{i,1})==string(model1.rxnNames{j,1})
            rxnPair(i,:)=[i,j];
            break;
        end
    end
end

% cmp genes diff
for i=1:a
    for j=1:b
        if string(model2.genes{i,1})==string(model1.genes{j,1})
            genePair(i,:)=[i,j];
            break;
        end
    end
end


% get diff mets and rxns info
diff_mets=find(metPair(:,1)==0);
diff_rxns=find(rxnPair(:,1)==0);
%diff_genes=find(genePair(:,1)==0);
num_diff_mets=size(diff_mets,1);
num_diff_rxns=size(diff_rxns,1);
%num_diff_genes=size(diff_genes,1);
same_mets=find(metPair(:,1));
same_rxns=find(rxnPair(:,1));
%same_genes=find(genePair(:,1));
num_same_mets=size(same_mets,1);
num_same_rxns=size(same_rxns,1);
%num_same_genes=size(same_genes,1);


% get diff gene info
genesid=findDiffGene(model2,diff_rxns);
diff_genes=genesid;
num_diff_genes=size(diff_genes,1);
same_genes=find(genePair(:,1));
num_same_genes=size(same_genes,1);


% construct new model
universal_model=model1;
diff_matrix=model2.S(:,diff_rxns);

% add left matrix
rightMatrix=zeros(x,num_diff_rxns);
for i=1:num_same_mets
    id=same_mets(i);
    rightMatrix(metPair(id,2),:)=diff_matrix(metPair(id,1),:);
end

% under zero matrix
underMatrix=zeros(num_diff_mets,y);

% right under diff met diff rxn matrix
rightUnderMatrix=model2.S(diff_mets,diff_rxns);

% S matrix of universal model
universal_model.S=[model1.S,rightMatrix;underMatrix,rightUnderMatrix];

% update all quant info
universal_model.lb=[model1.lb;model2.lb(diff_rxns)];
universal_model.ub=[model1.ub;model2.ub(diff_rxns)];
universal_model.b=[model1.b;zeros(num_diff_mets,1)];
universal_model.c=[model1.c;zeros(num_diff_rxns,1)];
universal_model.rev=[model1.rev;model2.rev(diff_rxns)];

% update all mets info
universal_model.mets=[model1.mets;model2.mets(diff_mets)];
universal_model.metNames=[model1.metNames;model2.metNames(diff_mets)];
universal_model.metFormulas=[model1.metFormulas;model2.metFormulas(diff_mets)];
if isfield(model1,'metCharge') && isfield(model2,'metCharge')
    universal_model.metCharge=[model1.metCharge;model2.metCharge(diff_mets)];
end

% update all rxns info
universal_model.grRules=[model1.grRules;model2.grRules(diff_rxns)];
universal_model.rxns=[model1.rxns;model2.rxns(diff_rxns)];
universal_model.rxnNames=[model1.rxnNames;model2.rxnNames(diff_rxns)];
universal_model.subSystems=[model1.subSystems;model2.subSystems(diff_rxns)];

% update gene info
universal_model.genes=model1.genes;
universal_model.genes=[universal_model.genes;model2.genes(diff_genes)];


% update rxnGeneMat info
southMatrix=sparse(num_diff_rxns,b);
diffRxnMat=model2.rxnGeneMat(diff_rxns,:);
for i=1:num_same_genes
    id=same_genes(i);
    southMatrix(:,genePair(id,2))=diffRxnMat(:,genePair(id,1));
end

eastMatrix=sparse(y,num_diff_genes);
diffGeneMat=model2.rxnGeneMat(:,diff_genes);
for i=1:num_same_rxns
    id=same_rxns;
    eastMatrix(rxnPair(id,2),:)=diffGeneMat(rxnPair(id,1),:);
end

southEastMatrix=model2.rxnGeneMat(diff_rxns,diff_genes);
universal_model.rxnGeneMat=[model1.rxnGeneMat,eastMatrix;southMatrix,southEastMatrix];

% update description
universal_model.description=char("United models. Core: "+string(model1.description) ...
    +", Addition: "+string(model2.description));

% end function
end