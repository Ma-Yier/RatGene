function [objLoc] = locateGene(nRxn,nGpr,nAux,indRea,objGene,model)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

% get name length
geneNameLen=size(objGene,2);

% if name is aux.. real.. or 
if geneNameLen>2
    
    % get first 3 alpha
    geneNameHead=objGene(1:3);
    switch geneNameHead

        % real case 
        case 'rea'
            ind=indRea(str2double(objGene(5:geneNameLen)));
            objLoc=ind+nRxn;

        % aux case
        case 'aux'
            objLoc=str2double(objGene(4:geneNameLen))+nRxn+nGpr;

        % gene case
        otherwise
            objLoc=find(contains(model.genes,objGene))+nRxn+nGpr+nAux;

    end
    
else
    
    % gene case
    objLoc=find(contains(model.genes,objGene))+nRxn+nGpr+nAux;    
end

% end fucntion
end

