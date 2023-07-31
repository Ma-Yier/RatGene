function [objLoc] = locateGene(nRxn,nGpr,nAux,indRea,objGene,model)
%Find the location for a gene in the variable.
%
%function [objLoc] = locateGene(nRxn,nGpr,nAux,indRea,objGene,model)
%
%INPUTS
%   nRxn    Number of reactions in the model
%   nGpr    Number of reactions with GPR rules
%   nAux    Number of auxiliary clauses for boolean functions
%   indRea  Index of reactions that trim reactions having no grRules
%   objGene The name of object gene
%   model   The same struct type as the .mat file downloaded from BiGG
%
%OUTPUTS
%   objLoc    The location of the gene
%
%
% July 31, 2023    Ma Yier
%

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

