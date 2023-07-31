function [parInfo,nRxn,nGen,nAux,nRelation,nEqual,gprLabel] = maxParseGPR(model)
%Analyze the GPR rules and parse the GPR rules to linear relationships.
%
%function [parInfo,nRxn,nGen,nAux,nRelation,nEqual,gprLabel] = maxParseGPR(model)
%
%INPUTS
%   model    The same struct type as the .mat file downloaded from BiGG
%
%OUTPUTS
%   parInfo    Decomposed clauses
%   nRxn    Number of reactions in the model
%   nGen    Number of genes in the model
%   nAux    Number of auxiliary clauses for boolean functions
%   nRelation Number of single boolean functions
%   nEqual    Number of reactions with single genes
%   gprLabel  Indicate a reactions with a GPR rule or without a GPR rule
%
%
% July 31, 2023    Ma Yier
%

% data scale 
nRxn=size(model.rxns,1);
nGen=size(model.genes,1);
nAux=1;
parInfo=cell(nRxn,1);
nRelation=0;
nEqual=0;
gprLabel=ones(nRxn,1);

% get GPR rules, gene information, 

for i=1:nRxn
    gpr=model.grRules{i,1};
    
    % has a GPR rule
    if ~isempty(gpr)
        % not a single gene
        if contains(gpr,' ')
            gprStr=string(gpr);
            gprStr="(" + gprStr + ")";
            gprChar=char(gprStr);
            maxUnit=size(strfind(gprStr,"("),2);
            nRelation=nRelation+2*maxUnit;
            parInfoUnit=cell(maxUnit,3);
            unitLabel=1;

            % iterately parsing a GPR rule
            while unitLabel<=maxUnit
               lbracket=strfind(gprStr,"(");
               rbracket=strfind(gprStr,")");            
               rpoint=rbracket(1);
               lpoint=find(rpoint>lbracket);           
               lpoint=lbracket(lpoint(size(lpoint,2)));
               gprUnit=gprChar(lpoint+1:rpoint-1);
               gprStr=replaceBetween(gprStr,lpoint,rpoint,"aux"+string(nAux));
               gprChar=char(gprStr);
               
               % assign the one-to-one gene
               if unitLabel==maxUnit
                   parInfoUnit{unitLabel,1}=char("real"+string(i));
               else
                   parInfoUnit{unitLabel,1}=char("aux"+string(nAux));
                   nAux=nAux+1;
               end

               % deal with single and/or
               gprUnitElement=split(gprUnit,' ');
               gprUnitEleNum=size(gprUnitElement,1);
               parInfoUnit{unitLabel,2}=gprUnitElement{2,1};
               for j=2:uint8((gprUnitEleNum+1)/2)
                  gprUnitElement(j,:)=[];                  
               end
               parInfoUnit{unitLabel,3}=gprUnitElement;
               unitLabel=unitLabel+1;
            end
        % a single gene        
        else
            parInfoUnit=cell(1,3);
            parInfoUnit{1,1}=char("real"+string(i));
            parInfoUnit{1,2}=[];
            parInfoUnit{1,3}=gpr;
            nEqual=nEqual+1;
        end
        
    % has no GPR rule  
    else
        parInfoUnit=cell(1,3);
        gprLabel(i,1)=0;
    end
    parInfo{i,1}=parInfoUnit; 
end
nAux=nAux-1;

% end function
end

