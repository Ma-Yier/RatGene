function [lessMatrixLeft,lessMatrixRight,equalMatrixLeft,equalMatrixRight,nGpr,nAux,nGen,indGPR] = constructMatrix(model)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明


% obtain parsing GPR info
[parInfo,nRxn,nGen,nAux,nRelation,nEqual,gprLabel]=maxParseGPR(model);
%[parInfo,nRxn,nGen,nAux,nRelation,nEqual,gprLabel]=minParseGPR(model);

% get exist gene-rxn info 
indGPR=find(gprLabel);  % index of rxn has grRules
indRea=zeros(size(gprLabel,1),1);
nGpr=size(indGPR,1);  % number of rxn has grRules
for i=1:nGpr
    indRea(indGPR(i),1)=i; % rxn-index in list that trim rxns having no grRules
end


% init matrix
nMet=size(model.S,1);

% and/or complex relation
gprMatrixLeft=zeros(nRelation,nRxn+nGpr+nGen+nAux);
gprMatrixRight=zeros(nRelation,1);

% lb*y<=v<=ub*y
gprLabelDiag=diag(gprLabel);
grCorrelationRxn=gprLabelDiag(indGPR,:);
boundMatrixLeft=[grCorrelationRxn,-diag(model.ub(indGPR,:)),zeros(nGpr,nGen+nAux); ...,
                -grCorrelationRxn,diag(model.lb(indGPR,:)),zeros(nGpr,nGen+nAux)];
boundMatrixRight=zeros(2*nGpr,1);

% one gene control
gprEqualLeft=zeros(nEqual,nRxn+nGpr+nGen+nAux);
gprEqualRight=zeros(nEqual,1);

% S*v=0
fluxEqualLeft=[model.S,zeros(nMet,nGpr+nGen+nAux)];
fluxEqualRight=zeros(nMet,1);

% row point when add GPR relation
rowPoint1=1;
rowPoint2=1;

% construct GPR matrix
for i=1:size(parInfo,1)
    
   % get each reaction GPR parsing info 
   parInfoUnit=parInfo{i,1};
   
   % has GPR
   if ~isempty(parInfoUnit{1,1})
       maxUnit=size(parInfoUnit,1);
       
       % single gene control
       if (maxUnit==1)&&isempty(parInfoUnit{1,2})
           geneID=find(contains(model.genes,parInfoUnit{1,3}));
           gprEqualLeft(rowPoint1,nRxn+indRea(i))=1;
           gprEqualLeft(rowPoint1,nRxn+nGpr+nAux+geneID)=-1;
           rowPoint1=rowPoint1+1;
           %gprMatrixLeft(rowPoint,nRxn+i)=1;
           %gprMatrixLeft(rowPoint,2*nRxn+nAux+geneID)=-1;
           %rowPoint=rowPoint+1;
           %gprMatrixLeft(rowPoint,nRxn+i)=-1;
           %gprMatrixLeft(rowPoint,2*nRxn+nAux+geneID)=1;
           %rowPoint=rowPoint+1;      

       % multiple gene control
       else

           % iteratly get each auxilary GPR
           for j=1:maxUnit
               eachAuxGPR=parInfoUnit(j,:);
               logicRelation=eachAuxGPR{1,2};
               n=size(eachAuxGPR{1,3},1);
               eachAuxGPR_ctrlGene=eachAuxGPR{1,3};
               switch logicRelation
                   case 'and'

                       % obj gene
                       objLoc=locateGene(nRxn,nGpr,nAux,indRea,eachAuxGPR{1,1},model);
                       gprMatrixLeft(rowPoint2,objLoc)=-1;
                       gprMatrixLeft(rowPoint2+1,objLoc)=n;
                       gprMatrixRight(rowPoint2,1)=n-1;

                       % ctrl gene
                       for k=1:n
                           objLoc=locateGene(nRxn,nGpr,nAux,indRea,eachAuxGPR_ctrlGene{k,1},model);
                           gprMatrixLeft(rowPoint2,objLoc)=1;
                           gprMatrixLeft(rowPoint2+1,objLoc)=-1;
                       end
                       rowPoint2=rowPoint2+2;

                   case 'or'

                       % obj gene
                       objLoc=locateGene(nRxn,nGpr,nAux,indRea,eachAuxGPR{1,1},model);
                       gprMatrixLeft(rowPoint2,objLoc)=-n;
                       gprMatrixLeft(rowPoint2+1,objLoc)=1;

                       % ctrl gene
                       for k=1:n
                           objLoc=locateGene(nRxn,nGpr,nAux,indRea,eachAuxGPR_ctrlGene{k,1},model);
                           gprMatrixLeft(rowPoint2,objLoc)=1;
                           gprMatrixLeft(rowPoint2+1,objLoc)=-1;
                       end
                       rowPoint2=rowPoint2+2;

                   otherwise
                       error('parsing GPR error')
               end   
           end
       end    
   end
end

% construct matrix
lessMatrixLeft=[boundMatrixLeft;gprMatrixLeft];
lessMatrixRight=[boundMatrixRight;gprMatrixRight];
equalMatrixLeft=[fluxEqualLeft;gprEqualLeft];
equalMatrixRight=[fluxEqualRight;gprEqualRight];

%objFunction=[transpose(model.c),zeros(1,nRxn+nAux+nGen)];

% end funciton
end

