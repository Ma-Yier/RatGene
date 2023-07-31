function [preNumKnockouts,afterNumKnockouts,newKnockouts] = cutKnockouts(model,knockouts,id_target,id_biomass)
% model should have exchange_targetRxn

afterNumKnockouts=-1;
preNumKnockouts=-1;
newKnockouts=inf;
% verify gene or rxn knockouts
if numel(knockouts)==size(model.genes,1)
    %disp('gene knouckouts');
    mode=1;
elseif numel(knockouts)==size(model.rxns,1)
    disp('rxn knockouts');
    mode=0;
    return
else
    disp('neither rxn knockout nor gene knockout');
    return;
end

% gene knockouts 
if mode==1
    preNumKnockouts=numel(find(knockouts==0));
    var_rxn=verifyRatioGene(model,knockouts);
    reservRxns=find(var_rxn==0);
    delRxns=find(var_rxn==1);
    reservGenes=find(knockouts==1);
    delGenes=find(knockouts==0);
    
    % get biomassValue and targetValue
    model1=model;
    model1.lb(delRxns)=0;
    model1.ub(delRxns)=0;
    model1.lb(id_biomass)=0;
    model1.c(id_biomass)=1;
    model1.c(id_target)=0;
    OPTIONS.Display='off';
    [X,~,EXITFLAG]=cplexlp(-model1.c,[],[],model1.S,model1.b,model1.lb,model1.ub);
    if EXITFLAG==1
        biomassVal=X(id_biomass);
        targetVal=X(id_target);
    else
        disp('no feasible solution');
        return
    end
   
    % formate matrix
   [lessMatrixLeft,lessMatrixRight,equalMatrixLeft,equalMatrixRight,nGpr,nAux,nGenes,indGPR]=constructMatrix(model);
   nVar=size(lessMatrixLeft,2);
   nRxns=nVar-nAux-nGpr-nGenes;
   
   % set biomass and target lb
   biotarMatrixLeft=zeros(2,nVar);
   biotarMatrixRight=zeros(2,1);
   biotarMatrixLeft(1,id_biomass)=-1;
   biotarMatrixLeft(2,id_target)=-1;   
   biotarMatrixRight(1,1)=-biomassVal;
   biotarMatrixRight(2,1)=-targetVal;
   
   % set rxn deletion number constraint
   %nRxnsConstraintLeft=zeros(1,nVar);
   %nRxnsConstraintRight=zeros(1,1);
   %nRxnsConstraintRight(1,1)=numel(delRxns);
   
   
   % set rxn reserve/deletion number constraint
   reservRxnLeft=zeros(2,nVar);
   reservRxnRight=zeros(2,1);
   rxnsIndex=zeros(nRxns,1);
  
   for i=1:nGpr
       rxnsIndex(indGPR(i),1)=i;
   end
   for i=1:numel(delRxns)
       reservRxnLeft(1,nRxns+rxnsIndex(delRxns(i),1))=1;
   end
   for i=1:nGpr
       if reservRxnLeft(1,nRxns+i)==0
           reservRxnLeft(2,i)=1;
       end
   end
   reservRxnRight(2,1)=sum(reservRxnLeft(2,:));
   
   % set gene deletion number constraint
   nGenesConstraintLeft=zeros(1,nVar);
   nGenesConstraintRight=zeros(1,1);
   nGenesConstraintRight(1,1)=numel(delGenes);
   for i=1:nGenesConstraintRight(1,1)
      nGenesConstraintLeft(1,nRxns+nGpr+nAux+delGenes(i))=1; 
   end
   
   % set revserve gene constraint
   reservGeneLeft=zeros(1,nVar);
   reservGeneRight=zeros(1,1);
   reservGeneRight(1,1)=numel(reservGenes);
   for i=1:reservGeneRight(1,1)
       reservGeneLeft(1,nRxns+nGpr+nAux+reservGenes(i))=1;
   end
   
   % construct problem
   Aeq=[equalMatrixLeft;reservGeneLeft;reservRxnLeft];
   beq=[equalMatrixRight;reservGeneRight;reservRxnRight];
   A=[lessMatrixLeft;biotarMatrixLeft;nGenesConstraintLeft];
   b=[lessMatrixRight;biotarMatrixRight;nGenesConstraintRight];
   lb=[model.lb;zeros(nGpr+nAux+nGenes,1)];
   ub=[model.ub;ones(nGpr+nAux+nGenes,1)];
   c_gene=zeros(nGenes,1);
   c_gene(delGenes,1)=1;
   c=[zeros(nRxns+nGpr+nAux,1);c_gene];
   c_label=join(repmat("C",nRxns,1),"");
   b_label=join(repmat("B",nGpr+nAux+nGenes,1),"");
   ctype=char(c_label+b_label);
   
   
   [X1,~,EXITFLAG1]=cplexmilp(-c,A,b,Aeq,beq,[],[],[],lb,ub,ctype);
   if EXITFLAG1==1
       newKnockouts=X1(nRxns+nGpr+nAux+1:nRxns+nGpr+nAux+nGenes,1);
       [x_target,~]=verifyGeneKnock(model,newKnockouts,indGPR,id_biomass,id_target,1,0.05);
       if x_target>0
           afterNumKnockouts=numel(find(newKnockouts==0));
       else
           newKnockouts=knockouts;
           afterNumKnockouts=preNumKnockouts;
       end
   else
       newKnockouts=knockouts;
       afterNumKnockouts=preNumKnockouts;
   end
end


% end fucntion
end

