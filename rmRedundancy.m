function [preNumKnockouts,afterNumKnockouts,newKnockouts] =  ...,
    rmRedundancy(model,knockouts,id_target,id_biomass,split)
%Apply two-step approach to reduce the size of strategies.
%
%function [preNumKnockouts,afterNumKnockouts,newKnockouts] =  ...,
%   rmRedundancy(model,knockouts,id_target,id_biomass,split)
%
%INPUTS
%   model        The same struct type as the .mat file downloaded 
%                from BiGG
%   id_biomass   The id of biomass reaction
%   id_target    The id of the target met exchanget reaction
%   knockouts    The obtained strategies after validation
%   split        The index indicate the endpoint id of original 
%                genes and reactions
%
%OUTPUTS
%   preNumKnockouts    The size of strategy before reduction
%   afterNumKnockouts  The size of strategy after reduction
%   newKnockouts       The new strategy with a reduced size
%
%
%July 31, 2023    Ma Yier
%

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
    
    
    % assign model to model1 and model2 when addition problem
    if split(1)>0
        model1=model;
        model2=model;
        model1.genes=model.genes(1:split(2),:);
        model1.grRules=model.grRules(1:split(1),:);
        model1.rxns=model.rxns(1:split(1),:);
        model2.genes=model.genes(split(2)+1:length(model.genes),:);
        model2.grRules=model.grRules(1+split(1):length(model.grRules),:);
        model2.rxns=model.rxns(1+split(1):length(model.rxns),:);
    else
       model1=model; 
    end
    

    if split(1)>0
        endRxn=split(1);
        endGene=split(2);
        additions=knockouts(split(2)+1:length(knockouts));
        knockouts=knockouts(1:split(2));
    else
        endRxn=numel(model1.rxns);
        endGene=numel(model1.genes);
    end
    % parse grRule and get rxns deletion sets
    %[parInfo,nRxn,nGen,nAux,nRelation,nEqual,gprLabel]=maxParseGPR(model1);
    [parInfo,~,~,~,~,~,~]=maxParseGPR(model1);
    knockouts(knockouts>0.1,1)=1;
    var_x=verifyRatioGene(model1,knockouts);

    % get to form mip problem from each rxn
    rxnNum=0;
    b=[];
    splitRxn=[];
    splitGene=[];
    splitAux=[];
    equalRxn=[];
    equalGene=[];
    for i=1:endRxn
        if var_x(i)==0 || isempty(model1.grRules{i,1})
           continue; 
        end
        if isempty(model1.grRules{i,1})
            continue;
        end
        rxnNum=rxnNum+1;
        
        if ~isempty(find(strcmp(model1.grRules{i,1},model1.genes),1))
            equalRxn=blkdiag(equalRxn,1);
            eg=zeros(1,numel(knockouts));
            eg(1,strcmp(model1.grRules{i,1},model1.genes))=-1;
            equalGene=[equalGene;eg];
            continue;
        end 

        canDel=getCandidateDel(parInfo{i,1},model1.genes(1:endGene,1),knockouts);
        if sum(canDel)==0
           error('error del'); 
        end

        
        if size(canDel,1)==1

           % form intermediate matrix
           interMatGene=[canDel;-1*canDel];
           interMatRxn=[-1;sum(canDel)];
           interMatb=[sum(canDel)-1;0];

           % form split matrix
           b=[b;interMatb];
           splitGene=[splitGene;interMatGene];
           splitAux=[splitAux;zeros(2,size(splitAux,2))];
           splitRxn=blkdiag(splitRxn,interMatRxn);

        else
           [matGene,matAux,matb]=formBoolMat(canDel);
           nAux=size(matAux,2);

           % form inermediate matrix 
           interMatGene=[zeros(2,size(matGene,2));matGene];
           interMatb=[0;0;matb];
           interMatAux=[-ones(1,nAux);ones(1,nAux);matAux];
           interMatRxn=[1;-1*nAux;zeros(size(matGene,1),1)];

           % form split matrix
           b=[b;interMatb];
           splitGene=[splitGene;interMatGene];
           splitAux=blkdiag(splitAux,interMatAux);
           splitRxn=blkdiag(splitRxn,interMatRxn);

        end
    end
    
    % prepare mip problem
    splitGene2=splitGene(:,knockouts==0);
    equalGene2=equalGene(:,knockouts==0);
    a=size(splitRxn,2)+size(equalRxn,2);
    bb=size(splitAux,2);
    c=size(splitGene2,2);
    
    % f A B Aeq beq lb ub ctype
    Aeq=[equalRxn,zeros(size(equalRxn,1),size(splitRxn,2)+bb),equalGene2];
    beq=zeros(size(equalRxn,1),1);
    A=[zeros(size(splitRxn,1),size(equalRxn,2)),splitRxn,splitAux,splitGene2];
    f=[zeros(a+bb,1);ones(c,1)];
    lb=[ones(a,1);zeros(bb+c,1)];
    ub=ones(a+bb+c,1);
    bLabel=join(repmat("B",a+bb+c,1),"");
    ctype=char(bLabel);
    
    % mip
    options=cplexoptimset('cplex');
    options.mip.tolerances.integrality=10^(-12);
    x=cplexmilp(f,A,b,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
    
    % parse result
    newDel=x(a+bb+1:a+bb+c,1);
    newDel=newDel-1;
    newDel(newDel==-1,:)=1;
    del=knockouts-1;
    del(del==-1,:)=1;
    del=diag(del);
    del=del(:,knockouts==0);
    newDel=del*newDel;
    newKnockouts=knockouts+newDel;
    
    % verify
    var_x2=verifyRatioGene(model1,newKnockouts);
    index1=find(var_x==1);
    index2=find(var_x2==1);
    if numel(find((index1==index2)==0))>0
        newKnockouts=knockouts;
    end
      
    preNumKnockouts=numel(find(knockouts==0));
    afterNumKnockouts=numel(find(newKnockouts==0));
    
%%    
    % reduce size for addition sets
    if split(1)==0
        return;
    end
    
    afterNumAdditions=-1;
    preNumAdditions=-1;
    newAdditions=inf;
    [parInfo,~,~,~,~,~,~]=maxParseGPR(model2);
    additions(additions>0.1,1)=1;
    var_a=verifyRatioGene(model2,additions);

    % get to form mip problem from each rxn
    rxnNumA=0;
    b=[];
    splitRxn=[];
    splitGene=[];
    splitAux=[];
    equalRxn=[];
    equalGene=[];
    for i=1:length(model2.rxns)
        if var_a(i)==1 || isempty(model2.grRules{i,1})
           continue; 
        end
        if isempty(model2.grRules{i,1})
            continue;
        end
        rxnNumA=rxnNumA+1;
        
        if ~isempty(find(strcmp(model2.grRules{i,1},model2.genes),1))
            equalRxn=blkdiag(equalRxn,1);
            eg=zeros(1,numel(additions));
            eg(1,strcmp(model2.grRules{i,1},model2.genes))=-1;
            equalGene=[equalGene;eg];
            continue;
        end 
        canAdd=getCandidateAdd(parInfo{i,1},model2.genes,additions);
        if sum(canAdd)==0
           error('error add'); 
        end
        
        if size(canAdd,1)==1

           % form intermediate matrix
           interMatGene=[canAdd;-1*canAdd];
           interMatRxn=[-1;sum(canAdd)];
           interMatb=[sum(canAdd)-1;0];

           % form split matrix
           b=[b;interMatb];
           splitGene=[splitGene;interMatGene];
           splitAux=[splitAux;zeros(2,size(splitAux,2))];
           splitRxn=blkdiag(splitRxn,interMatRxn);

        else
           [matGene,matAux,matb]=formBoolMat(canAdd);
           nAux=size(matAux,2);

           % form inermediate matrix 
           interMatGene=[zeros(2,size(matGene,2));matGene];
           interMatb=[0;0;matb];
           interMatAux=[-ones(1,nAux);ones(1,nAux);matAux];
           interMatRxn=[1;-1*nAux;zeros(size(matGene,1),1)];

           % form split matrix
           b=[b;interMatb];
           splitGene=[splitGene;interMatGene];
           splitAux=blkdiag(splitAux,interMatAux);
           splitRxn=blkdiag(splitRxn,interMatRxn);

        end
    end
    
    % prepare mip problem
    splitGene2=splitGene(:,additions==1);
    equalGene2=equalGene(:,additions==1);
    a=size(splitRxn,2)+size(equalRxn,2);
    bb=size(splitAux,2);
    c=size(splitGene2,2);
    
    % f A B Aeq beq lb ub ctype
    Aeq=[equalRxn,zeros(size(equalRxn,1),size(splitRxn,2)+bb),equalGene2];
    beq=zeros(size(equalRxn,1),1);
    A=[zeros(size(splitRxn,1),size(equalRxn,2)),splitRxn,splitAux,splitGene2];
    f=[zeros(a+bb,1);ones(c,1)];
    lb=[ones(a,1);zeros(bb+c,1)];
    ub=ones(a+bb+c,1);
    bLabel=join(repmat("B",a+bb+c,1),"");
    ctype=char(bLabel);
    
    % mip
    options=cplexoptimset('cplex');
    options.mip.tolerances.integrality=10^(-12);
    x=cplexmilp(f,A,b,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
    
    % parse result
    newAdd=x(a+bb+1:a+bb+c,1);
    newAdd=newAdd-1;
    newAdd(newAdd==-1,:)=1;
    add=additions;
    add=diag(add);
    add=add(:,additions==1);
    newAdd=add*newAdd;
    newAdditions=additions-newAdd;
    
    % verify
    var_a2=verifyRatioGene(model2,newAdditions);
    index1=find(var_a==1);
    index2=find(var_a2==1);
    if numel(find((index1==index2)==0))>0
        newAdditions=additions;
    end
      
    preNumAddtions=numel(find(additions==1));
    afterNumAdditions=numel(find(newAdditions==1));
    preNumDel=preNumKnockouts;
    aferNumDel=afterNumKnockouts;
    preNumKnockouts=[];
    afterNumKnockouts=[];
    preNumKnockouts=[preNumDel;preNumAddtions];
    afterNumKnockouts=[aferNumDel;afterNumAdditions];
    newKnockouts=[newKnockouts;newAdditions];
    
    
end    

%%
function [canAdd]=getCandidateAdd(parInfo,genes,add)
    % get basic information
    nRel=size(parInfo,1);
    
    % get 'and' grRule
    if string(parInfo{nRel,2})=="or"
        parUnit=parInfo{nRel,3};
        iter=0;
        for i=1:numel(parUnit)
           g=parUnit{i,1};
           if string(g(1:3))=="aux"
               getID=find(strcmp(g,parInfo(:,1)));
               a=getCandidateAdd(parInfo(1:getID,:),genes,add);
               iter=iter+1;
               canDelCell{iter,1}=a; 
           elseif add(strcmp(g,genes))==1
               a=zeros(1,size(add,1));
               a(find(strcmp(g,genes)))=1;
               iter=iter+1;
               canDelCell{iter,1}=a;
           end
        end
        canAdd=[];
        len=numel(add);
        for i=1:iter
            a=canDelCell{i,1};
            
            % drop zero option
            if numel(find(a(1,:)==1))==0
               continue; 
            end
            if len>numel(find(a(1,:)==1))
                canAdd=a;
                len=numel(find(a(1,:)==1));
            elseif len==numel(find(a(1,:)==1))
                canAdd=[canAdd;a];
            end   
        end
        
        % check if empty
        if isempty(canAdd)
           canAdd=zeros(1,numel(add)); 
        end
        
    % get 'or' grRule
    elseif string(parInfo{nRel,2})=="and"
        parUnit=parInfo{nRel,3};
        iter=0;
        andSet=zeros(1,size(add,1));
        sign=0;
        for i=1:numel(parUnit)     
            g=parUnit{i,1};
            if sign==1
                break;
            end
            
            % sign==0
            if string(g(1:3))=="aux"
                getID=find(strcmp(g,parInfo(:,1)));
                a=getCandidateAdd(parInfo(1:getID,:),genes,add);
                if  numel(find(a(1,:)==1))==0
                    andSet=zeros(1,size(add,1));
                    sign=1;
                else
                    iter=iter+1;
                    canDelCell{iter,1}=a; 
                end
            elseif add(strcmp(g,genes),1)==1
                andSet(strcmp(g,genes))=1;
            elseif add(strcmp(g,genes),1)==0
                andSet=zeros(1,size(add,1));
                sign=1;
            end
        end
        
        % merge and set
        mode=1;
        if exist('canDelCell','var')==0
           canDelCell=cell(0,0); 
        end
            
        if sign==1
            canAdd=zeros(1,numel(add)); 
        else
            canAdd=mergeSet(andSet,canDelCell,iter,mode);
        end
    else
        canAdd=zeros(1,numel(add));
    end

end


%%
function [canDel]=getCandidateDel(parInfo,genes,del)
    % get basic information
    nRel=size(parInfo,1);
    
    % get 'and' grRule
    if string(parInfo{nRel,2})=="and"
        parUnit=parInfo{nRel,3};
        iter=0;
        for i=1:numel(parUnit)
           g=parUnit{i,1};
           if string(g(1:3))=="aux"
               getID=find(strcmp(g,parInfo(:,1)));
               a=getCandidateDel(parInfo(1:getID,:),genes,del);
               iter=iter+1;
               canDelCell{iter,1}=a; 
           elseif del(find(strcmp(g,genes)),1)==0
               a=zeros(1,size(del,1));
               a(find(strcmp(g,genes)))=1;
               iter=iter+1;
               canDelCell{iter,1}=a;
           end
        end
        canDel=[];
        len=numel(del);
        for i=1:iter
            a=canDelCell{i,1};
            
            % drop zero option
            if numel(find(a(1,:)==1))==0
               continue; 
            end
            if len>numel(find(a(1,:)==1))
                canDel=a;
                len=numel(find(a(1,:)==1));
            elseif len==numel(find(a(1,:)==1))
                canDel=[canDel;a];
            end   
        end
        
        % check if empty
        if isempty(canDel)
           canDel=zeros(1,numel(del)); 
        end
        
    % get 'or' grRule
    elseif string(parInfo{nRel,2})=="or"
        parUnit=parInfo{nRel,3};
        iter=0;
        andSet=zeros(1,size(del,1));
        sign=0;
        for i=1:numel(parUnit)     
            g=parUnit{i,1};
            if sign==1
                break;
            end
            
            % sign==0
            if string(g(1:3))=="aux"
                getID=find(strcmp(g,parInfo(:,1)));
                a=getCandidateDel(parInfo(1:getID,:),genes,del);
                if  numel(find(a(1,:)==1))==0
                    andSet=zeros(1,size(del,1));
                    sign=1;
                else
                    iter=iter+1;
                    canDelCell{iter,1}=a; 
                end
            elseif del(find(strcmp(g,genes)),1)==0
                andSet(find(strcmp(g,genes)))=1;
            elseif del(find(strcmp(g,genes)),1)==1
                andSet=zeros(1,size(del,1));
                sign=1;
            end
        end
        
        % merge and set
        mode=1;
        if exist('canDelCell','var')==0
           canDelCell=cell(0,0); 
        end
            
        if sign==1
            canDel=zeros(1,numel(del)); 
        else
            canDel=mergeSet(andSet,canDelCell,iter,mode);
        end
    else
        canDel=zeros(1,numel(del));
    end

end

%%
function [canDel]=mergeSet(andSet,canDelCell,iter,mode)
    if mode~=0 && mode ~=1
        error('mode only 1 or 0');
    end 

    if ~isempty(find(andSet==1))
        iter=iter+1;
        canDelCell{iter,1}=andSet;
    end        

    for i=2:iter
        if i==2
            set1=canDelCell{i-1,1};
            set2=canDelCell{i,1};
            canDel=crossSet(set1,set2,mode);
        else
            set1=canDel;
            set2=canDelCell{i,1};
            canDel=crossSet(set1,set2,mode);
        end
    end
    
    if mode==1 && iter>1
        len=numel(andSet);
        lenSet=zeros(size(canDel,1),1);
        for i=1:size(canDel,1)
            lenSet(i,1)=sum(canDel(i,:));
            if len>lenSet(i,1)
                len=lenSet(i,1); 
            end
        end
        canDel=canDel(lenSet==len,:);
    else
        canDel=canDelCell{1,1};
    end
    

end

function [finalSet]=crossSet(set1,set2,mode)
% mode=0 get shortest between two sets
% mode=1 not filter shortest set between two sets

    len=size(set1,2);
    finalSet=[];
    if mode==0
        for i=size(set1,1)
           for j=size(set2,1)
               a=set1(i,:)+set2(j,:);
               b=zeros(1,len);
               b(1,find(a>0))=1;
               if numel(find(b==1))<len
                   len=numel(find(b==1));
                   finalSet=b;
               elseif numel(find(b==1))==len
                   c=finalSet*b.';
                   if isempy(find(c==len))
                       finalSet=[finalSet;b];
                   else
                       continue;
                   end
               else
                   continue;
               end
           end
        end
    elseif mode==1
        lenSet=[];
        for i=size(set1,1)
           for j=size(set2,1)
               a=set1(i,:)+set2(j,:);
               b=zeros(1,len);
               b(1,find(a>0))=1;
               lenb=numel(find(b==1));
               cmpSet=finalSet(find(lenSet==lenb),:);
               if i==1 && j==1
                   finalSet=[finalSet;b];
                   lenSet=[lenSet;lenb];
                   continue;
               end
               
               if isempty(cmpSet)
                  cmpSet=0; 
               end
               
               if ~isempty(find(cmpSet*b.'==lenb))
                   continue;
               else
                   finalSet=[finalSet;b];
                   lenSet=[lenSet;lenb];
               end
           end
        end

    end

end

%%
function [matE,matC,matR]=formBoolMat(canDel)
    [numRela,numElem]=size(canDel);
    matE=zeros(2*numRela,numElem);
    matC=zeros(2*numRela,numRela);
    matR=zeros(2*numRela,1);
    for i=1:numRela
        matE(2*i-1,canDel(i,:)==1)=1;
        matE(2*i,canDel(i,:)==1)=-1;
        matC(2*i-1,i)=-1;
        matC(2*i,i)=sum(canDel(i,:));
        matR(2*i-1,1)=sum(canDel(i,:))-1;
        matR(2*i,1)=0;
    end
end










