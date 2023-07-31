function [x_target,alpha,knockout] = ratioGene(model,id_biomass,id_target,TMGR,max_loop,gap,numMultiStrat,timeLimit)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明

% load GPR rules
[lessMatrixLeft,lessMatrixRight,equalMatrixLeft,equalMatrixRight,nGpr,nAux,nGen,indGPR]=constructMatrix(model);


% init parameters
alpha=0;
n=size(model.S,2);
row_ratio=size(equalMatrixLeft,1)+1;
%y=size(lessMatrixLeft,2);

% construct A and b
A=lessMatrixLeft;
b=lessMatrixRight;
%A=[lessMatrixLeft,zeros(x,n); ...,
%   eye(n),zeros(n,y-n),-eye(n); ...,
%   -eye(n),zeros(n,y-n),-eye(n)];
%b=[lessMatrixRight;zeros(2*n,1)];

% construct Aeq and beq
ratioCons=zeros(1,n+nGpr+nAux+nGen);
ratioCons(1,id_target)=1;
ratioCons(1,id_biomass)=0;
Aeq=[equalMatrixLeft;ratioCons];
beq=[equalMatrixRight;0];
%Aeq=[equalMatrixLeft,zeros(size(equalMatrixLeft,1),n);zeros(1,y+n)];
%beq=[equalMatrixRight;0];

% set biomass o2 glc nonzero value
%Aeq(m+2,n+id_biomass)=1;
%Aeq(m+3,n+id_our)=1;
%Aeq(m+4,n+id_gur)=1;
%Aeq(m+5,n+id_target)=1;

% set lb and ub
lb=[model.lb;zeros(nGpr+nAux+nGen,1)];
ub=[model.ub;ones(nGpr+nAux+nGen,1)];

% set object function 
f=[zeros(n,1);TMGR*ones(nGpr,1);zeros(nAux+nGen,1)];
f(id_biomass)=-1;

% set milp labels-cplex/incons-gurobi
c_label=join(repmat("C",n,1),"");
b_label=join(repmat("B",nGpr+nAux+nGen,1),"");
ctype=char(c_label+b_label);
intcon=n+1:n+nGpr+nAux+nGen;

%options=cplexoptimset('cplex');
%options.mip.tolerances.integrality=10^(-12);
OPTIONS.Display='off';
OPTIONS.MaxTime=100;
OPTIONS.IntegerTolerance=10^(-9);
OPTIONS.NumSol=numMultiStrat;
OPTIONS.SearchMode=1;

% time limit
x_target=0;
knockout=[];
tStart=tic;
for i=1:max_loop    
    %ratio=ratio+gap;
    %lp.Aeq(colPoint,gid)=-ratio;
    
    % set ratio for each loop
    alpha=alpha+gap;
    Aeq(row_ratio,id_biomass)=-alpha;
    Aeq(row_ratio,id_target)=1;
    
    
    % solve milp
    %x=cplexmilp(f,A,b,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
    %[~,~,~,OUTPUT]=INTLINPROG(f,intcon,A,b,Aeq,beq,lb,ub,[],OPTIONS);
    %if ~isfield(OUTPUT,'pool')
    %   continue; 
    %end
    [~,pool]=CPLEXINTLINPROG(f,A,b,Aeq,beq,lb,ub,ctype);

    %[opt.x, opt.f, opt.stat, opt.output] = ...
    %    cplexmilp(lp.f, lp.A, lp.b, lp.Aeq, lp.beq,[],[],[],lp.lb, lp.ub,lp.ctype,[],options);

    % verify
    if ~isempty(pool)
        knockout=zeros(nGen,1);
        for j=1:size(pool,2)
            x=pool(:,j);
            if size(x,1)==n+nGpr+nAux+nGen
                [x_target,knockout_gene]=verifyGeneKnock(model,x(n+nGpr+nAux+1:n+nGpr+nAux+nGen,:),indGPR,id_biomass,id_target,model.lb(id_biomass));
                %x_target = verifyKnock_zero(model,x(1:n,:),id_biomass,id_target,1,0.05);
                if x_target>0
                    knockout(:,1)=knockout_gene;
                    return;
                else
                    if j==size(pool,2)
                        x_target=0;
                    end
                end
            else
                x_target=0;
            end 
        end
    end

    % time limit
    timerun=toc(tStart);
    if timerun>timeLimit
        break;
    end 
end



% end function
end



