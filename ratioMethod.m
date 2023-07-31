function [x_target,alpha,knockout] = ratioMethod(model,id_biomass,id_target,max_loop,gap,timeLimit)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明

alpha=0;
[m,n]=size(model.S);
knockout=[];
x_target=0;

tStart=tic;
for i=1:max_loop
    alpha=alpha+gap;
    ratio=zeros(1,size(model.rxns,1));
    ratio(1,id_biomass)=-alpha;
    ratio(1,id_target)=1;
    Aeq=[[model.S;ratio],zeros(m+1,n)];
    beq=[model.b;0];
    A=[eye(n),-eye(n);-eye(n),-eye(n)];
    b=zeros(2*n,1);
    c=[zeros(n,1);ones(n,1)];
    lb=[model.lb;zeros(n,1)];
    ub=[model.ub;999999*ones(n,1)];
    
    % Cplex
    %x=cplexlp(c,A,b,Aeq,beq,lb,ub);
    
    % Gurobi
    OPTIONS.Display='off';
    x=LINPROG(c,A,b,Aeq,beq,lb,ub,OPTIONS);
    
    if size(x,1)==2*n
        [x_target,tknockout] = verifyKnock(model,x(1:n,:),id_biomass,id_target,model.lb(id_biomass));
        %x_target = verifyKnock_zero(model,x(1:n,:),id_biomass,id_target,1,0.05);
        if x_target>0
            knockout=tknockout;
            break;
        end
    else
        x_target=0;
    end
    
    timesum=toc(tStart);
    if timesum>timeLimit
        break;
    end
end



% end function
end

