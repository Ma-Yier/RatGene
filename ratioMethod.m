function [x_target,alpha,knockout] = ratioMethod(model,id_biomass,id_target,max_loop,gap,timeLimit)
%Contruct the problem and solve to obtain reaction deletion or reaction deletion/addition
%strategies.
%
%function [x_target,alpha,knockout] = ratioMethod ...,
%         (model,id_biomass,id_target,max_loop,gap,timeLimit)
%
%INPUTS
%   model        The same struct type as the .mat file downloaded from BiGG
%   id_biomass   The id of biomass reaction
%   id_target    The id of the target met exchange reaction
%   max_loop     Maximum number of iterations
%   gap          Change of the value of alpha in each loop
%   timeLimit    Time limit for the computation
%
%OUTPUTS
%   x_target    The reaction rate of target reaction after modification
%   alpha       The value of alpha when the strategy is obtained
%   knockout    The obtained strategies after validation
%
%
% July 31, 2023    Ma Yier
%

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

