function [cplex,x] = CPLEXINTLINPROG(f,A,b,Aeq,beq,lb,ub,ctype)
%Form IBM CPLEX Cplex class for an MILP problem and solve it with pooling
%by IBM CPLEX solver.
%
%function [cplex,x] = CPLEXINTLINPROG(f,A,b,Aeq,beq,lb,ub,ctype)
%
%INPUTS
%   f    Objective function
%   A    Inequality left part
%   b    Inequality right part
%   Aeq  Equality left part
%   beq  Equality right part
%   lb   Lower bound of variables
%   ub   Upper bound of variables
%   ctype Type of variables - continues or integer or binary
%
%OUTPUTS
%   cplex    IBM CPLEX cplex class with results of solver computation
%   x        The solutions if have. [] for no solution after computation
%
%
% July 31, 2023    Ma Yier
%


model.obj=f;
model.lb=lb;
model.ub=ub;
model.A=[sparse(A);sparse(Aeq)];
model.lhs=[-inf(size(A,1),1);beq];
model.rhs=[b;beq];
model.ctype=ctype;
model.sense='minimize';
cplex=Cplex(model);
cplex.Model=model;
%cplex.Param=cplexoptimset('cplex');
% set params
cplex.Param.timelimit.Cur=100;
cplex.Param.mip.limits.populate.Cur=10;
cplex.DisplayFunc=[];
cplex.Param.mip.tolerances.integrality.Cur=10^(-12);

%no clonelog 
cplex.Param.output.clonelog.Cur=0;
%cplex.solve();
cplex.populate();

numsol=size(cplex.Solution.pool.solution,1);
if numsol>0
    x=zeros(size(cplex.Solution.pool.solution(1).x,1),numsol);
    for i=1:numsol
        x(:,i)=cplex.Solution.pool.solution(i).x;   
    end
else
    x=[];
end




% end function
end

