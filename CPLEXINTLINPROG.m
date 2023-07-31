function [cplex,x] = CPLEXINTLINPROG(f,A,b,Aeq,beq,lb,ub,ctype)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明


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

