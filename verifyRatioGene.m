function [geneKnock] = verifyRatioGene(model,geneValue)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

% get gr rules and genes
grRules=model.grRules;
genes=model.genes;
num_gene=size(geneValue,1);
num_rxn=size(grRules,1);
geneKnock=zeros(num_rxn,1);

% replace and/or with */+
grRules=strrep(grRules,'or','+');
grRules=strrep(grRules,'and','*');

% replace gene to the value
[genes1,index]=sortrows(genes,'descend'); % some gene has similar name like YER060W_A and YER060W, YER060W_A need to be repleced first

for i=1:num_gene
   grRules=strrep(grRules,genes1{i,1},num2str(geneValue(index(i),1))); 
end

for j=1:num_rxn
   gpr=grRules{j,1};
   if ~isempty(gpr)
      if eval(gpr)<0.9
          geneKnock(j,1)=1;
      end
   end
end


% end function
end
