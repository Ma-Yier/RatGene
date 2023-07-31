function [genesid] = findDiffGene(model,diff_rxns);
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

label=1;

for i=1:size(diff_rxns,1)
   grRule=model.grRules{diff_rxns(i)};
   if ~isempty(grRule)
       if contains(grRule,' ')
           grRule=replace(grRule,'(','');
           grRule=replace(grRule,')','');
           glist=split(grRule);
           gn=size(glist,1);
           for j=1:gn
              if strcmp(glist{j},'and') || strcmp(glist{j},'or')
                 continue; 
              else
                  genes{label,1}=glist{j};
                  label=label+1;
              end
           end
           
       else
           genes{label,1}=grRule;
           label=label+1;
       end
   end
end

point=1;
while point<=size(genes,1)
    pattern=genes{point,1};
    genesid(point,1)=find(strcmp(model.genes,pattern));
    repid=find(strcmp(genes,pattern));
    repid(1,:)=[];
    genes(repid,:)=[];
    point=point+1;  
end

genesid=sort(genesid);

% end function
end

