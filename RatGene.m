function [x_target,varargout] = RatGene(model,targetMet,varargin)
%RatGene adopts a ratio-based method to iteratively construct a series of
%mixed integer linear programming problems and solve them one by one to
%obtain candidate modification strategies for various networks until
%generating a solution that will satisfy the criterion of problem setting
%after the validation procedures. Optionally, the output results could be
%processed by a two-step approach to reduce the size of the modification
%strategy after the validation procedures and the new trimmed-size strategy
%will satisfy the criterion as well.
%
%function [x_target,varargout] = RatGene(model,targetMet,varargin)
%
%INPUTS
%   model      The same struct type as the .mat file downloaded from BiGG
%   targetMet  The target metabolite, should be a cell or a char
%
%OPTIONAL INPUTS
%   biomass    The biomass reaction in the model and it should be a cell
%              or a char (default: the growth reaction in the model)
%   carbon     The carbon source of the defined problem and it should be
%              a cell or a char (default: EX_glc__D_e)
%   oxygen     The input oxygen of the defined problem and it should be a
%              cell or a char (default: EX_o2_e)
%   LBbiomass  The lower thresold of the biomass reaction (default: 0.05)
%   LBcarbon   The lower threshold of the input carbon source exchange
%              reaction and it should be a negative value (default: -15)
%   LBoxygen   The lower threshold of the input oxygen source exchange
%              reaction and it should be a negative value (default: -15)
%   maxLoop    The maximum iterations assigned to RatGene ratio-based
%              procedures (default: 1000)
%   timeLimit  The maximum computation time for the method (default: inf)
%   pool       The number of solutions obtained from the IBM CPLEX solution
%              pool (default: 10)
%   type       The type of modification strategy, 1 for gene and 0 for
%              reaction (default: 1)
%   size       Whether reduce the size of the output strategy, 1 for ture
%              and 0 for false (default: 1)
%   addition   Another model with the same struct type as the .mat file
%              downloaded from BiGG to supply additional information if the
%              problem is defined as deletion/addition problem.
%
%OUTPUTS
%   x_target   The exchange reaction rate for the production of the target
%              metabolite under the condition of applying the output
%              modification strategy to the model.
%
%OPTIONAL OUTPUTS
%   Knockouts  A cell of the name list of knockout strategy indicates which
%              genes or reactions to be knocked out.
%   Additions  A cell of the name list of addition strategy indicates which
%              genes or reactions to be added up.
%
%
% August 1, 2023    Ma Yier
%


    
    % taeget id
    if iscell(targetMet)
        targetMet=targetMet{1,1};
    end
    if ~ischar(targetMet)
        error('targetMet should be either cell type or char type');
    end
    metID=find(strcmp(targetMet,model.mets));
    
    % set biomass id
    label=find(strcmp('biomass',varargin));
    if numel(label)==1
        biomassName=varargin{label(1)+1};
        if iscell(biomassName)
           biomassName=biomassName{1,1}; 
        end
        idSet=find(strcmp(biomassName,model.rxns));
        if isempty(idSet)
            error('no such biomass name in the input model');
        else
            id_biomass=idSet(1);
        end
    elseif numel(label)==0
        nonEmpty=find(~cellfun(@isempty,strfind(model.rxns,'BIOMASS')));
        if numel(nonEmpty)>1
            potential=model.rxns(nonEmpty);
            ni=find(~cellfun(@isempty,strfind(potential,'core')));
            id_biomass=nonEmpty(ni);
        elseif numel(nonEmpty)==1
            id_biomass=nonEmpty(1);
        else
            error('no biomass rxn in the model')
        end
        if id_biomass~=find(model.c)
            id_biomass=find(model.c); 
        end
    else
        error('biomass name error, please check the name');
    end
    
    
    
    % get carbon source id
    label=find(strcmp('carbon',varargin));
    if numel(label)==1
        carbonName=varargin{label(1)+1};
        if iscell(carbonName)
           carbonName=carbonName{1,1}; 
        end
        idSet=find(strcmp(carbonName,model.rxns));
        if isempty(idSet)
            error('no such carbon source name in the input model');
        else
            id_carbon=idSet(1);
        end
    elseif numel(label)==0
        nonEmpty=find(~cellfun(@isempty,strfind(model.rxns,'EX_glc__D_e')));
        if numel(nonEmpty)==1
            id_carbon=nonEmpty(1);
        else
            error('glucose name error in the model')
        end
    else
        error('carbon source name error, please check the name');
    end
    
    % get oxygen id
    label=find(strcmp('oxygen',varargin));
    if numel(label)==1
        oxygenName=varargin{label(1)+1};
        if iscell(oxygenName)
           oxygenName=oxygenName{1,1}; 
        end
        idSet=find(strcmp(oxygenName,model.rxns));
        if isempty(idSet)
            error('no such oxygen exchange name in the input model');
        else
            id_oxygen=idSet(1);
        end
    elseif numel(label)==0
        nonEmpty=find(~cellfun(@isempty,strfind(model.rxns,'EX_o2_e')));
        if numel(nonEmpty)==1
            id_oxygen=nonEmpty(1);
        else
            error('oxygen exchange name error in the model')
        end
    else
        error('oxygen name error, please check the name');
    end
    
    % set oxygen, carbon source, biomass threshold
    label=find(strcmp('LBbiomass',varargin));
    if ~isempty(label)
        LBbiomass=varargin{label(1)+1};
    else
        LBbiomass=0.05;
    end
    if LBbiomass<=0
       error('biomass lower bound should be posiive'); 
    end
    
    label=find(strcmp('LBcarbon',varargin));
    if ~isempty(label)
        LBcarbon=varargin{label(1)+1};
    else
        LBcarbon=-15;
    end
    
    label=find(strcmp('LBoxygen',varargin));
    if ~isempty(label)
        LBoxygen=varargin{label(1)+1};
    else
        LBoxygen=-15;
    end
    
    % set maximum loop and time limit
    label=find(strcmp('maxLoop',varargin));
    if ~isempty(label)
        maxLoop=varargin{label(1)+1};
    else
        maxLoop=1000;
    end
    
    label=find(strcmp('timeLimit',varargin));
    if ~isempty(label)
        timeLimit=varargin{label(1)+1};
    else
        timeLimit=inf;
    end
    
    % set pooling number
    label=find(strcmp('pool',varargin));
    if ~isempty(label)
        numMultiStrat=varargin{label(1)+1};
    else
        numMultiStrat=10;
    end
    
    % get type 1-gene/0-rxn
    label=find(strcmp('type',varargin));
    if ~isempty(label)
        type=varargin{label(1)+1};
    else
        type=1;
    end

    % get rmredundancy 1-cut/0-no cut
    label=find(strcmp('size',varargin));
    if ~isempty(label)
        cutSize=varargin{label(1)+1};
    else
        cutSize=1;
    end
    
    % get another input model
    label=find(strcmp('addition',varargin));
    if ~isempty(label)
        model2=varargin{label(1)+1};
        if ~isstruct(model2)
            error('second input model is not a sturct');
        end
        addition=1;
    else
        addition=0;
    end
         
    changeCobraSolver('ibm_cplex');
    % gene or rxn + addition        
    if type==1 && addition==0
        model1=model;
        model1.lb(id_biomass)=LBbiomass;
        model1.lb(id_carbon)=LBcarbon;
        model1.lb(id_oxygen)=LBoxygen;
        model1.ub(id_carbon)=0;
        model1.ub(id_oxygen)=0;
        [new_model,id_target,TMPR] = introExchange(model1,id_biomass,[id_carbon,id_oxygen],metID);
        if TMPR<=0
           error('theoretical maximum production rate is not positive'); 
        end
        opt=optimizeCbModel(model1);
        TMGR=opt.f;
        [x_target,~,knockout] = ratioGene(new_model,id_biomass,id_target,TMGR,maxLoop,TMPR/(maxLoop*LBbiomass),numMultiStrat,timeLimit);
        if cutSize==0
            newKnockout=knockout;
        elseif cutSize==1
            [~,~,newKnockout]=rmRedundancy(new_model,knockout,id_target,id_biomass,[0;0]);
        else
            error('cutsize can only be 1 or 0');
        end
        varargout{1}=new_model.genes(newKnockout==0);
    elseif type==1 && addition==1
        model1=uniteGeneModel(model,model2);
        model1.lb(id_biomass)=LBbiomass;
        model1.lb(id_carbon)=LBcarbon;
        model1.lb(id_oxygen)=LBoxygen;
        model1.ub(id_carbon)=0;
        model1.ub(id_oxygen)=0;
        [new_model,id_target,TMPR] = introExchange(model1,id_biomass,[id_carbon,id_oxygen],metID);
        if TMPR<=0
           error('theoretical maximum production rate is not positive'); 
        end
        opt=optimizeCbModel(model1);
        TMGR=opt.f;
        [x_target,~,knockout] = ratioGene(new_model,id_biomass,id_target,TMGR,maxLoop,TMPR/(maxLoop*LBbiomass),numMultiStrat,timeLimit);
        if cutSize==0
            newKnockout=knockout;
        elseif cutSize==1
            [~,~,newKnockout] = rmRedundancy(new_model,knockout,id_target,id_biomass,[length(model.rxns);length(model.genes)]);
        else
            error('cutsize can only be 1 or 0');
        end
        varargout{1}=model.genes(newKnockout(1:length(model.genes))==0);
        addGene=new_model.genes(length(model.genes)+1:length(new_model.genes));
        varargout{2}=addGene(newKnockout(length(model.genes)+1:length(new_model.genes))==1);
        
    % rxn del case
    elseif type==0
        if addition==1
            model1=uniteGeneModel(model,model2);
        end
        model1.lb(id_biomass)=LBbiomass;
        model1.lb(id_carbon)=LBcarbon;
        model1.lb(id_oxygen)=LBoxygen;
        model1.ub(id_carbon)=0;
        model1.ub(id_oxygen)=0;
        [new_model,id_target,TMPR] = introExchange(model1,id_biomass,[id_carbon,id_oxygen],metID);
        if TMPR<=0
           error('theoretical maximum production rate is not positive'); 
        end
        [~,~,newKnockout]=ratioMethod(new_model,id_biomass,id_target,maxLoop,TMPR/(maxLoop*LBbiomass),timeLimit);
        if addition==0
            varargout{1}=new_model.rxns(newKnockout==0);
        elseif addition==1
            varargout{1}=model.rxns(newKnockout(1:length(model.rxns))==0);
            if id_target==length(new_model.rxns)
                addRxn=new_model.rxns(length(model.rxns)+1:length(new_model.rxns)-1);
                varargout{2}=addRxn(newKnockout(length(model.rxns)+1:length(new_model.rxns)-1)==1);
            else
                addRxn=new_model.rxns(length(model.rxns)+1:length(new_model.rxns));
                varargout{2}=addRxn(newKnockout(length(model.rxns)+1:length(new_model.rxns))==1);
            end
        end
    else
        error('type can only be 1-gene or 0-reation');
    end
    

% end function
end

