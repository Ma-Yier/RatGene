# RatGene

Growth to production ratio-based design algorithms for constraint-based metabolic networks for growth-coupled production: RatGene


## Introduction

A ratio-based generalized approach that is capable of fulfilling a variety of modification criteria at both the gene and reaction levels. Additionally, RatGene also delivers the results with reduced-size strategies for the deletion or addition on both gene levels.


## Dependencies

+ MATLAB <=2021a
+ [IBM ILOG CPLEX](https://www.ibm.com/docs/en/icos/12.10.0?topic=SSSA5P_12.10.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_home.htm) <=12.10.0
+ [CobraToolBox](https://opencobra.github.io/cobratoolbox/stable/index.html)

## Instructions

RatGene adopts a ratio-based method to iteratively construct a series of mixed integer linear programming problems and solve them one by one to obtain candidate modification strategies for various networks until generating a solution that will satisfy the criterion of the problem setting after the validation procedures. Optionally, the output results could be processed by a two-step approach to reduce the size of the modification strategy after the validation procedures and the new trimmed-size strategy will satisfy the criterion as well.

### Function

```
[xTarget,varargout] = RatGene(model,targetMet,varargin)
```
#### INPUTS
   |Parameters| |
   |:---|:---|
   |*model*|      The same struct type as the .mat file downloaded from BiGG|  
   |*targetMet*|  The target metabolite, should be a cell or a char|  
   |biomass|    The biomass reaction in the model and it should be a cell or a char (default: the growth reaction in the model)|  
   |carbon|     The carbon source of the defined problem and it should be a cell or a char (default: EX_glc__D_e)|  
   |oxygen|     The input oxygen of the defined problem and it should be a cell or a char (default: EX_o2_e)|  
   |LBbiomass|  The lower thresold of the biomass reaction (default: 0.05)|  
   |LBcarbon|   The lower threshold of the input carbon source exchange reaction and it should be a negative value (default: -15)|  
   |LBoxygen|   The lower threshold of the input oxygen source exchange reaction and it should be a negative value (default: -15)|  
   |maxLoop|    The maximum iterations assigned to RatGene ratio-based procedures (default: 1000)|  
   |timeLimit|  The maximum computation time for the method (default: inf)|  
   |pool|       The number of solutions obtained from the IBM CPLEX         solution pool (default: 10)|        
   |type|       The type of modification strategy, 1 for gene and 0 for reaction (default: 1)|   
   |size|       Whether reduce the size of the output strategy, 1 for ture and 0 for false (default: 1)|   
   |addition|   Another model with the same struct type as the .mat file downloaded from BiGG to supply additional information if the problem is defined as deletion/addition problem.| 
 
NOTE: Parameters except for *model* and *targetMet* are optional inputs.    

#### OUTPUTS
   |Solutions| |
   |:---|:---|
   |*xTarget*|   The exchange reaction rate for the production of the target metabolite under the condition of applying the output modification strategy to the model.|  
   |Knockouts|  A cell of the name list of knockout strategy indicates which genes or reactions to be knocked out.|  
   |Additions|  A cell of the name list of addition strategy indicates which genes or reactions to be added up.|
  
NOTE: Solutions except for *xTarget* are optional outputs.
  
***
### Usage
+ `[~,knockouts] = RatGene(model,targetMet)` returns the deletion strategy of the `targetMet` in the `model`.
+ `[~,knockouts] = RatGene(model,targetMet,'biomass',biomass,'carbon',carbon,'oxygen',oxygen)` returns the deletion strategy of the `targetMet` in the `model` of which the biomass growth reaction, the carbon source exchange reaction and the oxygen source exchange reaction are assigned.
+ `[~,knockouts] = RatGene(model,targetMet,'LBbiomass',LBbiomass,'LBcarbon',LBcarbon,'LBoxygen',LBoxygen)` returns the deletion strategy of the `targetMet` in the `model` with the minimum threshold of three reactions assigned.
+ `[~,knockouts] = RatGene(model,targetMet,'maxLoop',maxLoop,'timeLimit',timeLimit,'pool',pool)` returns the deletion strategy of the `targetMet` in the `model` and the parameters limit the maximum number of loops, the maximum computation time, the number of pools, the type of modification strategyies and the excution of reduction approaches.
+ `[~,knockouts] = RatGene(model,targetMet,'type',type,'size',size)` returns the deletion strategy of the `targetMet` in the `model`. The parameter `type` resolves the gene modification strategy or the reaction modification strategy and `size` decides whether the function executes the reduction processes. Please note that the reduction processes are only available to gene modification strategies.
+ `[~,knockouts,additions] = RatGene(...,'addition',model2)` returns the deletion and addition strategy at the same time with various optional parameters as described above. `model2` is another model which supplies additional information for the addition strategy.
+ `[xTarget,...] = RatGene(...)` returns an extra solution which is the target reaction rate when the modification strategy is applied to the model with various optional parameters as described above.


***
### Running Examples
#### Example 1
The `test1.m` demostrates the process of computing a gene knockout strategy for the production of glycerol which is used as a solvent, sweetening agent, and also as medicine in yeast fermentation. The model used is [`iMM904`](http://bigg.ucsd.edu/models/iMM904) downloaded from [BiGG](http://bigg.ucsd.edu/) database. Run the following code in the command line of MATLAB:
```
test1
```
Then the results will be save in `results_test1.mat` and be printed as:
```
The TMGR is: 0.9662 mmol/gDW/h 
continue to RatGene...
The number of gene knockouts is: 528 
The target reaction rate with the strategy applied: 0.6583 mmol/gDW/h 
--------------------
The target reaction rate in worst case: 0.6583 mmol/gDW/h 
The target reaction rate in best case: 0.6583 mmol/gDW/h  
```
#### Example 2
The `test2.m` displays another instance of computing a gene modification strategy with both deletions and additions for the production of formate which is a potential building block for flavors and fragrances used in various industries by yeast fermentation. The two models used are [`iMM904`](http://bigg.ucsd.edu/models/iMM904) and[`iJR904`](http://bigg.ucsd.edu/models/iJR904). Run the following code in the command line of MATLAB:
```
test2
```
Then the results will be save in `results_test2.mat` and be printed as:
```
The metabolite is Formate  
continue to RatGene...
The number of gene knockouts is: 205 
The number of gene additions is: 351 
The target reaction rate with the strategy applied: 192.6837 mmol/gDW/h 
--------------------
The target reaction rate in worst case: 192.6837 mmol/gDW/h 
The target reaction rate in best case: 238.3193 mmol/gDW/h 
```
#### Example 3
The `test3.m` demostrates the process of integrating two models. The models used are the toy network models. Run the following code in the command line of MATLAB:
```
test3
```
Then the results will be save in `results_test3.mat` and be printed as:
```
------------------------------------------
Integrate core and edge model...
The id of biomass is 8 
The id of target is 9 
Solve the problem to find the strategy...
Target reaction rate: 1.00 
------------------------------------------
Deletion gene/(s) is/are: 
 gene1
 gene2
 gene3
 gene4
 gene6
 gene7
Addition gene/(s) is/are: 
 gene8
 gene9
------------------------------------------
Reduce size processing...
Reduce deletion size: 6-->5 
Reduce addition size: 2-->2 
------------------------------
Reduced deletion gene/(s) is/are: 
 gene1
 gene2
 gene4
 gene6
 gene7
Reduced addition gene/(s) is/are: 
 gene8
 gene9
```