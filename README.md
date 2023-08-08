# RatGene

Growth to production ratio-based design algorithms for constraint-based metabolic networks for growth-coupled production: RatGene


## Introduction

A ratio-based generalized approach that is capable of fulfilling a variety of modification criteria at both the gene and reaction levels. Additionally, RatGene also delivers the results with reduced-size strategies for the deletion or addition on both gene levels.


## Dependencies

+ MATLAB <=2021a
+ [IBM ILOG CPLEX](https://www.ibm.com/docs/en/icos/12.10.0?topic=SSSA5P_12.10.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_home.htm) <=12.10.0
+ [CobraToolBox](https://opencobra.github.io/cobratoolbox/stable/index.html)

## Instructions

RatGene adopts a ratio-based method to iteratively construct a series of mixed integer linear programming problems and solve them one by one to obtain candidate modification strategies for various networks until generating a solution that will satisfy the criterion of problem setting after the validation procedures. Optionally, the output results could be processed by a two-step approach to reduce the size of the modification strategy after the validation procedures and the new trimmed-size strategy will satisfy the criterion as well.

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

#### OUTPUTS
   |Solutions| |
   |:---|:---|
   |*xTarget*|   The exchange reaction rate for the production of the target metabolite under the condition of applying the output modification strategy to the model.|  
   |Knockouts|  A cell of the name list of knockout strategy indicates which genes or reactions to be knocked out.|  
   |Additions|  A cell of the name list of addition strategy indicates which genes or reactions to be added up.|
  
### Usage

### Running Examples
#### Example 1
The `test1.m` demostrates the process of computing a gene knockout strategy for the production of ethanol which is a valuable biofuel and industrial solvent in yeast fermentation. The model used is [`iMM904`](http://bigg.ucsd.edu/models/iMM904) downloaded from [BiGG](http://bigg.ucsd.edu/) database. Run the following code in the command line of MATLAB:
```
test1
```
And the results will be save in `results_test1.mat` and be printed as:
```
The TMGR is: 0.9662 mmol/gDW/h 
continue to RatGene...
The number of gene knockouts is: 455 
The target reaction rate with the strategy applied: 11.4718 mmol/gDW/h 
```
#### Example 2
The `test2.m` displays another instance of computing a gene modification strategy with both deletions and additions for the production of ethanol in yeast fermentation. The two models used are [`iMM904`](http://bigg.ucsd.edu/models/iMM904) and[`iJR904`](http://bigg.ucsd.edu/models/iJR904). Run the following code in the command line of MATLAB:
```
test2
```
And the results will be save in `results_test2.mat` and be printed as:
```
The TMGR is: 0.9662 mmol/gDW/h 
continue to RatGene...
The number of gene knockouts is:  
The number of gene additions is: 
The target reaction rate with the strategy applied:  mmol/gDW/h 
```