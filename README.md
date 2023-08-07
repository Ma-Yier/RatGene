# RatGene

Growth to production ratio-based design algorithms for constraint-based metabolic networks for growth-coupled production: RatGene


## Introduction
A ratio-based generalized approach that is capable of fulfilling a variety of modification criteria at both the gene and reaction levels. Additionally, RatGene also delivers the results with reduced-size strategies for the deletion or addition on both gene levels.


## Dependencies

+ MATLAB <=2021a
+ [IBM ILOG CPLEX](https://www.ibm.com/docs/en/icos/12.10.0?topic=SSSA5P_12.10.0/ilog.odms.studio.help/Optimization_Studio/topics/COS_home.htm) <=12.10.0
+ [CobraToolBox](https://opencobra.github.io/cobratoolbox/stable/index.html)

## Instructions

### Running Examples
The `test1.m` demostrates the process of computing a gene knockout strategy for the production of ethanol which is a valuable biofuel and industrial solvent in yeast fermentation. The model used is `iMM904` downloaded from [BiGG](http://bigg.ucsd.edu/models/iMM904) database. Run the following code in the command line of MATLAB:
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