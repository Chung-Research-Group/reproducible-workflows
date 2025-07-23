clc
format long

profile=parcluster('Processes');
profile.NumWorkers=20;
saveProfile(profile);

addpath('./Collections/NGPM -- A NSGA-II Program in Matlab v1.4/')

N = 30 ;
type = 'EconomicEvaluation' ;

Function = @(x) PSACycleSimulation_together_with_cost( x, type, N ) ; % Function to simulate the PSA cycle

options            = nsgaopt();                          % create default options structure
options.popsize    = 90;                                 % populaion size
options.maxGen     = 60;                                 % max generation

options.vartype    = [1, 1, 1, 1, 1, 1, 1, 1, 1] ;
options.outputfile = 'high_19_modified_Economic_30_data_99Pu_BFG.txt'         ;

options.numObj  = 1 ;                                    % number of objectives
options.numVar  = 9 ;                                    % number of design variables
options.numCons = 2 ;                                    % number of constraints
options.lb      = [1e5,  10, 0.01, 0.1, 0.01, 0.1e5, 10, 10, 2]  ;               % lower bound of x
options.ub      = [3e5, 500, 0.99, 1, 0.99, 1e5, 100, 100, 9]  ;               % upper bound of x
options.nameObj = {'Cost'} ;           % the objective names are showed in GUI window.
options.objfun  = Function                   ;           % objective function handle

options.useParallel = 'yes' ;                            % parallel computation is non-essential here
options.poolsize     = 20   ;                            % number of worker processes

result = nsga2(options)     ;                            % begin the optimization!
