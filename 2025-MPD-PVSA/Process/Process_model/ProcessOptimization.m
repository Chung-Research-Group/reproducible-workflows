clc
format long

profile=parcluster('Processes');
profile.NumWorkers=20;
saveProfile(profile);

addpath('./Collections/NGPM -- A NSGA-II Program in Matlab v1.4/')

N = 30 ;
type = 'ProcessEvaluation';

Function = @(x) PSACycleSimulation(x, type, N) ; % Function to simulate the PSA cycle

options         = nsgaopt() ;                            % create default options structure
options.popsize = 80        ;                            % populaion size
options.maxGen  = 80        ;                            % max generation

options.vartype    = [2, 2, 1, 1, 1, 2, 2, 2]         ;
options.outputfile = 'AFG_1_2comp_10_90_pu_rec_MPD_real_vec_organized_version.txt' ;

options.numObj  = 2 ;                                    % number of objectives
options.numVar  = 8 ;                                    % number of design variables
options.numCons = 3 ;                                    % number of constraints
options.lb      = [2e5,  10, 0.01, 0.1, 0.01, 1e4, 10, 10]  ;               % lower bound of x
options.ub      = [10e5, 500, 0.99, 0.8, 0.99, 2e5, 50, 50]  ;             % upper bound of x
options.nameObj = {'Obj1','Obj2'} ;               % the objective names are showed in GUI window.
options.objfun  = Function               ;               % objective function handle

options.useParallel = 'yes' ;                            % parallel computation non essential here
options.poolsize     = 20   ;                            % number of worker processes

result = nsga2(options)     ;                            % begin the optimization!
