k-Neighborhood Policy Updating for Influence Diagrams
=====================================================

This is a C++ implementation of the *k-neighborhood local search algorithms for selecting strategies in (limited memory) influence diagrams* described in

>    D.D. Maua and F.G. Cozman (2014), Speeding Up k-Neighborhood Local Search in Limited Memory Influence Diagrams. In _Proceedings of the Seventh European Conference on Probabilistic Graphical Models_, LNAI 8754, pp. 334-349.
>    D.D. Maua and F.G. Cozman (2015), _Fast local search methods for solving limited memory influence diagrams_. International Journal of Approximate Reasoning (in press).    
    
If you use this code in an academic work, please cite any of the the publications above.

As by-products, this package also implements the strategy selection algorithms described in

>    D.D. Maua, C.P. de Campos and M. Zaffalon (2012), _Solving Limited Memory Influence Diagrams_. Journal of Artificial Intelligence Research 44, pp. 97-140. 

as well as a variable elimination scheme for posterior probability computations in Bayesian networks represented in the UAI file format. With some minimal effort it is possible to use this code to implement the algorithms for MAP and credal network inferences described in

>   D.D. Maua and C.P. de Campos (2012), Anytime Marginal MAP Inference. In _Proceedings of the 28th International Conference on Machine Learning_, pp. 1471–1478.
>   D.D. Maua, C.P. de Campos and M. Zaffalon (2012), _Updating Credal Networks is Approximable in Polynomial Time_. International Journal of Approximate Reasoning 53(8), pp. 1183–1199.

LICENSE
-------
    
Copyright by Denis D. Maua (2014)

See the LICENSE file.

This package is released so that others can reproduce the experiments in the paper, and potentially use the algorithms in derivative work. You may contact-me with simple questions, but please do not expect to get real support from me.

INSTALLATION
------------

Type `make` from a command line in the project main  directory to compile. The executables will be stored in the directory `bin/`

USAGE
-----

Typing
   
  `bin/solve_limid`
   
from the project main directory will show usage on how to run kPU. This command reads limids in a special format described in next section.

FILE FORMAT
-----------

The diagrams can be specified in the following format inspired by UAI File Format:

  LIMID N M O SIZE1 ... SIZE(N+M) PA1 ... PA(N+M+O) CPT1 ... CPTN UTIL1 ... UTILO
  where:
  *  - LIMID is a file header string
  *  - N, M, and O denote, resp., the no. of chance, decision and value nodes.
  *  - SIZEi is an integer describing the cardinality of variable i (order is chance variables then decision variables)
  *  - PAi is NUMPA PA, where NUMPA is the number of parents of the variable and PA is a space-separated list of VAR_IDs for each of the parents. 
  *  - CPTi is NUMVALUES VALUE1 ... VALUE(NUMVALUES), with NUMVALUES numeric values describing the conditional probability values P(C_ID|PA_ID), arranged such that the least significant value (the one that increases the fastest) is C_ID and the most signficant value (the one that increases the slowest) is the value of the rightmost parent.
  *  - UTILi is NUMVALUES VALUE1 ... VALUE(NUMVALUES), with NUMVALUES numeric values specifying the utility function U(PA_ID) with least significant digit denoting the value of the leftmost parent.

 Linebreaks are discarded as white spaces, thus they can be used to improve readalibity. It is possible to have a C-style comment block (i.e., a block enclosed by /* and */) at the beginning (before LIMID appears); this will be simply ignored.
 
 See the examples in directory limids.
