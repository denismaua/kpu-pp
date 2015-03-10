k-Neighborhood Policy Updating for Influence Diagrams
=====================================================

This is a C++ implementation of the *k-neighborhood local search algorithms for selecting strategies in (limited memory) influence diagrams* described in

>    D.D. Maua and F.G. Cozman (2014), Speeding Up k-Neighborhood Local Search in Limited Memory Influence Diagrams. In _Proceedings of the Seventh European Conference on Probabilistic Graphical Models_, LNAI 8754, pp. 334-349.
    
If you use this code in an academic work, please cite the publication above.

As a by-product, this package also implements a variable elimination scheme for posterior probability computations in Bayesian networks.

LICENSE
-------
    
Copyright by Denis D. Maua (2014)

See the LICENSE file.

This package is released so that others can reproduce the experiments in the paper, and potentially use the algorithms in derivative work. You may contact-me with simple questions, but please do not expect to get real support from me.

INSTALLATION
------------

Go to directory src/ and type make to compile. The executables will be stored in the directory bin/

USAGE
-----

Typing
   
  `bin/mpu`
   
from the project main directory will run show usage on how to run kPU.

The chain diagrams used in the experiments in the paper are in the compressed file chain.tar.gz inside the directory limids.

The command above reads limids in a special format described next.

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

 Linebreaks are discarded as white spaces, thus they can be used to improve readalibity.
 
 See the examples in directory limids.
