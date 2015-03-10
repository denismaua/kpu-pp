// Copyright (c) 2014 Denis Maua
// All Rights Reserved.
//
// This file is part of MSP library
//
// MSP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MSP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with MSP.  If not, see <http://www.gnu.org/licenses/>.

/** Describes Limid class interface, which implements the data structure for a limited-memory influence diagram (limid). */

#ifndef MSP_LIMID_H
#define MSP_LIMID_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "variable.h"
#include "factor.h"
//#include "dag.h"

namespace msp {


  enum limid_solver_t {SPU_solver, kPU_solver, MPU_solver}; // define types of solvers available

  /** Limid interface */
  class Limid
  {

  private:
    unsigned N, M, O; /*< number of chance, decision and value variables, resp. */
    std::vector<Variable > _variables; /*< vector of (chance, decision and value) variables */
    std::vector<Factor > _factors; /*< vector of factors representing CPTs, (randomized) policies and utility functions */
    //std::unordered_map<Variable*, std::unordered_set<Variable* > > _dag; /*< DAG represented as a map of parent sets */
    //std::unordered_map<Variable*, std::unordered_set<Variable* > > _dag; /*< DAG represented as a map of parent sets */

  public:

    /** data structure to store a LIMID solution */
    struct solution_t
    {
      std::string method;
      unsigned num_iters;
      unsigned num_tables;
      double value;
    };
    

    /** Default constructor. */
    Limid () : N(0), M(0), O(0) { }

    unsigned num_chance_nodes() { return N; } /*< number of chance nodes */
    unsigned num_decision_nodes() { return M; } /*< number of decision nodes */
    unsigned num_value_nodes() { return O; } /*< number of value nodes */

    /** Find an elimination ordering for variables.
     *
     * @param verbosity a positive integer specifying how much output to produce (0 means no output).
     * @return a vector of variables such that the leftmost variable is to be eliminated first.
     */
    std::vector<Variable> find_order(unsigned verbosity=0);


    /** Solves the limid. 
     *
     * Finds a strategy by kneighborhood search and computes its expected value. 
     *
     * @param k an unsigned integer describing the size of the neighborhood in search
     * @param verbosity a positive integer specifying how much output to produce (0 means no output).
     * @return a solution structure containing the expected value of an optimal strategy and statistics of the used method.
     */
    solution_t solve(unsigned k, unsigned verbosity=0);

    /** Single Policy Updating (SPU) algorithm with Pareto-dominance pruning.
     *
     * Implements a 1-neighborhood local search on the space of policies, using set-valued variable elimination with Pareto-dominance pruning.
     *
     * @param order a vector of variables defining a variable elimination sequence.
     * @param iter an unsigned used to store the number of iterations until convergence.
     * @param numtables an unsigned that will be used to store that maximum number of tables generated during the procedure.
     * @param verbosity a positive integer specifying how much output to produce (0 means no output).
     * @return a double containing the expected value of a locally optimal strategy.
     */
    double SPU(std::vector<Variable>& order, unsigned& iter, unsigned& numtables, unsigned verbosity=0);
    double SPU2(std::vector<Variable>& order, unsigned& iter, unsigned& numtables, unsigned verbosity=0);

    /** k-Policy Updating (kPU) algorithm.
     *
     * Implements a k-neighborhood local search on the space of policies, using set-valued variable elimination with Pareto-dominance pruning.
     *
     * @param order a vector of variables defining a variable elimination sequence.
     * @param k the number of variables in the search neighborhood.
     * @param iter an unsigned used to store the number of iterations until convergence.
     * @param numtables an unsigned that will be used to store that maximum number of tables generated during the procedure.
     * @param verbosity a positive integer specifying how much output to produce (0 means no output).
     * @return a double containing the expected value of a locally optimal strategy.
     */
    double kPU(std::vector<Variable>& order, unsigned k, unsigned& iter, unsigned& numtables, unsigned verbosity=0);
    double kPU2(std::vector<Variable>& order, unsigned k, unsigned& iter, unsigned& numtables, unsigned verbosity=0);

    /** Multiple Policy Updating (MPU) algorithm.
     *
     * Computes the maximum expected utility by set-valued variable elimination with Pareto-dominance pruning.
     *
     * @param order a vector of variables defining a variable elimination sequence.
     * @param numtables an unsigned used to store the maximum number of tables generated during the procedure.
     * @param verbosity a positive integer specifying how much output to produce (0 means no output).
     * @return a double containing the expected value of an optimal strategy.
     */
    double MPU(std::vector<Variable>& order, unsigned& numtables, unsigned verbosity=0);
    double MPU2(std::vector<Variable>& order, unsigned& numtables, unsigned verbosity=0);

    /** Read limited-memory influence diagram model from file.
     *
     *  File format is as follows:
     *  LIMID N M O SIZE1 ... SIZE(N+M) PA1 ... PA(N+M+O) CPT1 ... CPTN UTIL1 ... UTILO
     *  where:
     *  - LIMID is a file header string
     *  - N, M, and O denote, resp., the no. of chance, decision and value nodes.
     *  - SIZEi is an integer describing the cardinality of variable i (order is chance variables then decision variables)
     *  - PAi is NUMPA PA, where NUMPA is the number of parents of the variable and PA is a space-separated list of IDs for each of the parents. 
     *  - CPTi is NUMVALUES VALUE1 ... VALUE(NUMVALUES), with NUMVALUES numeric values describing the conditional probability values P(C_ID|PA_ID), arranged such that the least significant value (the one that increases the fastest) is C_ID and the most signficant value (the one that increases the slowest) is the value of the rightmost parent.
     *  - UTILi is NUMVALUES VALUE1 ... VALUE(NUMVALUES), with NUMVALUES numeric values specifying the utility function U(PA_ID) with least significant digit denoting the value of the leftmost parent.
     *
     * Linebreaks are discarded as white spaces, thus they can be used to improve readalibity.
     */
    void load();

    void print()
    {
      for (auto v: _variables)
	std::cout << v << std::endl;
      for (auto f: _factors)
	std::cout << f << std::endl;
    }
  
  };
  
}

#endif
