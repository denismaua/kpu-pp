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

/** This file describes the interface of Graph class, which contains routines for manipulating domain graphs */

#ifndef MSP_GRAPH_H
#define MSP_GRAPH_H

#include "variable.h"
#include "factor.h"
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <cmath>

namespace msp {

  /** Undirected graph class.
   *
   * Undirected graph representation as list of parents
   * and commom functions for computing eliminatin orders.
   */

  class Graph 
  {

  protected:
    std::unordered_map<Variable*, std::unordered_set<Variable* > > _g; /*< map from nodes to adjancency sets */
    std::unordered_map<Variable*, bool> _p; /*< node flag (used for marking processed nodes internally) */
    std::vector<Variable* > _nodes; /*< node ordering */
    unsigned _u_width, _l_width; /*< treewidth estimates */
    double _u_wwidth, _l_wwidth; /*< weighted treewidth estimates */
    
  public:

    /** Class constructor.
     *
     * Creates the domain graph for a vector of factors.
     */
    Graph( std::vector<Factor >& factors );

    /** Class constructor.
     *
     * Creates the domain graph for a vector of cliques.
     */
    Graph( std::vector<std::vector<Variable* > >& cliques );


    /** Return treewidth estimates for the graph.
     *
     * Return lower and upper bounds on the treewidth of the graph. The graph is intializied with trivial bounds which are updated by call to other methods (e.g., find_order).
     *
     * @return a pair of unsigned integers.
     */
    std::pair<unsigned, unsigned> treewidth() { return std::make_pair(_l_width, _u_width); }

    /** Return weighted treewidth estimates for the graph.
     *
     * Return lower and upper bounds on the weighted treewidth of the graph. The weighted treewidth is the log_2 of the product of the variable cardinalities in the largest clique of an optimal chordalization of the graph. The graph is intializied with trivial bounds which are updated by call to other methods (e.g., find_order).
     *
     * @return a pair of doubles.
     */
    std::pair<double, double> w_treewidth() { return std::make_pair(_l_wwidth, _u_wwidth); }


    /** Min fill heuristic.
     *
     * Returns the next node and its score in ordering according to the min fill heuristic.
     * @return a pair containing the best node as first element and its score as second element
     */
    std::pair<Variable*, unsigned> min_fill();

    /** Min degree heuristic.
     *
     * Returns the next node and its score in ordering according to the min fill heuristic.
     * @return a pair containing the best node as first element and its score as second element
     */
    std::pair<Variable*, unsigned> min_degree();


    /** Triangulates the graph and find a suitable variable elimination sequence.
     *
     *  Uses a heuristic to triangulate the graph (i.e., make it chordal) and then
     *  finds a perfect elimination sequence for the resulting chordal graph.
     *
     */
    void triangulate();

    /** Returns a simplicial node.
     *
     * @return a pointer to a simplicial node or a NULL if no simplicial exists.
     */
    Variable* find_simplicial();


    /** Returns an ordered vector of Variables.
     *
     * @return a vector of variable objects ordered accordingly
     */
    std::vector<Variable > ordering();



    /** output factor content info. */
    std::ostream& print_edges( std::ostream &o );
    friend std::ostream &operator<<(std::ostream &o, Graph &g); 

  };

}


#endif
