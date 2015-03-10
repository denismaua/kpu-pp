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

/** Directed Acyclic Graph interface */

#ifndef MSP_DAG_H
#define MSP_DAG_H

#include <unordered_set>
#include <unordered_map>
#include <iostream>

#include "variable.h"
#include "graph.h"

namespace msp {


  /** Direct Acyclic Graph over Variables */
  class DAG {

  private:
    std::unordered_map<Variable*, std::unordered_set<Variable* > > _pa; /*< map to adjancency lists */
    std::unordered_map<Variable*, std::unordered_set<Variable* > > _ch; /*< map to adjancency lists */
    std::unordered_set<Variable* > _nodes; /*< vector of nodes */

  public:

    /** Return number of nodes */
    unsigned order() { return _nodes.size(); }

    /** Add node. 
     *
     * @param v a variable
     */
    void add_node(Variable& v) { _nodes.insert(&v); }

    /** Add arc.
     * @param from a variable
     * @param to a variable
    */
    void add_arc(Variable& from, Variable& to) { _pa[&to].insert(&from); _ch[&from].insert(&to); }

    /** Build and return moralization.
     *
     * The moralization of a DAG is the undirected graph obtained by connecting any two parents of a variable and dropping arc directions.
     *
     * @return a graph object containing the morazliation of the dag.
     */
    Graph moralization();

    /** Returns a topological ordering of variables. */
    std::vector<Variable> topsort();


    /** output factor content info. */
    std::ostream& print_arcs( std::ostream &o );
    friend std::ostream &operator<<(std::ostream &o, DAG &g); 

  };

}


#endif
