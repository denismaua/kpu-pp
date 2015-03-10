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

/** DAG implementation */

#include "dag.h"
#include <vector>
#include <queue>

namespace msp {


  /** Moralization */
  Graph DAG::moralization()
  {
    std::vector<std::vector<Variable* > > scopes;
    scopes.resize(_nodes.size());
    unsigned i=0;
    for (auto v: _pa) // for every node
      scopes[i++].assign(v.second.begin(),v.second.end());
    return Graph(scopes);
  }

  /** Topological sorting */
  std::vector<Variable> DAG::topsort()
  {
    std::vector<Variable> topo;
    topo.reserve(_nodes.size());
    // Kahn's algorithm
    std::unordered_map<Variable*, unsigned> d; // node's in-degree
    // Compute in-degrees and store root nodes (i.e., those of zero indegree)
    std::queue<Variable*> Q;
    for (auto& v: _nodes)
      { 
	d[v] = _pa[v].size();
	if (!_pa[v].size()) Q.push(v);
      }
    while (!Q.empty())
      {
	Variable* u = Q.front(); Q.pop();
	topo.push_back( *u );
	for (auto& v: _ch[u])
	  {
	    d[v]--;
	    if (!d[v]) Q.push(v);
	  }
      }
    if (topo.size() != _nodes.size()) throw "ERROR! Directed cycle detected!";
    return topo;
  }

  /** Auxiliar function for printing out the arcs. */
  std::ostream& DAG::print_arcs( std::ostream &o )
  {
    for (auto u: _nodes)
      {	
	o << u->name() << ": ";
	for (auto v: _pa[u])
	  o << v->name() << " ";
	o << ",\r\n";
      }
    return o;
  }

  /** Default printing. */
  std::ostream& operator<<(std::ostream &o, DAG &g) 
  {
    o << "DAG:{arcs:{\r\n";
    g.print_arcs(o);
    o << "}}";
    return o;
  }

}
