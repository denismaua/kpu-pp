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

/** This file implements the Graph class, which contains routines for manipulating domain graphs */

#include "graph.h"
#include <utility>

namespace msp {

  /** Constructor for given list of factors */
  Graph::Graph( std::vector<Factor >& factors ) : _l_width(0), _l_wwidth(0) // trivial lower bounds on treewidth
  {
    // turn factor scopes into cliques
    for (Factor f: factors)
      for (unsigned i=0; i<f.width()-1; i++)
	for (unsigned j=i+1; j<f.width(); j++)
	  { _g[f.var_at(i)].insert(f.var_at(j));  _g[f.var_at(j)].insert(f.var_at(i)); }

    _nodes.resize(_g.size());

    unsigned card=0;
    // unmark all nodes and find highest size of variable
    for (auto v: _g)
      {
	_p[v.first] = false;
    	if (v.first->size() > card) card = v.first->size();
      }
    _u_width  = _g.size(); // trivial upper bound on treewidth
    _u_wwidth = (_u_width+1)*log2(card); // trivial upper bound on weighted treewidth

    // // run min-degree heuristic to find lower bound on treewidth
    std::pair<Variable*, unsigned> r;
    for (unsigned i=0; i < _nodes.size(); i++)
      {
     	// select variable
    	r = min_degree();
    	// add it to order
    	_nodes[i] = r.first;
     	// mark variable
    	_p[r.first] = true;
     	// for detecting maximum variable cardinality
    	if (r.first->size() > card) card = r.first->size();

    	double degree = 0; double weight = log2(r.first->size());
    	// comput weigth
    	  for (auto x: _g[r.first])
    	    if (!_p[x])
    	      {
    		degree++; // node's degree in chordalization
    		weight += log2(x->size()); // node's weight in cordalization
    	      }
    	  if (degree>_l_width) _l_width=degree; // update lower bound on treewidth
    	  if (weight>_l_wwidth) _l_wwidth=weight; // update lower bound on weighted treewidth
    }
  }

  /** Constructor for given vector of variable cliques */
  Graph::Graph( std::vector<std::vector<Variable* > >& cliques ) : _l_width(0), _l_wwidth(0) // trivial lower bounds on treewidth
  {
    // turn factor scopes into cliques
    for (auto clique: cliques)
      for (auto u: clique)
	for (auto v: clique)
	  if (u != v) _g[u].insert(v);

    _nodes.resize(_g.size());

    unsigned card=0;
    // unmark all nodes and find highest size of variable
    for (auto v: _g)
      {
	_p[v.first] = false;
    	if (v.first->size() > card) card = v.first->size();
      }
    _u_width  = _g.size(); // trivial upper bound on treewidth
    _u_wwidth = (_u_width+1)*log2(card); // trivial upper bound on weighted treewidth

    // // run min-degree heuristic to find lower bound on treewidth
    std::pair<Variable*, unsigned> r;
    for (unsigned i=0; i < _nodes.size(); i++)
      {
     	// select variable
    	r = min_degree();
    	// add it to order
    	_nodes[i] = r.first;
     	// mark variable
    	_p[r.first] = true;
     	// for detecting maximum variable cardinality
    	if (r.first->size() > card) card = r.first->size();

    	double degree = 0; double weight = log2(r.first->size());
    	// comput weigth
    	  for (auto x: _g[r.first])
    	    if (!_p[x])
    	      {
    		degree++; // node's degree in chordalization
    		weight += log2(x->size()); // node's weight in cordalization
    	      }
    	  if (degree>_l_width) _l_width=degree; // update lower bound on treewidth
    	  if (weight>_l_wwidth) _l_wwidth=weight; // update lower bound on weighted treewidth
    }
  }

  /** Min fill heuristic.
   *
   * Chooses node with least number of unconnected neighbors and returns a pair (node, degree*), where degree* is the number of higher-ordered neighbors.
   *
   */
  std::pair<Variable*, unsigned> Graph::min_fill()
  {
    unsigned score = _g.size();
    std::pair<Variable* , unsigned> res(0, _g.size());
    for (auto var: _p)
      if (!var.second)
	{ 
	  unsigned this_score = 0; // score for this variable
	  // count the number of fill-in edges for each variable
	  for (auto x: _g[var.first])
	    if (!_p[x])
	      for (auto y: _g[var.first])
		if (*x < *y && !_p[y] && !_g[x].count(y)) // y is not a adjacent to x -- fill-in edge found!
		      this_score++;
	  if (this_score < res.second)
	    { // update incubent solution
	      res.second = this_score; res.first=var.first;
	    }
	}
    return res;
  }

  /** Min degree heuristic.
   *
   * Chooses node with least number of unconnected neighbors, and  returns a pair (node, degree*), where degree* is the number of higher-ordered neighbors.
   *
   */
  std::pair<Variable*, unsigned> Graph::min_degree()
  {
    unsigned score = _g.size();
    std::pair<Variable* , unsigned> res(NULL, _g.size());
    for (auto var: _p)
      if (!var.second)
	{ 
	  unsigned this_score = 0; // score for this variable
	  for (auto x: _g[var.first])
	    if (!_p[x]) 	  // count the number of higher-order neighbors
	      this_score++;
	  if (this_score < res.second)
	    { // update incubent solution
	      res.second = this_score; res.first=var.first;
	    }
	}
    return res;
  }

  /** Triangulates the graph and finds a suitable variable elimination sequence.
   *
   *  Uses a heuristic to triangulate the graph (i.e., make it chordal) and then
   *  finds a perfect elimination sequence for the resulting chordal graph.
   *
   */
  void Graph::triangulate()
  {
    std::pair<Variable*, unsigned> r;

    // unmark all nodes
    for (auto v: _p)
      _p[v.first] = false;

    unsigned i = 0, degree = 0;
    double weight =0, wwidth = 0; // node and order's weight, resp.
    unsigned width = 0; // width of elimination order 
    while (i < _nodes.size())
      {
	// select variable to eliminate
	r = min_fill();
	// add it to order
	_nodes[i++] = r.first;
	// mark variable
	_p[r.first] = true;
	degree = 0; weight = log2(r.first->size());
	// connect neighbors
	  for (auto x: _g[r.first])
	    if (!_p[x])
	      {
		degree++; // node's degree in chordalization
		weight += log2(x->size()); // node's weight in cordalization
		for (auto y: _g[r.first])
		  if (*x < *y && !_p[y])
		    { _g[x].insert(y); _g[y].insert(x); }
	      }
	  if (degree>width) width=degree; // update upper bound on treewidth
	  if (weight>wwidth) wwidth=weight; // update upper bound on weighted treewidth
      }
    if (width<_u_width) _u_width=width; // update upper bound on treewidth (if necessary)
    if (wwidth<_u_wwidth) _u_wwidth=wwidth; // update upper bound on treewidth (if necessary)
  }

  /** Returns a simplicial node.
   *
   * @return a pointer to a simplicial node or a NULL if no simplicial exists.
   */
  Variable* Graph::find_simplicial()
  {
    for (auto v: _nodes)
      {
	bool simplicial = true;
	for (auto x: _g[v])
	  if (!_p[x])
	    for (auto y: _g[v])
	      if (*x < *y && !_p[y] && !_g[x].count(y))
		{
		  simplicial = false; break;
		}
	if (simplicial) return v;
      }
    return NULL;
  }


  /** Returns an ordered vector of Variables 
   *
   * @return a vector of variable objects ordered accordingly
   */
  std::vector<Variable > Graph::ordering()
  {
    std::vector<Variable > order(_nodes.size());
    for (unsigned i=0; i < _nodes.size(); i++)
      order[i] = *(_nodes[i]);
    return order;
  }

  /** Auxiliar function for printing out the edges. */
  std::ostream& Graph::print_edges( std::ostream &o )
  {
    for (auto u: _g)
      {	
	o << (_p[u.first]?"*":" ") << u.first->name() << ": ";
	for (Variable* v: _g[u.first])
	  o << v->name() << " ";
	o << "," << std::endl;
      }
    return o;
  }
  /** Default printing. */
  std::ostream& operator<<(std::ostream &o, Graph &g) 
  {
    std::pair<unsigned,unsigned> tw = g.treewidth();
    std::pair<double,double> wtw = g.w_treewidth();
    o << "Graph(treewidth:[" << tw.first << "," << tw.second << "], w_treewidth:[" << wtw.first << "," << wtw.second << "], edges:{" << std::endl;
    g.print_edges(o);
    o << "})";
    return o;
  }

}
