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

/** Describes Limid class implementation */
#define VERBOSITY 0 // change this to 1 or higher to increase the verbosity of model loading (useful for debugging)
#define CHAIN 0

#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <random>

#include "limid.h"
#include "graph.h"
#include "inference.h"
#include "setfactor.h"
#include "valuation.h"

#include "operations.h"

namespace msp {

  std::vector<Variable> Limid::find_order(unsigned verbosity)
  {
    std::vector<Variable> order(_variables.size());
    if (CHAIN) 
      {
	// use reverse topological ordering
	// order[N+M+O-1] = _variables[N+M+O-1];
	// for (unsigned i=0,j=0,k=0; i<N+M+O; i++)
	//   {  
	//     if (i>=N+M) { order[i] = _variables[i]; }
	//     else if (i%2==0 && j<N) { order[i] = _variables[j]; ++j; }
	//     else { order[i] = _variables[N+k]; ++k; }
	//   }
	// eliminate value node then decision variables then chance nodes in reverse topological ordering
	order.assign(_variables.begin(),_variables.end());
	std::reverse(order.begin(),order.end());

      } else {  // otherwise find good ordering using min-fill

      // TO-DO: Implement ordering that combines reverse toposort and minfill
      // DAG d
      // for (auto& v: _variables) d.add_node(v);
      // for (auto f: _factors) 
      //  for (unsigned i=1; i<f.width(); ++i) d.add_arc(*f.var_at(i),f.var_at(0));
      // std::cout << d << std::endl
      // order.assign( d.topsort() );

      // first add value variables in given order
      for (unsigned i=0;i<O;++i)
	order[i] = _variables[N+M+i];
      // then decision variables in given order
      // for (unsigned i=0;i<M;++i)
      // 	order[O+i] = _variables[N+i];
      // at last, add chance variables in min-fill order

      // build domain graph
      Graph g(_factors);
      if (verbosity > 2) std::cout << g << std::endl;
      // find elimination sequence
      g.triangulate();
      if (verbosity > 2) std::cout << g << std::endl;
      else if (verbosity > 0) std::printf("treewidth: %d\t\tweighted_treewidth: %g\r\n\r\n", g.treewidth().second, g.w_treewidth().second);
      std::vector<Variable> g_order = g.ordering();

      // for (unsigned j=0, i=M+O; j<g_order.size(); ++j)
      // 	if (g_order[j].name()[0] == 'C') { order[i] = g_order[j]; ++i; }
      for (unsigned j=0, i=O; j<g_order.size(); ++j)
	if (g_order[j].name()[0] != 'V') { order[i] = g_order[j]; ++i; }

    }

    if (verbosity > 1)
      {
	std::cout << "elimination sequence: ";
	for (auto& v: order)
	  std::cout << v.name() << " ";
	std::cout << std::endl ;
      }

    return order;
  }

  /** Solve Limid using Single Policy Updating Algorithm (SPU) */
  double Limid::SPU(std::vector<Variable>& order, unsigned& it, unsigned& numtables, unsigned verbosity)
  {
    if (verbosity) std::cout << "SPU" << std::endl;
    // INITIALIZATION
    std::vector< SetFactor > sfactors;
    sfactors.reserve(_factors.size());
    for (unsigned i=0; i<N+M+O; ++i)
      sfactors.emplace_back(_factors[i]);
    sfactors.shrink_to_fit();

    std::vector<SetFactor> p_space; // space of policies for each decision variable
    p_space.reserve(M);
    for (unsigned i=N; i<N+M; ++i)
      {
	SetFactor p(_factors[i]);
	for (unsigned j=0; j<_variables[i].size(); j++)
	  {
	    p.set_label(j,std::to_string(j));
	    for (unsigned k=0; k<_variables[i].size(); k++)
	      if (k==j) p.set( j, k, 1.0);
	      else p.set( j, k, 0.0 );
	    if (j < _variables[i].size()-1) p.add_table(); 
	  }
	p_space.push_back(p);
      }
    p_space.shrink_to_fit();

    // SEARCH
    SetFactor::max_tables = 0;     // reset maximum number of tables in a valuation
    // evaluate initial strategy
    SetFactor r = set_variable_elimination(order, sfactors, verbosity);
    std::pair<double, std::pair<unsigned,unsigned> > x = r.max2();
    double E = x.first;
    
    double prevE=-INFINITY; 
    it=0;
    while (E>prevE) // iterate until convergence
      {
    	prevE=E;
	if (verbosity) std::cout << "it:"<< it << std::endl;
    	// do 1-neighborhood search
	for (unsigned i=N+M-1; i>=N; --i)
	  {
	    // find best policy for i-th decision variable 
	    sfactors[i] = p_space[i-N];
	    SetFactor r = set_variable_elimination(order, sfactors, verbosity);
	    std::pair<double, std::pair<unsigned,unsigned> > x = r.max2();
	    if (x.first > E)
	      { // found improving move
		std::string s = r.label(x.second.first);
		// update incumbent policy for this variable
		Factor p(_factors[i]);
		for (unsigned j=0; j<p.size(); ++j)
		  p.set( j, p_space[i-N].get( std::stoi(s) , j) );
		_factors[i] = p;
		sfactors[i] = SetFactor(p);
		E = x.first;
	      } else sfactors[i] = SetFactor(_factors[i]); // restore policy
	  }
	++it;
      }
    numtables = SetFactor::max_tables;
    
    return E;
  }

  /** Solve Limid using Single Policy Updating Algorithm (SPU) */
  double Limid::SPU2(std::vector<Variable>& order, unsigned& it, unsigned& numtables, unsigned verbosity)
  {
    if (verbosity) std::cout << "SPU2" << std::endl;
    
    // INITIALIZATION    
    std::vector<Valuation> vals;
    vals.reserve(N+M+O);
    // initialize chance and decision variable valuations
    for (unsigned i=0; i<N+M; ++i)
      {
	vals.emplace_back( _factors[i] );
	vals[i].add( std::vector<double> (_factors[i].begin(),_factors[i].end()), std::vector<double> (_factors[i].size(),0.0), std::string() );
      }

    // initialize value variable valuations
    for (unsigned i=N+M; i<N+M+O; ++i)
      {
	vals.emplace_back( _factors[i] );
	vals[i].add( std::vector<double> (_factors[i].size(),1.0), std::vector<double> (_factors[i].begin(),_factors[i].end()), std::string() );
      }
    vals.shrink_to_fit();

    // pre-compute space of policies for each decision variable
    std::vector<Valuation> p_space;
    p_space.reserve(M);
    for (unsigned i=N; i<N+M; ++i)
      {
	p_space.emplace_back( _factors[i] );
	for (unsigned j=0; j<_variables[i].size(); j++)
	  {
	    std::vector<double> left(_factors[i].size(),0.0);
	    left[j]=1.0;
	    p_space[i-N].add( std::move(left), std::vector<double> (_factors[i].size(),0.0), std::to_string(j) );
	  }
      }
    p_space.shrink_to_fit();

    // SEARCH

    // reset maximum number of tables in a valuation
    Valuation::max_tables = 0;
    // evaluate initial strategy
    Valuation r = variable_elimination(order, vals, verbosity);
    std::pair<double, std::pair<unsigned,unsigned> > x = r.maxr(); // find maximum in case there are more than one table
    double E=x.first;
    double prevE=-INFINITY; 
    it=0;
    while (E>prevE) // iterate until convergence
      {
    	prevE=E;
	if (verbosity) std::cout << "it:"<< it << std::endl;
    	// do 1-neighborhood search
	for (unsigned i=N+M-1; i>=N; --i)
	  {
	    // find best policy for i-th decision variable 
	    vals[i] = p_space[i-N];
	    Valuation r = variable_elimination(order, vals, verbosity);
	    std::pair<double, std::pair<unsigned,unsigned> > x = r.maxr(); // find maximum in case there are more than one table
	    if (x.first > E)
	      { // found improving strategy
		E = x.first;
		// update incumbent policy for this variable
		for (unsigned j=0; j < _factors[i].size(); ++j) _factors[i].set(j,0.0);
		_factors[i].set(std::stoi(r.label(x.second.first)),1.0);
	      }
	    // transform current valuation back into single-table valuation
	    vals[i].clear();
	    vals[i].add( std::vector<double> (_factors[i].begin(),_factors[i].end()), std::vector<double>(_factors[i].size(),0.0), std::string() );
	  }
	++it;
      }
    numtables=Valuation::max_tables;
    
    return E;
  }

  
  /** Solve Limid using k-Policy Updating Algorithm (k-Policy Local Search, kPU) */
  double Limid::kPU(std::vector<Variable>& order, unsigned k, unsigned& it, unsigned& numtables, unsigned verbosity)
  {    
    if (k>M) k=M;
    if (verbosity) std::cout << k << "PU" << std::endl;

    std::vector< SetFactor > sfactors;
    sfactors.reserve(_factors.size());
    for (unsigned i=0; i<N+M+O; ++i)
      sfactors.emplace_back(_factors[i]);
    sfactors.shrink_to_fit();

    std::vector<SetFactor> p_space; // space of policies for each decision variable
    p_space.reserve(M);
    for (unsigned i=N; i<N+M; ++i)
      {
	SetFactor p(_factors[i]);
	for (unsigned j=0; j<_variables[i].size(); j++)
	  {
	    p.set_label(j,std::string("V")+std::to_string(i)+"T"+std::to_string(j));
	    for (unsigned k=0; k<_variables[i].size(); k++)
	      if (k==j) p.set( j, k, 1.0);
	      else p.set( j, k, 0.0 );
	    if (j < _variables[i].size()-1) p.add_table(); 
	  }
	p_space.push_back(p);
      }
    p_space.shrink_to_fit();

    SetFactor::max_tables = 0;

    // evaluate intial strategy
    SetFactor r = set_variable_elimination(order, sfactors, verbosity);
    // interpret result to get locally optimal strategy
    std::pair<double, std::pair<unsigned,unsigned> > x = r.max2();
    double E = x.first;
    double prevE=-INFINITY;

    it=0;
    while (E>prevE) // iterate until convergence
      {
    	prevE=E;
	if (verbosity) std::cout << "it:"<< it << std::endl;
    	// do k-neighborhood search
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,M-1);
	// samples n subsets of {1,...,M} of cardinality k such that
	// the i-th sample contains element i
	for (unsigned i=0; i < M; ++i)	    
	  {
	    std::set<unsigned> ne;
	    ne.insert(i); 
	    while(ne.size() < k) ne.insert( distribution(generator) );
	  
	    for (unsigned j: ne) sfactors[N+j] = p_space[j];
	    if (verbosity > 1) { std::cout << "neighborhood: "; for (unsigned j: ne) std::cout << j << " "; std::cout << std::endl; } // print out search dimensions
	    // run set-valued variable eliminatoin
	    SetFactor r = set_variable_elimination(order, sfactors, verbosity);
	    // interpret result to get locally optimal strategy
	    std::pair<double, std::pair<unsigned,unsigned> > x = r.max2();
	    if (x.first > E)
	      { // found improving move
		std::string s = r.label(x.second.first);
		unsigned c=0;
		unsigned vnum,tnum;
		// parse labels to recover partial strategy
		while (c < s.size())
		  {
		    if (s[c] != 'V') { std::cerr << "weird thing happened here:" << s.substr(0,c) << '_' << s[c] << '_' << s.substr(c+1) << std::endl; }
		    for (vnum=0,c++;c<s.size();c++)
		      {
			if (s[c] < '0' || s[c] > '9') break;
			vnum *= 10; vnum += s[c]-'0';
		      }
		    if (s[c] != 'T') { std::cerr << "weird thing happened here:" << s.substr(0,c) << '_' << s[c] << '_' << s.substr(c+1) << std::endl; }
		    for (tnum=0,c++;c<s.size();c++)
		      {
			if (s[c] < '0' || s[c] > '9') break;
			tnum *= 10; tnum += s[c]-'0';
		      }
		    // update incumbent strategy
		    Factor p(_factors[vnum]);
		    for (unsigned j=0; j<p.size(); ++j)
		      p.set( j, p_space[vnum-N].get( tnum , j) );
		    _factors[vnum] = p;
		    sfactors[vnum] = SetFactor(p);
		  }
		E = x.first;
	      } 
	    else for (unsigned j: ne) sfactors[N+j] = SetFactor(_factors[N+j]); // restore policies
	  }

	// // search all n chooses k subsets (very inneficient)
	// // ne is a bitmask for which decision variables to include in the neighborhood.
	// // start with the last set in lexicographic order:	
	// std::string ne(k, 1); // K leading 1's
	// ne.resize(M, 0); // M-k trailing 0's
	// do
	//   {
	//   // variables in search neighborhood are associated to sets of all policies
	//     for (unsigned i=0; i<M; ++i)
	//       if (ne[i]) sfactors[N+i] = p_space[i];
	//     // run set-valued variable eliminatoin
	//     SetFactor r = set_variable_elimination(order, sfactors, verbosity);
	//     // interpret result to get locally optimal strategy
	//     std::pair<double, std::pair<unsigned,unsigned> > x = r.max2();
	//     if (x.first > E)
	//       { // found improving move
	// 	std::string s = r.label(x.second.first);
	// 	unsigned c=0;
	// 	unsigned vnum,tnum;
	// 	// parse labels to recover partial strategy
	// 	while (c < s.size())
	// 	  {
	// 	    if (s[c] != 'V') { std::cerr << "weird thing happened here:" << s.substr(0,c) << '_' << s[c] << '_' << s.substr(c+1) << std::endl; }
	// 	    for (vnum=0,c++;c<s.size();c++)
	// 	      {
	// 		if (s[c] < '0' || s[c] > '9') break;
	// 		vnum *= 10; vnum += s[c]-'0';
	// 	      }
	// 	    if (s[c] != 'T') { std::cerr << "weird thing happened here:" << s.substr(0,c) << '_' << s[c] << '_' << s.substr(c+1) << std::endl; }
	// 	    for (tnum=0,c++;c<s.size();c++)
	// 	      {
	// 		if (s[c] < '0' || s[c] > '9') break;
	// 		tnum *= 10; tnum += s[c]-'0';
	// 	      }
	// 	    // update incumbent strategy
	// 	    Factor p(_factors[vnum]);
	// 	    for (unsigned j=0; j<p.size(); ++j)
	// 	      p.set( j, p_space[vnum-N].get( tnum , j) );
	// 	    _factors[vnum] = p;
	// 	    sfactors[vnum] = SetFactor(p);
	// 	  }
	// 	E = x.first;
	//       } 
	//     else for (unsigned i=0; i<M; ++i) if (ne[i]) sfactors[N+i] = SetFactor(_factors[N+i]); // restore policies

	//   } while (std::prev_permutation(ne.begin(),ne.end()));
	
	++it;
      }
    numtables = SetFactor::max_tables;
    return E;
  }


  /** Solve Limid using k-Policy Updating Algorithm (k-Policy Local Search, kPU) */
  double Limid::kPU2(std::vector<Variable>& order, unsigned k, unsigned& it, unsigned& numtables, unsigned verbosity)
  {    
    if (k>M) k=M;
    if (verbosity) std::cout << k << "PU" << std::endl;

    // INITIALIZATION

    std::vector<Valuation> vals;
    vals.reserve(N+M+O);
    // initialize chance and decision variable valuations
    for (unsigned i=0; i<N+M; ++i)
      {
	vals.emplace_back( _factors[i] );
	vals[i].add( std::vector<double> (_factors[i].begin(),_factors[i].end()), std::vector<double> (_factors[i].size(),0.0), std::string() );
      }

    // initialize value variable valuations
    for (unsigned i=N+M; i<N+M+O; ++i)
      {
	vals.emplace_back( _factors[i] );
	vals[i].add( std::vector<double> (_factors[i].size(),1.0), std::vector<double> (_factors[i].begin(),_factors[i].end()), std::string() );
      }
    vals.shrink_to_fit();

    // pre-compute space of policies for each decision variable
    std::vector<Valuation> p_space;
    p_space.reserve(M);
    for (unsigned i=N; i<N+M; ++i)
      {
	p_space.emplace_back( _factors[i] );
	for (unsigned j=0; j<_variables[i].size(); j++)
	  {
	    std::vector<double> left(_factors[i].size(),0.0);
	    left[j]=1.0;
	    p_space[i-N].add( std::move(left), std::vector<double> (_factors[i].size(),0.0), std::string("V")+std::to_string(i)+"T"+std::to_string(j) );
	  }
      }
    p_space.shrink_to_fit();

    Valuation::max_tables = 0; // reset max. no. of tables in a valuation

    // SEARCH
    
    // evaluate initial strategy
    Valuation r = variable_elimination(order, vals, verbosity);
    // interpret result to get locally optimal strategy
    std::pair<double, std::pair<unsigned,unsigned> > x = r.maxr();
    double E = x.first;
    double prevE=-INFINITY;
    
    it=0;
    while (E>prevE) // iterate until convergence
      {
    	prevE=E;
	if (verbosity) std::cout << "it:"<< it << std::endl;
    	// do k-neighborhood stochastic search
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,M-1);
	// samples n subsets of {1,...,M} of cardinality k such that
	// the i-th sample contains element i
	for (unsigned i=0; i < M; ++i)
	  { 
	    std::set<unsigned> ne;
	    ne.insert(i); 
	    while(ne.size() < k) ne.insert( distribution(generator) );
	    for (unsigned j: ne) vals[N+j] = p_space[j];
	    if (verbosity) { std::cout << "neighborhood: "; for (unsigned j: ne) std::cout << j << " "; std::cout << std::endl; } // print out search dimensions
	    // run variable elimination
	    Valuation r = variable_elimination(order, vals, verbosity);
	    // interpret result to get locally optimal strategy
	    std::pair<double, std::pair<unsigned,unsigned> > x = r.maxr();
	    if (x.first > E)
	      { // found improving strategy
		E = x.first;
		// parse labels to recover partial strategy
		std::string s = r.label(x.second.first);
		unsigned c=0;
		unsigned vnum,tnum;
		while (c < s.size())
		  {
		    if (s[c] != 'V') { std::cerr << "weird thing happened here:" << s.substr(0,c) << '_' << s[c] << '_' << s.substr(c+1) << std::endl; }
		    for (vnum=0,c++;c<s.size();c++)
		      {
			if (s[c] < '0' || s[c] > '9') break;
			vnum *= 10; vnum += s[c]-'0';
		      }
		    if (s[c] != 'T') { std::cerr << "weird thing happened here:" << s.substr(0,c) << '_' << s[c] << '_' << s.substr(c+1) << std::endl; }
		    for (tnum=0,c++;c<s.size();c++)
		      {
			if (s[c] < '0' || s[c] > '9') break;
			tnum *= 10; tnum += s[c]-'0';
		      }
		    // update incumbent strategy
		    for (unsigned j=0; j<_factors[vnum].size(); ++j) _factors[vnum].set(j,0.0);
		    _factors[vnum].set( tnum, 1.0 );
		  }
	      } 
	    // transform valuations in the neighborhood back into single-table valuations
	    for (unsigned j: ne)
	      {
		vals[N+j].clear();
		vals[N+j].add( std::vector<double> (_factors[N+j].begin(),_factors[N+j].end()), std::vector<double>(_factors[N+j].size(),0.0), std::string() );
	      }
	  } 	 
	++it;
      }
    numtables = Valuation::max_tables;
    
    return E;
  }


  /** Solve Limid using Multiple Policy Updating Algorithm (MPU) */
  double Limid::MPU(std::vector<Variable>& order, unsigned& numtables, unsigned verbosity)
  {
    if (verbosity) std::cout << "MPU" << std::endl;
    std::vector< SetFactor > sfactors;
    sfactors.reserve(_factors.size());
    for (unsigned i=0; i<N; ++i)
      sfactors.emplace_back(_factors[i]);
    for (unsigned i=N; i<N+M; ++i)
      {
	SetFactor p(_factors[i]);
	for (unsigned j=0; j<_variables[i].size(); j++)
	  {
	    p.set_label(j,std::string("V")+std::to_string(i)+"T"+std::to_string(j));
	    for (unsigned k=0; k<_variables[i].size(); k++)
	      if (k==j) p.set( j, k, 1.0);
	      else p.set( j, k, 0.0 );
	    if (j < _variables[i].size()-1) p.add_table(); 
	  }
	sfactors.push_back(p);
      }
    for (unsigned i=N+M; i<N+M+O; ++i)
      sfactors.emplace_back(_factors[i]);
    sfactors.shrink_to_fit();

    // do variable elimination
    SetFactor::max_tables = 0;
    SetFactor f = set_variable_elimination(order, sfactors, verbosity);    
    double val=-INFINITY;
    // get maximum
    for (unsigned i=0; i<f.num_tables(); ++i)
      for (unsigned j=0; j<f.domain_size(); ++j)
	{
	  double v = f.get(i,j);
	  if (v > val) val=v;
	}
    numtables = SetFactor::max_tables;
    return val;

  }

   /** Solve Limid using Multiple Policy Updating Algorithm (MPU) */
  double Limid::MPU2(std::vector<Variable>& order, unsigned& numtables, unsigned verbosity)
  {

    // INITIALIZATION
    std::vector<Valuation> vals;
    vals.reserve(N+M+O);
    // chance variable valuations
    for (unsigned i=0; i<N; ++i)
      {
	vals.emplace_back( _factors[i] );
	vals[i].add( std::vector<double> (_factors[i].begin(),_factors[i].end()), std::vector<double> (_factors[i].size(),0.0), std::string() );
      }
    // decision variable valuations
    for (unsigned i=N; i<N+M; ++i)
      {
	vals.emplace_back( _factors[i] );
	for (unsigned j=0; j<_variables[i].size(); j++)
	  {
	    std::vector<double> left(_factors[i].size(),0.0);
	    left[j]=1.0;
	    vals[i].add( std::move(left), std::vector<double> (_factors[i].size(),0.0), std::string("V")+std::to_string(i)+"T"+std::to_string(j) );
	  }
      }
    // value variable valuations
    for (unsigned i=N+M; i<N+M+O; ++i)
      {
	vals.emplace_back( _factors[i] );
	vals[i].add( std::vector<double> (_factors[i].size(),1.0), std::vector<double> (_factors[i].begin(),_factors[i].end()), std::string() );
      }
    vals.shrink_to_fit();

    // SEARCH
    // do variable elimination
    Valuation::max_tables = 0; // reset max. no. of tables in a valuation
    Valuation r = variable_elimination(order, vals, verbosity);    
    // interpret result to get optimal strategy
    std::pair<double, std::pair<unsigned,unsigned> > x = r.maxr();
    // parse labels to recover strategy
    std::string s = r.label(x.second.first);
    unsigned c=0;
    unsigned vnum,tnum;
    while (c < s.size())
      {
	// parse deterministic policy
	for (vnum=0,c++;c<s.size();c++)
	  {
	    if (s[c] < '0' || s[c] > '9') break;
	    vnum *= 10; vnum += s[c]-'0';
	  }
	for (tnum=0,c++;c<s.size();c++)
	  {
	    if (s[c] < '0' || s[c] > '9') break;
	    tnum *= 10; tnum += s[c]-'0';
	  }
	// save policy
	for (unsigned j=0; j<_factors[vnum].size(); ++j) _factors[vnum].set(j,0.0);
	_factors[vnum].set( tnum, 1.0 );
      }
    numtables = Valuation::max_tables;
    return x.first;
  } 

  /** solve limid with given method */
  Limid::solution_t Limid::solve(unsigned k, unsigned verbosity)
  {
    solution_t sol;
    sol.method = std::to_string(k)+std::string("PU");
    unsigned iter,numtables;
    std::vector<Variable> order = find_order(verbosity);
    if (k==1) 
      if (O > 1)
	sol.value = SPU2(order,sol.num_iters,sol.num_tables,verbosity);
      else
	sol.value = SPU(order,sol.num_iters,sol.num_tables,verbosity);
    else if (k < M)
      if (O > 1)
	{
	  unsigned iters, num_tables;
	  sol.value = SPU2(order,iters,num_tables,verbosity);
	  sol.value = kPU2(order,k,sol.num_iters,sol.num_tables,verbosity);
	  sol.num_iters += iters;
	  if (num_tables > sol.num_tables) sol.num_tables = num_tables;
	}
      else
	{
	  unsigned iters, num_tables;
	  sol.value = SPU(order,iters,num_tables,verbosity);
	  sol.value = kPU(order,k,sol.num_iters,sol.num_tables,verbosity);
	  sol.num_iters += iters;
	  if (num_tables > sol.num_tables) sol.num_tables = num_tables;
	}
    else
      {
	sol.method = std::string("MPU");
	sol.num_iters = 1;
	if (O>1) sol.value = MPU2(order,sol.num_tables,verbosity);
	else sol.value = MPU(order,sol.num_tables,verbosity);
      }
    return sol;
  }

  /** Load model from file in LIMID format. Read input from stdin --
   * redirect cin to read from file.
   */
  void Limid::load()
  {
    std::string str;
    // assume cin is set accordingly
    // process header
    std::cin >> str;
    // skip possible C-style comment block at beginning
    if (str[0] == '/' && str[1] == '*') { while(str[0] != '*' or str[1] != '/') std::cin >> str; std::cin >> str; }
    if (str.compare("LIMID") != 0) 
      {
	std::cerr << "ERROR! Expected 'LIMID' file header, found: " << str << std::endl;
	throw "Limid::load: Expected 'LIMID' file header!";
      }
    // read number of variables
    std::cin >> N;        std::cin >> M;        std::cin >> O;    
    if (VERBOSITY > 1) {
      std::cout << "N:" << N << "\tM:" << M << "\tO:" << O << std::endl;
      std::cout << std::flush;
    }

    // allocate memory
    _variables.reserve(N+M+O);    _factors.reserve(N+M+O);

    // read variables' sizes and initialize them
    for (unsigned i = 0; i < N; ++i) 
      {
    	unsigned size;
    	std::string s;
    	std::cin >> size; // cardinality
    	s += 'C'; s += std::to_string( i );
	_variables.emplace_back( size, s );
      }

    for (unsigned i = 0; i < M; ++i) 
      {
    	unsigned size;
    	std::string s;
    	std::cin >> size; // cardinality
    	s += 'D'; s += std::to_string( i );
	_variables.emplace_back( size, s );
      }
    
    for (unsigned i = 0; i < O; ++i) 
      {
    	unsigned size;
    	std::string s;
    	s += 'V'; s += std::to_string( i );
	_variables.emplace_back( 1, s );
      }

    if (VERBOSITY > 1) {
      for (auto v: _variables)
	std::cout << v << std::endl;
      std::cout << std::flush;
    }

    // read parent sets 
    for (unsigned i = 0; i < N+M+O; ++i) 
      {
    	unsigned numpa;
    	std::vector<Variable* > scope;
    	scope.push_back( &(_variables[i]) );
	
    	std::cin >> numpa; // no. of parents
    	for (unsigned j=0; j < numpa; j++)
    	  {
    	    unsigned id;
	    std::cin >> id;
	    if (id >= N+M)
	      {
		std::cerr << "Limid::load Variable id is outside allowed range: " << id << " <> [0," << (N+M) << "]!" << std::endl;
		throw;
	      }
	    scope.push_back( &(_variables[id]) );
    	  }
	if (i >= N && i < N+M)
	  _factors.emplace_back( scope, 1.0/_variables[i].size() ); // initialize policies with uniform strategy
	else
	  _factors.emplace_back( scope );
	if (VERBOSITY > 1) {
	  std::cout << _factors[i] << std::endl << std::flush;
	}
      }
	

    // read CPTs
    for (unsigned i = 0; i < N; ++i) 
      {
    	double val;
    	unsigned j=0;
	std::cin >> j;
	  if (j != _factors[i].size()) 
	    {
	      std::cerr << "Limid::load: Expected " << _factors[i].size() << " values, found: " << j << "!" << std::endl;
	      throw;
	    }
    	for (j=0; j < _factors[i].size(); j++)
    	  { std::cin >> val; _factors[i].set( j, val ); }
      }

    // read Utilities
    for (unsigned i = N+M; i < N+M+O; ++i) 
      {
    	double val;
    	unsigned j=0;
	std::cin >> j;
	  if (j != _factors[i].size()) 
	    {
	      std::cerr << "Limid::load: Expected " << _factors[i].size() << " values, found: " << j << "!" << std::endl;
	      throw;
	    }
    	for (j=0; j < _factors[i].size(); j++)
    	  { std::cin >> val; _factors[i].set( j, val ); }
      }

    if (VERBOSITY > 1)
      for (auto f: _factors)
	std::cout << f << std::endl;
  }
  
}
