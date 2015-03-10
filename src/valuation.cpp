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

// This file implements the Valuation class

#include "valuation.h"
#include <vector>

namespace msp {


  unsigned Valuation::_num_instances = 0;
  unsigned Valuation::max_tables = 0;

  /** Constructor from factor. */
  Valuation::Valuation(const Factor& f)
  {
    std::vector<Variable* > scope(f.width());
    for (unsigned i=0; i < f.width(); i++)
      scope[i] = f.var_at(i);
    _domain = new Domain(scope);
  }

  /** Constructor from two valuations. */
  Valuation::Valuation(const Valuation& f1, const Valuation& f2)
  {
    _id = _num_instances++;
    _domain=union_of(f1._domain,f2._domain);
  }

  /** Constructor from two factors and a variables. */
  Valuation::Valuation(const Valuation& f1, const Valuation& f2, const Variable& var)
  {
    _id = _num_instances++;
    _domain=union_of(f1._domain,f2._domain,&var);
  }

  /** Add table filled with default value. */
  void Valuation::add()
  {
    _l_values.push_back( std::vector<double> (_domain->size(), MSP_DEFAULT_VALUE) );
    _r_values.push_back( std::vector<double> (_domain->size(), MSP_DEFAULT_VALUE) );
    _labels.push_back( std::string() );
    if ((_labels.size() % 100000)==0) std::cerr << "+ num_tables: " << _labels.size()/1000 << "k\tdim:" << _domain->size() << std::endl;
    if (_labels.size() > max_tables) max_tables = _labels.size();
  }

  /** Add given labeled paired tables. */
  void Valuation::add( std::vector<double>&& left, std::vector<double>&& right, std::string&& str)
  {
    if (left.size() != _domain->size() || right.size() != _domain->size()) throw "Valuation::add: Wrong table dimension!";
    _l_values.push_back( left );
    _r_values.push_back( right );
    _labels.push_back( str );
    if ((_labels.size() % 100000)==0) std::cerr << "+ num_tables: " << _labels.size()/1000 << "k\tdim:" << _domain->size() << std::endl;
    if (_labels.size() > max_tables) max_tables = _labels.size();
  }  
  
  /** Pruning. */
  void Valuation::reduce()
  {
    do_pruning();
    do_bucketing();
    //do_clustering();
  }

  /** Pareto dominance pruning */
  void Valuation::do_pruning()
  {
    if (num_tables() < MSP_MIN_NUM_TABLES) return;

    std::vector< std::vector<double> > lvalues, rvalues;
    std::vector< std::string > labels;
    bool dominated;
    for (unsigned t1=0; t1<_labels.size(); ++t1)
      { // check if table t1 is dominated by some other table
    	dominated = false; 
    	for (unsigned t2=0; t2<_labels.size(); ++t2)
    	  if (t1 != t2)
    	    { 
    	      unsigned i;
    	      // check whether t1 is strictly greater than t2 somewhere
    	      for (i=0; i<_domain->size(); ++i)
    		if (_l_values[t1][i] > _l_values[t2][i] || _r_values[t1][i] > _r_values[t2][i]) // then t2 does not dominate t1
    		  break;
    	      if (i == _domain->size()) 
    		{ // t2[i] >= t1[i] for all i
    		  // check whether t1 equals t2
    		  for (i=0; i<_domain->size(); ++i)
    		    if (_l_values[t1][i] != _l_values[t2][i] || _r_values[t1][i] != _r_values[t2][i])
    		      break;
    		  if (i < _domain->size() || t1 < t2) { // if t2 != t1 or t1 < t2 (the latter makes sure one of multiple copies is included)
    		    dominated = true; break;            // then t2 dominates t1
    		  }
    		}
    	    }
    	if (!dominated)
    	  { // add undominated table
    	    lvalues.push_back(std::vector<double> ( _l_values[t1].begin(), _l_values[t1].end()) );
    	    rvalues.push_back(std::vector<double> ( _r_values[t1].begin(), _r_values[t1].end()) );
    	    labels.push_back(std::string (_labels[t1]) );
    	  }
      }
    
    std::swap(_l_values, lvalues); std::swap(_r_values, rvalues);
    std::swap(_labels,labels);
  }

  /** Bucket pruning. */
  void Valuation::do_bucketing()
  {
    if (num_tables() < MSP_MAX_NUM_TABLES) return;

    unsigned scaling = std::pow(MSP_MAX_NUM_TABLES,1.0/_domain->size()/2.0);

    std::vector< std::vector<double> > lvalues, rvalues;
    std::vector< std::string > labels;
    bool hit;
    for (unsigned t1=0; t1<_labels.size(); ++t1)
      { // check if table t1 maps into same bucket as other table already included
    	hit = false; 
    	for (unsigned t2=0; t2<labels.size(); ++t2)
	  {
	    unsigned i;
	    // check whether transformation of t1 differs from that of t2 somewhere
	    for (i=0; i<_domain->size(); ++i)
	      if ((unsigned)(_l_values[t1][i]*scaling) != (unsigned)(lvalues[t2][i]*scaling) || (unsigned)(_r_values[t1][i]*scaling) != (unsigned)(rvalues[t2][i]*scaling))
		break;
	    if (i == _domain->size())
	      {
		hit = true;  break;
	      }
	  }
	if (!hit) { // add table
    	    lvalues.push_back(std::vector<double> ( _l_values[t1].begin(), _l_values[t1].end()) );
    	    rvalues.push_back(std::vector<double> ( _r_values[t1].begin(), _r_values[t1].end()) );
    	    labels.push_back(std::string (_labels[t1]) );
	    if (labels.size() >= MSP_MAX_NUM_TABLES) break;
	}
      }

    std::swap(_l_values, lvalues); std::swap(_r_values, rvalues);
    std::swap(_labels,labels);    
  }


  /** Bucket pruning. */
  void Valuation::do_clustering()
  {
    if (num_tables() < MSP_MAX_NUM_TABLES) return;

    // unsigned scaling = std::pow(MSP_MAX_NUM_TABLES,1.0/_domain->size()/2.0);

    // std::vector< std::vector<double> > lvalues, rvalues;
    // std::vector< std::string > labels;
    // bool hit;
    // for (unsigned t1=0; t1<_labels.size(); ++t1)
    //   { // check if table t1 maps into same bucket as other table already included
    // 	hit = false; 
    // 	for (unsigned t2=0; t2<labels.size(); ++t2)
    // 	  {
    // 	    unsigned i;
    // 	    // check whether transformation of t1 differs from that of t2 somewhere
    // 	    for (i=0; i<_domain->size(); ++i)
    // 	      if ((unsigned)(_l_values[t1][i]*scaling) != (unsigned)(lvalues[t2][i]*scaling) || (unsigned)(_r_values[t1][i]*scaling) != (unsigned)(rvalues[t2][i]*scaling))
    // 		break;
    // 	    if (i == _domain->size())
    // 	      {
    // 		hit = true;  break;
    // 	      }
    // 	  }
    // 	if (!hit) { // add table
    // 	    lvalues.push_back(std::vector<double> ( _l_values[t1].begin(), _l_values[t1].end()) );
    // 	    rvalues.push_back(std::vector<double> ( _r_values[t1].begin(), _r_values[t1].end()) );
    // 	    labels.push_back(std::string (_labels[t1]) );
    // 	}
    //   }

    // std::swap(_l_values, lvalues); std::swap(_r_values, rvalues);
    // std::swap(_labels,labels);

    // calculate pairwise distances
    std::vector< std::vector<double> > dist(_labels.size());
    for (unsigned t1=0; t1<_labels.size(); ++t1)
	dist[t1].resize(_labels.size());
    double maxdist = 0.0;
    for (unsigned t1=0; t1<_labels.size(); ++t1)
      {
	dist[t1].resize(_labels.size());
      for (unsigned t2=t1; t2<_labels.size(); ++t2)
	{
	  double d = 0;
	  for (unsigned i=0; i<_domain->size(); ++i)
	    d += std::pow(_l_values[t1][i]-_l_values[t2][i],2) + std::pow(_r_values[t1][i]-_r_values[t2][i],2);
	  d = std::sqrt(d);
	  dist[t1][t2] = dist[t2][t1] = d;
	  if (d > maxdist) maxdist = d;
	}
      }
	
    // randomly select k medoids
    std::vector< unsigned > med(MSP_MAX_NUM_TABLES);
    //std::vector< unsigned > cluster(_labels.size());
    for (unsigned i=0; i<MSP_MAX_NUM_TABLES; ++i) med[i]=i;
    // compute intial cost (sum of distances from each point to closest medoid)
    double cost=0.0;
    for (unsigned t=0; t<_labels.size(); ++t)
      {
	// find closest medoid to point t
	double mindist=maxdist;
	unsigned medoid = 0;
	for (unsigned k=0; k<MSP_MAX_NUM_TABLES; ++k)
	  if (dist[t][med[k]] < mindist) { mindist = dist[t][med[k]]; medoid=med[k]; }
	cost += dist[t][medoid];
	//cluster[t] = medoid;
      }
    // SEARCH
    for (unsigned it=0; it<10; ++it)
      {
	// replace each medoid with the point that causes the greatest reduction in cost
	for (unsigned k=0; k<MSP_MAX_NUM_TABLES; ++k)
	  {
	    unsigned oldmed = med[k];
	    double oldcost = cost;
	    for (unsigned t=0; t<_labels.size(); ++t)
	      {
		// compute new cost if medoid k is swaped with point t
		cost = 0.0; med[k] = t;
		for (unsigned s=0; s<_labels.size(); ++s)
		  {
		    // find closest medoid to point s
		    double mindist=maxdist;
		    unsigned medoid = 0;
		    for (unsigned j=0; j<MSP_MAX_NUM_TABLES; ++j)
		      if (dist[s][med[j]] < mindist) { mindist = dist[s][med[j]]; medoid=med[j]; }
		    cost += dist[s][medoid];
		  }
		// no improvement? backtrack!
		if (cost > oldcost) { med[k] = oldmed; cost = oldcost; }
	      }
	  }
      }
    // discard non-medoids
    std::vector< std::vector<double> > lvalues, rvalues;
    std::vector< std::string > labels;
    for (unsigned k: med)
      {
	lvalues.push_back(std::vector<double> ( _l_values[k].begin(), _l_values[k].end()) );
	rvalues.push_back(std::vector<double> ( _r_values[k].begin(), _r_values[k].end()) );
	labels.push_back(std::string (_labels[k]) );
      }

    std::swap(_l_values, lvalues); std::swap(_r_values, rvalues);
    std::swap(_labels,labels);    
    
  }  

  /** Default printing. */
  std::ostream& operator<<(std::ostream &o, Valuation &f) {
    o << "Valuation(" << *f._domain << ", tables:" << f.num_tables();
    if (f.size() <= MSP_DISP_MAX_VALS)
      {
	o << ", values:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  {
	    o << "([";
	    for (unsigned j=0; j < f.domain_size(); ++j) 
	      { o << f.l_value(i,j); if (j < f.domain_size()-1) o << ", "; }
	    o << "],[";
	    for (unsigned j=0; j < f.domain_size(); ++j) 
	      { o << f.r_value(i,j); if (j < f.domain_size()-1) o << ", "; }
	    o << "])";
	    if (i < f.num_tables()-1) o << ",";
	  }
	o << "}, labels:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  {
	    o << "[" << f.label(i) << "]";
	    if (i < f.num_tables()-1) o << ",";
	  }
	o << "}";
      }
    o << ")";
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const Valuation &f) {
    o << "Valuation(" << *f._domain << ", tables:" << f.num_tables();
    if (f.size() <= MSP_DISP_MAX_VALS)
      {
	o << ", values:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  {
	    o << "([";
	    for (unsigned j=0; j < f.domain_size(); ++j) 
	      { o << f.l_value(i,j); if (j < f.domain_size()-1) o << ", "; }
	    o << "],[";
	    for (unsigned j=0; j < f.domain_size(); ++j) 
	      { o << f.r_value(i,j); if (j < f.domain_size()-1) o << ", "; }
	    o << "])";
	    if (i < f.num_tables()-1) o << ",";
	  }
	o << "}, labels:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  {
	    o << "[" << f.label(i) << "]";
	    if (i < f.num_tables()-1) o << ",";
	  }
	o << "}";
      }
    o << ")";
    return o;
  }


  // OPERATIONS

  /** Valuation combination 
   *  @param f1 the first valuation
   *  @param f2 the second valuation
   *  @return the product of f1 and f2
   */
  Valuation product(const Valuation& f1, const Valuation& f2) {

    Valuation res(f1,f2); // creates valuation with union domain (and no tables)

    unsigned i, j=0, k=0, l;
    std::vector<unsigned> c(res.width(),0);
    Variable* v;

    // first cache indexing
    std::vector<unsigned> i1(res.domain_size()), i2(res.domain_size());

    for (i=0; i<res.domain_size(); ++i)
      {
	i1[i] = j; i2[i] = k;
    	for (l=0; l<res.width(); l++)
    	  {
    	    c[l]++;
    	    v = res.var_at(l);
    	    if (c[l] == v->size())
    	      {
    		c[l] = 0;
    		j -= (v->size()-1)*f1.offset( *v );
    		k -= (v->size()-1)*f2.offset( *v );
    	      }
    	    else
    	      {
    		j += f1.offset( *v ); k += f2.offset( *v );
    		break;
    	      }
    	  }
      }

    // now actually compute values
    for (unsigned r=0; r < f1.num_tables(); r++)
      for (unsigned s=0; s < f2.num_tables(); s++)
	{
	  std::vector<double> left(res.domain_size()), right(res.domain_size());	  	  
	  for (i=0; i<res.domain_size(); ++i)
	    {
	      left[i] = f1.l_value(r,i1[i])*f2.l_value(s,i2[i]); // left part 
	      right[i] = f1.l_value(r,i1[i])*f2.r_value(s,i2[i])+f2.l_value(s,i2[i])*f1.r_value(r,i1[i]); // right part
	    }
	  res.add( std::move(left), std::move(right), f1.label(r)+f2.label(s) );
	}
    
    return res;
    
  }  

  /** Joint valuation combination and variable elimination 
   *  @param f1 the first valuation
   *  @param f2 the second valuation
   *  @param var a pointer to a variable
   *  @return the elimination of var from the combination of f1 and f2
   */
  Valuation sum_product(const Valuation& f1, const Valuation& f2, const Variable& var) {

    Valuation res(f1,f2,var); // creates valuation with union domain setminus var (and no tables)

    unsigned i, j=0, k=0, l, m;
    std::vector<unsigned> c(res.width()+1,0);
    Variable* v;
    // first cache indexing
    std::vector<unsigned> i1(res.domain_size()*var.size()), i2(res.domain_size()*var.size());
    for (i=0; i<res.domain_size(); i++)
	for (m=0; m<var.size(); m++) 
	  {
	    i1[i*var.size()+m] = j; i2[i*var.size()+m] = k;
	    c[0]++;
	    if (c[0] == var.size())
	      {
		c[0] = 0;
		j -= (var.size()-1)*f1.offset( var );
		k -= (var.size()-1)*f2.offset( var );

		for (l=1; l<res.width()+1; l++)
		  {
		    c[l]++;
		    v = res.var_at(l-1);
		    if (c[l] == v->size())
		      {
			c[l] = 0;
			j -= (v->size()-1)*f1.offset( *v );
			k -= (v->size()-1)*f2.offset( *v );
		      }
		    else
		      {
			j += f1.offset( *v ); k += f2.offset( *v );
			break;
		      }
		  }

	      }
	    else
	      {
		j += f1.offset( var ); k += f2.offset( var );
	      }
	  }

    for (unsigned r=0; r < f1.num_tables(); r++)
      for (unsigned s=0; s < f2.num_tables(); s++)
    	{
	  std::vector<double> left(res.domain_size(), 0.0), right(res.domain_size(), 0.0);	  	  
	  
    	  for (i=0; i<res.domain_size(); ++i) // combination
    	    for (m=0; m<var.size(); ++m) // marginalization
	      {
		left[i] += f1.l_value(r,i1[i*var.size()+m])*f2.l_value(s,i2[i*var.size()+m]); // left part 
		right[i] += f1.l_value(r,i1[i*var.size()+m])*f2.r_value(s,i2[i*var.size()+m])+f2.l_value(s,i2[i*var.size()+m])*f1.r_value(r,i1[i*var.size()+m]); // right part
	      }
	  res.add( std::move(left), std::move(right), f1.label(r)+f2.label(s) );	  
    	}

    return res;

  }

  
}
