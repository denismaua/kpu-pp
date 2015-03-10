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

// This file implements the SetFactor class

#include "setfactor.h"
#include "trieset.h"
#include <vector>

namespace msp {


  unsigned SetFactor::_num_instances = 0;
  unsigned SetFactor::max_tables = 0;

  /** Constructor from factor */
  SetFactor::SetFactor(const Factor& f)
  {
    std::vector<Variable* > scope(f.width());
    for (unsigned i=0; i < f.width(); i++)
      scope[i] = f.var_at(i);
    _domain = new Domain(scope);
    add_table();
    for (unsigned i=0; i < f.size(); i++)
      _values[0][i] = f[i];
  }

  /** Constructor from two set-factors */
  SetFactor::SetFactor(const SetFactor& f1, const SetFactor& f2)
  {
    _id = _num_instances++;
    // get union of scopes
    _domain=union_of(f1._domain,f2._domain);
  }

  /** Constructor from two factors and a variables*/
  SetFactor::SetFactor(const SetFactor& f1, const SetFactor& f2, const Variable& var)
  {
    _id = _num_instances++;
    _domain=union_of(f1._domain,f2._domain,&var);
  }

  /** Add table filled with default value */
  void SetFactor::add_table()
  {
    _values.push_back( std::vector<double> (_domain->size(), MSP_DEFAULT_VALUE) );
    _labels.push_back( std::string() );
    if ((_values.size() % 100000)==0) std::cerr << "+ num_tables: " << _values.size()/1000 << "k\tdim:" << _domain->size() << std::endl;
    if (_values.size() > max_tables) max_tables = _values.size();
  }

  /** Pruning */
  void SetFactor::reduce()
  {
    pruning();
    bucketing();
  }


  /** Pareto dominance pruning. */
  void SetFactor::pruning()
  {
    if (_labels.size() < MSP_MIN_NUM_TABLES) return;

    std::vector< std::vector<double> > values;
    std::vector< std::string > labels;
    bool dominated;
    for (unsigned t1=0; t1<_values.size(); ++t1)
      { // check if table t1 is dominated by some other table
    	dominated = false; 
    	for (unsigned t2=0; t2<_values.size(); ++t2)
    	  if (t1 != t2)
    	    { 
    	      unsigned i;
    	      // check if t1 is strictly greater than t2 somewhere
    	      for (i=0; i<_domain->size(); ++i)
    		if (_values[t1][i] > _values[t2][i]) // t2 does not dominate t1
    		  break;
    	      if (i == _domain->size()) 
    		{ // t2[i] >= t1[i] for all i
    		  // check if t1 equals t2
    		  for (i=0; i<_domain->size(); ++i)
    		    if (_values[t1][i] != _values[t2][i])
    		      break;
    		  if (i < _domain->size() || t1 < t2) { 
    		    dominated = true; break; 
    		  } // t2 dominates t1
    		}
    	    }
    	if (!dominated)
    	  { // add undominated table
    	    values.push_back(std::vector<double> ( _values[t1].begin(),_values[t1].end()) );
    	    labels.push_back(std::string (_labels[t1]) );
    	  }
      }
    
    std::swap(_values, values);
    std::swap(_labels, labels);    
  }  

  /** Bucket pruning. */
  void SetFactor::bucketing()
  {
    if (num_tables() < MSP_MAX_NUM_TABLES) return;

    unsigned scaling = std::pow(MSP_MAX_NUM_TABLES,1.0/_domain->size());

    std::vector< std::vector<double> > values;
    std::vector< std::string > labels;


    TrieSet t(scaling);
    for (unsigned t1=0; t1<_labels.size(); ++t1) // check if table t1 maps into same backet as other table already included
    	if (!t.count(_values[t1])) {
    	    t.insert(_values[t1]);
    	    values.push_back( std::vector<double> ( _values[t1].begin(), _values[t1].end()) );
    	    labels.push_back( std::string (_labels[t1]) );
	    if (labels.size() >= MSP_MAX_NUM_TABLES) break; 
    	}     

    std::swap(_values, values);
    std::swap(_labels,labels);    
  }
  

  /** Default printing. */
  std::ostream& operator<<(std::ostream &o, SetFactor &f) {
    o << "SetFactor(tables:" << f.num_tables() << ", " << *f._domain;
    if (f.size() <= MSP_DISP_MAX_VALS)
      {
	o << ", values:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  {
	    o << "[";
	    for (unsigned j=0; j < f.domain_size(); ++j) 
	      { o << f.get(i,j); if (j < f.domain_size()-1) o << ", "; }
	    o << "]";
	  }
	o << "}, labels:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  o << "[" << f.label(i) << "]";
	o << "}";
      }
    o << ")";
    return o;
  }

  std::ostream& operator<<(std::ostream &o, const SetFactor &f) {
    o << "SetFactor(tables:" << f.num_tables() << ", " << *f._domain;
    if (f.size() <= MSP_DISP_MAX_VALS)
      {
	o << ", values:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  {
	    o << "[";
	    for (unsigned j=0; j < f.domain_size(); ++j) 
	      { o << f.get(i,j); if (j < f.domain_size()-1) o << ", "; }
	    o << "]";
	  }
	o << "}, labels:{";
	for (unsigned i=0; i < f.num_tables(); ++i)
	  o << "[" << f.label(i) << "]";
	o << "}";
      }
    o << ")";
    return o;
  }

}
