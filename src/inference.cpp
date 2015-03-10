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

/** Implements direct inference algorithms */

/* whether to use Pareto pruninig in set-valued variable elimination */

#include <forward_list>
#include <utility>
#include <iostream>

#include "constants.h"
#include "inference.h"
#include "operations.h"
#include "utils.h"

namespace msp {

  /** variable elimination algorithm.
   *
   * eliminates variables in the ordering given.
   *
   * @param variables a vector of variables
   * @param factors a forward list (singly linked list) of factors
   * @param verbosity an integer
   * @return a factor 
   */
  Factor variable_elimination(const std::vector<Variable >& variables, std::vector<Factor >& factors, int verbosity)
  {

    std::forward_list<Factor > flist(factors.begin(),factors.end()); // pool of factors
    std::forward_list<Factor > bucket;
    for (auto var: variables)
      {
	if (verbosity) std::cout << "-" << var << std::endl;
	// collect all factors with variable *v in their scope and remove them from factor list
	bucket.clear();
	unsigned b = 0; // bucket size
	for (std::forward_list<Factor >::const_iterator pf=flist.before_begin(), f=flist.begin(); f != flist.end(); pf=f, f++)
	  {
	    if ( f->in_scope( var ) ) { 
	      bucket.push_front( std::move( *f ) );
	      b++; f = pf;
	      flist.erase_after(pf);
	    }
	  }      
	// multiply all factors in bucket and eliminate variable *v
	if (b>0) {	    
	    Factor p(1.0);
	    std::forward_list<Factor >::const_iterator f = bucket.begin();
	    while (b>1)
	      {
		p = product(p,*f);
		f++; b--;
	      }
	    flist.push_front( sum_product(p,*f,var) );	    
	    if (verbosity > 1) { f = flist.begin(); std::cout << "+" << *f << std::endl; }	    
	  } 
      }    
    // generate result by multiplying all remaining factors in the pool
    Factor p(1.0); // product-identity factor
    for (std::forward_list<Factor >::const_iterator f = flist.begin(); f != flist.end(); f++)
    	p = product(p,*f);
    return p;
  }

  /** variable elimination algorithm for set-factors.
   *
   * eliminates variables in the ordering given.
   *
   * @param variables a vector of variables
   * @param factors a forward list (singly linked list) of set-factors
   * @param verbosity an integer
   * @return a set-factor 
   */
  SetFactor set_variable_elimination(const std::vector<Variable>& variables, std::vector<SetFactor>& factors, int verbosity)
  {

    std::forward_list<SetFactor> flist(factors.begin(),factors.end()); // pool of factors
    std::forward_list<SetFactor> bucket;
    for (auto var: variables)
      {
	if (verbosity) std::cout << "-" << var << std::endl;
	// collect all factors with variable *v in their scope and remove them from factor list
	bucket.clear();
	unsigned b = 0; // bucket size
	for (std::forward_list<SetFactor >::const_iterator pf=flist.before_begin(), f=flist.begin(); f != flist.end(); pf=f, f++)
	  {
	    if ( f->in_scope( var ) ) { 
	      bucket.push_front( std::move( *f ) );
	      b++; f = pf;
	      flist.erase_after(pf);
	    }
	  }      
	// multiply all factors in bucket and eliminate variable *v
	if (b>0) {	
	  SetFactor p( Factor (1.0));
	    std::forward_list<SetFactor >::const_iterator f = bucket.begin();
	    while (b>1)
	      {		
		p = product(p,*f);
		p.reduce();
		f++; b--;
	      }
	    p = sum_product(p,*f,var);
	    p.reduce();
	    flist.push_front( std::move( p ) );
	    if (verbosity) { 
	      f = flist.begin(); 
	      std::cout << "num_tables:" << f->num_tables() << "\tdomain_size:" << f->domain_size() << "\tlog_size:" << f->logsize() << std::endl;
	      if (verbosity > 1) std::cout << *f << std::endl;
	    }
	  } 
      }    
    // generate result by multiplying all remaining factors in the pool
    SetFactor p(Factor (1.0)); // product-identity set-factor
    for (std::forward_list<SetFactor >::const_iterator f = flist.begin(); f != flist.end(); f++)
      { p = product(p,*f); p.reduce(); }
    return p;
  }

  /** Variable Elimination Algorithm for valuations.
   *
   * Eliminates variables from combination of valuations in the given order.
   *
   * @param variables a vector of variables
   * @param factors a forward list (singly linked list) of valuations
   * @param verbosity an integer
   * @return a valuation
   */
  Valuation variable_elimination(const std::vector<Variable>& variables, std::vector<Valuation>& vals, int verbosity)
  {

    std::forward_list<Valuation> flist(vals.begin(),vals.end()); // pool of valuations
    std::forward_list<Valuation> bucket;
    for (auto var: variables)
      {
	if (verbosity) std::cout << "-" << var << std::endl;
	// collect all valuations with variable *v in their scope and remove them from factor list
	unsigned b = 0; // bucket size
	for (std::forward_list<Valuation>::const_iterator pf=flist.before_begin(), f=flist.begin(); f != flist.end(); pf=f, f++)
	  {
	    if ( f->in_scope( var ) ) { 
	      bucket.push_front( std::move( *f ) );
	      b++; f = pf;
	      flist.erase_after(pf);
	    }
	  }      
	// multiply all valuations in bucket and eliminate variable *v
	if (b>0) {	
	  Valuation p; // identity
	  std::forward_list<Valuation>::const_iterator f = bucket.begin();
	  while (b>1)
	    {		
	      p = product(p,*f);
	      p.reduce();
	      f++; b--;
	    }
	  p = sum_product(p,*f,var);
	  p.reduce();
	  flist.push_front( std::move( p ) );
	  if (verbosity) { 
	    f = flist.begin(); 
	    std::cout << "num_tables:" << f->num_tables() << "\tdomain_size:" << f->domain_size() << "\tlog_size:" << f->logsize() << std::endl;
	    std::cout << "mem usage: " <<  get_mem_usage()/1024 << "MB" << std::endl;
	    if (verbosity > 1) std::cout << *f << std::endl;
	  }
	  bucket.clear();
	}
      }    
    // generate result by multiplying all remaining valuations in the pool
    Valuation p; // product-identity valuation
    for (std::forward_list<Valuation>::const_iterator f = flist.begin(); f != flist.end(); f++)
      { p = product(p,*f); p.reduce(); }
    return p;
  }

}
