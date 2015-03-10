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

// This file implements operations on factors and variables

#include <vector>
#include "variable.h"
#include "factor.h"
#include "operations.h"
#include "constants.h"

namespace msp {

  /** Factor multiplication 
   *  @param f1 the first factor
   *  @param f2 the second factor
   *  @return the product of f1 and f2
   */
  Factor product(const Factor& f1, const Factor& f2) {

    Factor res(f1,f2); // creates factor with union domain

    unsigned i, j=0, k=0, l;
    std::vector<unsigned> c(res.width(),0);
    Variable* v;

    for (i=0; i<res.size(); i++)
      {
	res.set(i, f1[j] * f2[k]);
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

    return res;

  }


  /** Combined factor multiplication and variable elimination 
   *  @param f1 the first factor
   *  @param f2 the second factor
   *  @param var a pointer to a variable
   *  @return the elimination of var from the product of f1 and f2
   */
  Factor sum_product(const Factor& f1, const Factor& f2, const Variable& var) {

    Factor res(f1,f2,var); // creates factor with union domain setminus var

    unsigned i, j=0, k=0, l, m;
    std::vector<unsigned> c(res.width()+1,0);
    Variable* v;

    for (i=0; i<res.size(); i++)
      {
	for (m=0; m<var.size(); m++) 
	  {
	    res.set(i, res[i] + f1[j] * f2[k]);
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
      }

    return res;

  }

  /** Set-Factor multiplication 
   *  @param f1 the first set-factor
   *  @param f2 the second set-factor
   *  @return the product of f1 and f2
   */
  SetFactor product(const SetFactor& f1, const SetFactor& f2) {

    SetFactor res(f1,f2); // creates factor with union domain

    unsigned i, j=0, k=0, l;
    std::vector<unsigned> c(res.width(),0);
    Variable* v;

    // first cache indexing
    std::vector<unsigned> i1(res.domain_size()), i2(res.domain_size());

    for (i=0; i<res.domain_size(); i++)
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
	  res.add_table();
	  res.set_label( r*f2.num_tables()+s, f1.label(r)+f2.label(s) );
	  for (i=0; i<res.domain_size(); i++)
	    res.set( r*f2.num_tables()+s, i, f1.get(r,i1[i])*f2.get(s,i2[i])  );
	}

    return res;

  }

  /** Combined set-factor multiplication and variable elimination 
   *  @param f1 the first set-factor
   *  @param f2 the second set-factor
   *  @param var a pointer to a variable
   *  @return the elimination of var from the product of f1 and f2
   */
  SetFactor sum_product(const SetFactor& f1, const SetFactor& f2, const Variable& var) {

    SetFactor res(f1,f2,var); // creates factor with union domain setminus var

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
    	  res.add_table();
	  res.set_label( r*f2.num_tables()+s, f1.label(r)+f2.label(s) );
    	  for (i=0; i<res.domain_size(); i++)
    	    for (m=0; m<var.size(); m++) 
    	      res.set( r*f2.num_tables()+s, i, res.get(r*f2.num_tables()+s, i) + f1.get(r,i1[i*var.size()+m]) * f2.get(s,i2[i*var.size()+m]) );
    	}

    return res;

  }

  /** Returns a new set-factor containing only non-dominated tables */
  SetFactor pareto_pruning(const SetFactor& f)
  {
    if (f.num_tables() < MSP_MIN_NUM_TABLES) return f;

    std::vector<Variable* > scope(f.width());
    for (unsigned i=0; i < f.width(); i++)
      scope[i] = f.var_at(i);
    SetFactor p(scope);
    bool dominated;
    for (unsigned t1=0; t1<f.num_tables(); ++t1)
      { // check if table t1 is dominated by some other table
	dominated = false; 
	for (unsigned t2=0; t2<f.num_tables(); ++t2)
	  if (t1 != t2)
	    { 
	      unsigned i;
	      // check if t1 is strictly greater than t2 somewhere
	      for (i=0; i<f.domain_size(); ++i)
		if (f.get(t1,i) > f.get(t2,i)) // t2 does not dominate t1
		  break;
	      if (i == f.domain_size()) 
		{ // t2[i] >= t1[i] for all i
		  // check if t1 equals t2
		  for (i=0; i<f.domain_size(); ++i)
		    if (f.get(t1,i) != f.get(t2,i))
		      break;
		  if (i < f.domain_size() || t1 < t2) { 
		    dominated = true; break; 
		  } // t2 dominates t1
		}
	    }
	if (!dominated)
	  { // add undominated table
	    p.add_table();  
	    p.set_label(p.num_tables()-1, f.label(t1));
	    for (unsigned i=0; i<f.domain_size(); ++i) p.set( p.num_tables()-1, i, f.get(t1, i) );
	  }
      }
    if (!p.num_tables())
      {
	p.add_table();
	p.set_label(p.num_tables()-1, f.label(0));
	for (unsigned i=0; i < f.domain_size(); i++)
	  p.set( 0, i, f.get(0,i) );
      }
    return p;
  }

}
