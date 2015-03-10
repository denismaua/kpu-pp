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

/** Interface for the Valuation class */

#ifndef MSP_VALUATION_H
#define MSP_VALUATION_H

#include "constants.h"
#include "variable.h"
#include "domain.h"
#include "factor.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

namespace msp {

  /**  Class of valuations (bicriteria set-valued factors).
   *
   * A valuation is a collection of pairs of real-valued factors, all defined on the
   * same domain. Algebraic operations such as product, sum, and variable elimination
   * are defined on valuations.
  */
  class Valuation {

  private:
    unsigned _id; /*< factor id. */
    Domain* _domain; /*< factor domain. */
    std::vector< std::vector<double> > _l_values, _r_values; /*< values. */
    std::vector< std::string > _labels; /*< for retrieving the solution */
    std::string _name; /*< factor name (optional). */
    static unsigned _num_instances; /*< total number of factor objects instantiated. */
    
  public:
    static unsigned max_tables; /*< maximum number of tables in a pontential (for generating statistics). */

    /** Constructor for identity valuation. */
    Valuation() : _id( _num_instances++ ), _domain( new Domain() ), _l_values( 1, std::vector<double> (1,1.0) ), _r_values( 1, std::vector<double> (1,0.0) ), _labels(1, std::string () ) {} /*< unidimensional factor constructor */

    
    /** Constructor given the domain scope, that is, a vector of (distinct) variables */
    Valuation(std::vector<Variable* >& scope) :  _id( _num_instances++ ), _domain ( new Domain(scope) ) {}

    /** Constructor given a temporary domain scope, that is, a vector of (distinct) variables */
    Valuation(std::vector<Variable* >&& scope) :  _id( _num_instances++ ), _domain ( new Domain(scope) ) {}

    /* Constructs Valuation from factor. */
    Valuation(const Factor& f);
    
    /** Constructs a valuation whose domain is the union of two given valuation. */
    Valuation(const Valuation& f1, const Valuation& f2);

    /** Constructs a valuation whose domain is the union of two given valuations with a given variable removed. */
    Valuation(const Valuation& f1, const Valuation& f2, const Variable& v);

    /** Getters */
    unsigned id() const { return _id; } 
    unsigned size() const { return _domain->size()*_l_values.size(); }
    unsigned logsize() const { return log2(_domain->size())+log2(_l_values.size()); }
    unsigned num_tables() const { return _l_values.size(); }
    unsigned domain_size() const { return _domain->size(); }
    unsigned width() const { return _domain->width(); }

    /** Size of i-th variable in scope 
     * @param i an integer
     * @result the size of the i-th variable in the scope
     */
    unsigned size_of_var(unsigned i) { return _domain->size_of_var(i); }

    /** Get the domain offset for the i-th variable in the scope.
     * @param i an integer specifiy a variable index.
     * @return the offset by which values of the i-th variable are shifted when linearizing the array.
     */
    unsigned offset(unsigned i) const { return _domain->offset(i); }

    /** Get the domain offset for a variable.
     * @param v a variable
     * @return the offset by which values of the variable are shifted when linearizing the array or zero if the variable is not in the scope.
     */
    unsigned offset(const Variable& v) const { return _domain->offset(v); }

    /** Determines whether variable is in the scope.
     * @param v a variable
     * @return true if it is, false otherwise
     */
    bool in_scope(const Variable& v) const { return _domain->in_scope(v); }

    /** Gets variable at localtion i in the scope.
     * @param i an integer
     * @return a pointer to the variable object
     */
    Variable* var_at(unsigned i) const { return _domain->var_at(i); }

    /** Adds paired tables filled with default value and label. */
    void add();

    /** Adds given paired tables with given label. */
    void add( std::vector<double>&& left, std::vector<double>&& right, std::string&& str);

    /** Removes all tables. */
    void clear() { _l_values.clear(); _r_values.clear(); _labels.clear(); }
    
    /** set j-th left and right values of i-th table.
     * @param i an integer specifying a table.
     * @param j an integer specifying a position in the table.
     * @param l_value a double specifying the new left value.
     * @param r_value a double specifying the new right value.
     */
    void set(unsigned i, unsigned j, double l_value, double r_value) { 
      if (i < _l_values.size() && j < _domain->size()) { _l_values[i][j] = l_value; _r_values[i][j] = r_value; }
      else throw "Valuation::set(): Index out of range!";
    }    

    /** get j-th left value of i-th table.
     * @param i an integer specifying a table.
     * @param j an integer specifying a position in the table.
     * @return a double
     */
    double l_value(unsigned i, unsigned j) const { 
      if (i < _l_values.size() && j < _domain->size()) return _l_values[i][j];
      else throw "Valuation::l_value: Index out of range!";
    }

    /** get j-th right value of i-th table.
     * @param i an integer specifying a table.
     * @param j an integer specifying a position in the table.
     * @return a double
     */
    double r_value(unsigned i, unsigned j) const { 
      if (i < _r_values.size() && j < _domain->size()) return _r_values[i][j];
      else throw "Valuation::r_value: Index out of range!";
    }


    /** Returns the maximum l-value in the domain and its position (table, index).
     * @return a double representing the largest value
     */
    std::pair< double, std::pair<unsigned,unsigned> > maxl() const {
      double m=-INFINITY;
      std::pair<unsigned, unsigned> pos;
      for (unsigned i=0; i<_l_values.size(); i++)
	for (unsigned j=0; j<_domain->size(); j++)
	  if (_l_values[i][j] > m) { m = _l_values[i][j]; pos.first = i; pos.second=j; }
      return std::make_pair(m,pos);
    }    
    
    /** Returns the maximum r-value in the domain and its position (table, index).
     * @return a double representing the largest value
     */
    std::pair< double, std::pair<unsigned,unsigned> > maxr() const {
      double m=-INFINITY;
      std::pair<unsigned, unsigned> pos;
      for (unsigned i=0; i<_r_values.size(); i++)
	for (unsigned j=0; j<_domain->size(); j++)
	  if (_r_values[i][j] > m) { m = _r_values[i][j]; pos.first = i; pos.second=j; }
      return std::make_pair(m,pos);
    }    

    /** set the label of the i-th table.
     * @param i an integer specifying a table.
     * @param i an integer specifying a table.
     * @param s a string
     */
    void set_label(unsigned i, std::string s) { 
      if (i < _labels.size()) _labels[i].assign(s);
      else throw "Valuation::set_label: Index out of range!";
    }

    /** set the label of the i-th table.
     * @param i an integer specifying a table.
     * @param i an integer specifying a table.
     * @param s an array of chars
     */
    void set_label(unsigned i, char* s) { 
      if (i < _labels.size()) _labels[i].assign(s);
      else throw "Valuation::set_label: Index out of range!";
    }

    /** get the label of the i-th table.
     * @param i an integer specifying a table.
     * @param a string
     */
    std::string label(unsigned i) const { 
      if (i < _labels.size()) return _labels[i];
      else throw "Valuation::label: Index out of range!";
    }


    std::string name() const { return _name; }
    static unsigned get_num_instances() { return _num_instances; }


    /** Attempts to reduce number of tables. */
    void reduce();
    
    /** Removes Pareto-dominated tables.
     *
     * A vector v Pareto-dominates a different vector u if v[i] >= u[i] for all i.
     *
     */
    void do_pruning();


    /** Bucket pruning.
     *
     * Partitions the tables into MSP_MAX_NUM_TABLES sets and discards all but one table from each set. The sets are selected to be hyperrectangles of equal size.
     */
    void do_bucketing();


    /** Clustering.
     *
     * Partitions the tables into MSP_MAX_NUM_TABLES sets and discards all but one table from each set. Differently from do_bucketing, it rans a k-medoids algorithm to determine the partition.
     */
    void do_clustering();

    
    /** output factor content info. */
    friend std::ostream &operator<<(std::ostream &o, Valuation &v); 
    friend std::ostream &operator<<(std::ostream &o, const Valuation &v); 

  };

  /** Combines two valuations. 
   * @param v1 the first valuation.
   * @param v2 the second valuation.
   * @result a valuation containing the combination of v1 and v2.
   */
  Valuation product(const Valuation& v1, const Valuation& v2);

  /** Joint valuation combination and variable elimination 
   *  @param f1 the first valuation
   *  @param f2 the second valuation
   *  @param var a pointer to a variable
   *  @return the elimination of var from the combination of f1 and f2
   */
  Valuation sum_product(const Valuation& f1, const Valuation& f2, const Variable& var);
  
}

#endif
