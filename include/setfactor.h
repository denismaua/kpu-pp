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

/** Interface for the SetFactor class */

#ifndef MSP_SETFACTOR_H
#define MSP_SETFACTOR_H

#include "constants.h"
#include "variable.h"
#include "domain.h"
#include "factor.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

namespace msp {

  /**  Class of set-valued factors.
   *
   * A set-factor is a collectiong of real-valued factors, all with the
   * same domain. Algebraic operations such as product, sum, and
   * variable elimination are defined on set-factors.
  */
  class SetFactor {

  private:
    unsigned _id; /*< factor id. */
    Domain* _domain; /*< factor domain. */
    std::vector< std::vector<double> > _values; /*< factor image. */
    std::vector< std::string > _labels; /*< for retrieving the solution */
    std::string _name; /*< factor name (optional). */
    static unsigned _num_instances; /*< total number of factor objects instantiated. */
    
  public:
    static unsigned max_tables; /*< maximum number of tables in a pontential (for generating statistics). */

    /** Constructor for factor given the domain scope, that is, a vector of (distinct) variables */
    SetFactor(std::vector<Variable* >& scope) :  _id( _num_instances++ ), _domain ( new Domain(scope) ) {}

    /* Construct SetFactor from Factor */
    SetFactor(const Factor& f);
    
    /** Constructs a factor whose domain is the union of two given factors */
    SetFactor(const SetFactor& f1, const SetFactor& f2);

    /** Constructs a factor whose domain is the union of two given factors with a given variable removed */
    SetFactor(const SetFactor& f1, const SetFactor& f2, const Variable& v);

    /** Getters */
    unsigned id() const { return _id; } 
    unsigned size() const { return _domain->size()*_values.size(); }
    unsigned logsize() const { return log2(_domain->size())+log2(_values.size()); }
    unsigned num_tables() const { return _values.size(); }
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

    /** Get variable at localtion i in the scope.
     * @param i an integer
     * @return a pointer to the variable object
     */
    Variable* var_at(unsigned i) const { return (*_domain)[i]; }

    /** Add table filled with default value */
    void add_table();

    /** set j-th value of i-th table.
     * @param i an integer specifying a table.
     * @param j an integer specifying a position in the table.
     * @param value a double specifying the new value.
     */
    void set(unsigned i, unsigned j, double value) { 
      if (i < _values.size() && j < _domain->size()) _values[i][j] = value;
      else throw "ERROR! Index out of range.";
    }

    /** get j-th value of i-th table.
     * @param i an integer specifying a table.
     * @param j an integer specifying a position in the table.
     * @return a double
     */
    double get(unsigned i, unsigned j) const { 
      if (i < _values.size() && j < _domain->size()) return _values[i][j];
      else throw "ERROR! Index out of range";
    }

    /** Returns the maximum value in the domain (over all tables).
     * @return a double representing the largest value
     */
    double max() const {
      double m=-INFINITY;
      for (unsigned i=0; i<_values.size(); i++)
	for (unsigned j=0; j<_domain->size(); j++)
	  if (_values[i][j] > m) m = _values[i][j];
      return m;
    }

    /** Returns the minimum value in the domain (over all tables).
     * @return a double representing the largest value
     */
    double min() const {
      double m=INFINITY;
      for (unsigned i=0; i<_values.size(); i++)
	for (unsigned j=0; j<_domain->size(); j++)
	  if (_values[i][j] < m) m = _values[i][j];
      return m;
    }


    /** Returns the maximum value in the domain and its position (table, index).
     * @return a double representing the largest value
     */
    std::pair< double, std::pair<unsigned,unsigned> > max2() const {
      double m=-INFINITY;
      std::pair<unsigned, unsigned> pos;
      for (unsigned i=0; i<_values.size(); i++)
	for (unsigned j=0; j<_domain->size(); j++)
	  if (_values[i][j] > m) { m = _values[i][j]; pos.first = i; pos.second=j; }
      return std::make_pair(m,pos);
    }

    /** set the label of the i-th table.
     * @param i an integer specifying a table.
     * @param i an integer specifying a table.
     * @param s a string
     */
    void set_label(unsigned i, std::string s) { 
      if (i < _labels.size()) _labels[i].assign(s);
      else throw "ERROR! Index out of range";
    }

    /** set the label of the i-th table.
     * @param i an integer specifying a table.
     * @param i an integer specifying a table.
     * @param s an array of chars
     */
    void set_label(unsigned i, char* s) { 
      if (i < _labels.size()) _labels[i].assign(s);
      else throw "ERROR! Index out of range";
    }

    /** get the label of the i-th table.
     * @param i an integer specifying a table.
     * @param a string
     */
    std::string label(unsigned i) const { 
      if (i < _labels.size()) return _labels[i];
      else throw "ERROR! Index out of range";
    }


    std::string name() const { return _name; }
    static unsigned get_num_instances() { return _num_instances; }
    
    /** returns i-th value.
     * @param i an integer specifying a variable index.
     * @return the value at the i-th position of the linearized array.
     */
    /* double operator[](unsigned i) const {  */
    /*   if (i < _values.size()) return _values[i]; */
    /*   else throw "Index out of range."; */
    /* } */


    /** Attemps to reduce number of tables
     *
     * Apply a pruning technique (currently pruning by Pareto dominance or bucketing).
     *
     */
    void reduce();

    /** Removes Pareto-dominated tables.
     *
     * A vector v Pareto-dominates a different vector u if v[i] >= u[i] for all i.
     *
     */
    void pruning();

    
    /** Bucket pruning.
     *
     * Partitions the tables into MSP_MAX_NUM_TABLES buckets and discards all but one table from bucket.
     */
    void bucketing();    

    /** output factor content info. */
    friend std::ostream &operator<<(std::ostream &o, SetFactor &v); 
    friend std::ostream &operator<<(std::ostream &o, const SetFactor &v); 

  };

}

#endif
