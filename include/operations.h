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

// This file contains interfaces to operations on factors and variables

#ifndef MSP_OPERATIONS_H
#define MSP_OPERATIONS_H

#include "variable.h"
#include "factor.h"
#include "setfactor.h"

namespace msp {

  /** Multiplies two factors 
   * @param f1 the first factor
   * @param f2 the second factor
   * @result a factor containing the product of f1 and f2
   */
  Factor product(const Factor& f1, const Factor& f2);

  /** Combined factor multiplication and variable elimination 
   *  @param f1 the first factor
   *  @param f2 the second factor
   *  @param var a pointer to a variable
   *  @return the elimination of var from the product of f1 and f2
   */
  Factor sum_product(const Factor& f1, const Factor& f2, const Variable& var);

  /** Multiplies two set-factors 
   * @param f1 the first set-factor
   * @param f2 the second set-factor
   * @result a set-factor containing the product of f1 and f2
   */
  SetFactor product(const SetFactor& f1, const SetFactor& f2);

  /** Combined set-factor multiplication and variable elimination 
   *  @param f1 the first set-factor
   *  @param f2 the second set-factor
   *  @param var a pointer to a variable
   *  @return the elimination of var from the product of f1 and f2
   */
  SetFactor sum_product(const SetFactor& f1, const SetFactor& f2, const Variable& var);

  /** Prunes a set-factor by removing Pareto-dominated tables.
   *
   * A vector v Pareto-dominates a vector u if v[i] >= u[i] for all i
   * and v[i] > u[i] for some i.
   *
   * @param f a set-factor.
   * @return a set-factor with the pareto-maximal tables of f.
   */
  SetFactor pareto_pruning(const SetFactor& f);

}

#endif
