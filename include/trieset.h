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

/** Interface for the TrieSet class */

#ifndef MSP_TRIESET_H
#define MSP_TRIESET_H

#include "constants.h"
#include <unordered_map>
#include <vector>

namespace msp {
    
  
  /** Trie-based set data structure whose keys are vectors of integers. */
  class TrieSet
  {

  private:
    /** Defines a node item of a TrieSet. */
    class Node {
    public:
      std::unordered_map<unsigned, Node* > children; /*< labeled outgoing arcs. */
    };

    Node* _root; /*< root node. */
    unsigned _scaling; /*< for number to bucket conversion. */

  public:

    /** Default Constructor. */
  TrieSet(unsigned scaling): _scaling(scaling), _root(new Node) { }
    
    /** Searches for an element of given key. 
     *
     *  Returns 0 if no such element exists and 1 otherwise.
     */
    unsigned count(std::vector<double>& key);
    
    /** Inserts an element of a given key if it does not exist otherwise does nothing. */
    void insert(std::vector<double>& key);

  };
}

#endif
