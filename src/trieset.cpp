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

#include "trieset.h"

namespace msp {

  /** Searchs for element of a given key. */ 
  unsigned TrieSet::count(std::vector<double>& key)
  {
    
    Node* n = _root;
    for (double d: key)
      {
	unsigned s = d*_scaling;
	if (!n->children.count(s)) return 0;
	n = (Node*) n->children[s];
      }
    return 1;
  }

  /** Inserts an element of a given key if it does not exist otherwise does nothing. */
  void TrieSet::insert(std::vector<double>& key)
  {
    Node* n = _root;
    for (double d: key)
      {
	unsigned s = d*_scaling;
	if (!n->children.count(s)) n->children[s] = new Node();
	n = (Node*) n->children[s];
      }
  }
  
}
