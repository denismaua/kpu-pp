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
    
    Node n = _root;
    for (double d: key)
      {
	if (!n->children.count(d)) return 0;
	n = n->children[d];
      }
    return 1;
  }

  /** Inserts an element of a given key if it does not exist otherwise does nothing. */
  void TrieSet::insert(std::vector<double>& key)
  {
    for (double d: key)
      {
	if (!n->children.count(d)) n->children[d] = new Node();
	n = n->children[d];
      }
  }

  
algorithm insert(root : node, s : string, value : any):
    node = root
    i    = 0
    n    = length(s)
 
    while i < n:
        if node.child(s[i]) != nil:
            node = node.child(s[i])
            i = i + 1
        else:
            break
 
    (* append new nodes, if necessary *)
    while i < n:
        node.child(s[i]) = new node
        node = node.child(s[i])
        i = i + 1
 
    node.value = value
  
}
