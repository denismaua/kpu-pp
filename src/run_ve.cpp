// Copyright (c) 2014 Denis Maua
// All Rights Reserved.
//
// This file is part of the MSP library
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
// along with MLS.  If not, see <http://www.gnu.org/licenses/>.

/*@file showroom */
#define VERBOSITY 1 /*< define the level of output verbosity */

//#include <boost/timer/timer.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include "io.h"
#include "variable.h"
#include "factor.h"
#include "operations.h"
#include "graph.h"
#include "inference.h"
#include "limid.h"

using namespace std;
using namespace msp;

int main(int argc, char* argv[]) {

  //boost::timer::auto_cpu_timer t;


  // select input: if filename is passed as argument, read input from file
  // otherwise, read input from stdin
  ifstream input;
  streambuf* orig_cin = 0;
  if (argc >= 2) {
    input.open(argv[1]);
    if (!input) return 1;
    orig_cin = cin.rdbuf(input.rdbuf());
    cin.tie(0); // tied to cout by default
  }

  // parse model from cin (supposedly set appropriately)
  try {
    vector<Variable > variables;     
    vector<Factor > factors;
    load_uai_model(variables, factors); // load model
    if (VERBOSITY > 0) 
      { 
	cout << "Model loaded: " << variables.size() << " variables, " << factors.size() << " factors. " << endl;
	//t.stop(); t.report(); t.resume();
      }

    if (VERBOSITY > 1) 
      {
	for (auto v: variables) 
	  cout << v << endl;
	for (auto f: factors) 
	  cout << f << endl;
	cout << endl;
      }

    // build domain graph
    Graph g(factors);
    if (VERBOSITY > 2) cout << g << endl;
    // find elimination sequence
    g.triangulate();
    if (VERBOSITY > 2) cout << g << endl;
    else if (VERBOSITY > 0) printf("treewidth: %d\t\tweighted_treewidth: %g\r\n\r\n", g.treewidth().second, g.w_treewidth().second);
    
    // do variable elimination
    Factor f = variable_elimination(g.ordering(),factors, VERBOSITY);    
    cout << "=" << f << endl;

  }
  catch (...) {
    if (orig_cin) cin.rdbuf(orig_cin);
    throw;
  }

  return 0;
}
