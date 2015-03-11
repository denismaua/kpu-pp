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

/*@file limid utilities showroom */

#include <iostream>
#include <string>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include "limid.h"
#include "utils.h"
#include "valuation.h"
#include "colormod.h" // namespace Color

using namespace std;
using namespace msp;

int main(int argc, char* argv[]) {
  
  unsigned verbose=0; // defines feedback level (0 min, 3 max)
  Color::Modifier red(Color::FG_RED);
  Color::Modifier green(Color::FG_GREEN);
  Color::Modifier blue(Color::FG_BLUE);
  Color::Modifier def(Color::FG_DEFAULT);
  
  /* select input: if filename is passed as argument, read input from file
     otherwise, read input from stdin */
  ifstream input;
  streambuf* orig_cin = 0;
  if (argc < 3) {
    cerr << red << "Usage: "<< def << argv[0] << blue << " model k [verbosity]" << def << endl;
    cerr << "   -" <<  blue << " model" << def << " is an string specifying the filename of a LIMID format." << endl;
    cerr << "   -" << blue << " k" << def << " is an integer specifying the size of the neighborhood in local search" << endl;
    cerr << "         [use k=1 for SPU and k=n, n the number of decision variables for MPU]" << endl;
    cerr << "   -" << blue << " verbosity" << def << " is an optional integer specifying the level of feedback"<< endl;
    return 0;
  }
  if (argc > 3) verbose = atoi(argv[3]);
   
  // redirect file to cin (because load function reads from stdin)
  input.open(argv[1]);
  if (!input) { cerr << red << "Unable to open file. Aborting..." << def << endl; return 1; }
  orig_cin = cin.rdbuf(input.rdbuf());
  cin.tie(0); // tied to cout by default

  Limid L;
  L.load(); // read model from file or std_in
  if (verbose) 
    { 
      cout << "Model loaded: " << green << L.num_chance_nodes() << def << " chance variables, " << green << L.num_decision_nodes() << def << " decision nodes, and " << green << L.num_value_nodes() << def << " value nodes. " << endl;
      //t.stop(); t.report(); t.resume();
    }

  unsigned k = atoi(argv[2]);
  if (k>L.num_decision_nodes()) {
    cerr << red << "Search-neighborhood size exceeds number of decision variables, using maximum allowed instead..." << def << endl;
    k = L.num_decision_nodes();
  }

  double start, end; // time measurement
  start = get_utime();
  Limid::solution_t sol = L.solve(k,verbose); // solve limid
  end = get_utime();
  cout << blue;
  if (verbose) printf("\n%40s\t%20s\t%10s\t%10s\t%15s\t%15s\t%20s\t%20s\n", "filename", "method", "# chance", "# decision", "# iters", "# tables", "runtime", "value");
  cout << green;
  printf( "%40s\t%20s\t%10d\t%10d\t%15d\t%15d\t%20s\t%20s\n", argv[1], sol.method.c_str(), L.num_chance_nodes(), L.num_decision_nodes(), sol.num_iters, sol.num_tables, to_string(end-start).c_str(), to_string(sol.value).c_str() );
  cout << def;
  return 0;
}
