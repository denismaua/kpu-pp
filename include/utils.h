/** @file Utility functions */

#ifndef MSP_UTILS_H
#define MSP_UTILS_H

namespace msp {

  /** Return elapsed user cpu time spent since process has started. */
  double get_utime();
  
  /** Returns the maximum physsical meory used so far (in bytes). */
  unsigned get_mem_usage();

}

#endif
