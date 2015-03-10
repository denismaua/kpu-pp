/** Utility functions */

#include <sys/resource.h>
#include <sys/times.h>

#include "utils.h"

namespace msp {

  // returns cpu time usage
  double get_utime()
  {
    struct rusage rusage;
    if ( getrusage( RUSAGE_SELF, &rusage ) != -1 )
      return (double)rusage.ru_utime.tv_sec +
	(double)rusage.ru_utime.tv_usec / 1000000.0;
    else return -1;
  }


  // returns memory usage
  unsigned get_mem_usage()
  {
    /* maximum physical memory use in bytes*/
    struct rusage r_usage;
    getrusage(RUSAGE_SELF,&r_usage);
    return (unsigned)r_usage.ru_maxrss;

    // if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    // 	/* Linux ---------------------------------------------------- */
    // 	long rss = 0L;
    // 	FILE* fp = NULL;
    // 	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
    // 		return (size_t)0L;		/* Can't open? */
    // 	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    // 	{
    // 		fclose( fp );
    // 		return (size_t)0L;		/* Can't read? */
    // 	}
    // 	fclose( fp );
    // 	return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

    // #else
    // 	/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    // 	return (size_t)0L;			/* Unsupported. */
    // #endif
  } 

}
