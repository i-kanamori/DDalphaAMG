#include "main.h"

void show_version(const global_struct *g){

  if ( g->my_rank == 0 ) {
    printf("\nDDalphaAMG solver library\n");
    printf("version info:\n");
#ifdef _COMMIT_ID
    printf("  %s\n", _COMMIT_ID);
#endif
#ifdef _COMMIT_AUTHOR
    printf("  %s\n", _COMMIT_AUTHOR);
#endif
#ifdef _COMMIT_DATE
    printf("  %s\n", _COMMIT_DATE);
#endif
    printf("\n");

#ifdef PROFILING
    printf("optimization flags:\n");

#ifdef OPENMP
    printf("OPENMP: on\n");
#else
    printf("OPENMP: off\n");
#endif
#ifdef SSE
    printf("SSE: on\n");
#else
    printf("SSE: off\n");
#endif
#ifdef K_OPT
    printf("K_OPT: on\n");
#else
    printf("K_OPT: off\n");
#endif


#endif // profile
  }
}
