/*
 * Copyright (C) 2018, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori, Ken-Ichi Ishikawa.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

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
