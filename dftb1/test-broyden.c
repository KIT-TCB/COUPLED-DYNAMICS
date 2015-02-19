#include<stdio.h>
#include<stdlib.h>
#include<math.h>

static int inverse(int lastm1, double *b, double *work, long *ipiv)
{
  extern void dsytrf_(char *, long *, double *, long *, long *, double *, long *, long *);
  extern void dsytri_(char *, long *, double *, long *, long *, double *, long *);
  static long info, imatsz=IMATSZ_BROYDEN, size;
  static char uplo='U';
  size = (long) lastm1;

  dsytrf_(&uplo, &size, b, &imatsz, ipiv, work, &imatsz, &info);
  printf("dsytrf info: %ld\n", info);
  if (info) return info;
  dsytri_(&uplo, &size, b, &imatsz, ipiv, work, &info);
  printf("dsytri info: %ld\n", info);
  return info;
} 
