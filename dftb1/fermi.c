#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

void fermi(int ndim, double *ev, double *occ, double *efermi, int nelectrons)
{
  const double degtol=1.0e-4;
  double occdg, nel;
  int i, nef1, nef2, nup, ndown, nocc2, ndeg;

  nel = (double) nelectrons;

  if (nel > 2*ndim) {
    printf("Too many electrons: %f > %d\n", nel, 2*ndim);
    exit(-1);
  }

  // start defining occupation numbers and their derivatives
  for (i=0; i<ndim; i++)
    occ[i] = 0.0;

  if (nel < 1.0e-5)
    return;

  // find energy range for fermi energy
  if (nel > (double)(int)nel) {
    nef1 = (int) (nel+2)/2;
    nef2 = (int) (nel+2)/2;
  } else {
    nef1 = (int) (nel+1)/2;
    nef2 = (int) (nel+2)/2;
  }

  *efermi = 0.5 * (ev[nef1]  + ev[nef2]);
  nup = ndown = nef1;

  for ( ; nup < ndim; nup++)
    if (fabs(ev[nup]-*efermi) > degtol)
      break;

  for ( ; ndown > 0; ndown--)
    if (fabs(ev[ndown-1]-*efermi) > degtol)
      break;

  ndeg = nup - ndown;
  nocc2 = ndown;

  for (i=0; i<nocc2; i++)
    occ[i] = 2.0;

  if (ndeg == 0)
    return;

  // for T = 0, occupy orbitals as usually
  occdg = (nel - 2*nocc2) / ndeg;
  for (i=nocc2; i<nocc2+ndeg; i++)
    occ[i] = occdg;

  return;
}
