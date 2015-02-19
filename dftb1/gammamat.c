#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

void gammamatrix(int nat, double (*rat)[3], double **gammamat, double uhubb[DFTB_MAXTYPES], int *izp)
{
  int i, j;
  double r[3], norm;

  for (i=0; i<nat; i++)
    for (j=0; j<=i; j++) {
      r[0] = rat[i][0] - rat[j][0];
      r[1] = rat[i][1] - rat[j][1];
      r[2] = rat[i][2] - rat[j][2];
      norm = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
      // get value for gamma
      gammamat[i][j] = gam12(norm, uhubb[izp[i]], uhubb[izp[j]]);
    }
  return;
}

double gam12(double r, double uhub1, double uhub2)
{
  const double zero=1.e-4;
  double a1, a2, src, avg, rrc, fac, efac;

  a1 = 3.2*uhub1;
  a2 = 3.2*uhub2;

  if (a1+a2 < zero)
    return 0.0;

  src = 1.0 / (a1+a2);
  fac = a1*a2*src;
  avg = 1.6 * (fac + fac*fac*src);

  if (r < zero)
    return 0.3125*avg;
  else {
    rrc = 1.0/r;
    if (fabs(a1-a2) < 1.e-5) {
      fac = avg * r;
      efac = exp(-fac)/48.0;
      return (1.0 - (48.0+33*fac+fac*fac*(9.0+fac))*efac) * rrc;
    } else
      return rrc - gamsub(a1, a2, r, rrc) - gamsub(a2, a1, r, rrc);
  }
}

double gamsub(double a, double b, double r, double rrc)
{
  double drc, fac;

  drc= 1.0 / (a*a - b*b);
  fac = (b*b*b*b*b*b - 3*a*a*b*b*b*b) * drc*drc*drc * rrc;
  return exp(-a*r) * (0.5 * a * b*b*b*b * drc*drc - fac);
}
