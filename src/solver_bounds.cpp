#include "solver.hpp"

void solver::applyBoundaryConditions(void)
{
  for (int bnd=0; bnd<nBounds; bnd++) {

#pragma omp parallel for
    for (int i=0; i<nBndPts[bnd]; i++) {
      int iv = bndPts(bnd,i);

      switch (bcList[bnd]) {
      case (CHAR):
        break;

      case SUP_IN:
        U(iv,0) = params->rhoBound;
        U(iv,1) = params->rhoBound*params->uBound;
        U(iv,2) = params->rhoBound*params->vBound;
        if (params->nDims == 3)
          U(iv,3) = params->rhoBound*params->wBound;
        U(iv,params->nDims+1) = params->pBound/(params->gamma-1) + 0.5*U(iv,0)*params->Uinf*params->Uinf;

      case SUP_OUT:
        // Extrapolated, so do nothing
        break;

      case SLIP_WALL:
        Vec3 norm = bndNorm(bnd,i);
        double rho = U(iv,0);
        double u = U(iv,1)/rho;
        double v = U(iv,2)/rho;
        double w = 0;
        if (nDims == 3)
          w = U(iv,3)/rho;
        double vn = u*norm(0)+v*norm(1)+w*norm(2);
        // Reflect velocity at wall
        u -= 2*vn*norm(0);
        v -= 2*vn*norm(1);
        w -= 2*vn*norm(2);
        U(iv,1) = rho*u;
        U(iv,2) = rho*v;
        if (nDims == 3)
          U(iv,3) = rho*w;
        break;
      }

    }
  }
}
