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
        break;
      }

    }
  }
}
