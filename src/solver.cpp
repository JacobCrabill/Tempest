/*!
 * \file solver.cpp
 * \brief Class to store all solution data & apply FR operators
 *
 * \author - Jacob Crabill
 *           Aerospace Computing Laboratory (ACL)
 *           Aero/Astro Department. Stanford University
 *
 * \version 0.0.1
 *
 * Flux Reconstruction in C++ (Flurry++) Code
 * Copyright (C) 2014 Jacob Crabill.
 *
 */

#include "../include/solver.hpp"

#include <sstream>
#include <omp.h>

#include "input.hpp"
#include "geo.hpp"

solver::solver()
{
  //tg = NULL;
}

solver::~solver()
{
//  if (tg != NULL)
//    delete tg;
}

void solver::setup(input *params, geo *Geo)
{
  this->params = params;
  this->Geo = Geo;
  //this->tg = Geo->tg; // Geo will have initialized this already if needed

  params->time = 0.;

  if (params->equation == NAVIER_STOKES) {
    if (params->nDims == 2)
      params->nFields = 4;
    else if (params->nDims == 3)
      params->nFields = 5;
  }

  nDims = params->nDims;
  nFields = params->nFields;

  gridID = Geo->gridID;
  gridRank = Geo->gridRank;
  nprocPerGrid = Geo->nprocPerGrid;

  /* Setup the final dual-mesh data structures */
  setupDualMesh();

  if (params->restart)
    readRestartFile();
  else
    initializeSolution();

  if (params->meshType == OVERSET_MESH)
    setupOverset();

  /* Additional Setup */

  // Time advancement setup
  switch (params->timeType) {
    case 0:
      nRKSteps = 1;
      RKb = {1};
      break;
    case 4:
      nRKSteps = 4;
      RKa = {.5, .5, 1.};
      RKb = {1./6., 1./3., 1./3., 1./6.};
      break;
    default:
      FatalError("Time-Stepping type not supported.");
  }

  F.setup(nEdges,nDims,nFields);

  divF.setup(nRKSteps,nVerts,nFields);

  if (params->viscous) {
    gradU.setup(nVerts,nDims,nFields);
  }

//#ifndef _NO_MPI
//  finishMpiSetup();
//#endif
}

void solver::update(void)
{
  // For RK time-stepping, store the starting solution values
  if (nRKSteps>1)
    copyU_U0();

  /* Intermediate residuals for Runge-Kutta time integration */

  for (int step=0; step<nRKSteps-1; step++) {

    if (step == 0)
      params->rkTime = params->time;
    else
      params->rkTime = params->time + RKa[step-1]*params->dt;

    moveMesh(step);

    calcResidual(step);

    /* If in first stage, compute CFL-based timestep */
    if (step == 0 && params->dtType == 1) calcDt();  // -- NOT CONSISTENT WITH MOVING MESH SEQUENCE --

    timeStepA(step);

  }

  /* Final Runge-Kutta time advancement step */

  if (nRKSteps == 1) {
    params->rkTime = params->time;
    /* Calculate CFL-based timestep */
    if (params->dtType == 1) calcDt();
  }
  else {
    params->rkTime = params->time + params->dt;
  }

  moveMesh(nRKSteps-1);

  calcResidual(nRKSteps-1);

  // Reset solution to initial-stage values
  if (nRKSteps>1)
    copyU0_U();

  for (int step=0; step<nRKSteps; step++)
    timeStepB(step);

  params->time += params->dt;
}


void solver::calcResidual(int step)
{

#ifndef _NO_MPI
  doCommunication();
#endif

  if (params->viscous)
    calcGradU();

  calcInviscidFlux();

  applyBoundaryConditions();

#ifndef _NO_MPI
  calcInviscidFlux_mpi();
#endif

  if (params->viscous) {

    calcViscousFlux();

#ifndef _NO_MPI
    calcViscousFlux_mpi();
#endif

  }

  calcFluxDivergence(step);

}

void solver::setupDualMesh(void) {

  if (params->rank==0) cout << "Solver: Setting up dual mesh arrays" << endl;

  nVerts = Geo->nVerts;
  nEdges = Geo->nEdges;

  xv = Geo->xv;
  v2ne = Geo->v2nv;
  v2e = Geo->v2e;
  v2v = Geo->v2v;
  Ae = Geo->e2A;
  bcList = Geo->bcList;
  bndNorm = Geo->bndNorm;
  bndPts = Geo->bndPts;
  nBndPts = Geo->nBndPts;

  // Solution at vertices
  U.setup(nVerts,nFields);

  // Flux at faces
  F.setup(nEdges,nFields);

  // Dual-mesh face areas
  A.setup(nVerts,getMax(v2ne));
  for (int i=0; i<nVerts; i++) {
    for (int j=0; j<v2ne[i]; j++) {
      A(i,j) = Ae[v2e(i,j)];
    }
  }
}

void solver::finishMpiSetup(void)
{
//  if (params->rank==0) cout << "Solver: Setting up MPI face communicataions" << endl;
}


void solver::calcInviscidFlux(void)
{

}

void solver::doCommunication()
{

}

void solver::calcInviscidFlux_mpi()
{

}

void solver::calcViscousFlux(void)
{

}

void solver::calcViscousFlux_mpi()
{

}


void solver::calcFluxDivergence(int step)
{
#pragma omp parallel for
  for (int i=0; i<nVerts; i++) {
    for (int k=0; k<nFields; k++) divF(step,i,k) = 0;

    for (int j=0; j<v2ne[i]; j++) {
      // Get properly-scaled normal vector
      int iv2 = v2v(i,j);
      Vec3 norm = point(xv[iv2]) - point(xv[i]);
      norm.norm();
      norm *= A(i,j);

      for (int k=0; k<nFields; k++) {
        for (int dim=0; dim<nDims; dim++) {
          divF(step,i,k) += F(v2e(i,j),dim,k)*norm[dim];
        }
      }
    }
  }
}

void solver::calcGradU(void)
{

}

//void solver::calcEntropyErr(void)
//{

//}

void solver::moveMesh(int step)
{

}

void solver::calcDt(void)
{
//  double dt = INFINITY;

//#pragma omp parallel for reduction(min:dt)
//  for (uint i=0; i<eles.size(); i++) {
//    dt = min(dt, eles[i].calcDt());
//  }

//#ifndef _NO_MPI
//  double dtTmp = dt;
//  MPI_Allreduce(&dtTmp, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//#endif

//  params->dt = dt;
}


void solver::timeStepA(int step)
{
#pragma omp parallel for
  for (int i=0; i<nVerts; i++) {
    for (int j=0; j<nFields; j++) {
      U(i,j) = U0(i,j) - RKa[step]*params->dt*divF(step,i,j)/vol[i];
    }
  }
}


void solver::timeStepB(int step)
{
#pragma omp parallel for
  for (int i=0; i<nVerts; i++) {
    for (int j=0; j<nFields; j++) {
      U(i,j) -= RKb[step]*params->dt*divF(step,i,j)/vol[i];
    }
  }
}

void solver::copyU_U0(void)
{
  U0 = U;
}

void solver::copyU0_U(void)
{
  U = U0;
}

void solver::readRestartFile(void) {

  ifstream dataFile;
  dataFile.precision(15);

  // Get the file name & open the file
  char fileNameC[50];
  string fileName = params->dataFileName;
#ifndef _NO_MPI
  /* --- All processors write their solution to their own .vtu file --- */
  sprintf(fileNameC,"%s_%.09d/%s_%.09d_%d.vtu",&fileName[0],params->restartIter,&fileName[0],params->restartIter,params->rank);
#else
  sprintf(fileNameC,"%s_%.09d.vtu",&fileName[0],params->restartIter);
#endif

  if (params->rank==0) cout << "Solver: Restarting from " << fileNameC << endl;

  dataFile.open(fileNameC);

  if (!dataFile.is_open())
    FatalError("Cannont open restart file.");

  // Find the start of the UnstructuredData region
  bool found = false;
  string str;
  while (getline(dataFile,str)) {
    stringstream ss;
    ss.str(str);
    ss >> str;
    if (str.compare("<UnstructuredGrid>")==0) {
      found = true;
      break;
    }
  }

  if (!found)
    FatalError("Cannot fine UnstructuredData tag in restart file.");

  dataFile.close();

  if (params->rank==0) cout << "Solver: Done reading restart file." << endl;

  // Finish setting up MPI faces
  if (params->rank==0 && params->nproc>1) cout << "MPIFace: Setting up MPI face communications." << endl;

}

void solver::initializeSolution()
{
  if (params->rank==0) cout << "Solver: Initializing Solution... " << flush;

  for (uint i=0; i<nVerts; i++) {
    switch (params->icType) {
    case 0: // Uniform Flow
      U(i,0) = params->rhoBound;
      U(i,1) = params->rhoBound*params->uBound;
      U(i,2) = params->rhoBound*params->vBound;
      if (nDims == 3)
        U(i,3) = params->rhoBound*params->wBound;
      U(i,nDims+1) = params->pBound/(params->gamma-1) + 0.5*params->rhoBound*params->Uinf*params->Uinf;
      break;

    case 1: // Isentropic Vortex test case
      // Grab from HiFiLES or Flurry
      break;
    }
  }

  /* If running a moving-mesh case and using CFL-based time-stepping,
   * calc initial dt for grid velocity calculation */
  /* If running a moving-mesh case and using CFL-based time-stepping,
   * calc initial dt for grid velocity calculation */
  //if ( (params->motion!=0 || params->slipPenalty==1) && params->dtType == 1 ) {
  if (params->dtType == 1) {

    double dt = INFINITY;

    #pragma omp parallel for reduction(min:dt)
    for (int i=0; i<nVerts; i++) {
      double rho = U(i,0);
      double rhovMagSq = U(i,1)*U(i,1) + U(i,2)*U(i,2);
      if (nDims == 3)
        rhovMagSq += U(i,3)*U(i,3);
      double p = (params->gamma-1)*U(i,nDims+1) - 0.5*rhovMagSq/rho;
      double a = sqrt(max(params->gamma*p/rho,0));
      //dx = ??;
      dt = std::min(dt, params->CFL*dx/a);
    }

#ifndef _NO_MPI
  double dtTmp = dt;
  MPI_Allreduce(&dtTmp, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

    params->dt = dt;
  }

  if (params->rank == 0) cout << "done." << endl;
}

/* ---- Overset-Grid Functions ---- */

void solver::setupOverset(void)
{
  //setupOversetData();

  // Give TIOGA a pointer to this solver for access to callback functions
  //tg->setcallback(this);

  //Geo->updateOversetConnectivity();
}

void solver::callDataUpdateTIOGA(void)
{
  //tg->dataUpdate(params->nFields,U.getData(),0);
}
