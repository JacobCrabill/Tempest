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
#pragma once

#include <memory>
#include <map>
#include <set>
#include <vector>

class oper;

#include "global.hpp"

#include "geo.hpp"
#include "input.hpp"

//class tioga;

//#include "tioga.h"

class solver
{
friend class geo; // Probably only needed if I make eles, opers private?

public:
  /* === Member Variables === */
  //! Pointer to geometry object for mesh-related operations
  geo *Geo;

  //! Pointer to Tioga object for processing overset grids
  //tioga* tg;

  /* === Setup Functions === */

  solver();

  ~solver();

  //! Setup the solver with the given simulation parameters & geometry
  void setup(input *params, geo *Geo);

  //! Setup the FR operators for all ele types and polynomial orders which will be used in computation
  void setupOperators();

  //! Finish setting up all dual-mesh data structures
  void setupDualMesh();

  //! If restarting from data file, read data and setup eles & faces accordingly
  void readRestartFile();

  //! Finish setting up the MPI faces
  void finishMpiSetup(void);

  /* === Functions Related to Basic FR Process === */

  //! Apply the initial condition to all elements
  void initializeSolution();

  void update(void);

  //! Perform one full step of computation
  void calcResidual(int step);

  //! Calculate the stable time step limit based upon given CFL
  void calcDt(void);

  //! Advance solution in time - Generate intermediate RK stage
  void timeStepA(int step);

  //! Advance solution in time - Final RK stage [assemble intermediate stages]
  void timeStepB(int step);

  //! For RK time-stepping - store solution at time 'n'
  void copyU_U0(void);

  //! For RK time-stepping - recover solution at time 'n' before advancing to 'n+1'
  void copyU0_U(void);

  //! Calculate the inviscid flux at all edges
  void calcInviscidFlux(void);

  //! Apply boundary conditions to all boundary nodes
  void applyBoundaryConditions(void);

  //! Calculate the inviscid interface flux at all MPI edges
  void calcInviscidFlux_mpi(void);

  //! Have all MPI faces begin their communication
  void doCommunication(void);

  //! Calculate the gradient of the solution at the solution points
  void calcGradU(void);

  //! Calculate the viscous flux at all edges
  void calcViscousFlux(void);

  //! Calculate the viscous interface flux at all MPI edges
  void calcViscousFlux_mpi(void);

  //! Wrapper to calc divergence of flux (using one of two possible methods)
  void calcFluxDivergence(int step);

  //! Apply mesh motion
  void moveMesh(int step);


  /* === Functions for Shock Capturing & Filtering=== */

  //! Use concentration sensor + exponential modal filter to capture discontinuities
  void shockCapture(void);

  /* === Functions Related to Adaptation === */

  //! Calculate an entropy-adjoint-based error indicator
  void calcEntropyErr_spts();

  /* === Functions Related to Overset Grids === */

  //! Do initial preprocessing for overset grids
  void setupOverset();

  /*!
   * \brief Initialize overset-related data storage
   *
   * Allocates global storage for overset data. Re-call if adding elements,
   * changing polynomial orders, etc.
   */
  void setupOversetData();

  //! Perform the overset data interpolation using TIOGA high-order
  void callDataUpdateTIOGA();

//private:
  //! Pointer to the parameters object for the current solution
  input *params;

  //! Set of all element types present in current simulation
  set<int> eTypes;

  //! Lists of cells to apply various adaptation methods to
  vector<int> r_adapt_cells, h_adapt_cells, p_adapt_cells;

  int nRKSteps;

  int gridID, gridRank, nprocPerGrid;

  vector<double> RKa, RKb;

  /* --- Solution  Variables --- */

  int nDims;
  int nFields;
  int nVerts;
  int nEdges;

  matrix<double> U;    //! Global solution array for solver (at vertices)
  matrix<double> U0;   //! Global solution array at beginning of RK step
  matrix<double> A;    //! Area of dual-mesh faces around each vertex 
  vector<double> vol;  //! Volume of dual-mesh elements
  Array<double,3> F;   //! Global flux array for solver (at edges)
  matrix<double> Fn;   //! Global normal flux array for solver (at edges)
  Array<double,3> divF;  //! Divergence of flux (at vertices) (for each RK step)
  Array<double,3> gradU; //! Gradient of solution (at vertices)
  vector<double> waveSp; //! Wave speed at each edge

  matrix<double> tempFL, tempFR;

  matrix<double> xv;   //! Coordinates of each vertex

  matrix<int> normDir; //! Direction of normal flux for each each around each vertex (+1 or -1)
  vector<double> Ae;   //! Face area for all dual-mesh faces
  vector<int> v2ne;    //! Number of edges (or dual-mesh faces) for each vertex
  matrix<int> v2e;     //! Vertex to edge connectivity
  matrix<int> v2v;     //! Vertex to vertex connectivity

  int nBounds;
  vector<int> bcList;  //! Boundary condition for each boundary
  matrix<int> bndPts;  //! List of boundary nodes on each boundary
  vector<int> nBndPts;  //! List of boundary nodes on each boundary
  //vector<int> bcType;  //! Boundary condition for each boundary node
  Array<Vec3,2> bndNorm;  //! Unit normal at all boundary nodes

  /* ---- Overset Grid Variables / Functions ---- */


};
