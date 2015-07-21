/*!
 * \file matrix.hpp
 * \brief Header file for matrix class
 *
 * The matrix class is really just a handy wrapper for a vector<vector<T>>, with
 * some functions for setting up and multiplying
 *
 * Yes, I'm well aware that I really should just be using boost::multi_array
 * instead for loads of reasons, but this was more fun.  Don't judge me.
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

#include <array>
#include <iomanip>   // for setw, setprecision
#include <iostream>
#include <vector>

#include "error.hpp"

struct point;

#include "global.hpp"

#define INSERT_AT_END -1

using namespace std;

typedef unsigned int uint;

template <typename T, uint N>
class matrixBase
{
public:
  /* --- Constructors, destructors --- */
  //! Default Constructor
  matrixBase();

  //! Secondary Constructor with Size Allocation
  matrixBase(uint inDim0, uint inDim1=1, uint inDim2=1, uint inDim3=1);

  //! Copy Constructor
  matrixBase(const matrixBase<T,N>& inMatrix);

  //! Assignment
  matrixBase<T,N> operator=(const matrixBase<T,N>& inMatrix);

  /*! Get dim0 [number of rows] */
  uint getDim0(void) {return this->dims[0];}

  /*! Get dim1 [number of columns] */
  uint getDim1(void) {return this->dims[1];}

  /*! Get the size of the underlying data array (total number of matrix elements) */
  uint getSize(void) {return data.size();}

  /* --- Member Functions --- */

  void setup(uint inDim0, uint inDim1=1, uint inDim2=1, uint inDim3=1);

  /*! Prints the contents of the matrix to the console */
  void print(void);

  /* --- Data-Access Operators --- */

  /*! Returns a pointer to the first element of row inDim0 */
  T* operator[](int inDim0);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j=0, int k=0, int l=0);

  /*! Returns the .data() pointer of the underlying vector<T> data */
  T* getData();

  /* --- Search Operations --- */

  /* --- Member Variables --- */
  //uint dim0, dim1;
  uint nDims;
  array<uint,4> dims;  //! Dimensions of the matrix

  vector<T> data;
};

template <typename T, uint N>
class Array : public matrixBase<T,N> { };

template <typename T>
class matrix : public matrixBase<T,2> {
public:

  matrix();

  matrix(uint inDim0, uint inDim1);

  /*! Standard (i,j) access operator */
  T &operator()(int i, int j);

  void initializeToZero(void);

  void initializeToValue(T val);

  //! Insert a column at the end of the matrix
  void addCol(void);

  //! Insert a column at the end of the matrix
  void addCols(int nCols);

  /*! Insert a row into the matrix at location rowNum [zero-indexed], with the default being at the end */
  void insertRow(const vector<T> &vec, int rowNum = -1);

  void insertRowUnsized(const vector<T> &vec);

  void insertRow(T* vec, uint rowNum, uint length);

  void insertRowUnsized(T* vec, int length);

  //! Remove columns from the end of the matrix
  void removeCols(int nCols = 1);

  vector<T> getRow(uint row);

  //! Return a sub-matrix view of the matrix using the rows in ind
  matrix<T> getRows(vector<int> ind);

  vector<T> getCol(int col);

  /*! Find all unique 'rows' in a matrix */
  void unique(matrix<T> &out, vector<int> &iRow);
};
