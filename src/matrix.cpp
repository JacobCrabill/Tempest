/*!
 * \file matrix.cpp
 * \brief Class for simplified matrix storage & manipulation
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
#include "../include/matrix.hpp"

#include <set>

template<typename T, uint N>
matrixBase<T,N>::matrixBase()
{
  data.resize(0);
  dims = {{0,0,0,0}};
}

template<typename T, uint N>
matrixBase<T,N>::matrixBase(uint inDim0, uint inDim1, uint inDim2, uint inDim3)
{
  data.resize(inDim0*inDim1*inDim2*inDim3);
  dims = {{inDim0,inDim1,inDim2,inDim3}};
}

template<typename T>
matrix<T>::matrix()
{

}

template<typename T>
matrix<T>::matrix(uint inDim0, uint inDim1)
{
  this->data.resize(inDim0*inDim1);
  this->dims = {{inDim0,inDim1,1,1}};
}

template<typename T, uint N>
matrixBase<T,N>::matrixBase(const matrixBase<T,N> &inMatrix)
{
  data = inMatrix.data;
  dims = inMatrix.dims;
}

template<typename T, uint N>
matrixBase<T,N> matrixBase<T,N>::operator=(const matrixBase<T,N> &inMatrix)
{
  data = inMatrix.data;
  dims = inMatrix.dims;
  return *this;
}

template<typename T, uint N>
void matrixBase<T,N>::setup(uint inDim0, uint inDim1, uint inDim2, uint inDim3)
{
  dims = {{inDim0,inDim1,inDim2,inDim3}};
  data.resize(inDim0*inDim1*inDim2*inDim3);
}

template<typename T, uint N>
T* matrixBase<T,N>::operator[](int inRow)
{
  if (inRow < (int)this->dims[0] && inRow >= 0) {
    return &data[inRow*dims[1]];
  }
  else {
    FatalError("Attempted out-of-bounds access in matrix.");
  }
}

template<typename T, uint N>
T& matrixBase<T,N>::operator()(int i, int j, int k, int l)
{
  if (i<(int)this->dims[0] && i>=0 && j<(int)this->dims[1] && j>=0 &&
      k<(int)this->dims[2] && k>=0 && l<(int)this->dims[3] && l>= 0)
  {
    return data[l+dims[3]*(k+dims[2]*(j+dims[1]*i))];
  }
  else {
    cout << "i=" << i << ", dim0=" << dims[0] << ", j=" << j << ", dim1=" << dims[1] << ", ";
    cout << "k=" << k << ", dim2=" << dims[2] << ", l=" << l << ", dim3=" << dims[3] << endl;
    FatalError("Attempted out-of-bounds access in Array.");
  }
}

template<typename T>
T& matrix<T>::operator()(int i, int j)
{
  if (i<(int)this->dims[0] && i>=0 && j<(int)this->dims[1] && j>=0) {
    return this->data[j+this->dims[1]*i];
  }
  else {
    cout << "i=" << i << ", dim0=" << this->dims[0] << ", j=" << j << ", dim1=" << this->dims[1];
    FatalError("Attempted out-of-bounds access in matrix.");
  }
}


template<typename T>
void matrix<T>::initializeToZero(void)
{
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      for (uint k=0; k<this->dims[2]; k++)
        for (uint l=0; l<this->dims[3]; l++)
          this->data[l+this->dims[3]*(k+this->dims[2]*(j+this->dims[1]*i))] = 0;
}

template<typename T>
void matrix<T>::initializeToValue(T val)
{
  for (uint i=0; i<this->dims[0]; i++)
    for (uint j=0; j<this->dims[1]; j++)
      for (uint k=0; k<this->dims[2]; k++)
        for (uint l=0; l<this->dims[3]; l++)
          this->data[l+this->dims[3]*(k+this->dims[2]*(j+this->dims[1]*i))] = val;
}

//template<typename T>
//void matrix<T,1>::insertRow(const vector<T> &vec, int rowNum)
//{
//  data.insert(data.begin()+rowNum,vec.begin(),1);

//  if (dims[0] == 0) {
//    dims = {{1,1,1,1}};
//  }
//  else {
//    dims = {{dims[0]+1,1,1,1}};
//  }
//}
template<typename T, uint N>
void matrixBase<T,N>::insertRow(const vector<T> &vec, int rowNum)
{
  if (N!=2)
    FatalError("InsertRow only supported for 2D arrays.");

  if (this->dims[1]!= 0 && vec.size()!=this->dims[1])
    FatalError("Attempting to assign row of wrong size to matrix.");

  if (rowNum==INSERT_AT_END || rowNum==(int)this->dims[0]) {
    // Default action - add to end
    this->data.insert(this->data.end(),vec.begin(),vec.end());
  }else{
    // Insert at specified location
    this->data.insert(this->data.begin()+rowNum*this->dims[1],vec.begin(),vec.end());
  }

  if (this->dims[1]==0) {
    this->dims[1] = vec.size(); // This may not be needed (i.e. may never have dim1==0). need to verify how I set up dim0, dim1...
    this->dims[2] = 1;
    this->dims[3] = 1;
  }
  this->dims[0]++;
}

template<typename T>
void matrix<T>::insertRow(const vector<T> &vec, int rowNum)
{
  if (this->dims[1]!= 0 && vec.size()!=this->dims[1])
    FatalError("Attempting to assign row of wrong size to matrix.");

  if (rowNum==INSERT_AT_END || rowNum==(int)this->dims[0]) {
    // Default action - add to end
    this->data.insert(this->data.end(),vec.begin(),vec.end());
  }else{
    // Insert at specified location
    this->data.insert(this->data.begin()+rowNum*this->dims[1],vec.begin(),vec.end());
  }

  if (this->dims[1]==0) this->dims[1]=vec.size(); // This may not be needed (i.e. may never have dim1==0). need to verify how I set up dim0, dim1...
  this->dims[0]++;
}

template<typename T>
void matrix<T>::insertRow(T *vec, uint rowNum, uint length)
{
  if (this->dims[1]!=0 && length!=(int)this->dims[1])
    FatalError("Attempting to assign row of wrong size to matrix.");

  if (rowNum==INSERT_AT_END || rowNum==(int)this->dims[0]) {
    // Default action - add to end
    this->data.insert(this->data.end(),vec,vec+length);
  }else{
    // Insert at specified location
    this->data.insert(this->data.begin()+rowNum*this->dims[1],vec,vec+length);
  }

  if (this->dims[0]==0)
    this->dims = {{0,length,1,1}};

  this->dims[0]++;
}


template<typename T>
void matrix<T>::insertRowUnsized(const vector<T> &vec)
{
  // Add row to end, and resize matrix (add columns) if needed
  if (vec.size() > this->dims[1]) addCols(vec.size()-this->dims[1]);

  this->data.insert(this->data.end(),vec.begin(),vec.end());

  // If row too short, finish filling with 0's
  if (vec.size() < this->dims[1]) this->data.insert(this->data.end(),this->dims[1]-vec.size(),(T)0);

  this->dims[0]++;
}

template<typename T>
void matrix<T>::insertRowUnsized(T* vec, int length)
{
  // Add row to end, and resize matrix (add columns) if needed
  if (length > this->dims[1]) addCols(length-this->dims[1]);

  this->data.insert(this->data.end(),vec,vec+length);

  // If row too short, finish filling with 0's
  if (length < this->dims[1]) this->data.insert(this->data.end(),this->dims[1]-length,(T)0);

  this->dims[0]++;
}

template<typename T>
void matrix<T>::addCol(void)
{
  typename vector<T>::iterator it;
  for (uint row=0; row<this->dims[0]; row++) {
    it = this->data.begin() + (row+1)*(this->dims[1]+1) - 1;
    this->data.insert(it,1,(T)0);
  }
  this->dims[1]++;
}

template<typename T>
void matrix<T>::addCols(int nCols)
{
  typename vector<T>::iterator it;
  for (uint row=0; row<this->dims[0]; row++) {
    it = this->data.begin() + (row+1)*(this->dims[1]+nCols) - nCols;
    this->data.insert(it,nCols,(T)0);
  }
  this->dims[1] += nCols;
}

//template<typename T>
//void matrix<T,1>::addCol(void)
//{
//  FatalError("addCol not supported for 1D matrices.");
//}


//template<typename T>
//void matrix<T,1>::addCols(int)
//{
//  FatalError("addCols not supported for 1D matrices.");
//}


template<typename T>
void matrix<T>::removeCols(int nCols)
{
  if (nCols == 0) return;

  typename vector<T>::iterator it;
  for (uint row=this->dims[0]; row>0; row--) {
    it = this->data.begin() + row*this->dims[1];
    this->data.erase(it-nCols,it);
  }
  this->dims[1] -= nCols;
}

template<typename T>
vector<T> matrix<T>::getRow(uint row)
{
  vector<T> out;
  out.assign(&(this->data[row*this->dims[1]]),&(this->data[row*this->dims[1]])+this->dims[1]);
  return out;
}

template<typename T>
matrix<T> matrix<T>::getRows(vector<int> ind)
{
  matrix<T> out;
  for (auto& i:ind) out.insertRow(&(this->data[i*this->dims[1]]),-1,this->dims[1]);
  return out;
}

template<typename T>
vector<T> matrix<T>::getCol(int col)
{
  vector<T> out;
  for (uint i=0; i<this->dims[0]; i++) out.push_back(this->data[i*this->dims[1]+col]);
  return  out;
}

template<typename T, uint N>
void matrixBase<T,N>::print()
{
//  for (uint i=0; i<dims[0]; i++) {
//    for (uint j=0; j<dims[1]; j++) {
//      cout << left << setw(10) << setprecision(8) << data[i*dims[1]+j] << " ";
//    }
//    cout << endl;
//  }
}

template<typename T>
void matrix<T>::unique(matrix<T>& out, vector<int> &iRow)
{
  out.setup(0,0);

  // Setup vector for rows of final matrix that map to each row of initial matrix
  iRow.resize(this->dims[0]);
  iRow.assign(this->dims[0],-1);

  /* --- For each row in the matrix, compare to all
     previous rows to get first unique occurence --- */
  typename vector<T>::iterator itI, itJ;
  for (uint i=0; i<this->dims[0]; i++) {
    itI = this->data.begin() + i*this->dims[1];
    for (uint j=0; j<i; j++) {
      itJ = this->data.begin() + j*this->dims[1];
      if (equal(itI,itI+this->dims[1],itJ)) {
        iRow[i] = iRow[j];
        break;
      }
    }

    // If no prior occurance found, put in 'out' matrix
    if (iRow[i]==-1) {
      out.insertRow(&(this->data[i*this->dims[1]]),-1,this->dims[1]);
      iRow[i] = out.getDim0() - 1;
    }
  }
}

template<typename T, uint N>
T *matrixBase<T,N>::getData(void)
{
  return data.data();
}


// Fix for compiler to know which template types will be needed later (and therefore must now be compiled):
template class matrixBase<int,1>;
template class matrixBase<int,2>;
template class matrixBase<int,3>;
template class matrixBase<int,4>;

template class matrixBase<double,1>;
template class matrixBase<double,2>;
template class matrixBase<double,3>;
template class matrixBase<double,4>;

template class matrixBase<double*,1>;
template class matrixBase<double*,2>;
template class matrixBase<double*,3>;
template class matrixBase<double*,4>;

template class matrixBase<point,1>;
template class matrixBase<point,2>;
template class matrixBase<point,3>;
template class matrixBase<point,4>;

template class matrixBase<set<int>,1>;
template class matrixBase<set<int>,2>;
template class matrixBase<set<int>,3>;
template class matrixBase<set<int>,4>;

template class matrix<int>;
template class matrix<double>;
