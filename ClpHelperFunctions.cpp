// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.

/*
    Note (JJF) I have added some operations on arrays even though they may
    duplicate CoinDenseVector.  I think the use of templates was a mistake
    as I don't think inline generic code can take as much advantage of
    parallelism or machine architectures or memory hierarchies.

*/
#include <cfloat>
#include <cmath>
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"
double 
maximumAbsElement(const double * region, int size)
{
  int i;
  double maxValue=0.0;
  for (i=0;i<size;i++) 
    maxValue = std::max(maxValue,fabs(region[i]));
  return maxValue;
}
void 
setElements(double * region, int size, double value)
{
  int i;
  for (i=0;i<size;i++) 
    region[i]=value;
}
void 
multiplyAdd(const double * region1, int size, double multiplier1,
		 double * region2, double multiplier2)
{
  int i;
  if (multiplier1==1.0) {
    if (multiplier2==1.0) {
      for (i=0;i<size;i++) 
	region2[i] = region1[i] + region2[i];
    } else if (multiplier2==-1.0) {
      for (i=0;i<size;i++) 
	region2[i] = region1[i] - region2[i];
    } else if (multiplier2==0.0) {
      for (i=0;i<size;i++) 
	region2[i] = region1[i] ;
    } else {
      for (i=0;i<size;i++) 
	region2[i] = region1[i] + multiplier2*region2[i];
    }
  } else if (multiplier1==-1.0) {
    if (multiplier2==1.0) {
      for (i=0;i<size;i++) 
	region2[i] = -region1[i] + region2[i];
    } else if (multiplier2==-1.0) {
      for (i=0;i<size;i++) 
	region2[i] = -region1[i] - region2[i];
    } else if (multiplier2==0.0) {
      for (i=0;i<size;i++) 
	region2[i] = -region1[i] ;
    } else {
      for (i=0;i<size;i++) 
	region2[i] = -region1[i] + multiplier2*region2[i];
    }
  } else if (multiplier1==0.0) {
    if (multiplier2==1.0) {
      // nothing to do
    } else if (multiplier2==-1.0) {
      for (i=0;i<size;i++) 
	region2[i] =  -region2[i];
    } else if (multiplier2==0.0) {
      for (i=0;i<size;i++) 
	region2[i] =  0.0;
    } else {
      for (i=0;i<size;i++) 
	region2[i] =  multiplier2*region2[i];
    }
  } else {
    if (multiplier2==1.0) {
      for (i=0;i<size;i++) 
	region2[i] = multiplier1*region1[i] + region2[i];
    } else if (multiplier2==-1.0) {
      for (i=0;i<size;i++) 
	region2[i] = multiplier1*region1[i] - region2[i];
    } else if (multiplier2==0.0) {
      for (i=0;i<size;i++) 
	region2[i] = multiplier1*region1[i] ;
    } else {
      for (i=0;i<size;i++) 
	region2[i] = multiplier1*region1[i] + multiplier2*region2[i];
    }
  }
}
double 
innerProduct(const double * region1, int size, const double * region2)
{
  int i;
  double value=0.0;
  for (i=0;i<size;i++)
    value += region1[i]*region2[i];
  return value;
}
void 
getNorms(const double * region, int size, double & norm1, double & norm2)
{
  norm1 = 0.0;
  norm2 = 0.0;
  int i;
  for (i=0;i<size;i++) {
    norm2 += region[i]*region[i];
    norm1 = std::max(norm1,fabs(region[i]));
  }
}

