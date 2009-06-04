/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
int main (int argc, const char *argv[])
{
  ClpSimplex  model;
  int status;
  if (argc<2)
    status=model.readMps("../../Data/Sample/p0033.mps");
  else
    status=model.readMps(argv[1]);
  if (!status) {
    model.primal();
  }
  return 0;
}    
