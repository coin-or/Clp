// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolvePsdebug_H
#define PresolvePsdebug_H
void presolve_hincol_ok(const int *mcstrt, const int *hincol,
	       const int *hinrow,
	       const int *hrow, int ncols);

void presolve_links_ok(presolvehlink *link, int *starts, int *lengths, int n);


void presolve_no_zeros(const int *mcstrt, const double *colels, const int *hincol, int ncols);
#endif
