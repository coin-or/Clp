// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveIsolated_H
#define PresolveIsolated_H
class isolated_constraint_action : public PresolveAction {
  double rlo_;
  double rup_;
  int row_;
  int ninrow_;
  const int *rowcols_;
  const double *rowels_;
  const double *costs_;

  isolated_constraint_action(double rlo,
			     double rup,
			     int row,
			     int ninrow,
			     const int *rowcols,
			     const double *rowels,
			     const double *costs,
			     const PresolveAction *next) :
    PresolveAction(next),
    rlo_(rlo), rup_(rup), row_(row), ninrow_(ninrow),
    rowcols_(rowcols), rowels_(rowels), costs_(costs) {}

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix * prob,
					 int row,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;
};



#endif
