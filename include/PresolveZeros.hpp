// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveZeros_H
#define PresolveZeros_H

#define	DROP_ZERO	8

class drop_zero_coefficients_action : public PresolveAction {

  const int nzeros_;
  const dropped_zero *const zeros_;

  drop_zero_coefficients_action(int nzeros,
				const dropped_zero *zeros,
				const PresolveAction *next) :
    PresolveAction(next),
    nzeros_(nzeros), zeros_(zeros)
{}

 public:
  const char *name() const { return ("drop_zero_coefficients_action"); }

  static const PresolveAction *presolve(PresolveMatrix *prob,
					 int *checkcols,
					 int ncheckcols,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  ~drop_zero_coefficients_action() { deleteAction(zeros_,dropped_zero*); }
};

const PresolveAction *drop_zero_coefficients(PresolveMatrix *prob,
					      const PresolveAction *next);

#endif
