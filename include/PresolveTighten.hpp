// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveTighten_H
#define PresolveTighten_H
// This action has no separate class;
// instead, it decides which columns can be made fixed
// and calls make_fixed_action::presolve.
const PresolveAction *tighten_zero_cost(PresolveMatrix *prob,
					 const PresolveAction *next);

#define	DO_TIGHTEN	30

class do_tighten_action : public PresolveAction {
  struct action {
    int col;
    int nrows;
    int direction;	// just for assertions
    int *rows;
    double *lbound;
    double *ubound;
  };

  const int nactions_;
  const action *const actions_;

  do_tighten_action(int nactions,
		      const action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions) {}

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix *prob,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;
};
#endif


