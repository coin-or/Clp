// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveForcing_H
#define PresolveForcing_H
#define	IMPLIED_BOUND	7

class forcing_constraint_action : public PresolveAction {
  struct action {
    int row;
    int nlo;
    int nup;
    const int *rowcols;
    const double *bounds;
  };

  const int nactions_;
  const action *const actions_;

  forcing_constraint_action(int nactions,
		      const action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions) {};

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix * prob,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  ~forcing_constraint_action() { deleteAction(actions_,action*); }
};




// HACK - from doubleton.cpp
void compact_rep(double *elems, int *indices, CoinBigIndex *starts, const int *lengths, int n,
		 const presolvehlink *link);

#endif
