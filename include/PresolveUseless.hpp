// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveUseless_H
#define PresolveUseless_H
#define	USELESS		20

class useless_constraint_action : public PresolveAction {
  struct action {
    double rlo;
    double rup;
    int row;
    int ninrow;
    const int *rowcols;
    const double *rowels;
  };

  const int nactions_;
  const action *const actions_;

  useless_constraint_action(int nactions,
		      const action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions) {};

 public:
  const char *name() const;

  // These rows are asserted to be useless,
  // that is, given a solution the row activity
  // must be in range.
  static const PresolveAction *presolve(PresolveMatrix * prob,
					 const int *useless_rows,
					 int nuseless_rows,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;
};


#endif
