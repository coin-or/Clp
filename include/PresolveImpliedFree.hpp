// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveImpliedFree_H
#define PresolveInpliedFree_H
#define	IMPLIED_FREE	9

class implied_free_action : public PresolveAction {
  struct action {
    int row, col;
    double clo, cup;
    double rlo, rup;
    int ninrow;
    const double *rowels;
    const int *rowcols;
    const double *costs;
  };

  const int nactions_;
  const action *const actions_;

  implied_free_action(int nactions,
		      const action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions) {}

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix * prob,
					 const PresolveAction *next,
					int & fillLevel);

  void postsolve(PostsolveMatrix *prob) const;

  ~implied_free_action() { delete[]actions_; }
};

#endif
