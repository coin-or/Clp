// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveDupcol_H
#define PresolveDupcol_H
#define	DUPCOL	10

class dupcol_action : public PresolveAction {
  struct action {
    double thislo;
    double thisup;
    double lastlo;
    double lastup;
    int ithis;
    int ilast;

    int nincol;
    int *colrows;
    double *colels;
  };

  const int nactions_;
  const action *const actions_;

  dupcol_action():PresolveAction(NULL),nactions_(0),actions_(NULL) {};
  dupcol_action(int nactions,
		const action *actions,
		const PresolveAction *next)/* :
    nactions_(nactions), actions_(actions),
    PresolveAction(next) {}*/;

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix *prob,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  //~dupcol_action() { delete[]actions_; }
};


class duprow_action : public PresolveAction {
  struct action {
    int row;
    double lbound;
    double ubound;
  };

  const int nactions_;
  const action *const actions_;

  duprow_action():PresolveAction(NULL),nactions_(0),actions_(NULL) {};
  duprow_action(int nactions,
		      const action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions) {}

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix *prob,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  //~duprow_action() { delete[]actions_; }
};

#endif

