// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveSubst_H
#define PresolveSubst_H
#define	SUBST_ROW	21


class subst_constraint_action : public PresolveAction {
public:
  struct action {
    int col;
    int rowx;
    int rowy;

    int nincol;
    int *rows;
    double *rlos;
    double *rups;

    double *coeffxs;
    
    int *ninrowxs;
    /*const*/ int *rowcolsxs;
    /*const*/ double *rowelsxs;

    const double *costsx;
  };

  const int nactions_;
  const action *const actions_;

  subst_constraint_action(int nactions,
			  action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions) {};

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix * prob,
					 int *implied_free,
					 const PresolveAction *next,
					int & fill_level);
  static const PresolveAction *presolveX(PresolveMatrix * prob,
				  const PresolveAction *next,
				  int fillLevel);

  void postsolve(PostsolveMatrix *prob) const;

  ~subst_constraint_action();
};





/*static*/ void implied_bounds(const double *els,
			   const double *clo, const double *cup,
			   const int *hcol,
			   CoinBigIndex krs, CoinBigIndex kre,
			   double *maxupp, double *maxdownp,
			   int jcol,
			   double rlo, double rup,
			   double *iclb, double *icub);
#endif
