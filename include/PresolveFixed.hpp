// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveFixed_H
#define PresolveFixed_H
#define	FIXED_VARIABLE	1

class remove_fixed_action : public PresolveAction {
 public:
  struct action {
    int col;
    int nincol;

    double sol;
    int *colrows;
    double *colels;
  };

  int nactions_;
  const action *actions_;

 private:
  remove_fixed_action(int nactions,
		      const action *actions,
		      const PresolveAction *next);

 public:
  const char *name() const;

  static const remove_fixed_action *presolve(PresolveMatrix *prob,
					 int *fcols,
					 int nfcols,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  ~remove_fixed_action();
};


const PresolveAction *remove_fixed(PresolveMatrix *prob,
				    const PresolveAction *next);



class make_fixed_action : public PresolveAction {
  struct action {
    double bound;
  };

  int nactions_;
  const action *actions_;

  const bool fix_to_lower_;
  const remove_fixed_action *faction_;

  make_fixed_action(int nactions,
		    const action *actions,
		    bool fix_to_lower,
		    const remove_fixed_action *faction,
		    const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions),
    fix_to_lower_(fix_to_lower),
    faction_(faction)
{}

 public:
  const char *name() const;

  static const PresolveAction *presolve(PresolveMatrix *prob,
					 int *fcols,
					 int hfcols,
					 bool fix_to_lower,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  ~make_fixed_action() { delete[]actions_; delete faction_;};
};


const PresolveAction *make_fixed(PresolveMatrix *prob,
				    const PresolveAction *next);
#endif
