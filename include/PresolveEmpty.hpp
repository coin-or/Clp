// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveEmpty_H
#define PresolveEmpty_H
// Drop all empty rows/cols from the problem.
// 
// This should only be done once, after all other presolving actions have
// been done.  

const int DROP_ROW = 3;
const int DROP_COL = 4;

class drop_empty_cols_action : public PresolveAction {
private:
  const int nactions_;

  struct action {
    int jcol;
    double clo;
    double cup;
    double cost;
    double sol;
  };
  const action *const actions_;

  drop_empty_cols_action(int nactions,
			 const action *const actions,
			 const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), 
    actions_(actions)
  {}

 public:
  const char *name() const { return ("drop_empty_cols_action"); }

  static const PresolveAction *presolve(PresolveMatrix *,
					 int *ecols,
					 int necols,
					 const PresolveAction*);

  static const PresolveAction *presolve(PresolveMatrix *prob,
					 const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  ~drop_empty_cols_action() { deleteAction(actions_); }
};



class drop_empty_rows_action : public PresolveAction {
private:
  struct action {
    double rlo;
    double rup;
    int row;
    int fill_row;	// which row was moved into position row to fill it
  };

  const int nactions_;
  const action *const actions_;

  drop_empty_rows_action(int nactions,
			 const action *actions,
			 const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions)
{}

 public:
  const char *name() const { return ("drop_empty_rows_action"); }

  static const PresolveAction *presolve(PresolveMatrix *prob,
					    const PresolveAction *next);

  void postsolve(PostsolveMatrix *prob) const;

  ~drop_empty_rows_action() { deleteAction(actions_); }
};
#endif

