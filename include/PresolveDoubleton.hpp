// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveDoubleton_H
#define PresolveDoubleton_H

#define	DOUBLETON	5

class doubleton_action : public PresolveAction {
 public:
  struct action {
    int icolx;
    int row;

    double clox;
    double cupx;
    double costx;
    
    int icoly;
    double cloy;
    double cupy;
    double costy;

    double rlo;
    double rup;

    double coeffx;
    double coeffy;

    int ncolx;
    double *colx;
    int *indx;

    int ncoly;
    double *coly;
    int *indy;
  };

  const int nactions_;
  const action *const actions_;

 private:
  doubleton_action(int nactions,
		      const action *actions,
		      const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions), actions_(actions)
{}

 public:
  const char *name() const { return ("doubleton_action"); }

  static const PresolveAction *presolve(PresolveMatrix *,
					 const PresolveAction *next);
  
  void postsolve(PostsolveMatrix *prob) const;

  ~doubleton_action();
};
#endif


