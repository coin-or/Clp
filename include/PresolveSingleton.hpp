// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveSingleton_H
#define PresolveSingleton_H
#define	SLACK_DOUBLETON	2

const int MAX_SLACK_DOUBLETONS	= 1000;

class slack_doubleton_action : public PresolveAction {
  struct action {
    double clo;
    double cup;

    double rlo;
    double rup;

    double coeff;

    int col;
    int row;
  };

  const int nactions_;
  const action *const actions_;

  slack_doubleton_action(int nactions,
			 const action *actions,
			 const PresolveAction *next) :
    PresolveAction(next),
    nactions_(nactions),
    actions_(actions)
{}

 public:
  const char *name() const { return ("slack_doubleton_action"); }

  // notFinished is set if action array filled up
  static const PresolveAction *presolve(PresolveMatrix *,
					   const PresolveAction *next,
					bool &notFinished);

  void postsolve(PostsolveMatrix *prob) const;


  ~slack_doubleton_action() { deleteAction(actions_,action*); }
};
#endif
