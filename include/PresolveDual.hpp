
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#ifndef PresolveDual_H
#define PresolveDual_H
class remove_dual_action : public PresolveAction {
 public:
  remove_dual_action(int nactions,
		     //const action *actions,
		      const PresolveAction *next);
  static const PresolveAction *presolve(PresolveMatrix *prob,
					 const PresolveAction *next);
};
#endif


