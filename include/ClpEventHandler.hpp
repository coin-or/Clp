// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpEventHandler_H
#define ClpEventHandler_H

#include "ClpSimplex.hpp"

/** Base class for Clp event handling
    
This is just here to allow for event handling
    
*/

class ClpEventHandler  {
  
public:
  /** enums for what sort of event
  */
  enum Event {
    endOfIteration = 100, // used to set secondary status
    endOfFactorization,
    endOfValuesPass
  };
  /**@name Virtual methods that the derived classes must provide */
  //@{
  /** This can do whatever it likes.  If return code -1 then carries on
      if >=0 sets status and exits
  */
  virtual int event(Event whichEvent);
  //@}
  
  
  /**@name Constructors, destructor */

  //@{
  /** Default constructor. */
  ClpEventHandler(ClpSimplex * model = NULL);
  /** Destructor */
  virtual ~ClpEventHandler();
  // Copy
  ClpEventHandler(const ClpEventHandler&);
  // Assignment
  ClpEventHandler& operator=(const ClpEventHandler&);
  /// Clone
  virtual ClpEventHandler * clone() const;

  //@}
  
  /**@name Sets/gets */

  //@{
  /** set model. */
  void setSimplex(ClpSimplex * model);
  /// Get model
  inline ClpSimplex * simplex() const
  { return model_;};
  //@}
  
  
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer to simplex
  ClpSimplex * model_;
  //@}
};
#endif
