// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef ClpEventHandler_H
#define ClpEventHandler_H

#include "ClpSimplex.hpp"
/** Base class for Clp event handling
    
This is just here to allow for event handling.  By event I mean a Clp event
e.g. end of values pass.

One use would be to let a user handle a system event e.g. Control-C.  This could be done
by deriving a class MyEventHandler which knows about such events.  If one occurs 
MyEventHandler::event() could clear event status and return 3 (stopped).

Clp would then return to user code.

As it is called every iteration this should be fine grained enough.

User can derive and construct from CbcModel  - not pretty
    
*/

class ClpEventHandler  {
  
public:
  /** enums for what sort of event.

      These will also be returned in ClpModel::secondaryStatus() as int
  */
  enum Event {
    endOfIteration = 100, // used to set secondary status
    endOfFactorization,
    endOfValuesPass,
    node, // for Cbc
    treeStatus, // for Cbc
    solution, // for Cbc
    theta // hit in parametrics
  };
  /**@name Virtual method that the derived classe should provide.
   The base class instance does nothing and as event() is only useful method
   it would not be very useful NOT providing one!
  */
  //@{
  /** This can do whatever it likes.  If return code -1 then carries on
      if 0 sets ClpModel::status() to 5 (stopped by event) and will return to user.
      At present if <-1 carries on and if >0 acts as if 0 - this may change
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
  { return model_;}
  //@}
  
  
protected:
  /**@name Data members
     The data members are protected to allow access for derived classes. */
  //@{
  /// Pointer to simplex
  ClpSimplex * model_;
  //@}
};
/** Base class for Clp disaster handling
    
This is here to allow for disaster handling.  By disaster I mean that Clp
would otherwise give up

*/

class ClpDisasterHandler  {
  
public:
  /**@name Virtual methods that the derived classe should provide.
  */
  //@{
  /// Into simplex
  virtual void intoSimplex()=0;
  /// Checks if disaster
  virtual bool check() const = 0;
  /// saves information for next attempt
  virtual void saveInfo() =0;
  //@}
  
  
  /**@name Constructors, destructor */

  //@{
  /** Default constructor. */
  ClpDisasterHandler(ClpSimplex * model = NULL);
  /** Destructor */
  virtual ~ClpDisasterHandler();
  // Copy
  ClpDisasterHandler(const ClpDisasterHandler&);
  // Assignment
  ClpDisasterHandler& operator=(const ClpDisasterHandler&);
  /// Clone
  virtual ClpDisasterHandler * clone() const =0;

  //@}
  
  /**@name Sets/gets */

  //@{
  /** set model. */
  void setSimplex(ClpSimplex * model);
  /// Get model
  inline ClpSimplex * simplex() const
  { return model_;}
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
