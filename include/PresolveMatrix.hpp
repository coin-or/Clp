// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef PresolveMatrix_H
#define PresolveMatrix_H

#include "ClpSimplex.hpp"

#include <cmath>
#include <cfloat>
using std::min;
using std::max;


// OSL had a fixed zero tolerance; we still use that here.
const double ZTOLDP      = 1e-12;

double *presolve_duparray(const double *d, int n, int n2);
double *presolve_duparray(const double *d, int n);
int *presolve_duparray(const int *d, int n, int n2);
int *presolve_duparray(const int *d, int n);
char *presolve_duparray(const char *d, int n, int n2);
char *presolve_duparray(const char *d, int n);

void presolve_delete_from_row(int row, int col /* thing to delete */,
		     const int *mrstrt,
		     int *hinrow, int *hcol, double *dels);

int presolve_find_row(int row, int kcs, int kce, const int *hrow);
int presolve_find_row1(int row, int kcs, int kce, const int *hrow);

//#define	DEBUG_PRESOLVE	1

#ifdef	DEBUG_PRESOLVE
inline void	DIE(char *s)	{ std::cout<<s; abort(); }
#else
  inline void	DIE(char *s)	{}
#endif

#if	DEBUG_PRESOLVE
#define	PRESOLVE_STMT(s)	s
#define PRESOLVEASSERT(x)	((x) ? 1 : ((std::cerr<< "FAILED ASSERTION at line "<< __LINE__ << ":  " #x "\n"), abort(), 0))
#else
#define PRESOLVEASSERT(x)
#define	PRESOLVE_STMT(s)
#endif

struct dropped_zero {
  int row;
  int col;
};

inline int ALIGN(int n, int m)	{ return (((n + m - 1) / m) * m); }
inline int ALIGN_DOUBLE(int n)	{ return ALIGN(n,sizeof(double)); }

class PostsolveMatrix;

// Note 77
// "Members and bases are constructed in ordre of declation
//  in the class and destroyed in the reverse order."  C++PL 3d Ed. p. 307
//
// That's why I put integer members (such as ncols) before the array members;
// I like to use those integer values during initialization.
// NOT ANYMORE

//
// This is the abstract base class of all presolve routines.
// The object just stores the information needed to postsolve;
// this information is not expected to be changed, so the fields are
// all const.
// It is expected that subclasses will declare static functions that
// attempt to perform the presolve; if it succeeds, the function creates a
// new presolve object and returns it, otherwise it returns 0.
// It is expected that these static functions will be the only things
// that can create new presolve_action objects;
// this is expressed by making the constructor(s) private.
// Every subclass must define a postsolve routine.  This gets two records,
// one that contains information that is also used in presolving (prob)
// and one with information that is only used in postsolving (prob2).
// (suggestions for better names welcome).
// Finally, there is a single postsolve driver (the friend function)
// that just calls the postsolve member function of all the postsolve objects
// in order.

// note that since the only field in a presolve_action is const,
// anything one can do with a variable declared:
// 	PresolveAction*
// can also be done with a variable declared:
// 	const PresolveAction*
//
// It is expected that all derived subclasses of PresolveAction also
// have this property.
class PresolveAction {
 public:
  // Exceptions are inefficient, particularly with g++.
  // Even with xlC, the use of exceptions adds a long prologue to a routine.
  // Therefore, rather than use throw directly in the routine,
  // I use it in a stub routine.
  static void throwCoinError(const char *error, const char *ps_routine);

  // The only thing the base class does is point to the next
  // presolve transform in the list.
  const PresolveAction *next;
  
  PresolveAction(const PresolveAction *next) : next(next) {}

  // A name for debug printing.
  // It is expected that the name is not stored in the transform itself.
  virtual const char *name() const = 0;

  // postsolve this particular presolve action
  virtual void postsolve(PostsolveMatrix *prob) const = 0;

  virtual ~PresolveAction() {}
};

// This collects all the information about the problem that is needed
// both in presolve and postsolve.
class PrePostsolveMatrix {
 public:
  /// enums for status of various sorts (matches CoinWarmStartBasis)
  enum Status {
    isFree = 0x00,
    basic = 0x01,
    atUpperBound = 0x02,
    atLowerBound = 0x03,
    superBasic = 0x04
  };
  double *sol_;
  double *rowduals_;
  double *acts_;

  double *rcosts_;
  unsigned char *colstat_;
  unsigned char *rowstat_;


  // Original objective offset
  double originalOffset_;
  // Message handler
   CoinMessageHandler *  handler_; 
   /// Messages
   CoinMessages messages_; 

   inline CoinMessageHandler * messageHandler() const 
  { return handler_; }
   /// Return messages
   inline CoinMessages messages() const 
  { return messages_; }
  // colrep
  int ncols_;
  const int ncols0_;

  int nelems_;

  int *mcstrt_;
  int *hincol_;
  int *hrow_;
  double *colels_;

  double *cost_;

  double *clo_;
  double *cup_;
  double *rlo_;
  double *rup_;

  // Original column numbers
  int * originalColumn_;

  const double ztolzb_;
  const double ztoldj_;

  double maxmin_;

  int whichpass_;	// mostly for debugging

  PrePostsolveMatrix(const ClpSimplex& si,
			int ncols_,
			int nrows_,
			int nelems_);

  ~PrePostsolveMatrix();
  // Status stuff
  
  inline void setRowStatus(int sequence, Status status)
  {
    unsigned char & st_byte = rowstat_[sequence];
    st_byte &= ~7;
    st_byte |= status;
  };
  inline Status getRowStatus(int sequence) const
  {return static_cast<Status> (rowstat_[sequence]&7);};
  inline bool rowIsBasic(int sequence) const
  {return (static_cast<Status> (rowstat_[sequence]&7)==basic);};
  inline void setColumnStatus(int sequence, Status status)
  {
    unsigned char & st_byte = colstat_[sequence];
    st_byte &= ~7;
    st_byte |= status;
  };
  inline Status getColumnStatus(int sequence) const
  {return static_cast<Status> (colstat_[sequence]&7);};
  inline bool columnIsBasic(int sequence) const
  {return (static_cast<Status> (colstat_[sequence]&7)==basic);};
  /// Sets status (non -basic ) using value
  void setRowStatusUsingValue(int iRow);
  void setColumnStatusUsingValue(int iColumn);

};





/*
 * Currently, the matrix is represented the same way an CoinPackedMatrix is.
 * Occasionally columns increase in size.
 * In order to check whether there is enough room for the column
 * where it sits, I wanted to know what the next column (in memory order)
 * in the representation was.
 * To do that, I maintain a linked list of columns; the "pre" and "suc"
 * fields give the previous and next columns, in memory order (that is,
 * the column whose mcstrt entry is next smaller or larger).
 * The same thing for the row representation.
 *
 * This is all likely to change, but I'm leaving it as it is for now.
 */
//  static const int	NO_LINK	= -66666666;
#define NO_LINK -66666666
#define	PRESOLVE_INF DBL_MAX

class presolvehlink {
public:
  int pre, suc;
};
  
static inline void PRESOLVE_REMOVE_LINK(presolvehlink *link, int i)
{ 
  int ipre = link[i].pre;
  int isuc = link[i].suc;
  if (ipre >= 0) {
    link[ipre].suc = isuc;
  }
  if (isuc >= 0) {
    link[isuc].pre = ipre;
  }
  link[i].pre = NO_LINK, link[i].suc = NO_LINK;
}

// inserts i after pos
static inline void PRESOLVE_INSERT_LINK(presolvehlink *link, int i, int pos)
{
  int isuc = link[pos].suc;
  link[pos].suc = i;
  link[i].pre = pos;
  if (isuc >= 0) {
    link[isuc].pre = i;
  }
  link[i].suc = isuc;
}

// rename i to j
// that is, position j should be unused, and i will take its place
// should be equivalent to:
//   int pre = link[i].pre;
//   PRESOLVE_REMOVE_LINK(link, i);
//   PRESOLVE_INSERT_LINK(link, j, pre);
// if pre is not NO_LINK (otherwise -- ??)
static inline void PRESOLVE_MOVE_LINK(presolvehlink *link, int i, int j)
{ 
  int ipre = link[i].pre;
  int isuc = link[i].suc;
  if (ipre >= 0) {
    link[ipre].suc = j;
  }
  if (isuc >= 0) {
    link[isuc].pre = j;
  }
  link[i].pre = NO_LINK, link[i].suc = NO_LINK;
}



// this really should never happen.
// it will if there isn't enough space to postsolve the matrix.
// see the note below.
static void check_free_list(int free_list)
{
  if (free_list < 0) {
    printf("RAN OUT OF LINKS!!\n");
    abort();
  }
}


  // This collects all the information about the problem that is only
  // needed during presolve.
class PresolveMatrix : public PrePostsolveMatrix {
 public:

  // Crude linked lists, modelled after the linked lists used in OSL factorization.
  presolvehlink *clink_;
  presolvehlink *rlink_;

  double dobias_;

  // rowrep
  int nrows_;	// note 77
  int *mrstrt_;
  int *hinrow_;
  double *rowels_;
  int *hcol_;

  char *integerType_;
  // bounds can be moved by this to stay feasible
  double feasibilityTolerance_;
  // Output status 0=feasible, 1 infeasible, 2 unbounded
  int status_;
  // Should use templates ?
  // Rows
  // Bits to say if row changed
  unsigned int * rowChanged_;
  // Input list
  int * rowsToDo_;
  int numberRowsToDo_;
  // Output list
  int * nextRowsToDo_;
  int numberNextRowsToDo_;

  inline bool rowChanged(int i) const {
    return (rowChanged_[i>>5]>>(i&31))!=0;
  }
  inline void setRowChanged(int i) {
    unsigned int & value = rowChanged_[i>>5];
    int bit = i&31;
    value |= (1<<bit);
  }
  inline void addRow(int i) {
    unsigned int & value = rowChanged_[i>>5];
    int bit = i&31;
    if ((value&(1<<bit))==0) {
      value |= (1<<bit);
      nextRowsToDo_[numberNextRowsToDo_++] = i;
    }
  }
  inline void unsetRowChanged(int i) {
    unsigned int & value = rowChanged_[i>>5];
    int bit = i&31;
    value &= ~(1<<bit);
  }

  // Columns
  // Bits to say if column changed
  unsigned int * colChanged_;
  // Input list
  int * colsToDo_;
  int numberColsToDo_;
  // Output list
  int * nextColsToDo_;
  int numberNextColsToDo_;

  inline bool colChanged(int i) const {
    return (colChanged_[i>>5]>>(i&31))!=0;
  }
  inline void setColChanged(int i) {
    unsigned int & value = colChanged_[i>>5];
    int bit = i&31;
    value |= (1<<bit);
  }
  inline void addCol(int i) {
    unsigned int & value = colChanged_[i>>5];
    int bit = i&31;
    if ((value&(1<<bit))==0) {
      value |= (1<<bit);
      nextColsToDo_[numberNextColsToDo_++] = i;
    }
  }
  inline void unsetColChanged(int i) {
    unsigned int & value = colChanged_[i>>5];
    int bit = i&31;
    value &= ~(1<<bit);
  }
  PresolveMatrix(int ncols0,
		    double maxmin,
		    // end prepost members

		    ClpSimplex &si,

		    // rowrep
		    int nrows,
		    int nelems,
		 bool doStatus);

  ~PresolveMatrix();

  void consistent(bool testvals = true);

  inline void change_bias(double change_amount);

  void update_model(ClpSimplex& si,
			    int nrows0,
			    int ncols0,
			    int nelems0);
};


  // This collects all the information about the problem that is needed
  // only in postsolve.
class PostsolveMatrix : public PrePostsolveMatrix {
 public:


  int free_list_;
  int maxlink_;
  int *link_;

  // debug
  char *cdone_;
  char *rdone_;
  int nrows_;

  // needed for presolve_empty
  int nrows0_;

  PostsolveMatrix(ClpSimplex& si,

		   int ncols0,
		   int nrows0,
		   int nelems0,
		     
		   double maxmin_,
		   // end prepost members

		   double *sol,
		   double *acts,

		   unsigned char *colstat,
		   unsigned char *rowstat);

  ~PostsolveMatrix();

  // debugging
  void check_nbasic();

  ////virtual void postsolve(const PresolveAction *paction);
};

void PresolveMatrix::change_bias(double change_amount)
{
  dobias_ += change_amount;
#if DEBUG_PRESOLVE
  assert(fabs(change_amount)<1.0e50);
#endif
  if (change_amount)
    PRESOLVE_STMT(printf("changing bias by %g to %g\n",
		    change_amount, dobias_));  
}

// useful functions
inline void swap(int &x, int &y) { int temp = x; x = y; y = temp; }
inline void swap(double &x, double &y) { double temp = x; x = y; y = temp; }

inline void swap(int *&x, int *&y) { int *temp = x; x = y; y = temp; }
inline void swap(double *&x, double *&y) { double *temp = x; x = y; y = temp; }
// This returns a non const array filled with input from scalar
// or actual array
template <class T> inline T*
copyOfArray( const T * array, const int size, T value)
{
  T * arrayNew = new T[size];
  if (array) {
    memcpy(arrayNew,array,size*sizeof(T));
  } else {
    int i;
    for (i=0;i<size;i++) 
      arrayNew[i] = value;
  }
  return arrayNew;
}

// This returns a non const array filled with actual array (or NULL)
template <class T> inline T*
copyOfArray( const T * array, const int size)
{
  if (array) {
    T * arrayNew = new T[size];
    memcpy(arrayNew,array,size*sizeof(T));
    return arrayNew;
  } else {
    return NULL;
  }
}
#define	PRESOLVEFINITE(n)	(-PRESOLVE_INF < (n) && (n) < PRESOLVE_INF)
int presolve_find_row2(int irow, int ks, int nc, const int *hrow, const int *link);
void presolve_make_memlists(int *starts, int *lengths,
		       presolvehlink *link,
		       int n);
int presolve_find_row3(int irow, int ks, int nc, const int *hrow, const int *link);
void presolve_delete_from_row(int row, int col /* thing to delete */,
		     const int *mrstrt,
			      int *hinrow, int *hcol, double *dels);
void presolve_delete_from_row2(int row, int col /* thing to delete */,
		      int *mrstrt,
			       int *hinrow, int *hcol, double *dels, int *link, int *free_listp);
#endif
