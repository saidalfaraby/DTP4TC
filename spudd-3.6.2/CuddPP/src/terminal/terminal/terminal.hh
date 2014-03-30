#ifndef _TERMINAL
#define _TERMINAL

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define absv(x) (((x)<0) ? -(x) : (x))

class Terminal{

public:
  //field: a single value
  double value;


  //Constructors for my class
  Terminal();
  Terminal(double);
  Terminal(const Terminal &);
  //destructor for my class
  ~Terminal();

  Terminal operator+(Terminal&);
  Terminal operator-(Terminal&);
  Terminal operator/(Terminal&);
  Terminal operator*(Terminal&);
  Terminal operator*(double);
  bool operator>=(Terminal&);
  bool operator<=(Terminal&);
  bool operator>(Terminal&);
  bool operator<(Terminal&);
  bool operator!=(Terminal&); 
  bool operator<(int);
  Terminal operator-();
  void operator=(Terminal);
  Terminal operator/(double);
  bool operator==(Terminal&);

  Terminal absVal();

  void set(double value);
  
  /* We need a get_val() method for each class that we create.....using the first element of
     any structure seems like it would work....this is used mostly in stuff like Hamming distance
     which we don't really use */
  const double get_val() {return value;}
  
  /*We don't need these guys anymore */
  //const double get_max() {return value;}
  const double get_min() {return value;}
  
  /* We need a hashCode() method for each class.  It has to return a double that will be unique
     for every possible value of the data structure that we are using.
     For terminal, returning the value is fine since we are using doubles.  For the pair class, 
     returning min + PI*max seems to guarrantee that we will have a unique hash code....for the 
     vector class, we have to use all the values and make some combination that will uniquely 
     determines each value */
  
  const double hashCode() {return value;}
  
  /* We need to stringafy all the value so that we can print them all at once....*/
  
  char* toString();

  // we need these two methods to allow reading and writing dds to file
  // reads from file pointer fp into character array buf according to
  // what the toString() method would have writte
  int readString(FILE *fp, char *buf);

  // parses the char array buf from readString to 'decode' all the internal state for the object
  void parseString(char *buf);
  

};
//ostream& operator <<(ostream & outs, const Terminal& source);

Terminal ceill(const Terminal& source);


/* We need this for every class */
typedef Terminal* PtrTerminal;
typedef Terminal Terminal;

#endif

