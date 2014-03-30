#ifndef _VECTOR
#define _VECTOR

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define vabs(x) (((x)<0) ? -(x) : (x))

class Vector {

public:
  
  double *data;
  int n;
  // 0,1,2,3 = TERMINAL, GO_DOWN, GO_UP, STAY_SAME
  int type;
  

  //Constructors for vector class
  Vector();
  Vector(int);
  Vector(double);
  //  Vector(double,double);
  Vector(int, double []);
  Vector(int, double);
  Vector(int, int);
  Vector(const Vector&);
  //destructor for my class
  ~Vector();

  Vector operator+(Vector&);
  Vector operator-(Vector&);
  Vector operator/(Vector&);
  Vector operator*(Vector&);
  Vector operator*(double);
  bool operator>=(Vector&);
  bool operator<=(Vector&);
  bool operator>(Vector&);
  bool operator<(Vector&);
  bool operator!=(Vector&); 
  bool operator<(int);
  Vector operator-();
  void operator=(Vector);
  Vector operator/(double);
  bool operator==(Vector&);
  
  Vector absVal();
  void span(Vector&); 
  void set_max(double max) {set(max);}
  void set_min(double min) { set(min); }
  void set_both(double min, double max) { set((min+max)/2.0); }
  void set(double val);
  void setAssign(double val[]);
  void setType(int typ);

  Vector maximize(Vector&);
  Vector minimize(Vector&);

  Vector threshold(Vector&);

  // for now, I'm just making this up. Seems like it should
  // return a double[n] array? 

  const double get_max() {return data[n];}
  const double get_min() {return data[0];}

  const double get_val() {  return data[0];}
  
  const double hashCode();
  
  /* We need to stringafy all the value so that we can print them all at once....*/
  
  char* toString();
  double value;


};
//ostream& operator <<(ostream & outs, const Vector& source);

// used for rounding off
Vector ceill(const Vector& source);

typedef Vector* PtrTerminal;
typedef Vector Terminal;

#endif

