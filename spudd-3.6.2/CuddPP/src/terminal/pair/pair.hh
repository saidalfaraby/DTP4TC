#ifndef _PAIR
#define _PAIR

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#define absv(x) (((x)<0) ? -(x) : (x))
#define PI 3.1415926535897
class Pair{

public:
  //fields: minimum and maximum
  double min, max;

  //Constructors for my class
  Pair();
  Pair(double,double);
  Pair(double);
  Pair(const Pair&);
  //destructor for my class
  ~Pair();

  Pair operator+(Pair&);
  Pair operator-(Pair&);
  Pair operator/(Pair&);
  Pair operator*(Pair&);
  Pair operator*(double);
  bool operator>=(Pair&);
  bool operator<=(Pair&);
  bool operator>(Pair&);
  bool operator<(Pair&);
  bool operator!=(Pair&); 
  bool operator<(int);
  Pair operator-();
  void operator=(Pair);
  Pair operator/(double);
  bool operator==(Pair&);

  Pair absVal();

  void span(Pair&); 
  void set_max(double max);
  void set_min(double min);
  void set_both(double min, double max);
  void set(double minmax);
  const double get_max() {return max;}
  const double get_min() {return min;}
  const double get_val() {return min + max;}
  const double hashCode();
  char *toString();
  int readString(FILE *, char *buf);
  void parseString(char *buf);
};
//ostream& operator <<(ostream & outs, const Pair& source);

Pair ceill(const Pair& source);

typedef Pair* PtrTerminal;
typedef Pair Terminal;

#endif

