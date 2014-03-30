#include "vector.hh"

//Constructors
Vector::Vector() {
  n=7;
  type = 0;
  data = new double[n];
}
Vector::Vector(int nv) {
  n = nv;
  type = 0;
  data = new double[n];
  set(0.0);
}
Vector::Vector(double val) {
  n = 7;
  type = 0;
  data = new double[n];
  set(val);
}
Vector::Vector(int nv, double val){
  n = nv;
  type = 0;
  data = new double[n];
  set(val);
}

Vector::Vector(int nv, double vals[])
{
  n= nv;
  type = 0;
  data = new double[n];
  setAssign(vals);
}
// to create a vector of type GO_DOWN, GO_UP or STAY_SAME
// use this constructor
Vector::Vector(int nv, int typ) 
{
  n=nv;
  type = typ;
  data = new double[n];
  set(0.0);
}

Vector::Vector(const Vector& gvector)
{
  n = gvector.n;
  type = gvector.type;
  data = new double[n];
  setAssign(gvector.data);
}

//destructor
Vector::~Vector(){
  delete data;
}

void Vector::set(double val){
  for (int i=0; i<n; i++)
    data[i] = val;
}
void Vector::setAssign(double val[]){
  for (int i=0; i<n; i++)
    data[i] = val[i];
}
void Vector::setType(int typ) {
  type = typ;
}
// I think this is null now...
void Vector::span(Vector &gvector) {
}

Vector Vector::operator+(Vector& gvector)
{
  // this check should be done for all functions in this class?
  // only operator* can deal with types other than 0
  if (!type && !gvector.type) {
    Vector newVector(n);
    
    for (int i=0; i<n; i++) 
      newVector.data[i] = data[i]+gvector.data[i];
    
    return newVector;
  } else {
    // error!
    exit(0);
  }
}

void Vector::operator=(Vector gvector){
  type = gvector.type;
  n = gvector.n;
  data = new double[n];
  setAssign(gvector.data);
}

bool Vector::operator==(Vector& gvector){
  
  int i=0;
  if (type == gvector.type) 
    while (i < n && data[i] == gvector.data[i]) 
      i++;
  return (i==n);
}

Vector Vector::operator-(Vector& gvector)
{  
  Vector newVector(n);
 
  for (int i=0; i<n; i++) 
    newVector.data[i] = data[i]-gvector.data[i];

  return newVector;
}

Vector Vector::operator*(Vector& gvector)
{
  Vector newVector(n);
  
  // type = {0,1,2,3} = {TERMINAL, GO_DOWN, GO_UP, STAY_SAME};
  // newVector is always of type 0

  switch (type) {
  case 0: 
    switch (gvector.type) {
    case 0:
      for (int i=0; i<n; i++) 
	newVector.data[i] = data[i]*gvector.data[i];
      break;
    case 1:
      newVector.data[0] = data[0];
      for (int i=1; i<n; i++) 
	newVector.data[i] = data[i-1];
      break;
    case 2:
      for (int i=0; i<n-1; i++) 
	newVector.data[i] = data[i+1];
      newVector.data[n-1] = data[n-1];
      break;
    case 3:
      for (int i=0; i<n; i++) 
	newVector.data[i] = data[i];      
      break;
    }
    break;
  case 1:
    if (!gvector.type) {
      newVector.data[0] = gvector.data[0];
      for (int i=1; i<n; i++) 
	newVector.data[i] = gvector.data[i-1];
    } else {
      exit(0);
    }
    break;
  case 2:
    if (!gvector.type) {
      for (int i=0; i<n-1; i++) 
	newVector.data[i] = gvector.data[i+1];
      newVector.data[n-1] = gvector.data[n-1];
    } else {
      exit(0);
    }
    break;
  case 3:
    if (!gvector.type) {
      for (int i=0; i<n; i++) 
	newVector.data[i] = gvector.data[i];      
    } else {
      exit(0);
    }
    break;
  }
  return newVector;
}


Vector Vector::operator *(double number)
{
  
  Vector newVector(n);

  for (int i=0; i<n; i++) 
    newVector.data[i] = data[i]*number;

  return newVector;
}

Vector Vector::operator /(Vector& gvector){
  
  Vector newVector(n);

  for (int i=0; i<n; i++) 
    newVector.data[i] = data[i]/gvector.data[i];

  return newVector;
}

bool Vector::operator!=(Vector& gvector)
{    
  int i=0;
  if (type == gvector.type)
    while (i < n && data[i] == gvector.data[i]) 
      i++;
  return !(i==n);
}

bool Vector::operator<(Vector& gvector)
{  
  int i=0;
  while (i < n && data[i] < gvector.data[i]) 
    i++;
  return (i==n);
  
}

bool  Vector::operator<=(Vector& gvector)
{
  
  int i=0;
  while (i < n && data[i] <= gvector.data[i]) 
    i++;
  return (i==n);
}

bool Vector::operator>(Vector& gvector){
  int i=0;
  while (i < n && data[i] > gvector.data[i]) 
    i++;
  return (i==n);  
  
}

bool Vector::operator>=(Vector& gvector){
  int i=0;
  while (i < n && data[i] >= gvector.data[i]) 
    i++;
  return (i==n);  
}

bool Vector::operator<(int number)
{
  int i=0;
  while (i < n && data[i] < number) 
    i++;
  return (i==n);  

}

Vector Vector::operator-(){

  Vector newVector(n);

  for (int i=0; i<n; i++) 
    newVector.data[i] = -data[i];

  return newVector;
}

Vector Vector::operator/(double number){

  Vector newVector(n);

  for (int i=0; i<n; i++) 
    newVector.data[i] = data[i]/number;

  return newVector;

}
Vector Vector::maximize(Vector& gvector) {
  Vector newVector(n);
  for (int i=0; i<n; i++) {
    if (data[i] >= gvector.data[i]) 
      newVector.data[i] = data[i];
    else
      newVector.data[i] = gvector.data[i];
  }
  return newVector;
}

Vector Vector::minimize(Vector& gvector) {
  Vector newVector(n);
  for (int i=0; i<n; i++) {
    if (data[i] <= gvector.data[i]) 
      newVector.data[i] = data[i];
    else
      newVector.data[i] = gvector.data[i];
  }
  return newVector;
}
Vector Vector::threshold(Vector& gvector) {
  Vector newVector(n);
  for (int i=0; i<n; i++) {
    if (data[i] >= gvector.data[i]) 
      newVector.data[i] = 1.0;
    else
      newVector.data[i] = 0.0;
  }
  return newVector;
}

Vector Vector::absVal() {
  Vector newVector(n);
  for (int i=0; i<n; i++)
    newVector.data[i] = vabs(data[i]);
  return newVector;
}

const double Vector::hashCode() {
  // I guess for now, I'll try using {data[0]+sum_{i=1}^{n-1} PI*prime[i]*data[i]}
  // where prime[i] is the i^th in a series of prime numbers {0,1,3,5,7} (since 
  // must also include the type in here - 
  // need n+2 prime numbers n = number of floors
  // the Sears Tower (Chicago) has 110 floors, as do the World Trade Center Buildings in NY.
  // No buildings with over this amount are scheduled to be built in the near future,
  // so 200 should be good enough.

  int i;
  double prime[202] = {0,1,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229};

  double sum = data[0];
  for (i=1; i<n; i++) {
    sum += M_PI*prime[i]*data[i];
  }
  sum += M_PI*prime[i++]*type;
  sum += M_PI*prime[i]*n;
  return sum;
}

char* Vector::toString(){

  char *newstring;
  char tempstring[32];
  newstring = (char *)malloc(sizeof(char)*256);
  sprintf(newstring,"%g \\n",data[0]);
  for (int i=1; i<n; i++) {
    sprintf(tempstring,"%g \\n",data[i]);
    strcat(newstring,tempstring);
  }
  
  return newstring;
}

Vector ceill(const Vector& source){
  
  Vector newVector(source.n);

  for (int i=0; i<source.n; i++) 
    newVector.data[i] = ceil(source.data[i]);

  return newVector;
}
