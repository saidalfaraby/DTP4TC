#include "pair.hh"

//Constructors
Pair::Pair(){
  min = 0.0;
  max = 0.0;
}
Pair::Pair(double mimax){
  min = max = mimax;
}
Pair::Pair(double mi, double ma){

  min = mi;
  max = ma;

}

Pair::Pair(const Pair& gpair){

  min = gpair.min;
  max = gpair.max;
}
//destructors?????
Pair::~Pair(){
}

void  Pair::set_max(double ma) {
  max = ma;
}

void Pair::set_min(double mi) {
  min = mi;
}

void Pair::set_both(double mi, double ma) {
  min = mi;
  max = ma;
}

void Pair::set(double minmax){
  min = minmax;
  max = minmax;
}

void Pair::span(Pair &gpair) {
  if (gpair.min < min)
    min = gpair.min;
  if (gpair.max > max)
    max = gpair.max; 
}

Pair Pair::operator+(Pair& gpair){
  
  Pair newPair;
  
  newPair.min =  min + gpair.min;
  newPair.max = max + gpair.max;

  return newPair;
}
void Pair::operator=(Pair gpair){
  
  min = gpair.min;
  max = gpair.max;

  //return *this; 
}
bool Pair::operator==(Pair& gpair){

  //if(min == gpair.min && max == gpair.max)
  if(fabs(min-gpair.min)<1.0e-12 && fabs(max-gpair.max)<1.0e-12)
    return true;
  else 
    return false;
  
}
const double Pair::hashCode() {
  return (1.0e-12*(floor(1.0e+12*(min + PI*max))));
}
Pair Pair::operator-(Pair& gpair){
  
  Pair newPair;
 
  newPair.min =  min - gpair.min;
  newPair.max = max - gpair.max;

  return newPair;
}

Pair Pair::operator *(Pair& gpair){
  
  Pair newPair;
 
  newPair.min =  min * gpair.min;
  newPair.max = max * gpair.max;
  
  //  fprintf(stderr,"%20.14f %20.14f %20.14f %20.14f %20.14f %20.14f\n",min,max,gpair.min,gpair.max,newPair.min,newPair.max);

  return newPair;
}


Pair Pair::operator *(double number){
  
  Pair newPair;
 
  newPair.min =  min * number;
  newPair.max = max * number;

  return newPair;
}

Pair Pair::operator /(Pair& gpair){
  
  Pair newPair;
  
  newPair.min = min/gpair.min;
  newPair.max = max/gpair.max;

  return newPair;
}

bool Pair::operator!=(Pair& gpair){
  
  //if( min == gpair.min && max == gpair.max)
  if(fabs(min-gpair.min)<1.0e-12 && fabs(max-gpair.max)<1.0e-12)
    return false;
  else
    return true;
}

bool Pair::operator<(Pair& gpair){
  
  if( min < gpair.min && max < gpair.max)
    return true;
  else 
    return false;
  
}

bool  Pair::operator<=(Pair& gpair){
  
  if( min <= gpair.min && max <= gpair.max)
    return true;
  else 
    return false;
  
}

bool Pair::operator>(Pair& gpair){
  
  if( min > gpair.min && max > gpair.max)
    return true;
  else 
    return false;
  
}

bool Pair::operator>=(Pair& gpair){
  
  if( min >= gpair.min && max >= gpair.max)
    return true;
  else 
    return false;
  
}
/*
bool Pair::operator<(Pair& gpair){  
  return( min < gpair.min);
}

bool  Pair::operator<=(Pair& gpair){
  return(min <= gpair.min);
}

bool Pair::operator>(Pair& gpair){
  return( max > gpair.max);
}

bool Pair::operator>=(Pair& gpair){
  return(max >= gpair.max);
}

*/
bool Pair::operator<(int number){

  if(max < (double) number)
    return true;
  else
    return false;
}

Pair Pair::operator-(){

  Pair newPair;

  newPair.min = -max;
  newPair.max = - min;
  
  return newPair;
}

Pair Pair::operator/(double number){

  Pair newPair;
  newPair.min = number/max;
  newPair.max = number/min;

  return newPair;

}

Pair Pair::absVal() {
  Pair newPair;
  newPair.min = absv(min);
  newPair.max = absv(max);
  return newPair;
}

char* Pair::toString(){

  char *newstring;
  newstring = (char *)malloc(sizeof(char)*256); 
  sprintf(newstring,"%20.14f %20.14f",min,max);
  
  return newstring;
}

int Pair::readString(FILE *fp, char * buf) {
  char buf1[100];
  int retval;
  retval = fscanf(fp,"%s %s\n",buf,buf1);
  strcat(buf," ");
  strcat(buf,buf1);
  return retval;
}
void Pair::parseString(char *buf) {
  double x1,x2;
  sscanf(buf,"%lf %lf",&x1,&x2);
  min = x1;
  max = x2;
} 






Pair ceill(const Pair& source){
  
  Pair newPair;

  newPair.min = ceil(source.min);
  newPair.max = ceil(source.max);

  return newPair;
}
