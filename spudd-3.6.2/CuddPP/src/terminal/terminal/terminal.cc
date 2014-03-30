#include "terminal.hh"


//Constructors
Terminal::Terminal(){
  value = 0.0;
}

Terminal::Terminal(double val){
  value = val;
}

Terminal::Terminal(const Terminal& term){
 value = term.value;
}

//destructors?????
Terminal::~Terminal(){
}

void Terminal::set(double val){
  value = val;
}

Terminal Terminal::operator+(Terminal& term){
  
  Terminal newTerminal;
  
  newTerminal.value =  value + term.value;
  return newTerminal;
}
void Terminal::operator=(Terminal term){
  
  value = term.value;
}
bool Terminal::operator==(Terminal& term){
  if(value == term.value)
    return true;
  else 
    return false;
}

Terminal Terminal::operator-(Terminal& term){
  
  Terminal newTerminal;
  
  newTerminal.value =  value - term.value;
  return newTerminal;
}

Terminal Terminal::operator *(Terminal& term){
  
  Terminal newTerminal;
 
  newTerminal.value =  value * term.value;

  return newTerminal;
}


Terminal Terminal::operator *(double number){
  
  Terminal newTerminal;
 
  newTerminal.value =  value * number;

  return newTerminal;
}

Terminal Terminal::operator /(Terminal& term){
  
  Terminal newTerminal;
  
  newTerminal.value = value/term.value;

  return newTerminal;
}

bool Terminal::operator!=(Terminal& term){
  
  if( value == term.value)
    return false;
  else
    return true;
}

bool Terminal::operator<(Terminal& term){
  
  if( value < term.value)
    return true;
  else 
    return false;
  
}

bool  Terminal::operator<=(Terminal& term){
  
  if( value <= term.value)
    return true;
  else 
    return false;
  
}

bool Terminal::operator>(Terminal& term){
  
  if( value > term.value)
    return true;
  else 
    return false;
  
}

bool Terminal::operator>=(Terminal& term){
  
  if( value >= term.value)
    return true;
  else 
    return false;
  
}

bool Terminal::operator<(int number){

  if(value < (double) number)
    return true;
  else
    return false;
}

Terminal Terminal::operator-(){

  Terminal newTerminal;

  newTerminal.value = -value;
  
  return newTerminal;
}

Terminal Terminal::operator/(double number){

  Terminal newTerminal;
  newTerminal.value = value/number;

  return newTerminal;

}

Terminal Terminal::absVal() {
  Terminal newTerminal;
  newTerminal.value = absv(value);
  return newTerminal;
}

char* Terminal::toString(){

  char *newstring;
  newstring = (char *)malloc(sizeof(char)*256); 
  sprintf(newstring,"%g",value);
  
  return newstring;
}


int Terminal::readString(FILE *fp, char * buf) {
  int retval;
  retval = fscanf(fp,"%s",buf);
  return retval;
}
void Terminal::parseString(char *buf) {
  double x1;
  sscanf(buf,"%lf",&x1);
  value = x1;
} 


Terminal ceill(const Terminal& source){
  
  Terminal newTerminal;

  newTerminal.value = ceil(source.value);

  return newTerminal;
}


