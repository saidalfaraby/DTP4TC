#include "pspudd.hh"
// these are in the homer project
#include "net2.h"
#include "networkSettings.h"

Client *spudd_c;
Server *spudd_s;
#define SUPERVISOR_HOST_TMP "sungod.cs.ubc.ca"
// for connecting to supervisor
int connectToSupervisor     ( void ) {
  //start a server and wait for a connection
  spudd_s = new Server (SPUDD_S,SPUDD_S+PORT_RANGE);
  if ( spudd_s ){
    spudd_s->startServer ( );
    if ( spudd_s->waitForClient ( ) ) {
      //Connect to the server
      sleep(1);
      spudd_c = new Client(SPUDD_C,SPUDD_C+PORT_RANGE,SUPERVISOR_HOST_TMP );
      if ( spudd_c && spudd_c->startClient ( ) ) 
	return true;
    }
  }
  return false;
}



void getAV(DdManager *dd, DdNode *val, DdNode *act, Pair & dval, Pair & aval, int *varvals) {
  // query the  policy action and value with the current values of varvals
  // build an array varass of 1s and 0s over the variables v
  // with assignments corresponding to the varvals (values for the original variables  ov)
  // then call *(Cudd_V(Cudd_Eval(dd,act, varass)))
  int nvars = numvars;
  int novars = numorigvars;
  rnum *v = vars;
  onum *ov = orig_vars;
  int *varass = new int[2*nvars];
  int i,j,nbv,tmp,nbvals;
  for (i=0; i<nvars*2; i++)
    varass[i] = 0;

  for (i=0; i<novars; i++) {
    nbv = ov[i].nbvars;
    nbvals = int(pow(2.0,nbv));
    tmp = nbvals-varvals[i]-1;
    for (j=nbv-1; j>=0; j--) {
      varass[Cudd_NodeReadIndex(v[ov[i].var1index+j].add_var)] = tmp%2;
      tmp = tmp/2;
    }
  }
  fprintf(stderr,"varass is ");
  for (i=0; i<nvars*2; i++)
    fprintf(stderr,"%d ",varass[i]);
  fprintf(stderr,"\n");

  dval = *(Cudd_V(Cudd_Eval(dd,val,varass)));
  aval = *(Cudd_V(Cudd_Eval(dd,act,varass)));
  
  delete [] varass;
}

void policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars) {
  Pair dval, aval;

  // connect to client and server
  if (!connectToSupervisor()) {
    fprintf(stderr,"Error - could not connect to supervisor\n");
    exit(0);
  }

  unsigned char * data = new unsigned char[numovars];
  int *varvals = new int[numovars];
  int bestact;
  char namestr[1024];
  while (true) {
    if (spudd_c->isDataAvailableToRead()) {
      int size = numovars*sizeof(int);
      spudd_c->net_get_data (varvals, size);
      
      /*
      fprintf(stderr,"received ");
      for (int i=0; i<numovars; i++)
	fprintf(stderr,"%d ",varvals[i]);
      */

      getAV(gb,act,val,aval,dval,varvals);
      
      // send back the action
      bestact = (int) (aval.get_min()-1);

      /*
      strcpy(namestr,"");
      aconvert(namestr,actionnames,aval.get_min());
      fprintf(stderr,"   sending %d which is %s\n",bestact,namestr);
      */

      spudd_s->net_write_data(&bestact,sizeof(int));
    } else {
      usleep(1000);
    }
  }
}
