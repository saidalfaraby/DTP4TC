#include "net2.h"
#include "networkSettings.h"

int main(int argc, char **argv) {
  // fire up server
  Server *spudd_s;
  Client *spudd_c;
  // connect client
  fprintf(stderr,"connecting client...");
  spudd_c = new Client(SPUDD_S,SPUDD_S+PORT_RANGE,"sungod.cs.ubc.ca");
  if (!spudd_c)
    exit(0);
  if (!spudd_c->startClient ( ))
    exit(0);
  spudd_s = new Server(SPUDD_C,SPUDD_C+PORT_RANGE);
  if (!spudd_s) {
    exit(0);
  }
  spudd_s->startServer ();
  if (!spudd_s->waitForClient())
    exit(0);
  
  
  int numovars = 5;
  int *varvals = new int[numovars];
  int bestact,i;
  // initial state
  for (i=0; i<numovars; i++) 
    varvals[i] = 1;
  varvals[0] = 0;
  varvals[4] = 3;

  
  fprintf(stderr,"sending ");
  for (i=0; i<numovars; i++) 
    fprintf(stderr,"%d ",varvals[i]);

  // write this to the spudd engine
  spudd_s->net_write_data(varvals,numovars*sizeof(int));
  
  // read the response
  spudd_c->net_get_data(&bestact,sizeof(int));
  fprintf(stderr,"       optimal action is %d\n",bestact);

  // change the state a bit
  varvals[4] = 0;

  // write this to the spudd engine
  spudd_s->net_write_data(varvals,numovars*sizeof(int));
  
  // read the response
  spudd_c->net_get_data(&bestact,sizeof(int));
  fprintf(stderr,"       optimal action is %d\n",bestact);


  // change the state a bit
  varvals[0] = 1;

  // write this to the spudd engine
  spudd_s->net_write_data(varvals,numovars*sizeof(int));
  
  // read the response
  spudd_c->net_get_data(&bestact,sizeof(int));
  fprintf(stderr,"       optimal action is %d\n",bestact);

  // change the state a bit
  varvals[3] = 0;

  // write this to the spudd engine
  spudd_s->net_write_data(varvals,numovars*sizeof(int));
  
  // read the response
  spudd_c->net_get_data(&bestact,sizeof(int));
  fprintf(stderr,"       optimal action is %d\n",bestact);
      
}
