#include "net2.h"
#include "networkSettings.h"
#define MDP_S  5400
#define MDP_C  5410
// where the server is running
#define MDP_HOST "neuron.ai.toronto.edu"  //128.100.3.14"

Client *mdp_c;
Server *mdp_s;
#define TRYSERVER 1

int connectToMDP() {
  mdp_c = new Client(MDP_C,MDP_C+PORT_RANGE,MDP_HOST );
  if ( mdp_c && mdp_c->startClient ( ) ) {
    //start a server and wait for a connection
    fprintf(stderr,"starting a new server\n");
    mdp_s = new Server (MDP_S,MDP_S+PORT_RANGE);
    if ( mdp_s ) {
      mdp_s->startServer();
      fprintf(stderr,"waiting for client to connect\n");
      return (mdp_s->waitForClient ( ));
    }
  }
  return false;
}

int main(int argc, char *argv[])
{
  int i,j,k;
  if (!connectToMDP()) {
    fprintf(stderr,"could not connect\n");
    exit(0);
  }
  int numovars,numvars;
  numovars = 6;
  numvars = 6;
  int numstates = 64;
  int * currState = new int[numovars];
  for (j=0; j<numovars; j++)
    currState[j] = 0;
  int bestact = 0;
  char response;
  int val,var;
  // send state
  k=0;
  //while (k<numstates) {
    fprintf(stderr,"sending  ");
    for (int i=0; i<numovars; i++)
      fprintf(stderr,"%d ",currState[i]);
    mdp_s->net_write_data(currState, numovars*sizeof(int));
    // wait for action
    while (!mdp_c->isDataAvailableToRead());
    mdp_c->net_get_data(&bestact, sizeof(int));
    
    fprintf(stderr,"best action in state is %d\n",bestact);
    j=0;
    while (j<numovars && currState[j] == 1) {
      currState[j] = 0;
      j++;
    }
    if (j<numovars)
      currState[j] = 1;
    k++;
      
    //}  
  delete mdp_s;
  delete mdp_c;
}

