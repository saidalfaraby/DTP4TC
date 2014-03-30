#include "mex.h"
#include "MDP.h"

/* Input arguments */
#define FILENAME_IN   0
#define STATE_IN   1

/* Output arguments */
#define OPTACT_OUT    0
#define MDP_S  5400
#define MDP_C  5510
#define MDP_HOST "128.100.3.14"
Client *mdp_c;
Server *mdp_s;

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


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const int* sizes; 
  double *optact, *doubleArray;

  /* retrieve filename */
  int filenameLength = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  char *filename = (char*)mxCalloc(filenameLength, sizeof(char));
  mxGetString(prhs[FILENAME_IN], filename, filenameLength);

  /* retrieve state */
  doubleArray = mxGetPr(prhs[STATE_IN]);
  sizes = mxGetDimensions(prhs[STATE_IN]);
  int *state = new int[sizes[1]];
  for (int i=0; i<sizes[1]; i++) {
    state[i] = (int)doubleArray[i];
    printf("%i\n",state[i]);
  }

  /* return action */
  plhs[OPTACT_OUT] = mxCreateDoubleMatrix(1,1,mxREAL);
  optact = mxGetPr(plhs[OPTACT_OUT]);
  
  // connect to server
  if (!connectToMDP()) {
    fprintf(stderr,"could not connect\n");
    return;
  }

  mdp_s->net_write_data(state, sizes[1]*sizeof(int));
  // wait for action
  while (!mdp_c->isDataAvailableToRead());
  mdp_c->net_get_data(optact, sizeof(int));
  
  fprintf(stderr,"best action in state is %d\n",*optact);

  delete mdp_s;
  delete mdp_c;


  /* call fn */
}







