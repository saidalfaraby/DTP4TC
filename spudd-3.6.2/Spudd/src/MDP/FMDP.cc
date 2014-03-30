#include "FMDP.h"

int connectToGame    ( Client **game_c, Server **game_s ) {
  // start new client connecting to server
  fprintf(stderr,"attempting to connect to Game server\n");
  *game_c = new Client(GAME_C,GAME_C+PORT_RANGE,GAME_HOST);
  (*game_c)->startClient ( );
 

  //connect to the server  
  fprintf(stderr,"started client - starting server\n");

  //start a server and wait for a connection
  *game_s = new Server (GAME_S,GAME_S+PORT_RANGE);
  if(*game_s) {
    fprintf(stderr,"starting server - waiting for client\n");
    (*game_s)->startServer ( );
    if ((*game_s)->waitForClient ( ) ) 
      fprintf(stderr,"started server and connected client\n");
    return true;
  }
  return false;
}
int connectToCMG   ( Client **cmg_c, Server **cmg_s ) {
  // start new client connecting to server
  fprintf(stderr,"attempting to connect to Cmg server\n");
  *cmg_c = new Client(CMG_C,CMG_C+PORT_RANGE,CMG_HOST);
  (*cmg_c)->startClient ( );
 

  //connect to the server  
  fprintf(stderr,"started client - starting server\n");

  //start a server and wait for a connection
  *cmg_s = new Server (CMG_S,CMG_S+PORT_RANGE);
  if(*cmg_s) {
    fprintf(stderr,"starting server - waiting for client\n");
    (*cmg_s)->startServer ( );
    if ((*cmg_s)->waitForClient ( ) ) 
      fprintf(stderr,"started server and connected client\n");
    return true;
  }
  return false;
}
// for connecting to face
int connectToFace    (Server **face_s) {
  //start a server and wait for a connection
  *face_s = new Server (FACE_S,FACE_S+PORT_RANGE);
  if ( *face_s ){
    (*face_s)->startServer ( );
    if ( (*face_s)->waitForClient ( ) )
      return true;
  }
  return false;
}
int FMDP::getnumvars() {
  return numvars;
}

void FMDP::policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars) {
  Pair dval, aval;
  int numgamevars;
  int agent_face, user_face, agent_action, user_action, actFlag, temp;
  Client *game_c, *cmg_c;
  Server *game_s, *face_s, *cmg_s;

  // connect to client and server
  if (!connectToGame(&game_c,&game_s)) {
    fprintf(stderr,"Error - could not connect to game\n");
    exit(0);
  }
  // read back number of variables provided by game
  while (!game_c->isDataAvailableToRead());
  game_c->net_get_data(&numgamevars, sizeof(int));

  // connect to CMG
  if (!connectToCMG(&cmg_c,&cmg_s)) {
    fprintf(stderr,"Error - could not connect to cmg\n");
    exit(0);
  }    

  // connect to Face
  fprintf(stderr,"connecting to face\n");
  if (!connectToFace(&face_s)) {
    fprintf(stderr,"Error = could not connect to face\n");
    exit(0);
  }
  // send number of faces (agent_faces) expected

  int *gamevarvals = new int[numgamevars];
  int *gamevarindices = new int[numgamevars];
  //gamevarvals[0] = payout
  //gamevarvals[1] = turn
  //gamevarvals[2] = user_face (supplied here)
  //gamevarvals[3] = bestcard
  //gamevarvals[4] = Av (user_action)
  gamevarindices[0] = 0;
  gamevarindices[1] = 4;
  // 2 is user_face - supplied by FMDP
  gamevarindices[2] = 3;
  gamevarindices[3] = 1;
  
  int *varvals = new int[numovars];
  int *privatevarvals = new int[numovars-numgamevars];
  int i,j;
  int bestact;
  char namestr[1024];
  while (true) {
    if (game_c->isDataAvailableToRead()) {
      int size = numgamevars*sizeof(int);
      game_c->net_get_data (gamevarvals, size);
      
      fprintf(stderr,"received state from GAME: ");
      for (i=0; i<numgamevars; i++) 
	fprintf(stderr,"%d ",gamevarvals[i]);
      fprintf(stderr,"\n");
      // add to varvals other private info (here faces)
      j=0;
      for (i=0; i<numovars; i++)
	varvals[i] = 0;
      for (i=0; i<numgamevars; i++) 
	varvals[gamevarindices[i]] = gamevarvals[i];
      
      // gather reward? 
      
      // consult F policy --> agent_face
      fprintf(stderr,"consulting face policy in state \n");
      for (i=0; i<numovars; i++)
	fprintf(stderr,"%d ",varvals[i]);
      fprintf(stderr,"\n");
      
      getAV(gb,act,val,aval,dval,varvals,vars,numvars,orig_vars,numovars);
      if ((temp = pickAction((int) (aval.get_min())))<0) {
	fprintf(stderr,"error in pickAction \n");
	exit(0);
      }
      
      // extract F from A,F (Which is temp)
      agent_face = temp%4;

      fprintf(stderr,"temp was %d performing face %d in state \n",temp,agent_face);
      for (i=0; i<numovars; i++)
	fprintf(stderr,"%d ",varvals[i]);
      fprintf(stderr,"\n");

      // perform the facial display - then return to neutral
      face_s->net_write_data(&agent_face,sizeof(int));
      // send neutral 
      i=5;
      face_s->net_write_data(&i,sizeof(int));
      
      //wait for user acted flag
      while (!game_c->isDataAvailableToRead());
      game_c->net_get_data(&actFlag,sizeof(int));

      // segment and classify F --> gives user_face
      // at prent the state is just whose turn it is
      int state=gamevarvals[1];
      fprintf(stderr,"writing state %d to CMG\n",state);
      cmg_s->net_write_data(&state,sizeof(int));
      while (!cmg_c->isDataAvailableToRead());
      cmg_c->net_get_data(&user_face,sizeof(int));
      fprintf(stderr,"read user_face %d from CMG\n",user_face);
      // consult A policy --> gives agent_action
      // after updating varvals
      varvals[2] = user_face;

      fprintf(stderr,"consutling action policy in state:\n");
      for (i=0; i<numovars; i++)
	fprintf(stderr,"%d ",varvals[i]);
      fprintf(stderr,"\n");

      getAV(gb,act,val,aval,dval,varvals,vars,numvars,orig_vars,numorigvars);

      temp = pickAction((int) (aval.get_min()));

      // extract A from A,F (Which is temp)
      // this hack is only temporary (for this game only)
      agent_action = temp/4;
 
      strcpy(namestr,"");
      aconvert(namestr,actionnames,aval.get_min());
      fprintf(stderr,"performing action %d which is %s in state \n",agent_action,namestr);
      for (i=0; i<numovars; i++)
	fprintf(stderr,"%d ",varvals[i]);
      fprintf(stderr,"\n");
      
      // perform the action
      game_s->net_write_data(&agent_action,sizeof(int));
      
      // wait for response
      while (!game_c->isDataAvailableToRead());
      game_c->net_get_data(&temp,sizeof(int));
      agent_action = temp/4;
      while (!game_c->isDataAvailableToRead());
      game_c->net_get_data(&temp,sizeof(int));
      user_action = temp;

      fprintf(stderr,"received actions agent: %d user: %d\n",agent_action,user_action);

      // possibly do some more state updates?
    } else {
      usleep(1000);
    }
  }
}
