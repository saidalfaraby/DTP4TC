#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>

#include <sys/poll.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <string.h>


#ifndef _NET2_H_
#define _NET2_H_

#define MAX_PORT_TRY 100

//default range of ports to start a server
#define DEFAULT_START_PORT 5000
#define DEFAULT_END_PORT   5100

class Server { 

 private:
  int ns;
  int start_port;
  int end_port;

  //from initServer local to class members
  struct sockaddr_in sin, sfrom;
  int sfrom_len;
  int s;

 public:
  Server ( void );
  Server ( unsigned int start_port);
  Server ( unsigned int start_port, unsigned int end_port );

  ~Server ( );

 public:
  int net_write_data (void *data, int size);
  int net_get_data (void *data, int size);
  int startServer ( void );
  int waitForClient ( void );
  int  isDataAvailableToRead ( void );

 private:
  int initServer (int* ns, unsigned int port);
  int initServer ( int* ns);
  int getServer ( int*ns, unsigned int start_port, unsigned int end_port );
  int getServer(int *ns, unsigned int port);  
  void closeServer ( int ns);
  int net_write_data (int ns, void *data, int size);
  
};


class Client {

 private:
  int ns;
  int start_port;
  int end_port;
  char* serverName;
  
 public:
  Client ( void );
  Client ( unsigned int start_port, char* name );
  Client ( unsigned int start_port, unsigned int end_port, char* name );

  ~Client ( void );

 public:
  int net_write_data (void *data, int size);
  int net_get_data (void *data, int size);
  int  isDataAvailableToRead ( void );
  int startClient ( void );


 private:
  int  initClient (int* ns,unsigned int port,char *serverName);
  int  getClient(int *ns, unsigned int port, char *serverName);
  int  getClient ( int* ns, unsigned int start_port, unsigned int end_port, 
		   char* serverName ); 
  void closeClient (int ns);
};



#endif
