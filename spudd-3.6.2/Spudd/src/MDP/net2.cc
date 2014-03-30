#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include "net2.h"
//#include "slm_msg.h"

#define MYPORT 5100


Server::Server ( void ){
  ns = -1;
  this->start_port = DEFAULT_START_PORT;
  this->end_port = DEFAULT_END_PORT;
}

Server::Server ( unsigned int start_port ){
  ns = -1;
  this->start_port = start_port;
  this->end_port = start_port; //it only tries at the port given
}


Server::Server ( unsigned int start_port, unsigned int end_port ){
  ns = -1;
  this->start_port = start_port;
  this->end_port = end_port;
}

Server::~Server ( ) {
  closeServer ( ns );
}

int Server::startServer ( void ) {
  return getServer ( &ns, start_port, end_port );
}

int Server::getServer ( int*ns, 
			unsigned int start_port, 
			unsigned int end_port ) {

  for ( int i = start_port; i<=end_port; i++ ){
    if ( getServer ( ns, i ) != 0 ) {
      return true;
    }
  }
  return false;
}

int Server::getServer(int *ns, unsigned int port) {

  if ( !initServer ( ns, port ) )
    return 0;
  else 
    return port;
}

//static int ns;
// It starts a server. The server waits for a single
// client to connect.
int Server::initServer ( int* ns, unsigned int port){
  //struct sockaddr_in sin, sfrom;
  // int sfrom_len;
  // int s;
  sin.sin_family = AF_INET;
  sin.sin_addr.s_addr = INADDR_ANY;
  sin.sin_port = htons(port);
  
  s = socket(AF_INET, SOCK_STREAM, 0);
  if ( bind(s, (struct sockaddr *) &sin, sizeof (sin)) ) {
    printf ("could not bind\n");
    return false;
  }
  printf ("bound to socket %d\n", sin.sin_port);
  
  listen(s, 5);
  
  return true;
}

int Server::waitForClient ( void ) {
  sfrom_len = sizeof (sfrom);
  printf ("\nwaiting for connection...\n");
  ns = accept(s, (struct sockaddr *)&sfrom, (socklen_t *)&sfrom_len);
  if ( ns < 0) {
    printf ("failed\n");
  }
  printf ("connected on %d\n", ns);  
  return ns;
}

// It starts a server. The server waits for a single
// client to connect.
int Server::initServer ( int* ns ){
  struct sockaddr_in sin, sfrom;
  int sfrom_len;
  int s;
  sin.sin_family = AF_INET;
  sin.sin_addr.s_addr = INADDR_ANY;
  sin.sin_port = htons(MYPORT);
  
  s = socket(AF_INET, SOCK_STREAM, 0);
  if ( bind(s, (struct sockaddr *) &sin, sizeof (sin)) ) {
    printf ("could not bind\n");
    return false;
  }
  printf ("bound to socket %d\n", sin.sin_port);
  
  listen(s, 5);
  
  sfrom_len = sizeof (sfrom);
  printf ("\nwaiting for connection...\n");
  *ns = accept(s, (struct sockaddr *)&sfrom, (socklen_t *)&sfrom_len);
  if ( ns < 0) {
    printf ("failed\n");
    return false;
  }
  printf ("connected on %d\n", *ns);
  return true;
}

int Server::net_get_data (void *data, int size) {
  int size_read = 0; 
  int total_read = 0 ;
  int total_remain = size;
  do { // loops until it reads all of the data or an error occurs.
    size_read    = read (ns, data, (size_t)total_remain);
    if ( size_read == -1 ){ // error
      printf("read has returned an error!\n");
      return false;
    }
    total_read   += size_read;
    total_remain = size-total_read;
    //printf("Read %d, Remain %d out of %d\n", 
    //total_read, total_remain, size );
    data = ((char*)data) + size_read;
  } while ( total_read < size );
  return true;
}

int Server::isDataAvailableToRead ( void ) {
  
  int ret =  0;
  struct pollfd res;
 
  res.fd = ns;
  res.events = 0;
  res.events |= POLLIN;

  ret = poll ( &res, 1, 10 /* ms */ );
  if ( ret == 0 ){
    //printf("No events for %d\n", ns );
  }
  else if ( ret == -1 ) {
    printf("Error calling poll on %d\n", ns);
  }
  else {
    //printf("Event for %d: ", ns);
    if ( res.revents | POLLIN ){
      //printf("POLLIN\n");
    }
    else if ( res.revents | POLLOUT ){
      //printf("POLLOUT\n");
    }
    else if ( res.revents | POLLPRI ){
      // printf("POLLPRI\n");
    }
    else if ( res.revents | POLLERR ){
      //printf("POLLERR\n");
    }
    else if ( res.revents | POLLHUP ){
      //printf("POLLHUP\n");
    }
    else if ( res.revents | POLLNVAL ){
      //printf("POLLVAL\n");
    }
  }  
  return ret;
}

int Server::net_write_data (void *data, int size) {
  return net_write_data ( ns, data, size );  
}

// Used to write data to the socket
int Server::net_write_data (int ns, void *data, int size){
  int size_wrote = 0; 
  int total_wrote = 0;
  int total_remain = size;
  do {

    //    size_wrote = send ( ns, 
    //			(const void *)data,(size_t) size, 
    //			MSG_NOSIGNAL | MSG_DONTWAIT);
    size_wrote = write (ns, (const void *)data,(size_t)total_remain); 
    
    if ( size_wrote == -1 ) {

      printf("send has return an error!\n");
      return false;
    }
    total_wrote += size_wrote;
    total_remain = size-total_wrote;
    data = ((char*)data) + size_wrote;
  } while ( total_wrote < size );
  /*
    if ( (size_wrote = send (ns, (const void *)data,(size_t) size, 
    MSG_NOSIGNAL | MSG_DONTWAIT) )!= size ) {   
    printf ("Error writing data to socket %d \n", size_wrote);
    return false;
    }
  */
  return true;
  //printf ( "  -- Wrote %d bytes of data to %d -- \n", size_wrote, ns );
}

// Closes the server.
void Server::closeServer ( int ns ){
  close ( ns );
}

/********************************************************************/
/********************************************************************/
/************                                              **********/
/************                   CLIENT                     **********/
/************                                              **********/
/********************************************************************/
/********************************************************************/
Client::Client ( void ) {
  ns =  -1;
  this->start_port = DEFAULT_START_PORT;
  this->end_port   = DEFAULT_END_PORT;
  serverName = NULL;
}
Client::Client ( unsigned int start_port, char* name ){
  ns = -1;
  this->start_port = start_port;
  this->end_port = start_port;
  serverName = new char[strlen(name)];
  strcpy ( serverName, name );
}

Client::Client ( unsigned int start_port, unsigned int end_port, char* name ){

  ns = -1;
  this->start_port = start_port;
  this->end_port   = end_port;
  serverName = new char[strlen(name)];
  strcpy ( serverName, name );
}

Client::~Client ( void ){
  closeClient ( ns );
}


int Client::startClient ( void ) {
  return getClient ( &ns, start_port, end_port, serverName);
}

int Client::getClient ( int* ns, 
			unsigned int start_port, 
			unsigned int end_port, 
			char* serverName ) {
  for ( int i = end_port; i >= start_port; i-- ) {
    if ( getClient  ( ns, i, serverName ) != 0 )
      return true;
  }
  return false;

}

int Client::getClient(int *ns, unsigned int port, char *serverName){
  int i=0;

  if ( !initClient( ns, port, serverName )) 
    return 0;
  else
    return port;
}
// It starts a client. The client connects to a server
// at the specified port.
int Client::initClient (int* ns, unsigned int port, char *serverName){
  hostent *hp;
  struct sockaddr_in server;
  
  hp = gethostbyname(serverName);
  if (hp == NULL) {
    fprintf(stderr, "rlogin: %s: unknown host\n", serverName);
    return 0;
  }
  //  printf ("found %s\n", hp->h_name);
  
  bzero((char *)&server, sizeof (server));
  bcopy(hp->h_addr, (char *) &server.sin_addr, hp->h_length);
  server.sin_family = hp->h_addrtype;
  server.sin_port = htons( port );
  
  *ns = socket(hp->h_addrtype, SOCK_STREAM, 0);
  if (*ns < 0) {
    //perror("rlogin: socket");
    return 0;
  }
  
  if (connect(*ns, (struct sockaddr *)&server, sizeof (server)) < 0) {
    //perror("rlogin: connect");
    return 0;
  }

  printf ("connected to port %d\n", port);
  return 1;
}
// Closes the client
void Client::closeClient ( int ns ){
  close ( ns );
}

int Client::net_write_data (void *data, int size) {
  int size_wrote; 
  //  if ( (size_wrote = write (ns, (const void *)data,(size_t) size)) != size ) {   
  if ( (size_wrote = send (ns, (const void *)data,(size_t) size, 
			   MSG_NOSIGNAL | MSG_DONTWAIT) )!= size ) {   
    printf ("Error writing data to socket %d\n", size_wrote);
    return false;
  }
  return true;  
}

// Used to read data from a socket.
int Client::net_get_data ( void *data, int size){
    
  int size_read = 0; 
  int total_read = 0 ;
  int total_remain = size;
  do { // loops until it reads all of the data or an error occurs.
    size_read    = read (ns, data, (size_t)total_remain);
    if ( size_read == -1 ){ // error
      printf("read has returned an error!\n");
      return false;
    }
    total_read   += size_read;
    total_remain = size-total_read;
    data = ((char*)data) + size_read;
    //printf("Read %d, Remain %d out of %d\n", 
    //total_read, total_remain, size );
  } while ( total_read < size );
  return true;
}

int Client::isDataAvailableToRead ( void ) {
  
  int ret =  0;
  struct pollfd res;
 
  res.fd = ns;
  res.events = 0;
  res.events |= POLLIN;

  ret = poll ( &res, 1, 10 /* ms */ );
  if ( ret == 0 ){
    //printf("No events for %d\n", ns );
  }
  else if ( ret == -1 ) {
    printf("Error calling poll on %d\n", ns);
  }
  else {
    //printf("Event for %d: ", ns);
    if ( res.revents | POLLIN ){
      //printf("POLLIN\n");
    }
    else if ( res.revents | POLLOUT ){
      //printf("POLLOUT\n");
    }
    else if ( res.revents | POLLPRI ){
      // printf("POLLPRI\n");
    }
    else if ( res.revents | POLLERR ){
      //printf("POLLERR\n");
    }
    else if ( res.revents | POLLHUP ){
      //printf("POLLHUP\n");
    }
    else if ( res.revents | POLLNVAL ){
      //printf("POLLVAL\n");
    }
  }  
  return ret;
}
