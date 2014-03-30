#include "MDP.h"

int findString(char **arrayOfStrings, int numStrings, char *string);
//#define HOME
#ifdef HOME
#define GAME_HOST  "198.162.38.54"
#define FACE_HOST "198.162.38.54"
#define CMG_HOST "198.162.38.54"
#else
#define GAME_HOST  "sungod.cs.ubc.ca"
#define FACE_HOST "sungod.cs.ubc.ca"
#define CMG_HOST "sungod.cs.ubc.ca"
#endif
#define GAME_C 5100
#define GAME_S 5300
#define CMG_C 5800
#define CMG_S 5700
#define FACE_S 5400
class FMDP : public MDP {
 private:
  // new policyServer
  int getnumvars();
  void policyServe(DdManager *gb, DdNode *act, DdNode *val, int numovars);
 public:
  FMDP() : MDP() {}
  FMDP(char *infile, double badd=BIGADD) :
    MDP(infile,badd) {}

};
