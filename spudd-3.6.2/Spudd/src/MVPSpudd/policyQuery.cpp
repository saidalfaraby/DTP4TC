#include "pspudd.hh"
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Menu.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/fl_ask.H>
#include <GL/glut.h>

#define MAXNUMPOLICIES 20

DdNode **act;
DdNode **val;
int numpolicies;
DdManager *dd;
int *query;
int window_width, window_height, field_width,choice_width;
int label_width,choice_height,menubar_height,char_width;
int max_chars;
Fl_Window *window;
Fl_Button *q_button, **ask;
Fl_Output **action_out, **value_out;
Fl_Text_Display **policyName;
Fl_Text_Buffer **policyNameBuf;
Fl_Choice **var_choice;
Fl_Menu_Bar *menu_bar;
int ***varindex, *askindex;
int *varvals;
int nvars, novars;
onum *ov;
rnum *v;
void quit_callback(Fl_Widget *w, void *data);
void var_callback(Fl_Widget *w, void *data);
void ask_callback(Fl_Widget *w, void *data);
void load_callback(Fl_Widget *w, void *data);
void getActionAndValue(Pair & dval, Pair & aval, int);
char *policyFile;
Fl_Menu_Item baritems[] = {
         {"&File", 0, 0, 0, FL_SUBMENU},
	 {"&Load",0,load_callback,(void *)1},
         {"&Quit", 0, quit_callback, (void *)2},
       {0},
{0}};

void pQuery(DdManager *ddman, DdNode *value, DdNode *action, char * pfile, onum *ovars, rnum *vars, int numvars, int numovars) { 
  int i,j;
  nvars = numvars;
  novars = numovars;
  ov = ovars;
  v = vars;
  dd = ddman;
  
  act = new DdNode*[MAXNUMPOLICIES];
  val = new DdNode*[MAXNUMPOLICIES];
  action_out = new Fl_Output*[MAXNUMPOLICIES];
  value_out = new Fl_Output*[MAXNUMPOLICIES];
  policyName= new Fl_Text_Display*[MAXNUMPOLICIES];
  policyNameBuf = new Fl_Text_Buffer*[MAXNUMPOLICIES];
  ask = new Fl_Button*[MAXNUMPOLICIES];
  askindex = new int[MAXNUMPOLICIES];
  numpolicies = 1;
  
  act[0] = action;
  Cudd_Ref(act[0]);
  val[0] = value;
  Cudd_Ref(val[0]);

  query = new int[numovars];
  choice_width = 80;
  label_width = 80;
  choice_height = 20;
  menubar_height = 30;
  max_chars = 10;
  field_width = choice_width+label_width;
  window_width = 2*(field_width)+300;
  window_height = menubar_height+max(int(ceil(novars/2.0)),5)*(choice_height+10);
  window = new Fl_Window(window_width,window_height);
  int k,R,G,B;
  R=125; G=125; B=255;
  window->color( fl_color_cube(R*FL_NUM_RED/256,G*FL_NUM_GREEN/256, B*FL_NUM_BLUE/256)); 

  var_choice = new Fl_Choice*[numovars];
  char *buf[numovars];
  char **buff[numovars];
  varindex = new int**[numovars];
  varvals = new int[numovars];
  int xpos,ypos;
  R=0; G=12; B=12;
  for (i=0; i<numovars; i++) {
    varindex[i] = new int*[ovars[i].nvals];
    varvals[i] = 0;
    buf[i] = new char[max_chars];
    buff[i] = new char*[ovars[i].nvals];
    for (j=0; j<max_chars-1; j++) 
      buf[i][j] = ov[i].name[j];
    buf[i][j] = '\0';

    xpos = i%2; ypos = i/2;
    var_choice[i] = new Fl_Choice(xpos*field_width+label_width,
				  menubar_height+ypos*(choice_height+10),
				  choice_width,choice_height,buf[i]);
    var_choice[i]->labelcolor(fl_color_cube(R*FL_NUM_RED/256,G*FL_NUM_GREEN/256, B*FL_NUM_BLUE/256)); 
    for (j=0; j<ovars[i].nvals; j++) {
      varindex[i][j] = new int[2];
      varindex[i][j][0] = i;
      varindex[i][j][1] = j;
      buff[i][j] = new char[max_chars];
      for (k=0; k<max_chars-1; k++) 
	buff[i][j][k] = ovars[i].valname[j][k];
      buff[i][j][k] = '\0';
      var_choice[i]->add(buff[i][j],"",&var_callback,varindex[i][j]);
    }
    var_choice[i]->value(0);
  }
  policyName[0] =new Fl_Text_Display(2*field_width+120,menubar_height,150,40);
  policyNameBuf[0] = new Fl_Text_Buffer(10);
  policyNameBuf[0]->append(pfile);
  policyName[0]->buffer(policyNameBuf[0]);
  
  ask[0] = new Fl_Button(2*field_width+120,menubar_height+40,150,30,"query");
  askindex[0] = 0;
  ask[0]->callback(ask_callback, askindex);

  action_out[0] = new Fl_Output(2*field_width+120,menubar_height+70,150,40,strdup("action"));
  value_out[0] = new Fl_Output(2*field_width+120,menubar_height+110,150,40,strdup("value"));

  menu_bar = new Fl_Menu_Bar(0,0,window_width,menubar_height-10);
  menu_bar->menu(baritems);

  window->end();
  window->show();
  Fl::run();
  
  // clean up
  for (i=0; i<numpolicies; i++) {
    Cudd_RecursiveDeref(dd,act[i]);
    Cudd_RecursiveDeref(dd,val[i]);
  }
  delete [] var_choice;
  for (i=0; i<numovars; i++) {
    for (j=0; j<ovars[i].nvals; j++) {
      delete [] varindex[i][j];
      delete [] buff[i][j];
    }
    delete [] buff[i];
    delete [] buf[i];
    delete [] varindex[i];
  }
  delete [] varvals;
  delete [] varindex;
}
void var_callback(Fl_Widget *w, void *data) {
  int *tmp = ((int *) data);
  //fprintf(stderr,"clicked var_callback data is %d %d\n",tmp[0],tmp[1]);
  varvals[tmp[0]] = tmp[1];
  //for (int i=0; i<novars;i++)
  //  fprintf(stderr,"variable %s has value %s\n",ov[i].name,ov[i].valname[varvals[i]]);
}
void load_callback(Fl_Widget *w, void *data) {

  //xfprintf(stderr,"load callback called\n");
  policyFile = fl_file_chooser("choose policy value file",NULL,NULL,1);
  //  fprintf(stderr,"file chosen: %s\n",policyFile);
  if (!policyFile) {
    policyFile = NULL;
    return;
  }
  
  DdNode **UPList= (DdNode **)malloc(2*(sizeof(DdNode *)));

  //load from this file into dd: act[numpolicies] and val[numpolicies]
  //load a new file into the gui...
  int numread = readDualOptimal(dd,policyFile,&UPList);
  if (numread == 0) {
    fl_message("Whoops! - that file didn't correspond to the problem you are working on");
    return;
  } else if (numpolicies+1 >= MAXNUMPOLICIES) {
    fl_message("Whoops! - too many policies");
    return;
  } 
    
  act[numpolicies] = UPList[0];
  Cudd_Ref(act[numpolicies]);
  
  val[numpolicies] = UPList[1];
  Cudd_Ref(val[numpolicies]);
  free(UPList);
  
  // resize window
  window_width += 150;
  window->resize(window->x(),window->y(),window_width,window_height);
  menu_bar->resize(menu_bar->x(),menu_bar->y(),window_width,menubar_height-10);
  
  
  // add a new action_out[numpolicies] and value_out[numpolicies] and ask[numpolicies] field
  policyName[numpolicies] =new Fl_Text_Display(2*field_width+120+numpolicies*150,menubar_height,150,40);
  policyNameBuf[numpolicies] =new Fl_Text_Buffer(10);
  policyNameBuf[numpolicies]->append(policyFile);
  policyName[numpolicies]->buffer(policyNameBuf[numpolicies]);
  
  ask[numpolicies] = new Fl_Button(2*field_width+120+numpolicies*150,menubar_height+40,150,30,"query");
  askindex[numpolicies] = numpolicies;
  ask[numpolicies]->callback(ask_callback, askindex+numpolicies);
  window->add(ask[numpolicies]);
  // position of this? 
  action_out[numpolicies] = new Fl_Output(2*field_width+120+numpolicies*150,menubar_height+70,150,40);
  value_out[numpolicies] = new Fl_Output(2*field_width+120+numpolicies*150,menubar_height+110,150,40);

  window->add(action_out[numpolicies]);
  window->add(value_out[numpolicies]);
  window->add(policyName[numpolicies]);
  // increment numpolicies
  numpolicies++;


  window->redraw();
}

void ask_callback(Fl_Widget *w, void *data) {
  Pair dval, aval;
  int polnum = *((int *) data);
  getActionAndValue(dval,aval,polnum);
  // put these two in the output fields
  char namestr[1024];
  strcpy(namestr,"");
  aconvert(namestr,actionnames,aval.get_min(),"/");
  action_out[polnum]->value(namestr);
  value_out[polnum]->value(dval.toString());
}
void quit_callback(Fl_Widget *w, void *data) {
  window->hide();
}
void getActionAndValue(Pair & dval, Pair & aval, int whichPolicy) {
  // query the whichPolicy policy action and value with the current values of varvals
  // fprintf(stderr,"querying policy and action diagrams %d\n",whichPolicy);
  // build an array varass of 1s and 0s over the variables v
  // with assignments corresponding to the varvals (values for the original variables  ov)
  // then call *(Cudd_V(Cudd_Eval(dd,act, varass)))
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

  dval = *(Cudd_V(Cudd_Eval(dd,val[whichPolicy],varass)));
  aval = *(Cudd_V(Cudd_Eval(dd,act[whichPolicy],varass)));
  
  delete [] varass;
}
