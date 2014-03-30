#include "MDP.h"

// new version of dumpdot to handle multi-variate variables
// does not handle within-time-step dependencies!!! - or primed variables!
// see dumpdotp for that!
DdNode **branches;
char **branchnodenames;
int numbranches, numexpbranches;
int DumpDot(DdManager *dd, DdNode *add, rnum *vars, int numvars, onum *orig_vars, int numorigvars, char ** lnames, FILE *fp) {
  int i;
  // first, we must apply an ordering in which all the orig_vars's 
  //corresponding binary variables are  together
  // so we just re-apply the original ordering
  // we'll also save the initial ordering so we can re-apply it after the dump.
  int *list = new int[2*numvars];
  int *orig_list = new int[2*numvars];
  for (i=0; i<numvars*2; i++) {
    // this is the original ordering
    list[i] = i;
    orig_list[i] = Cudd_ReadInvPerm(dd,i);
  }
  /*
  fprintf(stderr,"ordering passed into dumpdot is ");
  for (i=0; i<numvars*2; i++) 
    fprintf(stderr,"%d ",orig_list[i]);
  fprintf(stderr,"\nnew ordering is                 ");
  for (i=0; i<numvars*2; i++) 
    fprintf(stderr,"%d ",list[i]);
  fprintf(stderr,"\n");
  */
  int res = Cudd_ShuffleHeap(dd,list);
  

  // write out the header and global attributes
  fprintf(fp,"digraph \"DD\" {\n");
  fprintf(fp,"size = \"7.5,10\"\nratio=1.0;\ncenter = true;\nedge [dir = none];\n");

   // write the links to a different file name (temporary)
  // and the nodes to the original file
  int howmanynodes(0);
  int dagCount = Cudd_DagSize(add);
  int dagLeafCount = Cudd_CountLeaves(add);

  numbranches  = 0;

  numexpbranches = dagCount;
  branches = new DdNode*[numexpbranches];
  branchnodenames = new char*[numexpbranches];
  FILE *links_fp = fopen("/tmp/spudd.dat","w");
  DumpDoth(gbm,add,vars,numvars,orig_vars,numorigvars, 0,lnames,fp,links_fp,&howmanynodes);
  fclose(links_fp);
  //now transfer all of links_fp over to fp
  links_fp = fopen("/tmp/spudd.dat","r");
  char tmp;
  while (!feof(links_fp)) {
    tmp = fgetc(links_fp);
    //?? why??
    if (!feof(links_fp))
      fputc(tmp,fp);
  }
  fprintf(fp,"}\n");
  fclose(links_fp);
  for (i=0; i<numbranches; i++)
    Cudd_RecursiveDeref(gbm,branches[i]);

  // reapply the original ordering
  res = Cudd_ShuffleHeap(dd,orig_list);

  delete [] branches;
  delete [] list;
  delete [] orig_list;
  return(1);
}
void addOne(int *phase, int n) {
  int i(0);
  int carry(1);
  while (i < n && carry == 1) {
    phase[i] = (phase[i]+carry)%2;
    carry = (phase[i]+carry)/2;
    i++;
  }  
  return;
}

// this version only works with unprimed variables in the ADD
char* DumpDoth(DdManager *dd, DdNode *add, rnum *vars, int numvars, onum *orig_vars, int numorigvars,  int rovar, char ** lnames, FILE *fp, FILE *lfp, int *hmn) {
  // descend the add through the root_ovar's binary variables and recursively
  // write out the nodes for all sub-adds
  DdNode *temp, *branch, *newadd;
  char *nodename, *newnodename;
  char namestr[1024];
  // make up a new name for this node
  // and write it to fp
  nodename = new char[256];
  newnodename = new char[256];
  sprintf(newnodename,"a%d",*hmn);
  (*hmn)++;


  // we're at a leaf
  if (Cudd_IsConstant(add)) {
    //write the node to fp
    fprintf(fp,"{ rank = same; node [shape=box, style=filled, color=goldenrod];\"%s\" ",newnodename);

    if (lnames == NULL) {
      fprintf(fp," [label = \"%s \"];}\n",(*Cudd_V(add)).toString());
    } else {
      strcpy(namestr,"");
      aconvert(namestr,lnames,(*Cudd_V(add)).get_min());
      fprintf(fp,"[label = \"%s \"];}\n",namestr);
    }
    return newnodename;
  } else {
    // move down to the rovar which is at the root of this ADD
    while (Cudd_NodeReadIndex(vars[orig_vars[rovar].var1index].add_var) < Cudd_NodeReadIndex(add) &&
	   Cudd_NodeReadIndex(vars[orig_vars[rovar].var1index+orig_vars[rovar].nbvars-1].add_var) < Cudd_NodeReadIndex(add))
      rovar++;
    
    //

    fprintf(fp,"{ rank = same; node [shape=ellipse, style=filled, color=cornflowerblue];\"%s\" [label=\"%s\"];}\n",
	    newnodename,orig_vars[rovar].name);
  }
  
  onum root_ovar = orig_vars[rovar];
  int next_rovar_index;
  onum next_rovar;
  if (rovar+1 < numorigvars) {
    next_rovar_index = orig_vars[rovar+1].var1index;
  } else {
    next_rovar_index = -1;
  }
  next_rovar_index *= 2;
  int nbv(root_ovar.nbvars),tmp;
  int nbvals = int(pow(2.0,nbv));
  int *phase = new int[nbv];

  DdNode **arrayofvars = new DdNode*[nbv];

  DdNode *subadds[MAXVALS];
  temp = add;
  Cudd_Ref(temp);
  
  int coval =0,i;
  // go over all the branches at this node. 
  // we want to figure out the add rooted at each branch of the original
  // variable and then write the parent pointing to that branch. 
  while (coval<root_ovar.nvals) {
    
    for (i=0; i<nbv; i++)
      phase[i] = 0;
    tmp = nbvals-coval-1;
    i=nbv-1;
    // get the branch in binary
    while (tmp > 0 && i >=0) {
      phase[i--] = tmp%2;
      tmp = tmp/2;
    }
    i=0;
    int shit;
    // descend the add to find that branch. However ,not all the variables will be present in the add in general.
    // therefore, we have to skip all the ones which are not there until we find one, and then skip
    // all the ones missing at the end. Do this by checking to see if the either
    // (a) the add is a constant node (then we've gone far enough)
    // (b) the index of the node is in the next original variable's domain (gone far enough)
    while (i < nbv && ( !Cudd_IsConstant(temp) && (next_rovar_index < 0 || Cudd_NodeReadIndex(temp) < next_rovar_index))) {
      if (Cudd_NodeReadIndex(temp) == 2*(orig_vars[rovar].var1index+i)+1) {
	if (phase[i]) 
	  branch = Cudd_Then(temp);
	else
	  branch = Cudd_Else(temp);
	Cudd_Ref(branch);
	Cudd_RecursiveDeref(gbm,temp);
	temp = branch;
      }
      i++;
    }
    // now check if its the same as one we've done so far
    bool different(true);
    for (i=0; i<numbranches && different; i++) 
      different = (temp != branches[i]);
    if (different) {
      // recursive call - writes this branch to lfp
      nodename = DumpDoth(dd,temp,vars,numvars,orig_vars,numorigvars,rovar+1,lnames,fp,lfp,hmn);
      branches[numbranches] = temp;
      Cudd_Ref(branches[numbranches]);
      branchnodenames[numbranches] = nodename;
      numbranches++;
      // may need more space than expected so check for that
      // this would work better with a list
      if (numbranches >= numexpbranches) {
	// need more space
	// save old branches
	DdNode ** newbranches = new DdNode *[numbranches*2];
	for (i=0; i<numbranches; i++) 
	  newbranches[i] = branches[i];

	// delete old
	delete [] branches;
	// re-allocoate
	numexpbranches = 2*numbranches;
	branches = new DdNode *[numexpbranches];
	for (i=0; i<numbranches; i++)
	  branches[i] = newbranches[i];
	delete [] newbranches;
      }
    } else {
      nodename = branchnodenames[i-1];
    }

    //now write our current node pointing to this one
    fprintf(lfp,"\"%s\" -> \"%s\" [label = \"%s\"];\n",newnodename,nodename,root_ovar.valname[coval]);
    coval++;
    // start back at top
    Cudd_RecursiveDeref(gbm,temp);
    temp = add;
    Cudd_Ref(temp);
  }
  return newnodename;
  delete [] phase;
}
//writes a constant node to fp
void writeConstantNode(FILE *fp, char *nname, DdNode *add, char **lnames) {
  char namestr[1024];

  fprintf(fp,"{ rank = same; node [shape=box, style=filled, color=goldenrod];\"%s\" ",nname);
  
  if (lnames == NULL) {
    fprintf(fp," [label = \"%s \"];}\n",(*Cudd_V(add)).toString());
  } else {
    strcpy(namestr,"");
    aconvert(namestr,lnames,(*Cudd_V(add)).get_min());
    fprintf(fp,"[label = \"%s \"];}\n",namestr);
  }
}

int DumpDot_p(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, char ** lnames, FILE *fp) {
  int i;
  // first, we must apply an ordering in which all the orig_vars's 
  //corresponding binary variables are  together
  // so we just re-apply the original ordering
  // we'll also save the initial ordering so we can re-apply it after the dump.
  int *list = new int[2*numvars];
  int *orig_list = new int[2*numvars];
  for (i=0; i<numvars*2; i++) {
    // this is the original ordering
    list[i] = i;
    orig_list[i] = Cudd_ReadInvPerm(dd,i);
  }
  int res = Cudd_ShuffleHeap(dd,list);
  

  // write out the header and global attributes
  fprintf(fp,"digraph \"DD\" {\n");
  fprintf(fp,"size = \"7.5,10\"\nratio=0.5;\ncenter = true;\nedge [dir = none];\n");

   // write the links to a different file name (temporary)
  // and the nodes to the original file
  int howmanynodes(0);
  int dagCount = Cudd_DagSize(add);
  int dagLeafCount = Cudd_CountLeaves(add);

  numbranches  = 0;

  numexpbranches = dagCount;
  branches = new DdNode*[numexpbranches];
  branchnodenames = new char*[numexpbranches];
  FILE *links_fp = fopen("/tmp/spudd.dat","w");
  DumpDoth_p(gbm,add,vars,pvars,numvars,orig_vars,numorigvars, 0,lnames,fp,links_fp,&howmanynodes);
  fclose(links_fp);
  //now transfer all of links_fp over to fp
  links_fp = fopen("/tmp/spudd.dat","r");
  char tmp;
  while (!feof(links_fp)) {
    tmp = fgetc(links_fp);
    //?? why??
    if (!feof(links_fp))
      fputc(tmp,fp);
  }
  fprintf(fp,"}\n");
  fclose(links_fp);
  for (i=0; i<numbranches; i++)
    Cudd_RecursiveDeref(gbm,branches[i]);

  // reapply the original ordering
  res = Cudd_ShuffleHeap(dd,orig_list);
  delete [] branches;
  delete [] list;
  delete [] orig_list;
  return(1);
}

// this version to work with primed variables
char* DumpDoth_p(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars,  int rovar, char ** lnames, FILE *fp, FILE *lfp, int *hmn) {
  int nindex;
  bool isprime(false);
  // descend the add through the root_ovar's binary variables and recursively
  // write out the nodes for all sub-adds
  DdNode *temp, *branch, *newadd;
  char *nodename, *newnodename;
  // make up a new name for this node
  // and write it to fp
  nodename = new char[256];
  newnodename = new char[256];
  sprintf(newnodename,"a%d",*hmn);
  (*hmn)++;


  // we're at a leaf
  if (Cudd_IsConstant(add)) {
    writeConstantNode(fp,newnodename,add, lnames);
    return newnodename;
  } else {
    // recover the original variable's index
    rovar = 0;
    nindex = Cudd_NodeReadIndex(add);
    if (nindex%2 == 0) 
      isprime = true;
    if (isprime) 
      while (rovar < numorigvars && Cudd_NodeReadIndex(pvars[orig_vars[rovar].var1index].add_var) <= Cudd_NodeReadIndex(add))
	rovar++;
    else
      while (rovar < numorigvars && Cudd_NodeReadIndex(vars[orig_vars[rovar].var1index].add_var) <= Cudd_NodeReadIndex(add))
	rovar++;
    rovar--;

    fprintf(fp,"{ rank = same; node [shape=ellipse, style=filled, color=cornflowerblue];\"%s\" [label=\"%s",newnodename,orig_vars[rovar].name);      
    // if its primed - add that
    if (isprime)
      fprintf(fp,"'");
    fprintf(fp,"\"];}\n");
  }
  
  onum root_ovar = orig_vars[rovar];
  int nbv(root_ovar.nbvars),tmp;
  int nbvals = int(pow(2.0,nbv));

  temp = add;
  Cudd_Ref(temp);
  
  int coval(0),i;
  // go over all the branches at this node. 
  // we want to figure out the add rooted at each branch of the original
  // variable and then write the parent pointing to that branch. 
  while (coval<root_ovar.nvals) {
    // build a cube for this branch and restrict the add to that cube
    if (isprime) 
      branch = restrictVal(dd,temp,pvars,numvars,orig_vars,numorigvars,rovar,coval);
    else 
      branch = restrictVal(dd,temp,vars,numvars,orig_vars,numorigvars,rovar,coval);
    Cudd_Ref(branch);
    Cudd_RecursiveDeref(gbm,temp);
    temp = branch;

    // now check if its the same as one we've done so far
    bool different(true);
    for (i=0; i<numbranches && different; i++) 
      different = (temp != branches[i]);
    if (different) {

      // recursive call - writes this branch to lfp
      nodename = DumpDoth_p(dd,temp,vars,pvars,numvars,orig_vars,numorigvars,rovar+1,lnames,fp,lfp,hmn);

      // and saves the branch just written
      branches[numbranches] = temp;
      Cudd_Ref(branches[numbranches]);
      branchnodenames[numbranches] = nodename;
      numbranches++;
      // may need more space than expected so check for that
      // this would work better with a list
      if (numbranches >= numexpbranches) {
	// need more space
	// save old branches
	DdNode ** newbranches = new DdNode *[numbranches*2];
	for (i=0; i<numbranches; i++) 
	  newbranches[i] = branches[i];

	// delete old
	delete [] branches;
	// re-allocoate
	numexpbranches = 2*numbranches;
	branches = new DdNode *[numexpbranches];
	for (i=0; i<numbranches; i++)
	  branches[i] = newbranches[i];
	delete [] newbranches;
      }
    } else {
      nodename = branchnodenames[i-1];
    }

    //now write our current node pointing to this one
    fprintf(lfp,"\"%s\" -> \"%s\" [label = \"%s\"];\n",newnodename,nodename,root_ovar.valname[coval]);
    coval++;
    // start back at top
    Cudd_RecursiveDeref(gbm,temp);
    temp = add;
    Cudd_Ref(temp);
  }
  return newnodename;
}
// checks if add is in the list already
// returns position if found or -1
int findAddInList(rnum ***addlist, int & numadds, DdNode *add)
{
  int i;
  bool different(true);
  rnum **parentaddlist = *addlist;
  for (i=0; i<numadds && different; i++) 
    different = (add != parentaddlist[i]->add_var);
  int foundat(-1);
  if (!different)
    foundat = i;
  return foundat;
}

// stores add in addlist - returns index of new entry if it 
// wasn't there, or existing index if it was
// grows the list if necessary
int storeAddInList(rnum ***addlist, int & numadds, int & numexpadds, rnum *add)
{
  int i;
  bool different(true);
  rnum **parentaddlist = *addlist;
  for (i=0; i<numadds && different; i++) 
    different = (add->add_var != parentaddlist[i]->add_var);
  if (different) {
    // add new entry
    parentaddlist[numadds] = new rnum;
    parentaddlist[numadds]->add_var = add->add_var;
    Cudd_Ref(parentaddlist[numadds]->add_var);
    parentaddlist[numadds]->name = strdup(add->name);
    numadds++;
    //fprintf(stderr,"numadds at %d\n",numadds);
    if (numadds >= numexpadds) {
      // need more space
      // save old branches
      rnum ** newaddlist = new rnum *[numadds*2];
      for (i=0; i<numadds; i++) {
	newaddlist[i] = new rnum;
	newaddlist[i]->add_var = parentaddlist[i]->add_var;
	newaddlist[i]->name = strdup(parentaddlist[i]->name);
	Cudd_Ref(newaddlist[i]->add_var);
	Cudd_RecursiveDeref(gbm,parentaddlist[i]->add_var);
	delete parentaddlist[i];
      }
      delete [] parentaddlist;
      
      // re-allocoate
      numexpadds = 2*numadds;
      parentaddlist = new rnum *[numexpadds];
      for (i=0; i<numadds; i++) {
	parentaddlist[i] = new rnum;
	parentaddlist[i]->add_var = newaddlist[i]->add_var;
	parentaddlist[i]->name = strdup(newaddlist[i]->name);
	delete newaddlist[i];
      }
      delete [] newaddlist;
    }
    return numadds;
  } else {
    return i;
  }
}
// function to print out an add in SPUDD notation - new version with dds
int printDdNode(DdManager *dd, rnum ***addlist, int & numadds, int & numexpadds, 
		DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, FILE *fp) {
  
  int final_node = printDdNode(dd,add,addlist,numadds,numexpadds,vars,pvars,numvars,orig_vars,numorigvars);
  return final_node;
}
// prints a whole array of dds
void printDdNode(DdManager *dd, DdNode **add, int numdds, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, FILE *fp)
{
  int i,j;
  int numexpadds = 10000;
  int numadds = 0;
  rnum **addlist = new rnum*[numexpadds];
  int *v = new int[numdds];
  for (i=0; i<numdds; i++) {
    v[i] = printDdNode(dd,add[i],&addlist,numadds,numexpadds,vars,pvars,numvars,orig_vars,numorigvars);
  }
  fprintf(stderr,"numadds found %d\n",numadds);
  for (i=0; i<numadds; i++) {
    // see if its one of the adds
    for (j=0; j<numdds && v[j] != i+1; j++) ;
    if (j < numdds) 
      fprintf(fp,"dd ADD%d\n",j);
    else
      fprintf(fp,"dd dd%d\n",i);
    fprintf(fp,"%s",addlist[i]->name);
    fprintf(fp,"\nenddd\n");
  }

  for (int i=0; i<numadds; i++) 
    delete addlist[i];
  delete [] addlist;
  delete [] v;

}

int printDdNode(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars, FILE *fp)
{
  int numexpadds = 10000;
  int numadds = 0;
  rnum **addlist = new rnum*[numexpadds];
  int v = printDdNode(dd,add,&addlist,numadds,numexpadds,vars,pvars,numvars,orig_vars,numorigvars);

  fprintf(stderr,"numadds found %d\n",numadds);
  for (int i=0; i<numadds; i++) {
    fprintf(fp,"dd dd%d\n",i);
    fprintf(fp,"%s",addlist[i]->name);
    fprintf(fp,"\nenddd\n");
  }
      
  for (int i=0; i<numadds; i++) 
    delete addlist[i];
  delete [] addlist;
  return v;
}

// helper function to print out an add in original scheme tree 
int printDdNode(DdManager *dd, DdNode *add, rnum ***addlist, int & numadds, int & numexpadds, 
		rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars) {
  int nindex;
  bool isprime(false);
  // descend the add through the root_ovar's binary variables and recursively
  // write out the nodes for all sub-adds
  DdNode *temp, *branch, *newadd;
  int rovar;
  int whereisit;
  // we're at a leaf
  if (Cudd_IsConstant(add)) {
    rnum newleaf;
    newleaf.name = new char[256];
    newleaf.add_var = add;
    Cudd_Ref(newleaf.add_var);
    sprintf(newleaf.name,"(%f)",(*Cudd_V(add)).get_min());
    whereisit = storeAddInList(addlist,numadds,numexpadds,&newleaf);
    Cudd_RecursiveDeref(gbm,newleaf.add_var);
    delete [] newleaf.name;
    return whereisit;
  } else {
    // first check if this ddNode already exists in the list
    // if it does - then return where
    if ((whereisit = findAddInList(addlist,numadds,add)) >=0) {
      //fprintf(stderr,"found the %d dd \n",whereisit);
      return whereisit;
    }
    
    // recover the original variable's index
    rovar = 0;
    nindex = Cudd_NodeReadIndex(add);
    if (nindex%2 == 0) 
      isprime = true;
    if (isprime) 
      while (rovar < numorigvars && Cudd_NodeReadIndex(pvars[orig_vars[rovar].var1index].add_var) <= Cudd_NodeReadIndex(add))
	rovar++;
    else
      while (rovar < numorigvars && Cudd_NodeReadIndex(vars[orig_vars[rovar].var1index].add_var) <= Cudd_NodeReadIndex(add))
	rovar++;
    rovar--;
  }
  
  onum root_ovar = orig_vars[rovar];
  int nbv(root_ovar.nbvars),tmp;
  int nbvals = int(pow(2.0,nbv));

  temp = add;
  Cudd_Ref(temp);
  
  int coval(0),i;

  // get all branches
  int *branches = new int[root_ovar.nvals];
  bool same(true);
  int currentval(0);
  while (coval<root_ovar.nvals) {
    // build a cube for this branch and restrict the add to that cube
    // don't even really have to do this - before checking that the branch doesn't
    // already exist
    if (isprime) 
      branch = restrictVal(dd,add,pvars,numvars,orig_vars,numorigvars,rovar,coval);
    else 
      branch = restrictVal(dd,add,vars,numvars,orig_vars,numorigvars,rovar,coval);
    Cudd_Ref(branch);
    // recursivecall
    // definitely don't do this here! First check if the branch exists already...
    branches[coval] = printDdNode(dd, branch, addlist, numadds, numexpadds, vars, pvars, numvars, orig_vars, numorigvars);
    Cudd_RecursiveDeref(gbm,branch);
    if (coval == 0) {
      currentval = branches[coval];
    } else {
      same = same && (currentval == branches[coval]);
    }
    coval++;
  }
  
  // all branches are the same
  int bval;
  if (same) {
    bval = branches[coval-1];
  } else{ 
    rnum newbranch;
    newbranch.add_var = add;
    Cudd_Ref(newbranch.add_var);
    newbranch.name = new char[256];
    sprintf(newbranch.name," (%s",orig_vars[rovar].name);
    // if its primed - add that
    if (isprime) 
      sprintf(newbranch.name,"%s'",newbranch.name);
    // add the children
    for (i=0; i<root_ovar.nvals; i++) 
      sprintf(newbranch.name,"%s (%s (dd%d))",newbranch.name,orig_vars[rovar].valname[i],branches[i]-1);
    sprintf(newbranch.name,"%s)",newbranch.name);
    bval = storeAddInList(addlist,numadds,numexpadds,&newbranch);
    Cudd_RecursiveDeref(gbm,newbranch.add_var);
    delete [] newbranch.name;
  }
  delete [] branches;
  return bval;
}
/*

// function to print out an add in original scheme tree 
void oldprintDdNode(DdManager *dd, DdNode *add, rnum *vars, rnum *pvars, int numvars, onum *orig_vars, int numorigvars,  FILE *fp, char *tabstop) {
  int nindex;
  bool isprime(false);
  // descend the add through the root_ovar's binary variables and recursively
  // write out the nodes for all sub-adds
  DdNode *temp, *branch, *newadd;
  int rovar;
  char newtabstop[1024];
  // we're at a leaf
  if (Cudd_IsConstant(add)) {
    fprintf(fp,"    %s(%f)\n",tabstop,(*Cudd_V(add)).get_min());
    return;
  } else {
    // recover the original variable's index
    rovar = 0;
    nindex = Cudd_NodeReadIndex(add);
    if (nindex%2 == 0) 
      isprime = true;
    if (isprime) 
      while (rovar < numorigvars && Cudd_NodeReadIndex(pvars[orig_vars[rovar].var1index].add_var) <= Cudd_NodeReadIndex(add))
	rovar++;
    else
      while (rovar < numorigvars && Cudd_NodeReadIndex(vars[orig_vars[rovar].var1index].add_var) <= Cudd_NodeReadIndex(add))
	rovar++;
    rovar--;
    //sprintf(newtabstop,"%s",tabstop);
#ifdef WHITESPACE
    sprintf(newtabstop,"%s    ",tabstop);
    for (int i=0; i<strlen(orig_vars[rovar].name); i++)
      strcat(newtabstop," ");
#else
    sprintf(newtabstop,"%s",tabstop);
#endif
    fprintf(fp,"%s(%s",tabstop,orig_vars[rovar].name);
    

    // if its primed - add that
    if (isprime) {
      fprintf(fp,"'");
      strcat(newtabstop," ");
    }
    fprintf(fp,"\n");
  }
  
  onum root_ovar = orig_vars[rovar];
  int nbv(root_ovar.nbvars),tmp;
  int nbvals = int(pow(2.0,nbv));

  temp = add;
  Cudd_Ref(temp);
  
  int coval(0),i;
  // go over all the branches at this node. 
  // we want to figure out the add rooted at each branch of the original
  // variable and then write the parent pointing to that branch. 
  while (coval<root_ovar.nvals) {
    // build a cube for this branch and restrict the add to that cube
    if (isprime) 
      branch = restrictVal(dd,temp,pvars,numvars,orig_vars,numorigvars,rovar,coval);
    else 
      branch = restrictVal(dd,temp,vars,numvars,orig_vars,numorigvars,rovar,coval);
    Cudd_Ref(branch);
    Cudd_RecursiveDeref(gbm,temp);
    temp = branch;
    
    //    if (coval > 0)
    fprintf(fp,"%s",newtabstop);
    fprintf(fp,"(%s",root_ovar.valname[coval]);
#ifdef WHITESPACE
    fprintf(fp,"\n");
#endif

    // recursive call - writes this branch to lfp
    printDdNode(dd,temp,vars,pvars,numvars,orig_vars,numorigvars,fp,newtabstop);
    fprintf(fp,"%s)",newtabstop);
#ifdef WHITESPACE
    fprintf(fp,"\n");
#endif


    coval++;
    // start back at top
    Cudd_RecursiveDeref(gbm,temp);
    temp = add;
    Cudd_Ref(temp);
  }
  fprintf(fp,"%s)",tabstop);
#ifdef WHITESPACE
  fprintf(fp,"\n");
#endif
  
}
*/
