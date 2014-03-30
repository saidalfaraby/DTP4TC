#include "pspudd.hh"

// new version of dumpdot to handle multi-variate variables
// 
DdNode **branches;
char **branchnodenames;
int numbranches, numexpbranches;
int DumpDot(DdManager *dd, DdNode *add, rnum *vars, onum *orig_vars,   char ** lnames, FILE *fp) {
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
  FILE *nodes_fp = fopen("/tmp/spudd.dat","w");
  DumpDoth(gbm,add,orig_vars,0,lnames,fp,nodes_fp,&howmanynodes);
  fclose(nodes_fp);
  //now transfer all of nodes_fp over to fp
  nodes_fp = fopen("/tmp/spudd.dat","r");
  char tmp;
  while (!feof(nodes_fp)) {
    tmp = fgetc(nodes_fp);
    //?? why??
    if (!feof(nodes_fp))
      fputc(tmp,fp);
  }
  fprintf(fp,"}\n");
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

char* DumpDoth(DdManager *dd, DdNode *add, onum *orig_vars, int rovar, char ** lnames, FILE *fp, FILE *nfp, int *hmn) {
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
    //write the node to nfp
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
      // recursive call - writes this branch to nfp
      nodename = DumpDoth(dd,temp,orig_vars,rovar+1,lnames,fp,nfp,hmn);
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
    fprintf(nfp,"\"%s\" -> \"%s\" [label = \"%s\"];\n",newnodename,nodename,root_ovar.valname[coval]);
    coval++;
    // start back at top
    Cudd_RecursiveDeref(gbm,temp);
    temp = add;
    Cudd_Ref(temp);
  }
  return newnodename;
  delete [] phase;
}
