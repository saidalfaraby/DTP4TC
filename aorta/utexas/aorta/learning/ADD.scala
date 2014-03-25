package utexas.aorta.learning

import scala.collection.mutable.LinkedHashMap
import java.io.FileWriter
import scala.collection.mutable
import scala.collection.immutable

object GenericNode{
  var ID = 0
  def generateID :Int = {
    ID += 1
    return ID
  }
}
abstract class GenericNode(label : String){
  var parent : Node = null
  var edge : String = null
  def getLabel = label
  var ID = GenericNode.generateID
  override def toString = getLabel
}

class Data(pv : immutable.Map[String, String], dv : String){
  var parentVal = pv
  var decisionVal = dv
}

class Node(label : String) extends GenericNode(label){
  var children = LinkedHashMap[String, GenericNode]()
  val score : Double = 0.0
  var untestedParents = 0 //just make sure that the root will have nParent - 1, and -1 for 1-level deeper node
  def getChild(edgeName : String) = children(edgeName)
  def addChild(edgeName : String, child : GenericNode) = children.+=(edgeName -> child)
}

class Leaf(label : String, decisionVals : List[String]) extends GenericNode(label){
  def this() {this("", List())}
  var splitCandidates = mutable.ListBuffer[String]()
  val listData = mutable.ListBuffer[Data]()
  val counter = mutable.LinkedHashMap[String, Int]()
  for (k <- decisionVals){
    counter.update(k, 0)
  }
  def addData(data : Data){
    listData.+=(data)
    counter.update(data.decisionVal, counter(data.decisionVal)+1)
  }
  def removeData(data : Data){
    listData.-=(data)
    counter.update(data.decisionVal, counter(data.decisionVal)-1)
  }
  def getValuesString : String = {
    var s = ""
    for (v <- counter.valuesIterator){s+=(v+" ")}
    s
  }
}

/**
 * Decision node is a node whose probability of its values are going to be predicted
 * @param decisionNodeName - label of decision node
 * @param decisionVals - all possible values/states of decision node
 * @param parentsMap - contain parents labels as key and all possible value of each parent as value
 */
class ADD (decisionNodeName : String, decisionVals : List[String], parentsMap : Map[String, List[String]], typeTree : String){
  //Assume the order of internalNode is fix, and the first one is the root
  var root : GenericNode = null
  var decisionNodeVals = decisionVals
  var refCandidates = mutable.ListBuffer[Leaf]()
  def getParents = parentsMap.keys
  def getName = decisionNodeName
  val allLeaves = mutable.ListBuffer[Leaf]() //register all leaves to minimize computation of scoring function
  var MDLthreshold = 0.0
  
  var typeOfTree = typeTree  //"COMPACT" //"FULL"
 
  var N = 0//sample size so far
  
  init()
  
  def getScore = score
  
  def init(){
    root = new Leaf(decisionNodeName, decisionNodeVals)
    refCandidates.+=(root.asInstanceOf[Leaf])
    root.asInstanceOf[Leaf].splitCandidates = parentsMap.keys.to[mutable.ListBuffer]
  }
  
  def addData(data : Pair[immutable.Map[String,String],String]){
    N = N+1
    root.asInstanceOf[Leaf].addData(new Data(data._1, data._2))
  }
  
  //now we need to know all possible value of each internal node beforehand
  def addSplit(Y : String,l : Leaf):Node = {
    refCandidates.-=(l)
    val t = new Node(Y)
    try{
    	l.parent.addChild(l.edge, t)
    	t.parent = l.parent
    	t.edge = l.edge
    	t.untestedParents = t.parent.untestedParents - 1
    }catch {
      case e: Exception  => {
        root = t 
        t.untestedParents = parentsMap.keySet.size - 1
      }
    }
    
    val temp = mutable.HashMap[String, Leaf]()
//    println("Y ins: "+Y)
    for (v <- parentsMap(Y)){
      val newL = new Leaf(decisionNodeName, decisionNodeVals)
      temp.update(v, newL)
      newL.parent = t
      newL.edge = v
      t.addChild(v, newL)
      newL.splitCandidates = l.splitCandidates.diff(List(Y))
      //add each new leaf to allLeaves list
      allLeaves+=newL
      //if this is not the last parent available for splitting, then add this leaf to split candidates
      if (newL.splitCandidates.length > 0)      
        refCandidates.+=(newL)
    }
    //split the data from old leaf into the new leaves
    for (d <- l.listData){
      temp(d.parentVal(Y)).addData(d)
    }
    return t
  }
  
  //switch a parent of leaves with the previous leaf at the same position. Here we don't copy the data into the leaf
  //because it still contain the old data we need
  def switchNode(Y : Node, l : Leaf){
    try
    	Y.parent.addChild(Y.edge, l)
    catch{
      case e : Exception =>   root = l
    }
    Y.parent = null
    Y.edge = null
    for (child <- Y.children.valuesIterator){
      refCandidates.-=(child.asInstanceOf[Leaf])
      allLeaves.-=(child.asInstanceOf[Leaf])
    }
    Y.children = LinkedHashMap[String, GenericNode]()
    refCandidates.+=(l)
  }
  
  //remove a parent of leaves. Then we need to relocate the data.
  def removeSplit(Y : Node){
    val l = new Leaf(decisionNodeName, decisionNodeVals)
    allLeaves.+=(l)
    l.splitCandidates+=(Y.getLabel)
    Y.parent.addChild(Y.edge, l)
    refCandidates.+=(l)
    Y.parent = null
    Y.edge = null
    //combine the data below the leaves of Y
    for ( leaf <- Y.children.valuesIterator){
      allLeaves.-=(leaf.asInstanceOf[Leaf])
      for (d <- leaf.asInstanceOf[Leaf].listData)
        l.addData(d)
    }
  }
  
  def extendTree : Tuple3[Double, String, Leaf] ={
    //select a node from refinement candidates with a probability (uniform)
    var l : Leaf = null
    try{
      if (typeOfTree == "COMPACT")
    	l = scala.util.Random.shuffle(refCandidates).head
    	else if(typeOfTree == "FULL")
    	  l = refCandidates.head
    }catch{
      // ideally we already built the tree
      case e: Exception => {println("no body in refCandidates");return (Double.PositiveInfinity, null, null)}
    }
    var minScore = Double.PositiveInfinity
    var bestCandidate : String = null
    for (Y <- l.splitCandidates){
      var t = addSplit(Y, l)
      val score_ = score
      if (score_ < minScore){
        minScore = score_
        bestCandidate = Y
      }
      switchNode(t,l)
    }
    return (minScore, bestCandidate, l)
  }
  
  def score : Double ={
    def DLstruct(testNode : GenericNode) : Double ={
      if (testNode.isInstanceOf[Leaf])
        return 1.0
      else {
        var sum = 0.0
        for (v <- testNode.asInstanceOf[Node].children.valuesIterator){
          sum += DLstruct(v)
        }
        return 1+ (Math.log(testNode.asInstanceOf[Node].untestedParents + 1)/Math.log(2))+sum
      }
    }
    
    def DLparam():Double = {
      return 1./2*(decisionNodeVals.length-1)*allLeaves.length* (Math.log(N)/Math.log(2))
    }
    
    def DLdata() : Double = {
      var sum = 0.0
      println("N leaves : "+allLeaves.length)
      for (l <- allLeaves){
        val total = l.counter.valuesIterator.sum
        val nVal = l.counter.size
        for (v <- l.counter.valuesIterator){
        	val tmp = - ((v.toFloat+1)/N)*(Math.log((v.toFloat+1)/(total + nVal))/Math.log(2))
        	sum+= tmp
        }
      }
      return sum*N
        
    }
    
    val first_ = DLstruct(root)
    val second_ = DLparam
    val third_ = DLdata
    println("DLStr:" +first_ + " DLparam:" + second_ + " DLdata:"+third_)
    val total =first_ + second_ + third_ 
    if (total !=Double.PositiveInfinity && total !=Double.NegativeInfinity && !total.isNaN())
      return total
      else
        return 999999 //just to make sure if total =  NaN then return big number. But this shouldn't happen anymore. Will check later.
  }
  
  def learning(){
    var curr_score = Double.PositiveInfinity
    var differ = Double.PositiveInfinity
    var result : Tuple3[Double, String, Leaf] = (0.0, null, null)
    
    while (refCandidates.length > 0) {
      println("previously: " + curr_score)
      println("refCandidates " + refCandidates.toString)
      result = extendTree
      differ = curr_score - result._1
      var testCond = false
      if (typeOfTree == "COMPACT")
        testCond = curr_score - result._1 > MDLthreshold
        else
          testCond = true
        //      if (curr_score - result._1 > MDLthreshold){
      if (testCond){
        addSplit(result._2, result._3)
        curr_score = result._1
      } else 
        refCandidates.-=(result._3)
      println("after: " + curr_score)
      println("differ : " + differ)
    } 
  }
  
  /**
   * @param state - values/states of all internal nodes at time step n
   * @param value - value/state of decision node at time step n+1
   */
  def printTree(fileName : String) {
    //traverse breadth first
    var queue : Array[GenericNode] = Array(root)
    var dotRelation = ""
    var dotShape = ""
    var dotString = ""
    while (queue.length > 0){
      var curNode = queue(0)
      try{
    	  dotRelation+= "\""+curNode.parent.ID+"\" -> \""+curNode.ID+"\" [label = \""+curNode.edge+"\"];\n"
      } catch {
        case e : Exception => //println("root : "+root.getLabel+" curNode : "+curNode.getLabel)
      }
      
      //remove proceeded node
      queue = queue.tail
      //expand current Node
      var scoreVal = ""
      if (curNode==root){
        scoreVal += "score="+getScore+"\\nN="+N+"\\n"
      }
      if (!curNode.isInstanceOf[Leaf]){
        dotShape +="{ rank = same; node [shape=ellipse, style=filled, color=cornflowerblue];\""+curNode.ID+"\" [label=\""+scoreVal+curNode.getLabel+"\"];}\n"
        for (nextNode <- curNode.asInstanceOf[Node].children.valuesIterator) queue :+= nextNode
      } else {
        dotShape +="{ rank = same; node [shape=box, style=filled, color=goldenrod];\""+curNode.ID+"\" [label=\""+scoreVal+curNode.getLabel+"\\n"+curNode.asInstanceOf[Leaf].getValuesString+"\"];}\n"
      }
      //remaining-=1
    }
    println(dotShape)
    println(dotRelation)
    dotString = "digraph \"DD\" {\nsize = \"7.5,10\"\nratio=1.0;\ncenter = true;\nedge [dir = none];\n"
    dotString += dotShape+dotRelation+"\n}"
    val fw = new FileWriter(fileName)
    try {
    	fw.write(dotString)
    }
    finally fw.close()


  }
  
  
  /**
   * This method return a string of the ADD as a dat file format
   */
  def printToString() : String = {
    //traverse depth first
    var string = ""
    string += getName
    def traverseHelper(node : GenericNode){
      string += "\t("
      if (node.isInstanceOf[Node]){
        string += node.toString
        val children = node.asInstanceOf[Node].children
        var count=children.size
        for (edge <-children.keysIterator){
          string += "\t("
          string += edge
          traverseHelper(children.getOrElse(edge, null))
          string += ")"
          if (count>1){
            string += "\n"
          } 
          count -=1
        }
      } else if (node.isInstanceOf[Leaf]){
        if (node.asInstanceOf[Leaf].counter.size>0)
          string += node.asInstanceOf[Leaf].getValuesString
        else 
          string += node.asInstanceOf[Leaf].getLabel
      }
      string += ")"
    }
    traverseHelper(root)
    string += "\n"
    string
  }
}


/**
 * This class is for storing all transition probabilities for all actions and all
 * features for each actions
 */
class Model{
  
  import scala.util.Marshal
  import java.io._

  var actionADD : mutable.HashMap[String, mutable.ListBuffer[ADD]] = mutable.HashMap()
  val tolerance = 0.1
  val discount = 0.9
  
  var total_N = 0
  val N_thres = 1000
  private var cnt = 0
  private var flag = 1
  var typeOfTree =  "COMPACT"//"FULL" //
  val path = "previousDATA/prevdata.data"
  var list_Data = new mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]
  
  def addModel(action : String, decisionNode : String, parents : List[String], params : Map[String, List[String]]){
    val act = actionADD.getOrElse(action, {val v = new mutable.ListBuffer[ADD];actionADD.put(action, v);v})
    act+= new ADD(decisionNode, params(decisionNode), params.filterKeys(parents.toSet), typeOfTree)
  }
  
  def gather_data_per_ADD(action : String, prevState : mutable.HashMap[String, String], curState : mutable.HashMap[String, String], keep_gathering : Boolean){
	 // keep gathering data only for this episode
    // save them to disk after it finishes
       println("N Data : "+total_N)
       // if it is an episode that we want to build the trees
       //val action_to_check = "ArrayBuffer(1, 3, 4, 5)"
        //var cnt2 = 0
       if (total_N == N_thres && !keep_gathering){
         if(flag == 1){
//           list_Data.foreach(t => {
//        	   if (t._1 == action_to_check)
//        	     cnt2 += 1
//           })
           
           list_Data.foreach(t => {
        	var ADD = actionADD.get(t._1).get
	  		ADD.foreach(smt => {
	  			var filter_map = t._2.filterKeys(smt.getParents.toSet)
	  			var data_add_pair = Pair(filter_map.toMap, t._3.get(smt.getName).get)
	  			//println("data add pair : "+data_add_pair)
	  			smt.addData(data_add_pair)
	  		})
           })
             
        	 start_building()
        	 flag = 0
        	 cnt += 1
        	 println("Built: " + cnt)
        	 val chAction = actionADD.keys.head
        	 val chADD = actionADD(chAction)(16)
        	 println("Chosen Action : "+chAction.toString+" Chosen ADD : "+chADD.getName)
        	 printToDatFile
        	 //println(action_to_check+": "+ cnt2)
//        	 chADD.printTree(chAction+"__"+chADD.getName+".dot")
          }
	  	  
	  	  //total_N = 0
	  	}else if (total_N == N_thres){
	  	  println("Finished gathering data. Saving to disk...")
	  	  println("List size:" + list_Data.size)
	  	  //println(list_Data)
	  	   saveprevdata(path)
	  	   println("Done.")
	  	}
	  	else{
	  	  // load beforehand from disk the previous counts
	  	  if (total_N == 0){
	  	    println("Loading data gathered previously...")
	  	    println(list_Data.size)
	  	    try{
	  	    	list_Data = loadprevdata(path)
	  	    }catch{
	  	      case e : Exception =>{
	  	        println("First time gathering data. Using new hashmap...")
	  	        
	  	      }
	  	    }
	  	    println("Done.")
	  	    println("Loaded list size:" + list_Data.size)
	  	    //println(list_Data)
	  	    
	  	  }
	  	  if (flag ==1){
	  	    list_Data.append(Tuple3(action, prevState, curState))
//	  		var ADD = actionADD.get(action).get
//	  		ADD.foreach(smt => {
//	  			var filter_map = prevState.filterKeys(smt.getParents.toSet)
//	  			var data_add_pair = Pair(filter_map.toMap, curState.get(smt.getName).get)
//	  			//println("data add pair : "+data_add_pair)
//	  			smt.addData(data_add_pair)
//	  		})
	  	  }
	  		
	  	}
	  total_N += 1
	  	
	  	
  }
  
  // Save previous ADD's to disk
  def saveprevdata(path : String){
    val out = new FileOutputStream(path)
    out.write(Marshal.dump(list_Data))
    out.close
    
  }
  
  // Load previous ADD's and continue working
  def loadprevdata(path : String) : mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]] = {   
    val in = new FileInputStream(path)
    val bytes = Stream.continually(in.read).takeWhile(-1 !=).map(_.toByte).toArray
    Marshal.load[mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]](bytes)
    //new mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]
  }
  
  def start_building(){
    val total = actionADD.keys.size
    var i=0
    
    
    
    for (key <- actionADD.keys){
      println(i+" of "+total)
      i+=1
       actionADD(key).foreach(add => {
         println(key+" " + add.getName)
         if (add.getName=="TrafficSignal" && typeOfTree == "FULL"){
           //do nothing
         } else {
           add.learning
           add.printTree(key+"__"+add.getName+".dot")
         }
       })
    }
  } 
  
  
  def printToDatFile{
    val filename = "test.dat"
    var string = ""
    if (actionADD.keySet.size >0){
    //assume all actions have the same ADD, take the first action, find all possible
    //values for each features/variable
    val listADD = actionADD(actionADD.keySet.head)
    val actionMap = mutable.HashMap[String, String]()
    var count = 0
    actionADD.keysIterator.foreach(f=> {count = count+1;actionMap.+=(f-> ("action_"+count))})
    string += "(variables "
    for (add <- listADD){
      string += "("+add.getName+" "
      //string += add.decisionNodeVals.mkString(" ")
      for (a <- add.decisionNodeVals){
        if (actionMap.contains(a))
          string += actionMap(a)+" "
          else
            string += a+" "
      }
      string += ") "
    }
    string += ")\n"
    string += "unnormalised\n"
    for ((k,v)<- actionADD){
      string += "action "+actionMap(k)+"\n"
      for (add <- v){
        string += add.printToString
      }
      string += "endaction\n"
    }
    string += """reward [+ (South_in_seg0 (low	(5.0))
				(medium	(-1.0))
				(high (-10.0)))
			(South_out_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))
			(North_in_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))
			(North_out_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))
			(East_in_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))
			(East_out_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))
			(West_in_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))
			(West_out_seg0 (low (5.0))
				(medium (-1.0))
				(high (-10.0)))]
    		
    		"""

    string += "discount "+discount+"\n"
    string += "tolerance "+tolerance+"\n"
    }
    print(string)
    
    val fw = new FileWriter(filename)
    try {
    	fw.write(string)
    }
    finally fw.close()
    
  }
}
