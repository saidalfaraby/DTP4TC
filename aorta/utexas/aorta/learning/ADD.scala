package utexas.aorta.learning

import scala.collection.mutable.LinkedHashMap
import java.io.FileWriter
import scala.collection.mutable
import scala.collection.immutable

abstract class GenericNode(label : String){
  var parent : Node = null
  var edge : String = null
  def getLabel = label
}

class Data(pv : immutable.HashMap[String, String], dv : String){
  var parentVal = pv
  var decisionVal = dv
}

class Node(label : String) extends GenericNode(label){
  val children = LinkedHashMap[String, GenericNode]()
  val score : Double = 0.0
  var untestedParents : Int //just make sure that the root will have nParent - 1, and -1 for 1-level deeper node
  def getChild(edgeName : String) = children(edgeName)
  def addChild(edgeName : String, child : GenericNode) = children.+=(edgeName -> child)
  override def toString = label
}

class Leaf(label : String, decisionVals : List[String]) extends GenericNode(label){
  def this() {this("", List())}
  var splitCandidates = mutable.ListBuffer[String]()
  val listData = mutable.ListBuffer[Data]()
  val counter = mutable.HashMap[String, Int]()
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
  override def toString : String = {
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
class ADD (decisionNodeName : String, decisionVals : List[String], parentsMap : Map[String, List[String]]){
  //Assume the order of internalNode is fix, and the first one is the root
  var root : GenericNode = null
  var decisionNodeVals = decisionVals
  var refCandidates = mutable.ListBuffer[Leaf]()
  def getParents = parentsMap.keys
  def getName = decisionNodeName
  val allLeaves = mutable.ListBuffer[Leaf]() //register all leaves to minimize computation of scoring function
  var MDLthreshold = 0.0
  
  var N = 0//sample size so far
  
  def init(data : mutable.ListBuffer[(Pair[immutable.HashMap[String,String],String])]){
    root = new Leaf(decisionNodeName, decisionNodeVals)
    refCandidates.+=(root.asInstanceOf[Leaf])
    root.asInstanceOf[Leaf].splitCandidates = parentsMap.keys.to[mutable.ListBuffer]
    N = data.length
    for (d<- data){
      root.asInstanceOf[Leaf].addData(new Data(d._1, d._2))
    }
  }
  
  //now we need to know all possible value of each internal node beforehand
  def addSplit(Y : String,l : Leaf):Node = {
    refCandidates.-=(l)
    val t = new Node(Y)
    l.parent.addChild(l.edge, t)
    t.untestedParents = t.parent.untestedParents - 1
    val temp = mutable.HashMap[String, Leaf]()
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
    Y.parent.addChild(Y.edge, l)
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
    var l = scala.util.Random.shuffle(refCandidates).head
    var minScore = Double.PositiveInfinity
    var bestCandidate : String = null
    for (Y <- l.splitCandidates){
      var t = addSplit(Y, l)
      if (score < minScore){
        minScore = score
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
        return 1+ Math.log(testNode.asInstanceOf[Node].untestedParents)+sum
      }
    }
    
    def DLparam():Double = {
      return 1./2*(decisionNodeVals.length-1)*refCandidates.length* N
    }
    
    def DLdata() : Double = {
      var sum = 0.0
      for (l <- allLeaves){
        val total = l.counter.valuesIterator.sum
        for (v <- l.counter.valuesIterator)
        	sum+= - (v.toFloat/N)*Math.log(v/total)
      }
      return sum
    }
    
    return DLstruct(root)+DLparam+DLdata
  }
  
  def learning(){
    var score = Double.NegativeInfinity
    var differ = Double.PositiveInfinity
    var result : Tuple3[Double, String, Leaf] = (0.0, null, null)
    while (differ > MDLthreshold){
      result = extendTree
      if (result._1 - score > MDLthreshold){
        addSplit(result._2, result._3)
      }
      differ = result._1 - score
      score = result._1
    }
  }
  
  /**
   * @param state - values/states of all internal nodes at time step n
   * @param value - value/state of decision node at time step n+1
   */
  /*
  def update(state : Array[String], value : String){
    require(state.length == internalNode.length)
    def recur(parent : Node, edgeNames : Array[String]){

      //println(count)
      val edgeLabel = edgeNames.head

      if (edgeNames.length == 1) {
        //check if the current value is already registered in decisionNodeVal or not
        if (!decisionNodeVal.contains(value))
          decisionNodeVal+= value
        try { //get leaf
          val leaf = parent.getChild(edgeLabel).asInstanceOf[Leaf]
          leaf.hit(value)
        } catch { //leaf doesn't exist, create the new one
          case e:Exception => {
            val leaf = new Leaf(decisionNodeName)//create new leaf
            leaf.hit(value)
            parent.addChild(edgeLabel,leaf)//add to parent
          }
        }
      } else { //internal node
      var childNode : Node = null
        try { // get internal node
          childNode = parent.getChild(edgeLabel).asInstanceOf[Node]
        } catch { //internal node doesn't exist, create the new one
          case e:Exception => {
            val nodeName = internalNode(internalNode.length-edgeNames.length+1)
            childNode = new Node(nodeName)
            parent.addChild(edgeLabel,childNode)//add to parent
          }
        }
        recur(childNode.asInstanceOf[Node],edgeNames.tail)
      }
    }
    recur(root, state)
  }
  
  def printTree {
    //traverse breadth first
    var queue : Array[GenericNode] = Array(root)
    var remaining = 0
    while (queue.length > 0){
      if (remaining==0) {
        println()
        remaining=queue.length
      }
      var curNode = queue(0)
      //remove proceeded node
      queue = queue.tail
      print (curNode.toString+"  ")
      //expand current Node
      if (!curNode.isInstanceOf[Leaf]){
        for (nextNode <- curNode.asInstanceOf[Node].children.valuesIterator) queue :+= nextNode
      }
      remaining-=1
    }

  }
  
  /**
   * This method will print the ADD into a dot file
   * @param filename - filename where the data will be appended to
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
        if (node.asInstanceOf[Leaf].getValues.length>0){
          //string += node.toString
          decisionNodeVal.foreach(f => {
            if (node.asInstanceOf[Leaf].valuesMap.contains(f))
              string+= f+"="+node.asInstanceOf[Leaf].valuesMap(f)+" "})
        } else {
          string += node.asInstanceOf[Leaf].getLabel
        }
        
      }
      string += ")"
    }
    traverseHelper(root)
    string += "\n"
    string
    /*
    val fw = new FileWriter(filename, true)
    try {
    	fw.write(string)
    }
    finally fw.close() 
    * 
    */
  }
  * *
  */
}





/**
 * This class is for storing all transition probabilities for all actions and all
 * features for each actions
 */
class Model{
  val actionADD : mutable.HashMap[String, mutable.ListBuffer[ADD]] = mutable.HashMap()
  
  
  def addModel(action : String, decisionNode : String, parents : List[String], params : Map[String, List[String]]){
    val act = actionADD.getOrElse(action, {val v = new mutable.ListBuffer[ADD];actionADD.put(action, v);v})
    act+= new ADD(decisionNode, params(decisionNode), params.filterKeys(parents.toSet))
  }
  
  /*
  def update(action : String, prevState : HashMap[String, String], curState : HashMap[String,String]){
    val actSome = actionADD.get(action)
    if (actSome == None)
      throw new Exception("action not found, add new action using addModel")
    else{
      val act = actSome.get
      act.iterator.foreach(add => {
        val decisionVal = curState.getOrElse(add.getName, throw new Exception(add.getName+" is not in current State Map"))
        val parentVal = Array.fill[String](add.getParents.length)("")
        for (i <- 0 until add.getParents.length){
          parentVal(i) = curState.getOrElse(add.getParents(i), throw new Exception(add.getParents(i)+" is not in the previous State Map"))
        }
        add.update(parentVal, decisionVal)
      })
    }
  }
  */
  
  def printToDotFile{
    val filename = "test.dot"
    var string = ""
    if (actionADD.keySet.size >0){
    //assume all actions have the same ADD, take the first action, find all possible
    //values for each features/variable
    val listADD = actionADD(actionADD.keySet.head)
    string += "(variables "
    for (add <- listADD){
      string += "("+add.getName+" "
      string += add.decisionNodeVals.mkString(" ")
      string += ") "
    }
    string += ")\n"
    for ((k,v)<- actionADD){
      string += "action "+k+"\n"
      for (add <- v){
        //string += add.printToString
      }
      string += "endaction\n"
    }
    }
    print(string)
    /*
    val fw = new FileWriter(filename, true)
    try {
    	fw.write(string)
    }
    finally fw.close()
    */
  }
}
