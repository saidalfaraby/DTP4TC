package utexas.aorta.learning

import scala.collection.mutable.LinkedHashMap
import java.io.FileWriter
import scala.collection.mutable.HashMap
import scala.collection.mutable

abstract class GenericNode(label : String){
  var parent : Node = null
  var edge : String = null
  var splitCandidates : List[String]
  def getLabel = label
}
/*
class nullLeaf extends GenericNode{
  var parent : Node = null
  var edge : String = null
  var splitCandidates : List[String]
}*/

class Node(label : String) extends GenericNode(label){
  val children = LinkedHashMap[String, GenericNode]()
  val score : Double = 0.0
  def getChild(edgeName : String) = children(edgeName)
  def addChild(edgeName : String, child : GenericNode) = children.+=(edgeName -> child)
  override def toString = label
}

class Leaf(label : String) extends GenericNode(label){
  def this() {this("")}
  val valuesMap : LinkedHashMap[String, Int] = LinkedHashMap()
  //for (v : String <- values){ valuesMap.+=(v -> 0)  }
  def hit(valueName : String) : Unit = {
    val prevV = valuesMap.get(valueName)
    if (prevV != None)
      valuesMap.update(valueName,prevV.get+1)
    else 
      valuesMap.+=(valueName -> 1)
  }
  def getValue(key : String) = valuesMap(key)
  def getValues = valuesMap.valuesIterator
  override def toString : String = {
    var s = ""
    for (v <- valuesMap.valuesIterator){s+=(v+" ")}
    s
  }
}

/**
 * Decision node is a node whose probability of its values are going to be predicted
 * @param decisionNodeName - label of decision node
 * @param decisionNodeVal - all possible values/states of decision node
 * @param internalNode - in out case these are array of parents labels
 */
class ADD (decisionNodeName : String, internalNode : Array[String], lane_config: List[String], action_config: List[String]){
  //Assume the order of internalNode is fix, and the first one is the root
  var root = new Node(internalNode(0))
  var decisionNodeVal = lane_config//mutable.ListBuffer[String]()
  var refCandidates = mutable.ListBuffer[Leaf]()
  //def this(name: String, root: Node) { this(name, Array("")); this.root = root; }
  def getParents = internalNode
  def getDecisionValues = decisionNodeVal
  def getName = decisionNodeName
  val Pattern = "(*seg*)".r
  
  val parent_values = collection.immutable.HashMap("lane" -> lane_config, "action" -> action_config)
  
  
  var N = 0//sample size so far
  
  def get_parent_values(parent: String) : List[String] = {
    if (parent contains "seg"){
      return parent_values("lane")
    }else{
      return parent_values("action")
    }
  }
  
  //now we need to know all possible value of each internal node beforehand
  val domVar = mutable.HashMap[String, mutable.ListBuffer[String]]()
  def addSplit(Y : String,l : Leaf):Node = {
    refCandidates.-=(l)
    val t = new Node(Y)
    t.splitCandidates = l.splitCandidates
    l.parent.addChild(l.edge, t)
    for (v <- domVar(l.edge)){
      val newL = new Leaf()
      newL.parent = t
      newL.edge = v
      t.addChild(v, newL)
      newL.splitCandidates = l.splitCandidates.diff(List(Y))
      refCandidates.+=(newL)
    }
    return t
  }
  
  def removeSplit(Y : Node, l : Leaf){
    Y.parent.addChild(l.edge, l)
    refCandidates.+=(l)
    Y.parent = null
    Y.edge = null
  }
  
  def extend{
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
      removeSplit(t,l)
    }
    addSplit(bestCandidate, l)
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
        return 1+ Math.log(testNode.asInstanceOf[Node].splitCandidates.length)+sum
      }
    }
    
    def DLparam():Double = {
      return 1./2*(decisionNodeVal.length-1)*refCandidates.length* N
    }
    
    return DLstruct(root)+DLparam
  }
  
  def learning(){
    
  }
  
  /**
   * @param state - values/states of all internal nodes at time step n
   * @param value - value/state of decision node at time step n+1
   */
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
}





/**
 * This class is for storing all transition probabilities for all actions and all
 * features for each actions
 */
class Model{
  val actionADD : HashMap[String, mutable.ListBuffer[ADD]] = HashMap()
  
  
  def addModel(action : String, decisionNode : String, parents : Array[String], params : mutable.Map[String, List[String]]){
    val act = actionADD.getOrElse(action, {val v = new mutable.ListBuffer[ADD];actionADD.put(action, v);v})
    act+= new ADD(decisionNode, parents, params(decisionNode), params("Action"))
  }
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
      string += add.getDecisionValues.mkString(" ")
      string += ") "
    }
    string += ")\n"
    for ((k,v)<- actionADD){
      string += "action "+k+"\n"
      for (add <- v){
        string += add.printToString
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
