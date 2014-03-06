package utexas.aorta.learning

import scala.collection.mutable.LinkedHashMap
import java.io.FileWriter
import scala.collection.mutable.HashMap
import scala.collection.mutable

/**
 * Decision node is a node whose probability of its values are going to be predicted
 * @param decisionNodeName - label of decision node
 * @param decisionNodeVal - all possible values/states of decision node
 * @param internalNode - in out case these are array of parents labels
 */
class ADD (decisionNodeName : String, internalNode : Array[String]){
  //Assume the order of internalNode is fix, and the first one is the root
  var root = new Node(internalNode(0))
  var decisionNodeVal : mutable.ListBuffer[String] = mutable.ListBuffer()
  def this(name: String, root: Node) { this(name, Array("")); this.root = root; }
  def getParents = internalNode
  def getDecisionValues = decisionNodeVal
  def getName = decisionNodeName
  
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
        if (!decisionNodeVal.contains(value)){
          decisionNodeVal+= value
        }
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
        recur(childNode.asInstanceOf[Node],edgeNames.slice(1, edgeNames.length))
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
      queue = queue.slice(1,queue.length)
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
    //print(getName)
    string += getName
    def traverseHelper(node : GenericNode){
      //print("\t(")
      string += "\t("
      if (node.isInstanceOf[Node]){
        //print(node.toString)
        string += node.toString
        val children = node.asInstanceOf[Node].children
        var count=children.size
        for (edge <-children.keysIterator){
          //print("\t(")
          string += "\t("
          //print(edge)
          string += edge
          traverseHelper(children.getOrElse(edge, null))
          //print(")")
          string += ")"
          if (count>1){
            //println
            string += "\n"
          } 
          count -=1
        }
      } else if (node.isInstanceOf[Leaf]){
        if (node.asInstanceOf[Leaf].getValues.length>0){
          //print(node.toString)
          //string += node.toString
          decisionNodeVal.foreach(f => {
            if (node.asInstanceOf[Leaf].valuesMap.contains(f))
              string+= f+"="+node.asInstanceOf[Leaf].valuesMap(f)+" "})
        } else {
          //print (node.getLabel)
          string += node.getLabel
        }
        
      }
      //print(")")
      string += ")"
    }
    traverseHelper(root)
    //print("\n")
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

abstract class GenericNode(label : String){
  val parents : LinkedHashMap[String, GenericNode] = LinkedHashMap()
  def getLabel = label
  override def toString = label
  def addParent(edge : String, parent : Node) = parents+= (edge->parent)
  def getParent = parents
}

class Node(label : String) extends GenericNode(label){
  val children : LinkedHashMap[String, GenericNode] = LinkedHashMap()
  def getChild(edgeName : String) = children(edgeName)
  def addChild(edgeName : String, child : GenericNode) = children.+=(edgeName -> child)
}

class Leaf(label : String) extends GenericNode(label){
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
 * This class is for storing all transition probabilities for all actions and all
 * features for each actions
 */
class Model{
  val actionADD : HashMap[String, mutable.ListBuffer[ADD]] = HashMap()
  
  
  def addModel(action : String, decisionNode : String, parents : Array[String]){
    val act = actionADD.getOrElse(action, {val v = new mutable.ListBuffer[ADD];actionADD.put(action, v);v})
    act+= new ADD(decisionNode, parents)
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
