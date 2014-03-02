package learning

import scala.collection.mutable.LinkedHashMap

/**
 * Decision node is a node whose probability of its values are going to be predicted
 * @param internalNode - in out case these are array of parents labels
 * @param decisionNodeName - label of decision node
 * @param decisionNodeVal - all possible values/states of decision node
 */
class ADD (internalNode : Array[String], decisionNodeName : String, decisionNodeVal : Array[String]){
  //Assume the order of internalNode is fix, and the first one is the root
  val root = new Node(internalNode(0)) 
  
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
        try { //get leaf
          val leaf = parent.getChild(edgeLabel).asInstanceOf[Leaf]
          leaf.hit(value)
        } catch { //leaf doesn't exist, create the new one
          case e:Exception => {
            val leaf = new Leaf(decisionNodeName, decisionNodeVal)//create new leaf
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
  
  def printToData {
    //traverse depth first
    print(getName)
    def traverseHelper(node : GenericNode){
      print("\t(")
      if (node.isInstanceOf[Node]){
        print(node.toString)
        val children = node.asInstanceOf[Node].children
        var count=children.size
        for (edge <-children.keysIterator){
          print("\t(")
          print(edge)
          traverseHelper(children.getOrElse(edge, null))
          print(")")
          if (count>1) 
            println
          count -=1
        }
        
      } else if (node.isInstanceOf[Leaf]){
        print(node.toString)
      }
      
      print(")")
      
      
    }
    
    traverseHelper(root)
    /*
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
    * 
    */

  }
}

abstract class GenericNode(label : String){
  var parents : Array[Node] = Array()
  def getLabel = label
  override def toString = label
  def addParent(parent : Node) = parents :+= parent
}

class Node(label : String) extends GenericNode(label){
  val children : LinkedHashMap[String, GenericNode] = LinkedHashMap()
  def getChild(edgeName : String) = children(edgeName)
  def addChild(edgeName : String, child : GenericNode) = children.+=(edgeName -> child)
}

class Leaf(label : String, values : Array[String]) extends GenericNode(label){
  val valuesMap : LinkedHashMap[String, Int] = LinkedHashMap()
  for (v : String <- values){ valuesMap.+=(v -> 0)  }
  def hit(valueName : String) : Unit = {valuesMap.update(valueName,valuesMap(valueName)+1)}
  def getValue(key : String) = valuesMap(key)
  def getValues = valuesMap.valuesIterator
  override def toString : String = {
    var s = ""
    for (v <- valuesMap.valuesIterator){s+=(v+" ")}
    s
  }
}
