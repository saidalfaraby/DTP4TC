package learning

import scala.collection.mutable.LinkedHashMap

class ADD (internalNode : Array[String], decisionNodeName : String, decisionNodeVal : Array[String]){

  val root = new Node(internalNode(0))

  def update(state : Array[String], value : String){
    require(state.length == internalNode.length)
    def recur(parent : Node, edgeNames : Array[String]){

      //println(count)
      val edgeLabel = edgeNames.head

      if (edgeNames.length == 1) {
        try { //get leaf
          val leaf = parent.getChild(edgeLabel).asInstanceOf[Leaf]
          leaf.hit(value)
        } catch {
          case e:Exception => {
            val leaf = new Leaf(decisionNodeName, decisionNodeVal)//create new leaf
            leaf.hit(value)
            parent.addChild(edgeLabel,leaf)//add to parent
          }
        }
      } else { //internal node
      var childNode : Node = null
        try {
          childNode = parent.getChild(edgeLabel).asInstanceOf[Node]
        } catch {
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
      queue = queue.slice(1,queue.length)
      print (curNode.toString+"  ")
      //expand current Node
      if (!curNode.isInstanceOf[Leaf]){
        for (nextNode <- curNode.asInstanceOf[Node].children.valuesIterator) queue :+= nextNode
      }
      remaining-=1
    }

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