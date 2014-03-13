package utexas.aorta.learning
import scala.io.Source
import scala.collection.mutable
import scala.util.matching.Regex

/**
 * This class is an interface to working with result from SPUDD
 */
class SPUDD {
  /*
  dotToTree
  //Method for reading the resulting file and convert it to
  //variables and actions name. Or we can also use policy tree structure
  
  //convert actions name to object Turn
  def dotToTree{
    val nodesMap = new mutable.HashMap[String, GenericNode]
    val filename = "SPUDD-OPTpolicy.dot"
    for (line <- Source.fromFile(filename).getLines()) {
      println(line)
      if (line.startsWith("{")){
        val pattern = "\"[^=]*\"".r
        var m = (pattern findAllIn line).toList
        if (line.contains("ellipse")){ //internal node
          nodesMap += (m(0).replaceAll("\"", "") -> new Node(m(1).replaceAll("\"", "")))
        } else if (line.contains("box")){ //leaf node
          nodesMap += (m(0).replaceAll("\"", "") -> new Leaf(m(1).replaceAll("\"", "")))
        }
      } else if (line.contains("->")){ //edge
        val pattern = "\"[^-=]*\"".r
        var m = (pattern findAllIn line).toList
        val parent = nodesMap.get(m(0).replaceAll("\"", "")).get.asInstanceOf[Node] //should always be a Node type 
        val child = nodesMap.get(m(1).replaceAll("\"", "")).get //can be Leaf or Node
        val edge = m(2).replaceAll("\"", "")
        parent.addChild(edge, child)
        child.addParent(edge, parent)
      }
    }
    var root : Node = null
    for ((k,v) <- nodesMap){
      if (v.getParent.size==0){
        root = v.asInstanceOf[Node]
      }
    }
    var add = new ADD("Policy", root)
    add.printToString
  }
  
  */
  
  
}