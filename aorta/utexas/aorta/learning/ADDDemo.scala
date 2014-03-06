package utexas.aorta.learning

object ADDDemo {
  def main(args: Array[String]) {
    var x = new ADD("Aprime", Array("A","B","C"))
    x.update(Array("True","True","True"),"True")
    x.update(Array("True","True","False"),"True")
    x.update(Array("True","False","True"),"True")
    x.update(Array("True","False","False"),"True")
    x.update(Array("False","True","True"),"True")
    x.update(Array("False","True","False"),"True")
    x.update(Array("False","False","True"),"True")
    x.update(Array("False","False","False"),"True")
    x.printToString
    
    val y = new SPUDD
    
    //var y : String = _
    //println(y)

  }

}
