package learning

object ADDDemo {
  def main(args: Array[String]) {
    var x = new ADD(Array("A","B","C"), "Aprime", Array("True","False"))
    x.update(Array("True","True","True"),"True")
    x.update(Array("True","True","False"),"True")
    x.update(Array("True","False","True"),"True")
    x.update(Array("True","False","False"),"True")
    x.update(Array("False","True","True"),"True")
    x.update(Array("False","True","False"),"True")
    x.update(Array("False","False","True"),"True")
    x.update(Array("False","False","False"),"True")
    x.printToData
    
    //var y : String = _
    //println(y)

  }

}
