package learning

object ADDDemo {
  def main(args: Array[String]) {
    var x = new ADD(Array("A","B","C"), "Aprime", Array("True","False","Entah"))
    x.update(Array("True","False","True"),"True")
    x.update(Array("True","False","True"),"True")
    x.update(Array("True","False","True"),"Entah")
    x.update(Array("True","False","True"),"False")
    x.update(Array("True","False","True"),"Entah")
    x.printTree


  }

}