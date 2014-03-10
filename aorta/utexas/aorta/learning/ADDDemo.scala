package utexas.aorta.learning

import scala.collection.mutable

class test(t:String){
  val tt = t
  override def toString():String={return tt}
}
object ADDDemo {
  def main(args: Array[String]) {
    /*
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
    */
	val a=new test("1")
	val b=new test("2")
	val c=new test("3")
    val l = mutable.ListBuffer[test](a,b,c)
    l.-=(b)
    print(l.toString())
    
    val s = mutable.HashMap[String, Int]("s"->1,"a"->2)
    
    print(s.toString)
    s.+=("s"->5)
    print(s.toString)
    
    //val y = new SPUDD
    
    //var y : String = _
    //println(y)

  }

}
