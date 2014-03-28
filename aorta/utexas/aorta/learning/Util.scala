package utexas.aorta.learning

import scala.collection.mutable
import java.io.FileInputStream
import scala.util.Marshal
import java.io.FileOutputStream

class Util {
  
  val path = "previousDATA/prevdata8actions.data"
  var list_Data = new mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]
  var totalValueInData = mutable.HashMap[String, mutable.HashMap[String,mutable.LinkedHashMap[String, Int]]]()
  var config = new Config
  
  def showStat{
    val inOut = mutable.LinkedHashMap[String, mutable.LinkedHashMap[String,Int]]()
    inOut.update("In", new mutable.LinkedHashMap[String,Int])
    inOut.update("Out", new mutable.LinkedHashMap[String,Int])
    for ((key,value) <- totalValueInData){
      print("Action : "+key+"\n")
      for ((key2, value2) <- value){
        print(key2+"\t:")
        var temp = mutable.LinkedHashMap[String, Int]()
        if (key2.contains("in"))
          temp = inOut("In")
        else if (key2.contains("out"))
          temp = inOut("Out")
        for ((key3, value3) <- value2){
          temp.update(key3, temp.getOrElse(key3, 0)+value3)
          print(key3+" = "+value3+"\t")
        }
        print("\n")
      }
      print("\n")
    }
    
    println("In Out Statistics for All Actions")
    println("Outgoing Lane :")
    inOut("Out").foreach(f => print(f._1+" = "+f._2+"\t"))
    println
    println("Incoming Lane :")
    inOut("In").foreach(f => print(f._1+" = "+f._2+"\t"))
  }
  
  // Load previous ADD's and continue working
  def loadprevdata(path : String) : mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]] = {   
    val in = new FileInputStream(path)
    val bytes = Stream.continually(in.read).takeWhile(-1 !=).map(_.toByte).toArray
    val loadData =Marshal.load[mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]](bytes)
    loadData.foreach(l => {
      //compute simple statistics of the data
      val actionStat = totalValueInData.getOrElseUpdate(l._1, new mutable.HashMap[String, mutable.LinkedHashMap[String,Int]])
      for ((key,value) <- l._2){
        //val varName = actionStat.getOrElseUpdate(key, new mutable.HashMap[String, Int])
        var varName = mutable.LinkedHashMap[String,Int]()
        try {
          varName = actionStat(key)
        } catch {
          case e:Exception => {
            for (value <-config.parameters(key)){
              varName.update(value, 0)
            }
            actionStat.update(key, varName)
          }
        }
        varName.update(value, varName.getOrElse(value, 0)+1)
      }
    //--------------------------------------
    })
    //new mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]
    return loadData
  }
  
  // Save previous ADD's to disk
  def saveprevdata(path : String){
    println("Saving data")
    val out = new FileOutputStream(path)
    out.write(Marshal.dump(list_Data))
    out.close
    println("Finish saving data")
  }
  
  def mergeData(listFilename : List[String]){
    listFilename.foreach(f => {
      var temp = mutable.ListBuffer[Tuple3[String, mutable.HashMap[String, String], mutable.HashMap[String, String]]]()
      temp = loadprevdata(f)
      list_Data.++=(temp)
    })
  }
}

object Util {
  def main(args: Array[String]) {
    var util = new Util
    var listFilename = List("previousDATA/prevdata8actions_100.data", "previousDATA/prevdata8actions_100_300.data")
    util.mergeData(listFilename)
    util.saveprevdata("previousDATA/merged.data")
  }

}
