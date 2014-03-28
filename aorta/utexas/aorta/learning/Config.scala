package utexas.aorta.learning

import scala.collection.mutable

class Config {
    val location = List("South_in", "North_in", "West_in", "East_in", "South_out", "North_out", "West_out", "East_out")
    val lane_config = List("low", "medium", "high")
    val action_config = List("ArrayBuffer(5, 7, 8, 9)","ArrayBuffer(11, 4, 5, 6)","ArrayBuffer(10, 11, 12, 3)","ArrayBuffer(1, 2, 3, 7)","ArrayBuffer(1, 3, 4, 5)", "ArrayBuffer(10, 3, 5, 8)", "ArrayBuffer(11, 12, 7, 9)", "ArrayBuffer(11, 2, 6, 7)")
    val segments = 3
	var parameters = mutable.Map[String, List[String]]()
	val keep_gathering = false
	
	for (i <- 0 to location.length - 1){
	   for (j <- 0 to segments - 1){
	     parameters(location(i)+"_seg"+j) = lane_config
	   }
	}
    parameters("TrafficSignal") = action_config
    /*
    for (a <- action_config){
    	parameters(a) = action_config
    }
    * 
    */
}
