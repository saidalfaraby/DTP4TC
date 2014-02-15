package utexas.aorta.learning

import scala.annotation.tailrec
import scala.collection.immutable
import scala.collection.mutable
import java.io.FileWriter
import scala.io.Source
import utexas.aorta.map.{Graph, Edge, Vertex, Turn}
import utexas.aorta.sim.intersections.{Intersection, Phase, SystemWallets, Policy}
import utexas.aorta.sim.make.{Scenario, MkAgent, Factory, RouterType}
import utexas.aorta.sim.drivers.Agent
import utexas.aorta.sim.{Simulation, EV_AgentSpawned, EV_Reroute, EV_AgentQuit, EV_TurnFinished,
                         EV_IntersectionOutcome, EV_Transition, EV_Signal_Change}
import utexas.aorta.common.{Util, cfg, StateWriter, StateReader, Flags, AgentID, Publisher,
                            VertexID, EdgeID, RoadID, Timer}
import utexas.aorta.analysis.RerouteCountMonitor
import utexas.aorta.experiments.MetricInfo
import Array._
//import utexas.aorta.learning.DBNState

class DBN (sim: Simulation) {
  
	//var agents: immutable.SortedSet[Agent] = immutable.TreeSet.empty[Agent]
	var state = new DBNState(sim)
	//var CPT = ofDim[Double](2,2)
	val high_losers = 100
	val high_ratioau = 5
	val high_delay = 100
	private var count = 1
	var CPT_per_lane = new mutable.HashMap[String, mutable.HashMap[String, Double]]
	private var key = ""
	private var freq = 0.0
	var CPT = new mutable.HashMap[String, Double]()

	// init hashmap for the CPT's of each lane
	for(i <- 0 to state.discrete_traffic.length - 1)
		CPT_per_lane.put("Lane "+i, CPT.clone())
	
	// query for stats 
	val check_stats = new Thread(new Runnable {
		def run() {
		    var previous_traffic = state.discrete_traffic.clone()
			while(true){
				//println("Average_losers: " + state.avg_losers)
				//println("Ratio au: " + state.ratio_au)
				//println("Average Delay: " + state.avg_delay)
				print("Cars present: ")
				state.print_carsPresent()
				print("Lane traffic: ")
				state.print_traffic()
				updateCPT(previous_traffic)
				previous_traffic = state.discrete_traffic.clone()  
				state.reset_carsPresent()
				printCPT(CPT_per_lane, count)
				count += 1
				//printCPT(CPT)
				Thread.sleep(3000)
			}
		}
		
	})
	
	check_stats.start()
	
	
	def updateCPT(previous_traffic: Array[String]){
			for( i <- 0 to state.discrete_traffic.length -1){
				      //CPT = CPT_per_lane.getOrElseUpdate("Lane"+i, CPT.clone())
				      CPT = CPT_per_lane.get("Lane "+i).get
				      key = previous_traffic(i)+" "+state.discrete_traffic(i)
				      try{
				    		 freq = CPT.get(key).get.+(1.0)
				    		 CPT.put(key, freq)
				    	  }
				    	  catch{
				          	case e: Exception => {
				          		CPT.put(key, 1.0)
				          	}
				    	  }
				}
	}
	
	
	def printCPT(CPT: mutable.HashMap[String, mutable.HashMap[String, Double]], count: Int){
		var key = ""
		//var check = 0.0
		println()
		for(i <- 0 to CPT.keySet.size - 1){
		  key = "Lane "+i
		  for (j <- CPT.get(key).get.keySet){
		      println(key+" "+j+" "+CPT.get(key).get.get(j).get./(count))
		      //check += CPT.get(key).get.get(j).get./(count)
		  }
		  // simple validation to check if the probabilities sum to 1
		  //println("Valid: "+check)
		  //check = 0.0
		  println()
		}
	}
}