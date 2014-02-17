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
	var CPT_per_lane = mutable.Map[String, mutable.Map[String, Double]]()
	var CPT_total = mutable.Map[String, Double]().withDefaultValue(0.0)
	private var key = ""
	private var freq = 0.0
	var CPT = mutable.Map[String, Double]().withDefaultValue(0.0)
	
	var location = List("South", "North", "West", "East")
	// init hashmap for the CPT's of each lane
	for(i <- 0 to state.discrete_traffic.length - 1)
		CPT_per_lane.put("Lane "+location(i), CPT.clone())
	
	
	// query for stats 
	val check_stats = new Thread(new Runnable {
		def run() {
		    var previous_traffic = state.discrete_traffic.clone()
		    var previous_actions = state.actions_per_lane.clone()
			while(true){
			    
			    
			    
				//println("Average_losers: " + state.avg_losers)
				//println("Ratio au: " + state.ratio_au)
				//println("Average Delay: " + state.avg_delay)
				
			    print("Cars present: ")
				state.print_carsPresent()
				print("Lane traffic: ")
				state.print_traffic()
				
				//updateCPT(previous_traffic)
				updateCPT(previous_traffic, previous_actions)
				updateCPT_action(previous_actions)
				previous_traffic = state.discrete_traffic.clone()
				previous_actions = state.actions_per_lane.clone()
				
				
				//println(state.actions_per_lane)
				printCPT_ML(CPT_per_lane, count)
				printCPT_Dir(CPT_per_lane, count, count * state.discrete_traffic.size, count.toFloat/state.discrete_traffic.size)
				//printCPT_total(count*state.discrete_traffic.size)
				
				count += 1
				//printCPT(CPT)
				state.reset_carsPresent()
				Thread.sleep(1500)
			}
		}
		
	})
	
	check_stats.start()
	
	/*
	def updateCPT(previous_traffic: Array[String]){
			for( i <- 0 to state.discrete_traffic.length -1){
				      //CPT = CPT_per_lane.getOrElseUpdate("Lane"+i, CPT.clone())
				      CPT = CPT_per_lane.get("Lane "+location(i)).get
				      key = previous_traffic(i)+" "+state.discrete_traffic(i)
				      
				      
				      try{
				             
				    		 freq = CPT.get(key).get.+(1.0)
				    		 CPT.put(key, freq)
				    		 
				      }catch{
				         case e: Exception => {
				          		CPT.put(key, 1.0)
				          		
				          }
				      }
				      
				      // global probabilities (for smoothing)
				      try{
				             freq = CPT_total.get(key).get.+(1.0)
				    		 CPT_total.put(key, freq)
				      }catch{
				        case f: Exception => {
				          CPT_total.put(key, 1.0)
				        }
				      }
				}
	}
	
	*/

	
	def updateCPT_action(previous_actions: mutable.HashMap[String, List[Edge]]){
	   var key2 = ""
	     for( i<- 0 to state.discrete_traffic.length -1){
		       CPT = CPT_per_lane.get("Lane "+location(i)).get
		       key2 = "P("+state.actions_per_lane.getOrElse(location(i), List()).toString()+"_(t+1)|"+previous_actions.getOrElse(location(i), List()).toString()+"_(t),"+state.discrete_traffic(i)+"_(t+1))"
		       
		       CPT(key2) += 1.0
		       
		       CPT_total(key2) += 1.0
		    
		    }
	   
	}
	
	def updateCPT(previous_traffic: Array[String], previous_actions: mutable.HashMap[String, List[Edge]] ){
		      
			for( i <- 0 to state.discrete_traffic.length -1){
				      CPT = CPT_per_lane.get("Lane "+location(i)).get
				      key = "P("+state.discrete_traffic(i)+"_(t+1)|"+previous_traffic(i)+"_(t),"+previous_actions.getOrElse(location(i), List()).toString()+"_(t))"
				      
				      
				      CPT(key) += 1.0
				      
				      CPT_total(key) += 1.0

				      
				}
	}
	
	def printCPT_Dir(CPT: mutable.Map[String, mutable.Map[String, Double]], count: Int, count_t: Int, mu: Double){
		var key = ""
		var dir = 0.0
		var ml_estimate = 0.0
		//var check = 0.0
		println("Estimates with Dirichlet prior smoothing: ")
		for(i <- 0 to CPT.keySet.size - 1){
		  key = "Lane "+location(i)
		  println(key+":")
		  for (j <- CPT_total.keySet.toList.sorted){
		      try{
		        ml_estimate = CPT.get(key).get.get(j).get
		      }catch{
		        case e: Exception => {
		          ml_estimate = 0.0
		        }
		      }
		      
		      dir = (ml_estimate + (mu * CPT_total.get(j).get./(count_t)))/((count) + mu)
		      println(j+" "+dir)
		      //check += CPT.get(key).get.get(j).get./(count)
		  }
		  // simple validation to check if the probabilities sum to 1
		  //println("Valid: "+check)
		  //check = 0.0
		  println()
		}
	}
	
	def printCPT_ML(CPT: mutable.Map[String, mutable.Map[String, Double]], count: Int){
		var key = ""
		//var check = 0.0
		println("Maximum Likelihood estimates:")
		for(i <- 0 to CPT.keySet.size - 1){
		  key = "Lane "+location(i)
		  println(key+":")
		  for (j <- CPT.get(key).get.keySet.toList.sorted){
			  //println(j+" "+CPT.get(key).get.get(j).get)
		      println(j+": "+CPT.get(key).get.get(j).get./(count))
		      //check += CPT.get(key).get.get(j).get./(count)
		  }
		  // simple validation to check if the probabilities sum to 1
		  //println("Valid: "+check)
		  //check = 0.0
		  println()
		}
	}
	
	def printCPT_total(count: Int){
		//var check = 0.0
		println("Whole intersection estimates: ")
		  for (j <- CPT_total.keySet.toList.sorted){ 
		      println(j +": "+CPT_total.get(j).get./(count))
		  }
		  // simple validation to check if the probabilities sum to 1
		  //println("Valid: "+check)
		  //check = 0.0
		  println()
		}
}