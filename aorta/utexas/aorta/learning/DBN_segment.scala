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

class DBN_segment(sim: Simulation) {
  
  
	var state = new DBNState_segment(sim)
	val high_losers = 100
	val high_ratioau = 5
	val high_delay = 100
	private var count = 1
	var CPT_per_lane = mutable.Map[String, mutable.Map[String, Int]]()
	var CPT_total = mutable.Map[String, Int]().withDefaultValue(0)
	private var key = ""
	private var freq = 0.0
	var CPT = mutable.Map[String, Int]().withDefaultValue(0)
	
	var location = List("South", "North", "West", "East")
	// init hashmap for the CPT's of each lane
	for(i <- 0 to state.discrete_traffic.length - 1)
		CPT_per_lane.put("Lane "+location(i), CPT.clone())
	
	// put for the action
	CPT_per_lane.put("Action", CPT.clone())
	
	
	// query for stats 
	val check_stats = new Thread(new Runnable {
		def run() {
		    var previous_traffic = state.discrete_traffic.clone()
		    var previous_actions = state.actions.clone()
		    sim.listen(classOf[EV_Signal_Change], _ match{
		    	case f: EV_Signal_Change => if (f.greens.size >= 4){

				
			    print("Cars present: ")
				state.print_carsPresent()
				print("Lane traffic: ")
				state.print_traffic()
				
				//updateCPT(previous_traffic)
				updateCPT(previous_traffic, previous_actions)
				updateCPT_action(previous_traffic, previous_actions)
				previous_traffic = state.discrete_traffic.clone()
				previous_actions = state.actions.clone()
				
				
				//println(state.actions_per_lane)
				printCPT_ML(CPT_per_lane, count)
				//printCPT_Dir(CPT_per_lane, count, count * state.discrete_traffic.size, count.toFloat/state.discrete_traffic.size)
				//printCPT_total(count*state.discrete_traffic.size)
				
				count += 1
				//printCPT(CPT)
				state.reset_carsPresent()
				//Thread.sleep(1500)
			}
		    }
		)}
		//}
		
	})
	
	check_stats.start()

	
	def updateCPT_action(previous_traffic: Array[Array[String]], previous_actions: mutable.HashMap[String, Seq[String]]){
	     //println(previous_traffic.deep.toString())
		 var key2 = ""
	     key2 = "P("+state.actions.getOrElse("signals", "None").toString()+"_(t+1)|"+previous_actions.getOrElse("signals", "None").toString()+"_(t),"+previous_traffic.deep.toString()+"_(t))"
	     CPT_per_lane.get("Action").get(key2) += 1
	     CPT_total(key2) += 1
	   
	}
	
	def updateCPT(previous_traffic: Array[Array[String]], previous_actions: mutable.HashMap[String, Seq[String]] ){
		      
			for( i <- 0 to state.discrete_traffic.length -1){
			  CPT = CPT_per_lane.get("Lane "+location(i)).get
			  //println("Lane "+location(i)+":")
			  for(j <- 0 to state.discrete_traffic(i).length - 1){
				      
				      if (j != state.discrete_traffic(i).length - 1){
				          // segment j at t+1 -> segment j at t, segment j + 1 (previous segment) at t, action
				    	  key = "P("+state.discrete_traffic(i)(j)+"_seg"+j+"(t+1)|"+previous_traffic(i)(j)+"_seg"+j+"(t),"+previous_traffic(i)(j+1)+"_seg"+(j+1)+"(t),"+previous_actions.getOrElse("signals", "None").toString()+"_(t))"
				      }
				      else{
				          // segment j at t+1 -> segment j at t, action (last segment, only depends on itself)
				          key = "P("+state.discrete_traffic(i)(j)+"_seg"+j+"(t+1)|"+previous_traffic(i)(j)+"_seg"+j+"(t),"+previous_actions.getOrElse("signals", "None").toString()+"_(t))"
				      }
				      CPT(key) += 1
				      
				      CPT_total(key) += 1
				      
				      //println(key + "_lane: " + CPT(key))
				      //println(key +"_total: " + CPT_total(key))
				      
				}
			}
	}
	
	def printCPT_Dir(CPT: mutable.Map[String, mutable.Map[String, Int]], count: Int, count_t: Int, mu: Double){
		var key = ""
		var dir = 0.0
		var ml_estimate = 0.0
		//var check = 0.0
		println("Estimates with Dirichlet prior smoothing: ")
		for(i <- 0 to CPT.keySet.size - 1){
		  if (i < location.length)
			  key = "Lane "+location(i)
		  else
		    key = "Action"
		  println(key+":")
		  for (j <- CPT_total.keySet.toList.sorted){
		      try{
		        ml_estimate = CPT.get(key).get.get(j).get
		      }catch{
		        case e: Exception => {
		          ml_estimate = 0.0
		        }
		      }
		      
		      dir = (ml_estimate + (mu * CPT_total.get(j).get./(count_t))).toFloat/((count) + mu)
		      println(j+" "+dir)
		      //check += CPT.get(key).get.get(j).get./(count)
		  }
		  // simple validation to check if the probabilities sum to 1
		  //println("Valid: "+check)
		  //check = 0.0
		  println()
		}
	}
	
	def printCPT_ML(CPT: mutable.Map[String, mutable.Map[String, Int]], count: Int){
		var key = ""
		//var check = 0.0
		println("Maximum Likelihood estimates:")
		for(i <- 0 to CPT.keySet.size - 1){
		  if (i < location.length)
			  key = "Lane "+location(i)
		  else
		    key = "Action"
		  println(key+":")
		  for (j <- CPT.get(key).get.keySet.toList.sorted){
			  //println(j+" "+CPT.get(key).get.get(j).get)
		      println(j+": "+CPT.get(key).get.get(j).get.toFloat/(count))
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