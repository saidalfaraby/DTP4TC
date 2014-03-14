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
	var treeCPT = new Model()
	val high_losers = 100
	val high_ratioau = 5
	val high_delay = 100
	private var count = 1
	var CPT_per_lane = mutable.Map[String, mutable.Map[String, Int]]()
	var CPT_total = mutable.Map[String, Int]().withDefaultValue(0)
	private var key = ""
	private var freq = 0.0
	var CPT = mutable.Map[String, Int]().withDefaultValue(0)
	var isParentMapInitialized = false
	val config = new Config()
	
	

	var location = List("South_in", "North_in", "West_in", "East_in", "South_out", "North_out", "West_out", "East_out")
	// init hashmap for the CPT's of each lane
	for(i <- 0 to state.discrete_traffic.length - 1)
		CPT_per_lane.put("Lane "+location(i), CPT.clone())

	// put for the action
	CPT_per_lane.put("Action", CPT.clone())
	
	
	var parentMap = mutable.HashMap[String, mutable.ListBuffer[String]]()
	
	
	// query for stats
	val check_stats = new Thread(new Runnable {
		def run() {
		    var previous_traffic = state.discrete_traffic.map(_.clone)
		    var previous_actions = state.actions.clone()
		    var cnt = 0
		    sim.listen(classOf[EV_Signal_Change], _ match{
		    	case f: EV_Signal_Change => if (f.greens.size >= 4){

		    	if(cnt > 1){
		    		//print("Cars present: ")
		    		//state.print_carsPresent()
		    		//print("Lane traffic: ")
		    		//state.print_traffic()

		    		updateCPT(previous_traffic, previous_actions)
		    		updateCPT_action(previous_traffic, previous_actions)
		    		isParentMapInitialized = true
		    		println(parentMap.toString)
		    		//============================================
		    		//Initialize ADD for new action, or just update the ADD
		    		val actionName = previous_actions.get("signals").get.toString
		    		if (!treeCPT.actionADD.contains(actionName)){ //if the action is not initialized yet, than initialize first
		    		  //initialize new ADDs for this action
		    		  //for each variable in time t+1 effected by this action, initialize the ADD
		    		  for ((decisionNodeName,parents)<-parentMap){
		    		    treeCPT.addModel(actionName, decisionNodeName, parents.toList, config.parameters.toMap)
		    		  }
		    		  
		    		}
		    		//update the ADD
		    		val prevState = mutable.HashMap[String, String]()
		    		val curState = mutable.HashMap[String, String]()
		    		for( i <- 0 until state.discrete_traffic.length){
		    		  for(j <- 0 until state.discrete_traffic(i).length){
		    		    prevState.+=(location(i)+"_seg"+j -> previous_traffic(i)(j))
		    		    curState.+=(location(i)+"_seg"+j -> state.discrete_traffic(i)(j))
		    		  }
		    		}
		    		prevState.+=("TrafficSignal" -> previous_actions.get("signals").get.toString())
		    		curState.+=("TrafficSignal" -> state.actions.get("signals").get.toString())
		    		treeCPT.gather_data_per_ADD(actionName, prevState, curState)
		    		//treeCPT.update(actionName, prevState, curState)
		    		
		    		//============================================
		    		
		    		//println(state.discrete_traffic)
		    		//println(state.cars_present)
		    		previous_traffic = state.discrete_traffic.map(_.clone)
		    		previous_actions = state.actions.clone()


		    		//println(state.actions_per_lane)
		    		//printCPT_ML(CPT_per_lane, count)
		    		//printCPT_Dir(CPT_per_lane, count, count * state.discrete_traffic.size, count.toFloat/state.discrete_traffic.size)
		    		//printCPT_total(count*state.discrete_traffic.size)

		    		count += 1
		    		//printCPT(CPT)
		    		state.reset_carsPresent()
		    		//Thread.sleep(1500)
		    		
		    		
		    	}else{
		    	  previous_traffic = state.discrete_traffic.map(_.clone)
		    	  previous_actions = state.actions.clone()
		    	}
				cnt += 1
				//treeCPT.printToDotFile
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
	     if (!isParentMapInitialized){
	       val keys = parentMap.keySet.toList
	       val vals = mutable.ListBuffer("TrafficSignal") ++ keys
	       parentMap.+=("TrafficSignal" -> vals)
	     }
	       

	}

	def updateCPT(previous_traffic: Array[Array[String]], previous_actions: mutable.HashMap[String, Seq[String]] ){

			for( i <- 0 to state.discrete_traffic.length -1){
			  CPT = CPT_per_lane.get("Lane "+location(i)).get
			  //println("Lane "+location(i)+":")
			  for(j <- 0 to state.discrete_traffic(i).length - 1){
				  	// if we are modeling incoming lanes
				  	 if (i < state.discrete_traffic.length / 2){
				  		 if (j != state.discrete_traffic(i).length - 1){
				  			 // segment j at t+1 -> segment j at t, segment j + 1 (previous segment) at t, action
				  			 key = "P("+state.discrete_traffic(i)(j)+"_"+location(i)+"_seg"+j+"(t+1)|"+previous_traffic(i)(j) +
				  					 	"_"+location(i)+"_seg"+j+"(t),"+previous_traffic(i)(j+1)+"_"+location(i)+"_seg"+(j+1) +
				  					 		"(t),"+previous_actions.getOrElse("signals", "None").toString()+"_(t))"
				  			if (!isParentMapInitialized)
				  			  parentMap.+=(location(i)+"_seg"+j -> mutable.ListBuffer(location(i)+"_seg"+j, location(i)+"_seg"+(j+1),
				  			    "TrafficSignal"))
				  		 }
				  		 else{
				  			 // segment j at t+1 -> segment j at t, action (last segment, only depends on itself)
				  			 key = "P("+state.discrete_traffic(i)(j)+"_"+location(i)+"_seg"+j+"(t+1)|"+previous_traffic(i)(j)+
				  					 "_"+location(i)+"_seg"+j+"(t),"+previous_actions.getOrElse("signals", "None").toString()+"_(t))"
				  			if (!isParentMapInitialized)
				  				parentMap.+=(location(i)+"_seg"+j -> mutable.ListBuffer(location(i)+"_seg"+j, 
				  			    "TrafficSignal"))
				  		 }
				  	 }
				  	 // now we are modeling outgoing lanes
				  	 else{
				  	     // if we are modeling a segment that is not immediately after the intersection
				  		 if (j != 0){
				  		   key = "P("+state.discrete_traffic(i)(j)+"_"+location(i)+"_seg"+j+"(t+1)|"+previous_traffic(i)(j)+
				  				   "_"+location(i)+"_seg"+j+"(t),"+previous_traffic(i)(j-1)+"_"+location(i)+"_seg"+(j-1)+"(t),"+
				  				   		previous_actions.getOrElse("signals", "None").toString()+"_(t))"
				  			if (!isParentMapInitialized)
				  			  parentMap.+=(location(i)+"_seg"+j -> mutable.ListBuffer(location(i)+"_seg"+j, location(i)+"_seg"+(j-1), 
				  			    "TrafficSignal"))
				  		 }
				  		 // the segment immediately after the intersection, it has more dependencies
				  		 else{
				  		   var indices = Tuple3(0,0,0)
				  		   if (i == 7)
				  			 indices = Tuple3(0,1,2)
				  		   else if (i == 6)
				  		     indices = Tuple3(0,1,3)
				  		   else if (i == 5)
				  		     indices = Tuple3(0,2,3)
				  		   else if (i == 4)
				  		   	  indices = Tuple3(1,2,3)
				  		   	  
				  		   key = "P("+state.discrete_traffic(i)(j)+"_"+location(i)+"_seg"+j+"(t+1)|"+previous_traffic(i)(j)+
				  				   "_"+location(i)+"_seg"+j+"(t),"+previous_traffic(indices._1)(state.segments-1)+"_"+location(indices._1)+"_seg"+
				  				   (state.segments - 1)+"(t),"+ previous_traffic(indices._2)(state.segments-1)+"_"+location(indices._2)+"_seg"+
				  				   (state.segments - 1)+"(t)," + previous_traffic(indices._3)(state.segments-1)+"_"+location(indices._3)+"_seg"+
				  				   (state.segments - 1)+"(t)," + previous_actions.getOrElse("signals", "None").toString()+"_(t))"
				  			if (!isParentMapInitialized)
				  			  parentMap.+=(location(i)+"_seg"+j -> mutable.ListBuffer(location(i)+"_seg"+j, 
				  			    location(indices._1)+"_seg"+(state.segments - 1),location(indices._2)+"_seg"+(state.segments - 1), 
				  			    location(indices._3)+"_seg"+(state.segments - 1), "TrafficSignal"))
				  		 }
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
