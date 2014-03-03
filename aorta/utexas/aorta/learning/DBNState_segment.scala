package utexas.aorta.learning

import utexas.aorta.sim.{Simulation, EV_AgentSpawned, EV_Reroute, EV_AgentQuit, EV_TurnFinished, EV_Heartbeat,
                         EV_IntersectionOutcome, EV_IntersectionOutcomeWin, EV_IntersectionTickets, EV_Transition, EV_Signal_Change}
import utexas.aorta.sim.make.{Scenario, MkAgent, Factory, RouterType}
import utexas.aorta.map.{Turn, Edge}
import utexas.aorta.common.{Util, Timer, cfg, Flags}
import scala.collection.mutable


class DBNState_segment (sim: Simulation){
  
  // State Variables
  var avg_losers = 0.0
  var astar_count = 0
  var unrealizable_count = 0
  var ratio_au = 0.0
  var avg_delay = 0.0
  var thres_away = 80.0
  var uniq_roads =  mutable.Map[String, Double]()
  val segments = 3

  var cars_present = Array.ofDim[Int](4,3)

  var discrete_traffic = Array.ofDim[String](4,3)
  
  private var from = " "
  
  var actions = new mutable.HashMap[String, Seq[String]]

  
  /*var last_tick = sim.tick
    sim.listen(classOf[EV_Heartbeat], _ match { case e: EV_Heartbeat => {
      Util.log("[%.0fx] %s".format(e.tick - last_tick, e.describe))
      last_tick = e.tick
  }})
  */
    
  // Listen to events in the simulator and update the values accordingly
  sim.listen(classOf[EV_Reroute], _ match {
    case EV_Reroute(_, _, _, method, unrealizable, _) => method match {
      case x if x != RouterType.Fixed && x != RouterType.Unusable => {
        astar_count += 1
        if (unrealizable) {
          unrealizable_count += 1
          ratio_au = unrealizable_count.toFloat/astar_count
        }
      }
      case _ =>
    }
    case _ =>
  })
  
  sim.listen(classOf[EV_TurnFinished], _ match {
    // could be accept_delay
    case e: EV_TurnFinished => if(avg_delay == 0.0){avg_delay = e.total_delay} 
    						    else {avg_delay = (avg_delay + e.total_delay).toFloat / 2}
   })
    
    private var i = 0
    private val losers_list = Array(0, 0, 0, 0, 0)
    sim.listen(classOf[EV_IntersectionOutcome], _ match {
    	case EV_IntersectionOutcome(policy, losers) => 
    	  //println("Losers:" + losers.size)
    	  
    	  if (i < 5){
    	    losers_list(i) = losers.size
    	    i+=1//if(losers.size > 0){losers.foreach(loser => println(loser.intersection))} 
    	  }else{
    	    i = 0
    	    var k = 0
    	    var sum = 0
    	    while(k < losers_list.length){
    	      sum += losers_list(k)
    	      k += 1
    	    }
    	    avg_losers = sum.toFloat / 5
    	  }
    	  //println(losers_list(4))
    	  //println(avg_losers)
     })
     
     sim.listen(classOf[EV_IntersectionTickets], _ match{
       case EV_IntersectionTickets(policy, all) =>
         //println("All tickets: ")
         //all.foreach(ticket => println(ticket))
         //previous_traffic = discrete_traffic.clone()
         
         
         all.foreach(ticket =>{
        	 // parse the length of each road
        	 // will create segments according to it later
        	  if (!uniq_roads.contains(ticket.turn.from.road.id.toString())){
               uniq_roads.put(ticket.turn.from.road.id.toString(), ticket.turn.from.road.length)
             }
        	  
           val segment_len = uniq_roads.get(ticket.turn.from.road.id.toString()).get./(segments)
           
           if (ticket.turn.from.road.dir.toString() == "+" ){
             if (ticket.turn.from.road.name == "South"){
               
               for (x <- 1 to segments){
            	   if (ticket.dist_away <= x * segment_len){
            	     cars_present(0)(x-1) += 1
            	   }
               }
             }
             else if (ticket.turn.from.road.name == "West"){
               
                for (x <- 1 to segments){
            	   if (ticket.dist_away <= x * segment_len){
            	     cars_present(2)(x-1) += 1
            	   }
               }
               
             }
           }
           else if (ticket.turn.from.road.dir.toString() == "-"){ 
             if (ticket.turn.from.road.name == "North"){
                 for (x <- 1 to segments){
            	   if (ticket.dist_away <= x * segment_len){
            	     cars_present(1)(x-1) += 1
            	   }
               }

             }
             else if (ticket.turn.from.road.name == "East"){
                 for (x <- 1 to segments){
            	   if (ticket.dist_away <= x * segment_len){
            	     cars_present(3)(x-1) += 1
            	   }
                 }	

             }
           }
         })
         
        // reduce to three measurements (low, medium, high)
        for(i <- 0 to cars_present.length - 1){
          for (j <- 0 to cars_present(i).length - 1 ){
            if(cars_present(i)(j) <= 10)
              discrete_traffic(i)(j) = "low"
            else if(cars_present(i)(j) > 10 && cars_present(i)(j) <= 20)
             discrete_traffic(i)(j) = "medium"
           else if (cars_present(i)(j) > 20)
             discrete_traffic(i)(j) = "high"
          }
        } 
         
     })
    
     
     sim.listen(classOf[EV_Signal_Change], _ match{
       case f: EV_Signal_Change => if (f.greens.size >= 4){
         
          var action_list: Seq[String] = {
        	 val result = mutable.ArrayBuffer[String]()
        			 f.greens.foreach { green =>
        			 result += green.id.toString()
        	 	}
        	 result
         }
         
         var act = action_list.sortBy(_.toString())
         
         actions.put("signals", act)

       }
     })
     
     def reset_carsPresent(){
      
	  for(i <- 0 to cars_present.length - 1){
        for (j <- 0 to cars_present(i).length - 1){
        	cars_present(i)(j) = 0
        	discrete_traffic(i)(j) = "low"
        }
      } 
  	}
      def print_carsPresent(){
         for(i <- 0 to cars_present.length - 1){
           for (j <- 0 to cars_present(i).length - 1){
        	   print("("+i+","+j+"):" + cars_present(i)(j) + " ")
           }
         }
         println()
      }
      
      def print_traffic(){
         for(i <- 0 to discrete_traffic.length - 1){
           for (j <- 0  to discrete_traffic(i).length - 1){
        	   print("("+i + "," + j + "):" + discrete_traffic(i)(j)+" ")
           }
         }
         println()
      }
      


}