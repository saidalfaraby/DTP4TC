package utexas.aorta.learning

import utexas.aorta.sim.{Simulation, EV_AgentSpawned, EV_Reroute, EV_AgentQuit, EV_TurnFinished, EV_Heartbeat,
                         EV_IntersectionOutcome, EV_IntersectionOutcomeWin, EV_IntersectionTickets, EV_Transition, EV_Signal_Change}
import utexas.aorta.sim.make.{Scenario, MkAgent, Factory, RouterType}
import utexas.aorta.map.{Turn, Edge}
import utexas.aorta.common.{Util, Timer, cfg, Flags}
import scala.collection.mutable

class DBNState (sim: Simulation){
  
  // State Variables
  var avg_losers = 0.0
  var astar_count = 0
  var unrealizable_count = 0
  var ratio_au = 0.0
  var avg_delay = 0.0
  var thres_away = 80.0
  // 2*south, 2*north, 2* west, 2*east
  // south(+1, +0), north(-1, -0), west(+1, +0), east(-1, -0)
  //var cars_present = Array(0,0,0,0,0,0,0,0)
  // south, north, west, east
  var cars_present = Array(0,0,0,0)
  // low, medium, high for now
  //var discrete_traffic = Array("low","low","low","low","low","low","low","low") 
  var discrete_traffic = Array("low","low","low","low")
  //var previous_traffic = Array("low","low","low","low")
  
  // 16 binary values, according to traffic lights
  // south +1-> 
  //var traffic_lights = Array(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  private var from = " "
  
  var actions_per_lane = new mutable.HashMap[String, List[Edge]]
  //var previous_actions = new mutable.HashMap[String, List[Edge]]
  //var greens: Set[Turn]
  
  
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
           /*
           if (ticket.turn.from.road.dir.toString() == "+" && ticket.turn.from.lane_num.toString == "1" ){
             if (ticket.turn.from.road.name == "South"){
               if(ticket.dist_away < thres_away){
                   //println("+1 south:" + ticket.a.id)
            	   cars_present(0) += 1
               }
             }
             else if (ticket.turn.from.road.name == "West"){
               if(ticket.dist_away < thres_away){
                   //println("+1 west:" + ticket.a.id)
            	   cars_present(4) += 1
               }
             }
           }
           else if (ticket.turn.from.road.dir.toString() == "-" && ticket.turn.from.lane_num.toString == "1"){
             if (ticket.turn.from.road.name == "North"){
               if(ticket.dist_away < thres_away){
            	   //println("-1 north:" + ticket.a.id)
            	   cars_present(2) += 1
               }
             }
             else if (ticket.turn.from.road.name == "East"){
               if(ticket.dist_away < thres_away){
                   //println("-1 east:" + ticket.a.id)
            	   cars_present(6) += 1
               }
             }
           }*/
           if (ticket.turn.from.road.dir.toString() == "+" ){//&& ticket.turn.from.lane_num.toString == "0"){
             if (ticket.turn.from.road.name == "South"){
               if(ticket.dist_away < thres_away){
                   //println("+0 south:" + ticket.a.id)
            	   //cars_present(1) += 1
                     cars_present(0) += 1
               }
             }
             else if (ticket.turn.from.road.name == "West"){
               if(ticket.dist_away < thres_away){
                   //println("+0 west:" + ticket.a.id)
            	   //cars_present(5) += 1
            	   cars_present(2) += 1
               }
             }
           }
           else if (ticket.turn.from.road.dir.toString() == "-"){ //&& ticket.turn.from.lane_num.toString == "0"){
             if (ticket.turn.from.road.name == "North"){
               if(ticket.dist_away < thres_away){
                   //println("-0 north:" + ticket.a.id)
            	   //cars_present(3) += 1
            	   cars_present(1) += 1
               }
             }
             else if (ticket.turn.from.road.name == "East"){
               if(ticket.dist_away < thres_away){
                   //println("-0 east:" + ticket.a.id)
            	   //cars_present(7) += 1
            	   	cars_present(3) += 1
               }
             }
           }
           
           //println(ticket.turn.from.road.dir+" " +ticket.turn.from.lane_num+ " " + ticket.turn.from.road.name)//ticket.turn.from)//*ticket.turn.from.lane_num)
         	
         }
         )
         
        // reduce to three measurements (low, medium, high)
        for(i <- 0 to cars_present.length - 1){
           if(cars_present(i) >= 5 && cars_present(i) <= 10)
             discrete_traffic(i) = "medium"
           else if (cars_present(i) > 10)
             discrete_traffic(i) = "high"
         } 
         
     })
    
     
     //sim.listen(classOf[EV_IntersectionOutcomeWin], _ match{
     //  case EV_IntersectionOutcomeWin(policy, winners) =>
         //println("Winners:"  + winners.size)
     //})
     
     sim.listen(classOf[EV_Signal_Change], _ match{
       case f: EV_Signal_Change => if (f.greens.size >= 4){
         //println(f.greens.size)
         //previous_actions = actions_per_lane.clone()
         f.greens.foreach(green => {
             //println(green)
        	 //println(green.leads_to)
        	 actions_per_lane.put(green.from.road.name, green.leads_to)
         })
       }
     })
     
     def reset_carsPresent(){
      for(i <- 0 to cars_present.length - 1){
        cars_present(i) = 0
        discrete_traffic(i) = "low"
      } 
  	}
      def print_carsPresent(){
         for(i <- 0 to cars_present.length - 1){
        	print(" " + cars_present(i))
         }
         println()
      }
      
      def print_traffic(){
         for(i <- 0 to discrete_traffic.length - 1){
        	print(" " + discrete_traffic(i))
         }
         println()
      }
      

}