// AORTA is copyright (C) 2012 Dustin Carlino, Mike Depinet, and Piyush
// Khandelwal of UT Austin
// License: GNU GPL v2

package utexas.aorta.analysis

import utexas.aorta.sim.{Simulation, EV_Heartbeat, EV_AgentSpawned, EV_Reroute}
import utexas.aorta.sim.make.RouterType
import utexas.aorta.common.IO

// Dump a file showing how fast simulations are running
class SimSpeedMonitor(sim: Simulation, fn: String) {
  private val freq_ms = 1000

  // Don't have to explicitly close the file when simulation finishes
  // TODO fn library. denote what we're running on.
  private val file = new IO(None).output_file("sim_speed_" + fn)
  file.println("realtime_sec sim_tick")
  private val start = System.currentTimeMillis
  private var last_record = start

  sim.listen("sim_speed", _ match {
    case e: EV_Heartbeat => {
      val now = System.currentTimeMillis
      if (now - last_record >= freq_ms) {
        record((now - start) / 1000, e.tick)
        last_record = now
      }
    }
    case _ =>
  })

  private def record(realtime: Long, tick: Double) {
    file.println(s"$realtime $tick")
  }
}

// Count how many times agents are rerouting
class RerouteCountMonitor(sim: Simulation) {
  var astar_count = 0
  var unrealizable_count = 0
  def discretionary_count = astar_count - unrealizable_count

  sim.listen("reroute_count", _ match {
    case EV_Reroute(_, _, _, method, unrealizable, _) => method match {
      case x if x != RouterType.Fixed && x != RouterType.Unusable => {
        astar_count += 1
        if (unrealizable) {
          unrealizable_count += 1
        }
      }
      case _ =>
    }
    case _ =>
  })

  def reset() {
    astar_count = 0
    unrealizable_count = 0
  }
}
