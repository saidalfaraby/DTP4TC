// AORTA is copyright (C) 2012 Dustin Carlino, Mike Depinet, and Piyush
// Khandelwal of UT Austin
// License: GNU GPL v2

package utexas.aorta.sim

import utexas.aorta.ui.GUIDebugger
import utexas.aorta.analysis.{SimSpeedMonitor, ReplayReader, ReplayWriter}

import utexas.aorta.common.{Util, Timer, cfg, Flags}

object Headless {
  def main(args: Array[String]): Unit = {
    val sim = Util.process_args(args)
    // TODO move elsewhere?
    Flags.string("--benchmark") match {
      case Some(fn) => new SimSpeedMonitor(sim, fn)
      case None =>
    }
    Flags.string("--record") match {
      case Some(fn) => new ReplayWriter(sim, Util.writer(fn))
      case None =>
    }
    Flags.string("--replay") match {
      case Some(fn) => new ReplayReader(sim, Util.reader(fn))
      case None =>
    }
    val gui = new GUIDebugger(sim)

    // Print an update every second
    var last_tick = sim.tick
    sim.listen("headless", _ match {
      case e: EV_Heartbeat => {
        Util.log("[%.0fx] %s".format(e.tick - last_tick, e.describe))
        last_tick = e.tick
      }
      case _ =>
    })

    Util.log("Starting simulation with time-steps of " + cfg.dt_s + "s")
    val t = Timer("headless simulation")
    while (!sim.done) {
      sim.step()
    }
    t.stop
    sys.exit
  }
}
