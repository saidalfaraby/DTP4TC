// AORTA is copyright (C) 2012 Dustin Carlino, Mike Depinet, and Piyush
// Khandelwal of UT Austin
// License: GNU GPL v2

package utexas.aorta.sim.intersections

import scala.collection.mutable

import utexas.aorta.map.{Vertex, Turn, Edge}
import utexas.aorta.sim.{Simulation, EV_TurnApproved, EV_TurnStarted}
import utexas.aorta.sim.make.{IntersectionType, OrderingType, Factory}

import utexas.aorta.common.{Util, StateWriter, StateReader, TurnID}

// Reason about collisions from conflicting simultaneous turns.
class Intersection(val v: Vertex, policy_type: IntersectionType.Value,
                   val ordering_type: OrderingType.Value, sim: Simulation)
{
  val policy = Factory.make_policy(this, policy_type, ordering_type, sim)

  // Multiple agents can be on the same turn; the corresponding queue will
  // handle collisions. So in fact, we want to track which turns are active...
  // but we have to know how many are on the turn to know when nobody's
  // attempting it.
  val turns = mutable.HashMap[Turn, Int]() // TODO count type

  // For reporting stats
  private var sum_waiting_times = 0.0
  private var number_finished = 0

  override def toString = "Intersection(" + v + ")"

  def request_turn(ticket: Ticket) = {
    // Sanity check...
    Util.assert_eq(ticket.turn.vert, v)
    policy.request_turn(ticket)
  }

  def cancel_turn(ticket: Ticket) {
    Util.assert_eq(ticket.turn.vert, v)
    policy.cancel_turn(ticket)
  }

  def enter(ticket: Ticket) = {
    val t = ticket.turn
    if (!turns.contains(t)) {
      // We don't care until there are at least two... and this only changes when
      // we add a new one...
      // It's not really that much better to do the checking in-place and only
      // atomic-check when there could be a problem.
      if (turns.size == 1) {
        ticket.a.sim.active_intersections += this
      }

      turns(t) = 0
    }
    turns(t) += 1
    if (!policy.validate_entry(ticket)) {
      // Misleading error message. They may be going 0 speed, but the agent step
      // code hasn't finished moving them.
      Util.log(s"!!! ${ticket.a} illegally entered $this, going ${ticket.a.speed} m/s")
      Util.log("  Illegal entry was near " + t.from + " and " + t + " (vert " + t.from.to.id + ")")
      Util.log("  Origin lane length: " + t.from.length + "; time " + ticket.a.sim.tick)
      ticket.a.debug()
      policy.dump_info()
      sys.exit()
    }
    ticket.a.sim.publish(EV_TurnStarted(ticket))
  }

  def exit(ticket: Ticket) = {
    val t = ticket.turn
    turns(t) -= 1
    if (turns(t) == 0) {
      turns -= t
      // Potentially we now don't care...
      if (turns.size == 1) {
        ticket.a.sim.active_intersections -= this
      }
    }
    policy.handle_exit(ticket)
    number_finished += 1
    sum_waiting_times += ticket.how_long_waiting
  }

  def average_waiting_time =
    if (number_finished > 0)
      sum_waiting_times / number_finished
    else
      0.0
}

object Intersection {
  // TODO wheres this belong?
  def detect_gridlock(turn: Turn): Boolean = {
    var current = turn.from
    val seen = new mutable.HashSet[Edge]()
    while (current != null && !current.queue.slot_avail) {
      // A cycle!
      if (seen(current)) {
        //Util.log(s"Gridlock detected, involving: $seen")
        return true
      }
      seen += current

      // Where's the head of that stuck queue trying to go?
      current = current.queue.head match {
        // requires invariant that we don't grab a ticket till we're committed to that lane
        case Some(a) => a.get_ticket(current) match {
          case Some(ticket) if !ticket.is_approved => ticket.turn.to
          case _ => null
        }
        case None => null
      }
    }
    return false
  }
}

abstract class Policy(val intersection: Intersection) {
  //////////////////////////////////////////////////////////////////////////////
  // State

  // This will have a deterministic order.
  protected var request_queue: List[Ticket] = Nil
  protected val accepted = new mutable.TreeSet[Ticket]()

  // Agents could be added to this in any order, but they'll wind up in
  // request_queue in a deterministic order. This is transient state.
  private val new_requests = new mutable.TreeSet[Ticket]()

  //////////////////////////////////////////////////////////////////////////////
  // State

  def serialize(w: StateWriter) {
    Util.assert_eq(new_requests.isEmpty, true)
    w.int(request_queue.size)
    for (ticket <- request_queue) {
      w.int(ticket.a.id.int)
      w.int(ticket.turn.id.int)
    }
  }

  protected def unserialize(r: StateReader, sim: Simulation) {}

  //////////////////////////////////////////////////////////////////////////////
  // Actions

  // Agents inform intersections of their intention ONCE and receive a lease
  // eventually.
  def request_turn(ticket: Ticket) = {
    synchronized {
      new_requests += ticket
      // TODO do extra book-keeping to verify agents aren't double requesting?
    }
  }

  def cancel_turn(ticket: Ticket) {
    synchronized {
      // TODO assert not in new_requests
      // TODO assert it was in here!
      request_queue = request_queue.filter(t => t != ticket)
    }
  }

  def react_tick() {
    request_queue ++= new_requests
    new_requests.clear()
    react()
  }

  protected def accept(ticket: Ticket) {
    ticket.approve()
    accepted += ticket
    unqueue(ticket)
    ticket.a.sim.publish(EV_TurnApproved(ticket))
  }

  protected def unqueue(ticket: Ticket) {
    request_queue = request_queue.filter(_ != ticket)
  }

  def react(): Unit

  // This could be called for several reasons, so assume they could be queued or
  // accepted.
  def handle_exit(ticket: Ticket) {
    accepted -= ticket
    unqueue(ticket)
  }

  def unserialize_accepted(ticket: Ticket) {
    accepted += ticket
  }

  // Check for collisions
  def end_step(): Unit = {
    if (intersection.turns.size < 2) {
      return
    }

    // O(n^2 / 2)
    for (t1 <- intersection.turns.keys; t2 <- intersection.turns.keys if t1.id.int < t2.id.int) {
      if (t1.conflicts_with(t2)) {
        throw new Exception(s"Intersection collision: $t1 and $t2 conflict")
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Queries

  def policy_type(): IntersectionType.Value

  def dump_info() {
    Util.log(s"$intersection is a $policy_type")
    Util.log(s"Accepted: $accepted")
    Util.log(s"Queued: $request_queue")
  }
  def approveds_to(target: Edge) = accepted.filter(_.turn.to == target)
  def validate_entry(ticket: Ticket) = accepted.contains(ticket)
  def current_greens() = accepted.map(_.turn).toSet
  def queued_count = request_queue.size
}

object Policy {
  def unserialize(policy: Policy, r: StateReader, sim: Simulation) {
    val num_requests = r.int
    for (i <- Range(0, num_requests)) {
      policy.request_queue :+= find_ticket(r, sim)
    }
    policy.unserialize(r, sim)
  }

  def find_ticket(r: StateReader, sim: Simulation): Ticket =
    find_ticket(sim, r.int, new TurnID(r.int))
  def find_ticket(sim: Simulation, agent_id: Int, turn_id: TurnID): Ticket =
    sim.get_agent(agent_id).get.get_ticket(sim.graph.turns(turn_id)).get
}

// Simplest base-line ever.
class NeverGoPolicy(intersection: Intersection) extends Policy(intersection) {
  def react() {}
  def policy_type = IntersectionType.NeverGo
}
