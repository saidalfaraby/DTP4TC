// AORTA is copyright (C) 2012 Dustin Carlino, Mike Depinet, and Piyush
// Khandelwal of UT Austin
// License: GNU GPL v2

package utexas.aorta.tests

import utexas.aorta.sim.drivers.Kinematic
import utexas.aorta.common.{cfg, Physics}

object FollowingDistance {
  def main(args: Array[String]) {
    val speed_lim = 30
    // Perfect start
    var leader = Kinematic(100, speed_lim, speed_lim)
    var follower = Kinematic(100 - following_dist(speed_lim), speed_lim, speed_lim)

    val leader_choices = List(0, 0, 0, -2.7, 0, 0, 0, 2.7)
    for (leader_choice <- leader_choices) {
      //val follower_choice = follower.accel_to_follow(leader, leader.dist - follower.dist)
      //val follower_choice = fixed_worst_case(follower, leader)
      //val follower_choice = fixed_best_case(follower, leader)
      //val follower_choice = dynamic_best_case(follower, leader)
      val follower_choice = levin(follower, leader)

      leader = step(leader, leader_choice)
      follower = step(follower, follower_choice)
      println(s"leader ${descr(leader)} and follower ${descr(follower)}.")
      println(s"leader used acceleration $leader_choice, follower used $follower_choice")
      val dist = leader.dist - follower.dist
      println(s"  $dist separation, ${follower.speed} speed of follower")
      //println(s"  ideally ${following_dist(leader.speed)} separation (one tick off though?)")
      if (dist < cfg.follow_dist) {
        println("  *** follower is too close")
      }
      println("")
    }
  }

  private def step(k: Kinematic, want_accel: Double): Kinematic = {
    val accel = math.max(-cfg.max_accel, math.min(want_accel, cfg.max_accel))
    return Kinematic(
      k.dist + Physics.dist_at_constant_accel(accel, cfg.dt_s, k.speed),
      Physics.update_speed(k.speed, accel, cfg.dt_s), k.speed_limit
    )
  }

  private def descr(k: Kinematic) = s"K(dist = ${k.dist}, speed = ${k.speed})"

  private def fixed_worst_case(follower: Kinematic, leader: Kinematic): Double = {
    val their_least_dist = leader.dist + leader.min_next_dist
    val want_final_dist = their_least_dist - cfg.follow_dist
    println(s"  their least total $their_least_dist")
    println(s"  want to wind up at $want_final_dist")
    return Physics.accel_to_cover(want_final_dist - follower.dist, follower.speed)
  }

  private def fixed_best_case(follower: Kinematic, leader: Kinematic): Double = {
    val their_max_dist = leader.dist + leader.max_next_dist
    val want_final_dist = their_max_dist - cfg.follow_dist
    println(s"  their max total $their_max_dist")
    println(s"  want to wind up at $want_final_dist")
    return Physics.accel_to_cover(want_final_dist - follower.dist, follower.speed)
  }

  private def dynamic_best_case(follower: Kinematic, leader: Kinematic): Double = {
    val their_max_dist = leader.dist + leader.max_next_dist
    val want_final_dist = their_max_dist - following_dist(leader.speed)
    println(s"  their max total $their_max_dist")
    println(s"  want to wind up at $want_final_dist")
    return Physics.accel_to_cover(want_final_dist - follower.dist, follower.speed)
  }

  private def levin(follower: Kinematic, leader: Kinematic): Double = {
    val pt1 = leader.dist + (leader.speed * cfg.dt_s) + (0.5 * -cfg.max_accel * cfg.dt_s * cfg.dt_s)
    val pt2 = -follower.dist - (follower.speed * cfg.dt_s) - cfg.follow_dist
    val pt3 = 0.5 * cfg.dt_s * cfg.dt_s
    return (pt1 + pt2) / pt3
  }

  private def following_dist(speed: Double) = cfg.follow_dist + Physics.min_next_dist(speed)
}
