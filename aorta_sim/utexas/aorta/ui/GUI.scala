// AORTA is copyright (C) 2012 Dustin Carlino, Mike Depinet, and Piyush
// Khandelwal of UT Austin
// License: GNU GPL v2

package utexas.aorta.ui

import swing._  // TODO figure out exactly what
import java.awt.{Color, Component}
import swing.Dialog
import javax.swing.WindowConstants
import java.io.File

import utexas.aorta.sim.{Simulation, EV_Heartbeat}
import utexas.aorta.common.{Util, cfg}

object Status_Bar {
  val zoom       = new Label("1.0") // TODO from cfg
  val agents     = new Label("0 moved / 0 live / 0 ready")
  val time       = new Label("0.0")
  val sim_speed = new Label("Paused / 1x")

  val panel = new GridBagPanel {
    maximumSize = new Dimension(Int.MaxValue, 10)
    border = Swing.MatteBorder(5, 5, 5, 5, Color.BLACK)

    // TODO generate these?

    // all of this to prevent the rightmost 'At' column from spazzing out when
    // the text changes length
    // row 1: labels
    val c = new Constraints
    c.gridx = 0
    c.gridy = 0
    c.ipadx = 50
    layout(new Label("Zoom")) = c
    c.gridx = 1
    layout(new Label("Agents Active/Ready/Routing")) = c
    c.gridx = 2
    layout(new Label("Time")) = c
    c.gridx = 3
    layout(new Label("Sim Speed (Actual/Cap)")) = c

    // row 2: columns
    c.weightx = 0.0
    c.ipadx = 50
    c.gridx = 0
    c.gridy = 1
    layout(Status_Bar.zoom) = c
    c.gridx = 1
    layout(Status_Bar.agents) = c
    c.gridx = 2
    layout(Status_Bar.time) = c
    c.gridx = 3
    layout(Status_Bar.sim_speed) = c
  }
}

// TODO SwingApplication has a startup, quit, shutdown...
object GUI extends SimpleSwingApplication {
  val road_types = List(
    "null", "residential", "unclassified", "secondary",
    "motorway_link", "motorway", "trunk_link", "secondary_link", "primary_link",
    "tertiary", "primary", "service"
  )
  // null just because it's parametric from argv
  var canvas_2d: MapCanvas = null

  val helper = new BoxPanel(Orientation.Vertical) {
    border = Swing.MatteBorder(5, 5, 5, 5, Color.BLACK)
    yLayoutAlignment = java.awt.Component.TOP_ALIGNMENT
    // TODO These're fixed now, but the idea is to tie them to configuration and
    // add/remove some context-sensitively. And also organize better.

    // Simulation controls
    contents += new Label("p   pause/resume")
    contents += new Label("[   slow down time")
    contents += new Label("]   speed up time")
    contents += new Label("-   slown down time faster")
    contents += new Label("=   speed up time faster")

    // Actions
    contents += new Label("m   make new agent on current edge")
    contents += new Label("c   choose edge for pathfinding")
    contents += new Label("d   object-sensitive debug")
    contents += new Label("f   follow agent")
    contents += new Label("x   delete agent (may crash!)")

    // Polygons
    contents += new Label("Shift+Left click   draw a polygon")
    contents += new Label("Shift+s   begin/end agents in polygon")
    contents += new Label("Shift+p   change intersection policies")

    // View
    contents += new Label("r   reset view")
    contents += new Label("g   toggle greenflood colors")
    contents += new Label("CTRL   cycle through turns")
    contents += new Label("arrow keys pan")

    // TODO expand to fill the whole column, or otherwise work on aesthetics
    // TODO and option to hide the panel
  }

  private var headless = false
  var closed = false

  override def main(args: Array[String]) = {
    val sim = Util.process_args(args)
    canvas_2d = new MapCanvas(sim)
    // TODO doesnt start drawn correctly!
    canvas_2d.repaint
    super.main(args)
  }

  def launch_from_headless(canvas: MapCanvas) = {
    headless = true
    canvas_2d = canvas
    super.main(Array())
  }

  def top = new MainFrame {
    title = "AORTA"
    preferredSize = new Dimension(800, 600)
    peer.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE)

    override def closeOperation() {
      if (headless) {
        println("Closing GUI...")
        close
        closed = true
      } else {
        sys.exit
      }
    }
    
    menuBar = new MenuBar {
      contents += new Menu("File") {
        contents += new MenuItem(Action("Configuration") {
          popup_config
        })
        contents += new MenuItem(
          Action("Save scenario for later resimulation")
        {
          // TODO
          val fn = "resim_log"
          //Simulation.save_log(fn)
          Dialog.showMessage(message = "TODO Scenario saved in '" + fn + "'")
        })
        contents += new Separator
        contents += new MenuItem(Action("Quit") {
          sys.exit
        })
      }

      contents += new Menu("View") {
        contents += new Menu("Highlight type of road") {
          contents ++= road_types.map(t => new MenuItem(Action(t) {
            canvas_2d.handle_ev(EV_Param_Set("highlight", Some(t)))
          }))
        }
        contents += new MenuItem(Action("Clear all highlighting") {
          canvas_2d.handle_ev(EV_Param_Set("highlight", None))
        })
      }

      contents += new Menu("Query") {
        contents += new MenuItem(Action("Teleport to Edge") {
          canvas_2d.handle_ev(EV_Action("teleport-edge"))
        })
        contents += new MenuItem(Action("Teleport to Road") {
          canvas_2d.handle_ev(EV_Action("teleport-road"))
        })
        contents += new MenuItem(Action("Teleport to Agent") {
          canvas_2d.handle_ev(EV_Action("teleport-agent"))
        })
        contents += new MenuItem(Action("Teleport to Vertex") {
          canvas_2d.handle_ev(EV_Action("teleport-vertex"))
        })
        
        // TODO these are kind of toggleable...
        contents += new MenuItem(Action("Pathfind") {
          canvas_2d.handle_ev(EV_Action("pathfind"))
        })
        contents += new MenuItem(Action("Clear Route") {
          canvas_2d.handle_ev(EV_Action("clear-route"))
        })
      }

      contents += new Menu("Simulate") {
        contents += new MenuItem(Action("Spawn Army") {
          canvas_2d.handle_ev(EV_Action("spawn-army"))
        })
        contents += new MenuItem(Action("Play/Pause") {
          canvas_2d.handle_ev(EV_Action("toggle-running"))
        })
      }
    }

    // TODO toggle between helper and other stuff in right pane
    val main_content = canvas_2d

    contents = new BorderPanel {
      background = Color.LIGHT_GRAY
      border = Swing.MatteBorder(2, 2, 2, 2, Color.RED)

      add(Status_Bar.panel, BorderPanel.Position.North)
      add(main_content, BorderPanel.Position.Center)
    }
  }

  def popup_config() {
    // TODO tabbed pane by category?
  }
}

// TODO make pause work with whoever's calling us
class GUIDebugger(sim: Simulation) {
  // When this file exists, launch a GUI for sudden interactive watching.
  private val gui_signal = new File(".headless_gui")
  private var gui: Option[MapCanvas] = None

  sim.listen("gui-debugger", _ match {
    case e: EV_Heartbeat => {
      if (gui_signal.exists) {
        gui_signal.delete()
        gui match {
          case Some(ui) => {
            if (GUI.closed) {
              println("Resuming the GUI...")
              GUI.top.open()
              GUI.closed = false
            }
          }
          case None => {
            println("Launching the GUI...")
            gui = Some(new MapCanvas(sim, headless = true))
            GUI.launch_from_headless(gui.get)
          }
        }
      }
      gui match {
        case Some(ui) if !GUI.closed => ui.handle_ev(EV_Action("step"))
        case _ =>
      }
    }
    case _ =>
  })
}
