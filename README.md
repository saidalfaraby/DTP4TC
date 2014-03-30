========= INIT ==============
- We used Scala IDE for the project
    - Can be downloaded here: http://scala-ide.org/
- Egit plugin for eclipse can quickly clone the git repository where we have the eclipse project (easier than importing the project from disk)
    - Can be downloaded here: https://www.eclipse.org/egit/
    - Or also from the eclipse marketplace
        - Help -> eclipse marketplace -> search for egit -> install
    - Once installed then right click on projects -> Import -> Git -> Project from Git -> fill in the necessary info
    - Will create all the necessary dependencies
- If it complains about no scala compiler found then:
    - Right click on project -> Scala -> Add Scala Library to Build Path


========= CODE =============
- All of our contribution are in the utexas.aorta.learning package
    - DBNState_segment.scala
        - Contains the class that parses all the necessary statistics so as to determine the state
    - DBN_segment.scala
        - DBN for the global transition model
        - Updates the counts based on the structure which is predefined
        - Also initializes and gives data to the tree representation
    - Config.scala
        - Contains some necessary prespecified configurations for the tree
            - keepGathering: whether to only gather data on this episode or build the tree (Boolean)
            - Possible values for the variables
            - Possible values for the action variables
            - number of segments
            - name of locations for the lanes
    - Util.scala
        - Contains misc usefull methods
            - Show stats of the data
            - Merge data from different episodes
    - ADD.scala
        - Contains everything needed for constructing the tree
        - In the Model class inside
            - Specify scoring function (MDL, BIC, Entropy)
            - Specify type of tree (Full, Compact)

======== MISC STUFF ===========
- All the differen scenarios we created are under the aorta/scenarios folder
    - e.g. only_traffic_signals_one_500
        - Dummy map we used (one intersection)
        - 500 cars will be spawned
        - Traffic signals is the only thing that affects the traffic
- The necessary maps can be found at aorta/maps
- The data are stored in the folder aorta/previousDATA with the extension .data

======== EXECUTION ============
- In order to run the code we created some run configurations on eclipse depending on the type of execution
    - GUI, Launches a simulation with the GUI
        - project : aorta_sim
        - main class : utexas.aorta.ui.GUI
        - include system libraries and inhereted mains
        - arguments: The scenario or map
            - e.g. scenarios/only_traffic_signals_one_1000
    - Headless, Launches a simulation without the GUI
        - project : aorta_sim
        - main class : utexas.aorta.sim.Headless
        - include system libraries and inhereted mains
        - arguments : The scenario or map
            - e.g. scenarios/only_traffic_signals_one_1000
    - Utils, Can be used to merge data from different scenarios
        - project : aorta_sim
        - main class : utexas.aorta.learning.Util
        - include system libraries and inhereted mains
        - var listFilename = List("previousDATA/merged_newactions.data", "previousDATA/merged_old_actions.data")
            - Specify the filenames of the data that should be merged
        - util.saveprevdata("previousDATA/merged_all.data)
            - Call the function and specify the filename for the meurged data


