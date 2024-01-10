""" Contains default parameters used in the simulation """
import numpy as np

class Inputs():

    def __init__(self) -> None:
        

        # North-America power grid
        self.america = dict()
        self.america["num_nodes"]   = 16167
        self.america["num_samples"] = 50      # Phase Diagram
        self.america["max_tstep"]   = 2000
        self.america["init0"]       = 0.05
        self.america["r0_list"]     = np.linspace(0,2,100)
        self.america["r1_list"]     = np.linspace(0,1.1*(1/self.america["init0"] ),100)
        self.america["event"]       = "EARL"
        self.america["r0"] = 10     # Risk map
        self.america["r1"]= 0.3
        self.america["path_edgelist"] = "../Data/Processed/Topologies/america/powergrid_north_america.el"
        self.america["path_nodelist"] = "../Data/Processed/Topologies/america/america.nodelist"
        self.america["path_risk"]     = f"./Output_OAD/america_r0_{self.america['r0']}_r1_{self.america['r1']}_samples_10_maxtime_2000.dat"

        # Europa power grid
        self.europe = dict()
        self.europe["num_nodes"]    = 13844
        self.europe["num_samples"]  = 25
        self.europe["max_tstep"]    = 2000
        self.europe["init0"]       = 0.05
        self.europe["r0_list"]     = np.linspace(0,2,100)
        self.europe["r1_list"]     = np.linspace(0,1.1*(1/self.america["init0"] ),100)        
        self.europe["event"]        = "ciaran"
        self.europe["r0"] = 8
        self.europe["r1"]= 0.3
        self.europe["path_edgelist"] = "../Data/Processed/Topologies/europe/powergrid_europe.el"
        self.europe["path_nodelist"] = "../Data/Processed/Topologies/europe/europe.nodelist"
        self.europe["path_risk"]     = f"./Output_OAD/europe_r0_{self.europe['r0']}_r1_{self.europe['r1']}_samples_10_maxtime_2000.dat"

        # Airports network
        self.airports = dict()
        self.airports["num_nodes"]  = 3182
        self.airports["num_samples"]= 50
        self.airports["max_tstep"]  = 2000
        self.airports["init0"]       = 0.05
        self.airports["r0_list"]     = np.linspace(0,0.1,100)
        self.airports["r1_list"]     = np.linspace(0,3,100)        
        self.airports["event"]      = "earthquakes"
        self.airports["r0"] = 0.3
        self.airports["r1"]= 2
        self.airports["path_edgelist"] = "../Data/Processed/Topologies/airports/airports.edgelist"
        self.airports["path_nodelist"] = "../Data/Processed/Topologies/airports/airports.nodelist"
        self.airports["path_risk"]     = f"./Output_OAD/airports_r0_{self.airports['r0']}_r1_{self.airports['r1']}_samples_10_maxtime_2000.dat"

        # Random network
        self.random = dict()
        self.random["num_nodes"]  = 5000
        self.random["num_samples"]= 50
        self.random["max_tstep"]  = 2000
        self.random["prob_link"]  = 0.05
        self.random["init0"]      = 0.05

        # Full network
        self.full = dict()
        self.full["num_nodes"]  = 5000
        self.full["num_samples"]= 50
        self.full["max_tstep"]  = 2000
        self.full["init0"]      = 0.05