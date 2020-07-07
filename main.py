from simulation import  ElastodynamicSimulationParameters, Simulation, CommonSimulationParameters, DDDbParameters

simulation = Simulation(simulation_parameters=ElastodynamicSimulationParameters(),
                        common_simulation_parameters=CommonSimulationParameters())

# simulation = Simulation(simulation_parameters=DDDbParameters(),
#                         common_simulation_parameters=CommonSimulationParameters())

simulation.run()
