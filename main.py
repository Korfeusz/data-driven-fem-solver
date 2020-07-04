from simulation import  ElastodynamicSimulationParameters, Simulation, CommonSimulationParameters

simulation = Simulation(simulation_parameters=ElastodynamicSimulationParameters(),
                        common_simulation_parameters=CommonSimulationParameters())

simulation.run()
