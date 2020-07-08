from simulation import  ElastodynamicSimulationParameters, Simulation, CommonSimulationParameters, DDDbParameters
import logging
import fenics
# logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
fenics.set_log_level(logging.WARNING)
simulation = Simulation(simulation_parameters=ElastodynamicSimulationParameters(),
                        common_simulation_parameters=CommonSimulationParameters())
#
# simulation = Simulation(simulation_parameters=DDDbParameters(),
#                         common_simulation_parameters=CommonSimulationParameters())

simulation.run()
