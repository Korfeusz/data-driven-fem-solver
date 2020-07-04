import abc
from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation, FieldUpdates, ProblemForm
from space_definition import MeshCreator, BoundaryMarkers, DirichletBCCreator
from time_stepping_parameters import TimeSteppingParameters


class CommonSimulationParameters:
    @property
    def mesh_creator(self) -> MeshCreator:
        pass


    @property
    def boundary_markers(self) -> BoundaryMarkers:
        pass

    @property
    def bc_creator(self) -> DirichletBCCreator:
        pass

    @property
    def boundary_excitation(self) -> ExternalExcitation:
        pass


    @property
    def field_updates(self) -> FieldUpdates:
        pass

    @property
    def problem(self) -> ProblemForm:
        pass

    @property
    def alpha_params(self) -> GeneralizedAlphaParameters:
        pass

    @property
    def time_params(self) -> TimeSteppingParameters:
        pass
