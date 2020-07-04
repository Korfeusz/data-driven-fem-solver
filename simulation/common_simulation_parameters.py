import fenics

from generalized_alpha_parameters import GeneralizedAlphaParameters
from problem_definition import ExternalExcitation, FieldUpdates, ProblemForm, Traction, ElastodynamicsFieldUpdates, \
    ElastodynamicsForm
from space_definition import MeshCreator, BoundaryMarkers, DirichletBCCreator, UnitSquareMeshCreator, \
    SquareBoundaryDefinition, square_side_settings, SquareSideMarkerNames
from time_stepping_parameters import TimeSteppingParameters


class CommonSimulationParameters:
    def __init__(self):
        self._boundary_excitation = None
        self._alpha_params = None
        self._time_params = None
        self._constitutive_relation = None


    @property
    def mesh_creator(self) -> MeshCreator:
        return UnitSquareMeshCreator(horizontal_cell_no=5, vertical_cell_no=5)


    @property
    def boundary_markers(self) -> BoundaryMarkers:
        return BoundaryMarkers(tolerance=1e-14,
                        boundary_definition_constructor=SquareBoundaryDefinition,
                        boundary_definition_settings=square_side_settings,
                        boundary_names=SquareSideMarkerNames)

    @property
    def bc_creator(self) -> DirichletBCCreator:
        return DirichletBCCreator(boundary_names_for_condition=SquareSideMarkerNames.X0,
                           value=fenics.Constant((0, 0)))

    @property
    def boundary_excitation(self) -> ExternalExcitation:
        if self._boundary_excitation is None:
            self._boundary_excitation = Traction(initial_value=1.0,
                               cutoff_time=self.time_params.total_time / 5,
                               traction_boundary=SquareSideMarkerNames.X1.value)
        return self._boundary_excitation


    @property
    def field_updates(self) -> FieldUpdates:
        return ElastodynamicsFieldUpdates(generalized_alpha_parameters=self.alpha_params,
                                           time_stepping_parameters=self.time_params)


    def problem(self, constitutive_relation) -> ProblemForm:
        return ElastodynamicsForm(mass_density=1.0,
                                  eta_m=fenics.Constant(0.0), eta_k=fenics.Constant(0.0),
                                  constitutive_relation=constitutive_relation,
                                  generalized_alpha_parameters=self.alpha_params,
                                  delta_t=self.time_params.delta_t,
                                  f_ext=self.boundary_excitation.value)


    @property
    def alpha_params(self) -> GeneralizedAlphaParameters:
        if self._alpha_params is None:
            self._alpha_params = GeneralizedAlphaParameters(alpha_m=0.2, alpha_f=0.4)
        return self._alpha_params

    @property
    def time_params(self) -> TimeSteppingParameters:
        if self._time_params is None:
            self._time_params = TimeSteppingParameters(total_time=4.0, number_of_steps=50)
        return self._time_params
