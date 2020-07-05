from space_definition import Spaces
from .fields import Fields

class ElastodynamicsFields(Fields):
    def initialize(self, spaces: Spaces) -> None:
        self.vector_space = spaces.vector_space

    def __init__(self):
        super().__init__()
