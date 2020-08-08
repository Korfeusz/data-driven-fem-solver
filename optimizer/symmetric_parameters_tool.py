import  numpy as np

def make_parameters_symmetric(parameters: np.ndarray):
    b = [(i, v) for i, v in zip(range(1, len(parameters), 3), parameters[1::3].copy())]
    c = parameters.copy()
    for i, v in reversed(b):
        c = np.insert(c, i, v)
    return c