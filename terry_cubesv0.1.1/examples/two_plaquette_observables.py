"""Original two-square symmetry channels; this graph is not two joined cubes."""
import json
from terry_cubes import TwoPlaquetteModel

if __name__ == "__main__":
    result = TwoPlaquetteModel().compute_spectrum(max_twice=12)
    print(json.dumps(result.correlators(), indent=2, allow_nan=False))
