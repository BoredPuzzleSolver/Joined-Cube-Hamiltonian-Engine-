"""Small numerical cutoff comparison; all outputs remain uncertified."""
import json
from terry_cubes import JoinedCubesModel, cutoff_study

if __name__ == "__main__":
    result = cutoff_study(JoinedCubesModel(1), [6, 8, 10, 12], embedding_cutoff=18)
    print(json.dumps(result, indent=2, allow_nan=False))
