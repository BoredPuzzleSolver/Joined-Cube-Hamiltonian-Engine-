"""Fresh exact endpoint checks for the supported single-cube interval."""
import json
from terry_cubes import verify_certificate

if __name__ == "__main__":
    print(json.dumps(verify_certificate(cutoff=12), indent=2, allow_nan=False))
