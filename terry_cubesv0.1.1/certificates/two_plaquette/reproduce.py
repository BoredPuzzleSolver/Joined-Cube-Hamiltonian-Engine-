"""Run the packaged studies and compare with the supplied reference results."""
from pathlib import Path
import json
import math
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def compare(reference, actual, path='root'):
    if isinstance(reference, dict):
        assert reference.keys() == actual.keys(), path
        for key in reference:
            if key in ['elapsed_seconds', 'seconds', 'numerical_ordering_consistent',
                       'roundoff_limited']:
                continue
            compare(reference[key], actual[key], path+'.'+key)
    elif isinstance(reference, list):
        assert len(reference) == len(actual), path
        for index, (expected, observed) in enumerate(zip(reference, actual)):
            compare(expected, observed, f'{path}[{index}]')
    elif isinstance(reference, float):
        assert math.isclose(reference, actual, rel_tol=1e-10, abs_tol=1e-10), (
            path, reference, actual
        )
    else:
        assert reference == actual, (path, reference, actual)


def main():
    log = []
    for script in ['shared_su2_cutoff_study.py', 'verify_su2_cutoff_certificates.py',
                   'su2_direct_haar_validation.py', 'su2_bounds_audit_check.py']:
        run = subprocess.run([sys.executable, str(ROOT/'work'/script)],
                             cwd=ROOT, text=True, capture_output=True)
        log.append(f'## {script}\nExit code: {run.returncode}\n{run.stdout}\n{run.stderr}')
        if run.returncode:
            (ROOT/'outputs'/'Reproduction_Log.txt').write_text('\n'.join(log), encoding='utf-8')
            raise RuntimeError(f'{script} failed; see outputs/Reproduction_Log.txt')
        print(script+': passed', flush=True)
    compared = []
    for source in sorted((ROOT/'reference_results').glob('*.json')):
        expected = json.loads(source.read_text(encoding='utf-8'))
        observed = json.loads((ROOT/'outputs'/source.name).read_text(encoding='utf-8'))
        compare(expected, observed)
        compared.append(source.name)
    summary = {
        'all_scripts_passed': True, 'reference_results_matched': compared,
        'comparison': 'Exact strings, integers and booleans except roundoff diagnostic flags; floats at 1e-10 relative/absolute tolerance; elapsed times ignored',
        'certificate_note': 'Exact rational endpoints and inertia counts match exactly, independent of float comparison tolerance',
    }
    (ROOT/'outputs'/'Reproduction_Log.txt').write_text('\n'.join(log), encoding='utf-8')
    (ROOT/'outputs'/'Reproduction_Check.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print('All reference results reproduced; exact certificate endpoints and sign counts match.')


if __name__ == '__main__':
    main()
