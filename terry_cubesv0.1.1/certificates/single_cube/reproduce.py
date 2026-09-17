"""Run the delivered cube study and compare all supplied reference JSON."""
from pathlib import Path
import json
import math
import os
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def compare(expected,actual,path='root'):
    if isinstance(expected,dict):
        assert expected.keys()==actual.keys(),path
        for key in expected:
            if key.lower().endswith('seconds'):
                continue
            compare(expected[key],actual[key],path+'.'+key)
    elif isinstance(expected,list):
        assert len(expected)==len(actual),path
        for index,(a,b) in enumerate(zip(expected,actual)):
            compare(a,b,f'{path}[{index}]')
    elif isinstance(expected,float):
        assert math.isclose(expected,actual,rel_tol=1e-9,abs_tol=1e-9),(path,expected,actual)
    else:
        assert expected==actual,(path,expected,actual)


def main():
    scripts=['su2_cube_model.py','su2_cube_independent_check.py','su2_sparse_spectrum.py',
             'su2_cube_study.py','su2_cube_energy_study.py','su2_cube_bounds.py',
             'su2_cube_rational_certificate.py','su2_interval_inertia.py',
             'su2_cube_energy_certificate.py','su2_cube_certificate_audit_check.py',
             'su2_interval_inertia_independent_audit.py']
    # An additional independent interval-method audit is included when present.
    for optional in sorted((ROOT/'work').glob('su2_cube_interval*audit*.py')):
        if optional.name not in scripts:scripts.append(optional.name)
    log=[]
    for script in scripts:
        process=subprocess.run([sys.executable,str(ROOT/'work'/script)],cwd=ROOT,
                               env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),
                               capture_output=True,text=True)
        log.append(f'## {script}\nExit code: {process.returncode}\n{process.stdout}\n{process.stderr}')
        if process.returncode:
            (ROOT/'outputs'/'Reproduction_Log.txt').write_text('\n'.join(log),encoding='utf-8')
            raise RuntimeError(f'{script} failed; see outputs/Reproduction_Log.txt')
        print(script+': passed',flush=True)
    matched=[]
    for reference in sorted((ROOT/'reference_results').glob('*.json')):
        expected=json.loads(reference.read_text(encoding='utf-8'))
        actual=json.loads((ROOT/'outputs'/reference.name).read_text(encoding='utf-8'))
        compare(expected,actual)
        matched.append(reference.name)
    result={'all_scripts_passed':True,'reference_results_matched':matched,
            'comparison':'Exact fractions, integer sign counts, labels and booleans must match exactly; floats allow 1e-9 relative/absolute variation; timing fields ignored',
            'scope':'Reproduction of finite-cube computations and certificates, not a continuum Yang-Mills proof'}
    (ROOT/'outputs'/'Reproduction_Check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    (ROOT/'outputs'/'Reproduction_Log.txt').write_text('\n'.join(log),encoding='utf-8')
    print('All reference results reproduced; exact endpoints and sign counts match.',flush=True)


if __name__=='__main__':
    main()
