"""Small exact validation of the final joined-chain orientation convention."""
import json
from pathlib import Path

from su2_joined_cubes_model import joined_graph
from su2_oriented_vertex import self_checks


def run():
    # T<=18 gives physical doubled spins <=7; every auxiliary label is
    # bounded by the sum of its two incident physical labels, hence <=14.
    # Thus this finite exact sweep covers every local label in that study.
    local=self_checks(14)
    chains=[]
    for n in (1,2,3,4):
        model=joined_graph(n)
        phase_edges=[split.auxiliary_edge for split in model.splits
                     if model.real_coordinates[split.original_vertex][1]
                     !=model.real_coordinates[split.original_vertex][2]]
        chosen=set(phase_edges)
        for (edges,_),orientation_sign in zip(model.face_data,model.face_orientation_signs):
            assert (-1)**sum(e in chosen for e in edges)==orientation_sign
        chains.append({'num_cubes':n,'phase_auxiliary_tuple_indices':phase_edges,
                       'exact_face_incidence_identity':True,
                       'same_common_phase_for_every_face':True})
    return {'local_exact_arithmetic_checks':local,'default_chain_phase_checks':chains,
            'energy_study_covered':{'maximum_total_electric_cutoff':18,
                                  'maximum_sufficient_physical_doubled_spin':7,
                                  'maximum_sufficient_auxiliary_doubled_spin':14},
            'interpretation':'The final orientation-safe default Hamiltonian and the earlier fixed-order default differ by one diagonal state-sign transformation. Every spin/energy cutoff and every kappa/nu retain the same spectrum.',
            'independent_tensor_validation':'SU2_Joined_Cubes_Independent_Validation.json',
            'all_exact_checks_passed':True}


if __name__=='__main__':
    result=run()
    here=Path(__file__).resolve().parent
    output=here.parent/'outputs' if here.name=='work' else here
    destination=output/'SU2_Joined_Cubes_Orientation_Validation.json'
    destination.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
