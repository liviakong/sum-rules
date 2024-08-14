from sr_classes import *
from fractions import Fraction

def input_reps(text,state):
    input_text = input(text).split(',')
    reps = [float(Fraction(x)) for x in input_text] if input_text != [''] else []
    reps = [x for x in reps if x != 0]
    reps = [Multiplet(rep,state) for rep in reps]
    return reps

#print('Input representations as comma-separated lists (ex. 1/2,1).')
in_reps = input_reps('In state: ','in')
#h_rep = input_reps('Hamiltonian: ','h')
out_reps = input_reps('Out state: ','out')
reps = in_reps+out_reps

if sum([2*rep.spin for rep in reps]) % 2 != 0:
    raise ValueError('Invalid process. Please enter a system with non-zero amplitudes.')

spin_sys = System(reps)
lattice = Lattice(spin_sys)

#print("amp pairs",spin_sys.amp_pairs)
print("lattice",lattice.lattice)
print("sum rules",lattice.sum_rules)
sr = lattice.sum_rules
for i in range(len(sr)):
    print("sum rule count",i)
    print("breaking",sr[i].b)
    print("subspace",sr[i].subspace)
    #for j in range(len(sr[i].amp_pairs)):
    #    print("amp pair number",sr[i].amp_pairs[j].number)