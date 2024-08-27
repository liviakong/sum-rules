from sr_classes import *
from fractions import Fraction

def input_reps(text,state):
    input_text = input(text).split(',')
    reps = [float(Fraction(x)) for x in input_text] if input_text != [''] else []
    reps = [x for x in reps if x != 0]
    reps = [Multiplet(rep,state) for rep in reps]
    return reps

print('Input representations as comma-separated lists (ex. 1/2,1).')
in_reps = input_reps('In state: ','in')
h_rep = input_reps('Hamiltonian: ','h')
out_reps = input_reps('Out state: ','out')
reps = in_reps+h_rep+out_reps

if sum([2*rep.spin for rep in reps]) % 2 != 0:
    raise ValueError('Invalid process. Please enter a system with non-zero amplitudes.')

spin_sys = System(reps)
lattice = Lattice(spin_sys)

sr = lattice.sum_rules

def print_sr_dec():
    for i in range(len(sr)):
        print("Breaking order: "+str(sr[i].b))
        if sr[i].b % 2 == 0:
            amp_comb = "a"
        else:
            amp_comb = "s"
        str_lt = []
        for j in range(len(sr[i].amp_pairs)):
            str_lt.append(str(round(sr[i].amp_pairs[j].mu,5))+amp_comb+str(sr[i].amp_pairs[j].number))
        print(*str_lt,sep=' + ',end=' = 0\n')

def print_sr_node():
    for i in range(len(sr)):
        print("Breaking order: "+str(sr[i].b))
        if sr[i].b % 2 == 0:
            amp_comb = "a"
        else:
            amp_comb = "s"
        str_lt = []
        for j in range(len(sr[i].subspace)):
            str_lt.append(str(round(sr[i].amp_pairs[j].mu,5))+amp_comb+'['+str(sr[i].subspace[j])+']')
        print(*str_lt,sep=' + ',end=' = 0\n')

print_sr_dec()