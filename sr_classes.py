from itertools import product,permutations
from math import comb,prod,sqrt,factorial

# def conjugate(ntuple):
#     conj_ntuple = '0'*ntuple.count('1')+'1'*ntuple.count('0')
#     return conj_ntuple

class Multiplet:
    def __init__(self,spin,space):
        self.spin = spin
        self.space = space
        self.ntuples = self.form_ntuples()

    def form_ntuples(self):
        spin = self.spin
        min_m = spin*-1
        m_lt = []
        for i in range(int(2*spin)+1):
            m_lt.append(min_m+i)
        
        ntuples = []
        for m in m_lt:
            ntuple_str = '0'*int(spin-m)+'1'*int(spin+m)
            ntuples.append(ntuple_str)
        return ntuples

class Amplitude:
    def __init__(self,state,out_lt):
        self.gen_ntuple = state
        self.number = self.bin_to_dec()
        self.q = self.find_q(out_lt)
    
    def bin_to_dec(self):
        return int(''.join(self.gen_ntuple),2)
    
    def find_q(self,out_lt):
        out_str_lt = []
        for i in out_lt:
            out_str_lt.append(self.gen_ntuple[i])
        q = (-1)**(''.join(out_str_lt).count('0'))
        return q

    #def rearrange(self): also, need check equal to conjugate

    #def __eq__(self,other):
    #    return self.gen_ntuple == other.gen_ntuple
    
class AmplitudePair(Amplitude):
    def __init__(self,amp,out_lt):
        super().__init__(amp.gen_ntuple,out_lt)
        self.gen_coord = self.map_to_gen_coord() #can turn into helper function
        self.gen_coords = self.all_gen_coords()
        self.y_coord = self.map_to_y_coord() #can turn into helper function
        self.mu = self.calc_mu()
    
    def map_to_gen_coord(self):
        gen_coord = []
        for i,ntuple in enumerate(self.gen_ntuple[1:]):
            n_minus = ntuple.count('0')
            if n_minus != 0:
                gen_coord += [i+1]*n_minus
        return gen_coord
    
    def all_gen_coords(self):
        gen_coords = [list(i) for i in (permutations(self.gen_coord))]
        return gen_coords
    
    def map_to_y_coord(self):
        y_coord = []
        for ntuple in self.gen_ntuple[1:]:
            y_coord.append(ntuple.count('0'))
        return y_coord
    
    def calc_mu(self):
        mu_js = []
        for ntuple in self.gen_ntuple:
            y_j = ntuple.count('0')
            bin_coeff = comb(len(ntuple),y_j)
            mu_j = sqrt(bin_coeff)*factorial(y_j)
            mu_js.append(mu_j)
        mu = prod(mu_js)
        return mu

class System:
    def __init__(self,reps):
        self.reps = sorted(reps,key=lambda x: x.spin)
        self.out_lt = self.out_indices()
        self.amps = self.form_amplitudes() #necessary?
        self.n = self.find_n()
        self.p = self.find_p()
        self.amp_pairs = self.form_amplitude_pairs()
    
    def out_indices(self):
        out_lt = []
        for i,rep in enumerate(self.reps):
            if rep.space == 'out':
                out_lt.append(i)
        return out_lt

    def form_amplitudes(self):
        ntuples_lt = []
        for rep in self.reps:
            ntuples_lt.append(rep.ntuples)
        
        states = [list(i) for i in product(*ntuples_lt)]
        amplitudes = []
        for state in states:
            bin_str = ''.join(state)
            if (bin_str[0] == '0') and (bin_str.count('0') == bin_str.count('1')):
                amplitudes.append(Amplitude(state,self.out_lt))
        return amplitudes
    
    def find_n(self):
        n = sum([2*rep.spin for rep in self.reps])
        return n

    def find_p(self):
        spin_sum = 0
        for i in self.out_lt:
            spin_sum += self.reps[i].spin
        p = 2*spin_sum-self.n/2
        return p
    
    def form_amplitude_pairs(self):
        amplitude_pairs = []
        for amp in self.amps:
            amplitude_pairs.append(AmplitudePair(amp,self.out_lt))
        return amplitude_pairs

class SumRule:
    def __init__(self,subspace,b,l):
        self.b = b
        self.subspace = subspace
        self.amp_pairs = [l[x] for x in subspace]

def match_subspace(node1,node2,n):
    match_list = []
    for i in range(n):
        match_list.append(node1.split(',')[i] == node2.split(',')[i])
    match = all(match_list)
    return match

def sum_rules_from_lattice(l_copy,b,d,l):
    n = d-b
    subspace = []
    subspace.append(l_copy.pop(0))
    if len(l_copy) > 0:
        subspace += [node for node in l_copy if match_subspace(subspace[0],node,n)]
        l_copy = [node for node in l_copy if node not in subspace]
    sum_rule = SumRule(subspace,b,l)
    return l_copy,sum_rule

class Lattice:
    def __init__(self,system):
        self.d = int((system.n/2)-1)
        self.lattice = self.nodes(system)
        self.sum_rules = self.find_sum_rules()
    
    def nodes(self,system):
        lattice = {}
        for amp in system.amp_pairs:
            for coord in amp.gen_coords:
                lattice[','.join(str(x) for x in coord)] = amp
        return lattice

    def find_sum_rules(self):
        d = self.d
        l = self.lattice
        sum_rules = []
        for b in range(0,d+1):
            lattice_copy = [x for x in l]
            while len(lattice_copy) > 0:
                lattice_copy,sum_rule = sum_rules_from_lattice(lattice_copy,b,d,l)
                sum_rules.append(sum_rule)
        return sum_rules

'''
class SumRule:
    def __init__(self,subspace,b,n):
        self.b = b
        self.nodes = self.remove_dups(subspace)
        self.yp_coords = self.find_yp_coords()
        self.mb_factors = self.calc_mb_factors()
    
    def remove_dups(self,subspace):
        if self.b < 2:
            return subspace
        else:
            nodes = []
            for node in subspace:
                node_sorted = [x for x in node]
                node_sorted.sort()
                if node == node_sorted:
                    nodes.append(node)
            return nodes

    def find_yp_coords(self,n):
        yp_coords = []
        for node in self.nodes:
            node_spl = node[n:]
            yp_coord = []
            for i in range(len(node_spl)):
                yp_coord.append(node_spl.count(i))
        return yp_coords
'''