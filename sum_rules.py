#from sr_classes import *
#from fractions import Fraction
from itertools import product,permutations
from math import comb,prod,sqrt,factorial

class Multiplet:
    def __init__(self,spin,space):
        self.spin = spin
        self.space = space
        self.particles = self.form_particles()
        self.ntuples = self.form_ntuples()

    def form_particles(self):
        self

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
        self.number = self.bin_to_dec(self.gen_ntuple)
        self.q = self.find_q(out_lt)
        self.cg = [[1/2,1/2],[1/2,1/2],[1,1]]
    
    def bin_to_dec(self,st):
        return int(''.join(st),2)
    
    def find_q(self,out_lt):
        out_str_lt = []
        for i in out_lt:
            out_str_lt.append(self.gen_ntuple[i])
        q = (-1)**(''.join(out_str_lt).count('0'))
        return q

    def symmetrize(self):
        #nt_str0 = self.gen_ntuple[0]
        nt_str1 = self.gen_ntuple[1]
        #m0 = -1/2 if nt_str0 == '0' else 1/2
        m0 = -1/2 #can be made more general, doesn't need to be hardcoded for 1/2
        u0 = 1/2
        u1 = len(nt_str1)/2
        m1 = -u1+nt_str1.count('1')
        u = u0+u1
        m = m0+m1
        self.cg = [[u0,m0],[u1,m1],[u,m]]

        self.gen_ntuple[0:2] = [''.join(sorted(self.gen_ntuple[0:2]))]
        #self.gen_ntuple[0] = str(''.join(sorted(self.gen_ntuple[0])))
        #self.sort_irreps()
        self.number = self.bin_to_dec(self.gen_ntuple)

    def conjugate(self):
        conj_ntuple = []
        for irrep in self.gen_ntuple:
            conj_irrep = ''.join(['1' if s == '0' else '0' for s in irrep])
            conj_ntuple.append(conj_irrep)
        conj_ntuple = [''.join(sorted(irrep)) for irrep in conj_ntuple]
        return conj_ntuple
    
    def __eq__(self,other):
        return self.gen_ntuple == other.gen_ntuple

class AmplitudePair(Amplitude):
    def __init__(self,amp,out_lt):
        super().__init__(amp.gen_ntuple,out_lt)
        self.gen_coords = self.map_to_gen_coords()
        self.y_coord = self.map_to_y_coord() #can turn into helper function
        self.mu = self.calc_mu() #[a, b], mu = sqrt(a)*b
    
    def map_to_gen_coords(self):
        gen_coord = []
        for i,ntuple in enumerate(self.gen_ntuple[1:]):
            n_minus = ntuple.count('0')
            if n_minus != 0:
                gen_coord += [i+1]*n_minus
        
        gen_coords = [list(i) for i in (permutations(gen_coord))]
        return gen_coords
    
    def map_to_y_coord(self): #necessary?
        y_coord = []
        for ntuple in self.gen_ntuple[1:]:
            y_coord.append(ntuple.count('0'))
        return y_coord
    
    def calc_mu(self):
        mu_j_sqrts = []
        mu_j_facts = []
        for ntuple in self.gen_ntuple:
            y_j = ntuple.count('0')
            bin_coeff = comb(len(ntuple),y_j)
            mu_j_sqrts.append(bin_coeff)
            mu_j_facts.append(factorial(y_j))
        mu = [prod(mu_j_sqrts),prod(mu_j_facts)]
        return mu

class System:
    def __init__(self,reps):
        self.aux = False
        self.reps = self.sort_aux(reps)
        self.out_lt = self.out_indices()
        self.amps = self.form_amplitudes() #combine with amp_pairs?
        self.n = self.find_n()
        self.p = self.find_p()
        self.amp_pairs = self.form_amplitude_pairs()
    
    def sort_aux(self,reps):
        sorted_reps = sorted(reps,key=lambda x: x.spin)

        if 1/2 not in [rep.spin for rep in sorted_reps]:
            self.aux = True
            sorted_reps[0] = Multiplet(sorted_reps[0].spin-1/2,sorted_reps[0].space)
            sorted_reps.insert(0,Multiplet(1/2,"out"))
        
        return sorted_reps
    
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

    def extract_amplitudes(self):
        math_amps = [] #mathematica amplitudes
        for amp in self.amp_pairs:
            math_amps.append([amp.number, amp.gen_ntuple, amp.gen_coords[0],
                              amp.q, amp.mu, amp.cg])
        return math_amps
    
    def symmetrize(self):
        if self.aux:
            for amp in self.amp_pairs:
                amp.symmetrize()
        return self
    
    def find_dups(self):
        if len(self.amp_pairs) > 1:
            conj_amps = []
            for amp in self.amp_pairs:
                conj_amps.append(amp.conjugate())

            for i,amp in enumerate(self.amp_pairs[1:],start=1):
                for j,conj_amp in enumerate(conj_amps[:i]):
                    if amp.gen_ntuple == conj_amp:
                        amp.number = self.amp_pairs[j].number
        return self.extract_amplitudes()
    
    def find_self_conjugates(self):
        self_conj_amps = []
        for i,amp in enumerate(self.amp_pairs):
            if amp.gen_ntuple == amp.conjugate():
                self_conj_amps.append(i+1)
        return self_conj_amps

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
            sr_b = []
            lattice_copy = [x for x in l]
            while len(lattice_copy) > 0:
                lattice_copy,sum_rule = sum_rules_from_lattice(lattice_copy,b,d,l)
                sr_b.append(sum_rule)
            sum_rules.append(sr_b)
        return sum_rules
    
    '''
    def extract_sr(self,n):
        sr_flat = []
        i = 0
        while i < self.d+1:
            sr_b = []
            for sr in self.sum_rules:
                if sr.b == i:
                    amps = []
                    for amp in sr.amp_pairs:
                        amps.append(amp.number)
                    sr_b.append(amps)
            sr_flat.append(sr_b)
            i += 1
        
        math_sr = []
        for b in sr_flat:
            sr_b = []
            for sr in b:
                sr_row = [0]*(2**(n-1))
                for i in range(len(sr_row)):
                    sr_row[i] = sr.count(i)
                sr_b.append(sr_row)
            math_sr.append(sr_b)
        return math_sr
    '''

    def extract_sr(self,amps):
        srs = self.sum_rules
        sr_nums = []
        for b in srs:
            sr_b = []
            for sr in b:
                nums = [amp.number for amp in sr.amp_pairs]
                sr_b.append(nums)
            sr_nums.append(sr_b)
        
        math_sr = []
        for b in sr_nums:
            sr_b = []
            sr_unique = []
            for sr in b:
                if sr not in sr_unique:
                    sr_unique.append(sr)
                    sr_row = [sr.count(amp) for amp in amps]
                    sr_b.append(sr_row)
            math_sr.append(sr_b)
        return math_sr
        '''
        math_sr = []
        for b in self.sum_rules:
            sr_b = []
            sr_unique = []
            for sr in b:
                if sr not in sr_unique:
                    sr_unique.append(sr)
                    sr_row = [sr.amp_pairs.count(amp) for amp in amps]
                    sr_b.append(sr_row)
            math_sr.append(sr_b)
        return math_sr
        '''
        

def form_reps(inputs):
    sys_reps = []
    for i,input in enumerate(inputs):
        if i == 0:
            state = 'in'
        elif i == 1:
            state = 'h'
        else:
            state = 'out'
        
        #name = str(input[0])
        #particles = [str(p) for p in input[1]]

        reps = [rep for rep in input if rep != 0]
        reps = [Multiplet(rep,state) for rep in reps]
        sys_reps += reps

    n_doublets = sum([2*rep.spin for rep in sys_reps])

    if n_doublets % 2 != 0:
        raise ValueError('Invalid process. Number of would-be doublets is '+str(n_doublets)+'. Please enter a system with an even number of doublets.')
    
    return sys_reps

'''
def print_sr_dec():
    for i in range(len(sr)):
        print('Breaking order: '+str(sr[i].b))
        if sr[i].b % 2 == 0:
            amp_comb = 'a'
        else:
            amp_comb = 's'
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

'''