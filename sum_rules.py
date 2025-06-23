#from sr_classes import *
#from fractions import Fraction
from itertools import product,permutations
from math import comb,prod,sqrt,factorial

class Multiplet:
    '''
    Attributes:
        spin (flt): Total U-spin u of representation
        space (str): Either 'in', 'h', or 'out' to indicate vector space in which
                     the Multiplet lives
        particles (list): ############################################## complete
        nt_strs (list): Contains strings of '0's and '1's for each possible m
                        QM number, where -u <= m <= u
    '''
    def __init__(self,spin,space):
        self.spin = spin
        self.space = space # necessary to indicate in and h???
        self.particles = self.form_particles() # incomplete
        self.nt_strs = self.form_str()

    def form_particles(self): # complete, eventually replace lt of strings with lt of particles with respective m qm numbers
        self

    def form_str(self):
        spin = self.spin
        min_m = spin*-1
        m_lt = []
        for i in range(int(2*spin)+1):
            m_lt.append(min_m+i)
        
        nt_strs = []
        for m in m_lt:
            nt_str = '0'*int(spin-m)+'1'*int(spin+m)
            nt_strs.append(nt_str)
        return nt_strs

class Amplitude:
    '''
    Attributes:
        ntuple (list): ntuple of strings of '0's and '1's
        number (int): i index used to refer to a particular Amplitude
        q (int): (-1)^q_i factor (i.e. q = -1 or +1 only) accounting for movements
                 of reps between states
        cg (list): Contains Mathematica ClebschGordan function inputs in the form
                   [[u1,m1],[u2,m2],[u,m]]. Set to evaluate to 1 by default and
                   changes only during symmetrization.
    '''
    def __init__(self,ntuple,out_lt):
        self.ntuple = ntuple
        self.number = self.bin_to_dec(self.ntuple)
        self.q = self.find_q(out_lt)
        self.cg = [[1/2,1/2],[1/2,1/2],[1,1]]
    
    def bin_to_dec(self,st):
        return int(''.join(st),2)
    
    def find_q(self,out_lt):
        out_str_lt = []
        for i in out_lt:
            out_str_lt.append(self.ntuple[i])
        q = (-1)**(''.join(out_str_lt).count('0'))
        return q

    def symmetrize(self):
        nt_str1 = self.ntuple[1]
        m0 = -1/2
        u0 = 1/2 # function can be made more general, doesn't need to be hardcoded for u = 1/2
        u1 = len(nt_str1)/2
        m1 = -u1+nt_str1.count('1')
        u = u0+u1
        m = m0+m1
        self.cg = [[u0,m0],[u1,m1],[u,m]]

        self.ntuple[0:2] = [''.join(sorted(self.ntuple[0:2]))]
        self.number = self.bin_to_dec(self.ntuple)

    def conjugate(self):
        conj_ntuple = []
        for irrep in self.ntuple:
            conj_irrep = ''.join(['1' if s == '0' else '0' for s in irrep])
            conj_ntuple.append(conj_irrep)
        conj_ntuple = [''.join(sorted(irrep)) for irrep in conj_ntuple]
        return conj_ntuple
    
    def __eq__(self,other):
        return self.ntuple == other.ntuple

class AmplitudePair(Amplitude): #combine subclass with parent class?
    '''
    Attributes:
        gen_coords (list): Contains the node (list of ints) corresponding to the
                           AmplitudePair along with all of its permutations
        mu (list): Written as [a,b], where the actual mu factor is sqrt(a)*b
    '''
    def __init__(self,amp,out_lt):
        super().__init__(amp.ntuple,out_lt)
        self.gen_coords = self.map_to_gen_coords()
        self.mu = self.calc_mu()
    
    def map_to_gen_coords(self):
        gen_coord = []
        for i,nt_str in enumerate(self.ntuple[1:]):
            n_minus = nt_str.count('0')
            if n_minus != 0:
                gen_coord += [i+1]*n_minus
        
        gen_coords = [list(i) for i in (permutations(gen_coord))]
        return gen_coords
    
    def calc_mu(self):
        mu_j_sqrts = []
        mu_j_facts = []
        for nt_str in self.ntuple:
            y_j = nt_str.count('0')
            bin_coeff = comb(len(nt_str),y_j)
            mu_j_sqrts.append(bin_coeff)
            mu_j_facts.append(factorial(y_j))
        mu = [prod(mu_j_sqrts),prod(mu_j_facts)]
        return mu

class System:
    '''
    Attributes:
        aux (bool): True only if the System contains an auxiliary doublet
        reps (list): Contains Multiplets arranged from lowest to highest u. For 
                     Systems without doublets, contains an auxiliary u = 1/2 Multiplet.  ## turn into np array
        out_lt (list): Contains indices of 'out' space Multiplets in self.reps ### can prob just use an np array mask
        n (int): Number of doublets needed to construct the System ###### enforce intger?
        p (int): p factor whose parity is used to relate U-spin conjugate pairs #### enforce int?
        amp_pairs (list): Contains all AmplitudePairs (representing a-/s-type
                          amplitudes) for the System
    '''
    def __init__(self,reps):
        self.aux = False
        self.reps = self.sort_aux(reps)
        self.out_lt = self.out_indices()
        self.n = self.find_n()
        self.p = self.find_p()
        self.amp_pairs = self.form_amp_pairs()
    
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

    def find_n(self):
        n = sum([2*rep.spin for rep in self.reps])
        return n

    def find_p(self):
        spin_sum = 0
        for i in self.out_lt:
            spin_sum += self.reps[i].spin
        p = 2*spin_sum-self.n/2
        return p
    
    def form_amp_pairs(self): # only forms i amplitudes to represent a-/s-type amplitudes
        ntuples_lt = []
        for rep in self.reps:
            ntuples_lt.append(rep.nt_strs)
        
        ntuples = [list(i) for i in product(*ntuples_lt)]
        amps = []
        for ntuple in ntuples:
            bin_str = ''.join(ntuple)
            if (bin_str[0] == '0') and (bin_str.count('0') == bin_str.count('1')):
                amps.append(Amplitude(ntuple,self.out_lt))
        
        amp_pairs = []
        for amp in amps:
            amp_pairs.append(AmplitudePair(amp,self.out_lt))
        return amp_pairs
    
    def extract_amplitudes(self):
        math_amps = [] #mathematica amplitudes
        for amp in self.amp_pairs:
            math_amps.append([amp.number, amp.ntuple, amp.gen_coords[0],
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
                    if amp.ntuple == conj_amp:
                        amp.number = self.amp_pairs[j].number
        return self.extract_amplitudes()
    
    def find_self_conjugates(self):
        self_conj_amps = []
        for i,amp in enumerate(self.amp_pairs):
            if amp.ntuple == amp.conjugate():
                self_conj_amps.append(i+1)
        return self_conj_amps

class SumRule:
    '''
    Attributes:
        b (int): Dimension of subspace and breaking order of SumRule
        subspace (list): Contains nodes (strs)
        amp_pairs (list): Contains AmplitudePairs corresponding to nodes (can
                          include duplicates)
    '''
    def __init__(self,subspace,b,lattice):
        self.b = b
        self.subspace = subspace
        self.amp_pairs = [lattice[x] for x in subspace]

def match_subspace(node1,node2,n):
    '''
    checks if node1 and node2 are in the same b-dim subspace (i.e. that they share n = d-b coordinates)
    '''
    match_list = []
    for i in range(n):
        match_list.append(node1.split(',')[i] == node2.split(',')[i])
    match = all(match_list)
    return match

def srs_from_lattice(l_copy,b,d,l):
    '''
    forms b-dim subspace, then forms sumrule obj
    modifies l_copy to remove subspace
    '''
    n = d-b # all nodes that share d-b coordinates form a b-dim subspace, should prob use diff name
    subspace = []
    subspace.append(l_copy.pop(0))
    if len(l_copy) > 0:
        subspace += [node for node in l_copy if match_subspace(subspace[0],node,n)]
        l_copy = [node for node in l_copy if node not in subspace]
    sum_rule = SumRule(subspace,b,l)
    return l_copy,sum_rule

class Lattice:
    '''
    Attributes:
        d (int): Lattice dimension
        lattice (dict): Contains all non-zero nodes (strs) as keys. Values are
                        AmplitudePairs corresponding to nodes.
        sum_rules (list): Contains **********************
    '''
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
            lattice_copy = [x for x in l] # list of all nodes, could just use .keys()
            while len(lattice_copy) > 0:
                lattice_copy,sum_rule = srs_from_lattice(lattice_copy,b,d,l)
                sr_b.append(sum_rule)
            sum_rules.append(sr_b)
        return sum_rules

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

def form_reps(inputs):
    '''
    Inputs:
        inputs (list): Contains inReps, hRep, and outReps, each of which is a list
                       of the total U-spins (flts) in the incoming state, Hamiltonian,
                       and outgoing state.

    Returns:
        sys_reps (list): Contains Multiplets corresponding to each non-trivial
                         (u > 0) input representation.
    '''

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