import numpy as np
from math import comb,prod,factorial
from itertools import product
from fractions import Fraction

class Multiplet:
    '''
    Attributes:
        spin (flt): Total U-spin u of representation
        space (str): Either 'in', 'h', or 'out' to indicate vector space in which
                     the Multiplet lives
        nt_strs (list): Contains strings of '0's and '1's for each possible m
                        QM number, where -u <= m <= u
    '''

    def __init__(self,spin,space):
        self.spin = spin
        self.space = space # necessary to indicate in and h???
        self.nt_strs = self.form_str()

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
    
    def bin_to_dec(self,bin_str):
        return int(''.join(bin_str),2)
    
    def find_q(self,out_lt):
        out_str_lt = []
        for i in out_lt:
            out_str_lt.append(self.ntuple[i])
        q = (-1)**(''.join(out_str_lt).count('0'))
        return q

    def conjugate(self):
        conj_ntuple = []
        for irrep in self.ntuple:
            conj_irrep = ''.join(['1' if s == '0' else '0' for s in irrep])
            conj_ntuple.append(conj_irrep)
        conj_ntuple = [''.join(sorted(irrep)) for irrep in conj_ntuple]
        return conj_ntuple
    
    def __eq__(self,other):
        return self.ntuple == other.ntuple

class AmplitudePair(Amplitude):
    '''
    Attributes:
        gen_coord (list): Node (coords are ints) corresponding to the AmplitudePair
        y_coord (NumPy array): y coordinate (coords are ints) used to find mu
        mu (list): Written as [a,b], where the actual mu factor is sqrt(a)*b
    '''

    def __init__(self,amp,out_lt,inputs):
        super().__init__(amp.ntuple,out_lt)
        self.gen_coord = self.map_to_gen_coord()
        self.y_coord = self.map_to_y_coord()
        self.mu = self.calc_mu()
        self.indices = self.find_indices()
        self.processes = self.find_processes(inputs)
    
    def map_to_gen_coord(self):
        gen_coord = []
        for i,nt_str in enumerate(self.ntuple[1:]):
            n_minus = nt_str.count('0')
            if n_minus != 0:
                gen_coord += [i+1]*n_minus
        return gen_coord
    
    def map_to_y_coord(self):
        y_coord = []
        for nt_str in self.ntuple[1:]:
            y_coord.append(nt_str.count('0'))
        return np.array(y_coord)
    
    def calc_mu(self): # simplify with y_coord? -- might be able to delete y_coord entirely
        mu_j_sqrts = []
        mu_j_facts = []
        for nt_str in self.ntuple[1:]:
            y_j = nt_str.count('0')
            bin_coeff = comb(len(nt_str),y_j)
            mu_j_sqrts.append(bin_coeff)
            mu_j_facts.append(factorial(y_j))
        mu = [prod(mu_j_sqrts),prod(mu_j_facts)]
        return mu
    
    def find_indices(self):
        return [self.number,self.bin_to_dec(self.conjugate())]
    
    def particle_indices(self,inputs,ntuple):
        u_index_lt = np.array([[len(nt_str)/2,nt_str.count('1')+1] for nt_str in ntuple])

        part_indices = []
        for i,input in enumerate(inputs):
            state_indices = []
            for u in input:
                u_matches = np.where(u_index_lt[:,0] == u)[0]
                index = 1

                if np.size(u_matches) > 0:
                    j = u_matches[0]
                    index = u_index_lt[j,1]
                    u_index_lt = np.delete(u_index_lt,j,0)
                
                if i == 2:
                    index = 2*u+2-index # negates for in state and H
                
                state_indices.append(int(index))
            part_indices.append(state_indices)
        return part_indices
    
    def find_processes(self,inputs):
        if len(inputs) > 0:
            return [self.particle_indices(inputs,self.ntuple),
                    self.particle_indices(inputs,self.conjugate())]
        else:
            return "n/a"
    
    def symmetrize(self,inputs):
        nt_str1 = self.ntuple[1]
        m0 = -1/2
        u0 = 1/2
        u1 = len(nt_str1)/2
        m1 = -u1+nt_str1.count('1')
        u = u0+u1
        m = m0+m1
        self.cg = [[u0,m0],[u1,m1],[u,m]]

        self.ntuple[0:2] = [''.join(sorted(self.ntuple[0:2]))]
        self.number = self.bin_to_dec(self.ntuple)
        self.indices = self.find_indices()
        self.processes = self.find_processes(inputs)

class System:
    '''
    Attributes:
        aux (bool): True only if the System contains an auxiliary doublet
        reps (list): Contains Multiplets arranged from lowest to highest u. For 
                     Systems without doublets, contains an auxiliary u = 1/2 Multiplet.  ## turn into np array
        out_lt (list): Contains indices of 'out' space Multiplets in self.reps ### can prob just use an np array mask
        n (int): Number of doublets needed to construct the System
        p (int): p factor whose parity is used to relate U-spin conjugate pairs #### enforce int?
        amp_pairs (list): Contains all AmplitudePairs (representing a-/s-type
                          amplitudes) for the System
    '''

    def __init__(self,reps,inputs=[]):
        self.inputs = inputs
        self.aux = False
        self.reps = self.sort_add_aux(reps)
        self.out_lt = self.out_indices()
        self.n = self.find_n()
        self.p = self.find_p()
        self.amp_pairs = self.form_amp_pairs()
        self.self_conj_amps = self.find_self_conj_amps()
    
    def sort_add_aux(self,reps):
        sorted_reps = sorted(reps,key=lambda x: x.spin)

        if 1/2 not in [rep.spin for rep in sorted_reps]:
            self.aux = True
            space = sorted_reps[0].space
            sorted_reps[0] = Multiplet(sorted_reps[0].spin-1/2,space)
            sorted_reps.insert(0,Multiplet(1/2,space))
        
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
        p = int(2*spin_sum-self.n/2)
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
            amp_pairs.append(AmplitudePair(amp,self.out_lt,self.inputs))
        return amp_pairs
    
    def extract_amplitudes(self): # mathematica amplitudes
        math_amps = []
        for amp in self.amp_pairs:
            ntuple = [str(coord) for coord in amp.ntuple]
            ntuple = '('+(','.join(ntuple).replace('0','-').replace('1','+'))+')'
            node = [str(coord) for coord in amp.gen_coord]
            node = '('+(','.join(node))+')'
            p = (-1)**self.p

            math_amps.append([amp.processes, amp.indices, ntuple, node, amp.q, p,
                              amp.mu, amp.cg])
        return math_amps
    
    def symmetrize(self):
        for amp in self.amp_pairs:
            amp.symmetrize(self.inputs)
        self.aux = False
        return self
    
    def find_self_conj_amps(self):
        self_conj_amps = []
        for i,amp in enumerate(self.amp_pairs):
            if amp.ntuple == amp.conjugate():
                self_conj_amps.append(i)
        self.self_conj_amps = self_conj_amps
        return self_conj_amps
    
    def find_dup_pairs(self):
        dup_pairs = []

        if len(self.amp_pairs) > 1:
            conj_amps = []
            for amp in self.amp_pairs:
                conj_amps.append(amp.conjugate())

            for i,amp in enumerate(self.amp_pairs[1:],start=1):
                for j,conj_amp in enumerate(conj_amps[:i]):
                    if amp.ntuple == conj_amp:
                        amp.number = self.amp_pairs[j].number
                        dup_pairs.append([j,i])
        return dup_pairs

class SumRule: ############################## might be able to modify into a set of functions??
    '''
    Attributes:
        b (int): Dimension of subspace and breaking order of SumRule
        subspace (list): Contains nodes (strs)
        amp_pairs (list): Contains AmplitudePairs corresponding to nodes (can
                          include duplicates)
        M_values (list): Contains multiplicative factors (ints) encoding number
                         of times a node repeats in a SumRule. Ordered according
                         to System.amp_pairs.
    '''

    def __init__(self,b,subspace,lattice):
        #self.b = b
        self.subspace = self.remove_zeros(lattice,subspace)
        self.amp_pairs = [lattice[node] for node in self.subspace]
        self.M_values = self.unit(lattice) if (b == 0) else self.count_amps(lattice)
    
    def remove_zeros(self,lattice,subspace):
        new_subspace = []
        for node in subspace:
            if node in lattice.keys():
                new_subspace.append(node)
        return new_subspace

    def unit(self,lattice):
        all_amps = list(lattice.values()) # technically not ordered though...
        M_values = [0]*len(all_amps)
        i = all_amps.index(self.amp_pairs[0])
        M_values[i] = 1
        return M_values

    def count_amps(self,lattice):
        sr_amps = self.amp_pairs
        sys_amps = list(lattice.values()) # technically not ordered though...
        
        M_values = []
        for amp in sys_amps:
            M_values.append(sr_amps.count(amp))
        return M_values

class Lattice: ###### modify sum_rules attribute (list of np matrices) if using .M_values method
    '''
    Attributes:
        d (int): Lattice dimension
        l (int): Length of each dimension
        lattice (dict): Contains all non-zero nodes (strs) as keys. Values are
                        AmplitudePairs corresponding to nodes.
        sum_rules (list): Contains lists of SumRules by order of breaking
    '''

    def __init__(self,system):
        self.d = int((system.n/2)-1)
        self.l = int(len(system.reps)-1)
        self.lattice = self.nodes(system)
        self.sum_rules = self.find_sum_rules()
    
    def nodes(self,system):
        lattice = {}
        for amp in system.amp_pairs:
            lattice[','.join(str(x) for x in amp.gen_coord)] = amp
        return lattice

    def find_sum_rules(self):
        d = self.d
        l = self.l
        latt = self.lattice
        sum_rules = []
        
        # b = 0 case
        sum_rules.append(np.array([SumRule(0,[node],latt).M_values for node in latt]))

        # b >= 1 cases, full lattice in NumPy array form
        if d > 0:
            dtype = np.uint8 if l <= 255 else np.uint64
            np_latt = np.indices((l,)*d,dtype=dtype).reshape(d,-1).T + 1
            np_latt = np.sort(np_latt,axis=1)

            for b in range(1,d+1):
                np_latt = np_latt.reshape(l**(d-b),l**b,d) # each element is a list of coordinates within a shared b-dim subspace
                np_b_subspaces = np.unique(np_latt,axis=0)

                b_srs = []
                for subspace in np_b_subspaces:
                    subspace = [','.join(map(str,coord)) for coord in subspace]
                    b_sr = SumRule(b,subspace,latt).M_values
                    if any(M != 0 for M in b_sr):
                        b_srs.append(b_sr)
                sum_rules.append(np.array(b_srs))

        return sum_rules

def define_system(inputs,phys=False):
    '''
    Parameters:
        inputs (list): Contains lists of the total U-spins (flts) in the incoming
                       state, Hamiltonian, and outgoing state.

    Returns:
        system (System): System formed from non-trivial (u > 0) reps in input
    '''

    sys_reps = []
    states = ['in','h','out'] ###################### might not need in and h, just make a boolean for out
    for i,input in enumerate(inputs):
        reps = [rep for rep in input if rep != 0]
        reps = [Multiplet(rep,states[i]) for rep in reps]
        sys_reps += reps
    
    n_doublets = sum([2*rep.spin for rep in sys_reps])
    if n_doublets % 2 != 0:
        raise ValueError('Invalid process. Number of would-be doublets is '+str(n_doublets)+'. Please enter a system with an even number of doublets.')
    
    system = System(sys_reps,inputs) if phys else System(sys_reps)
    return system

def generate_srs(system):
    lattice = Lattice(system)
    sum_rules = lattice.sum_rules # M values for aux/original system, when multiplied by mu factors they become SRs

    if system.aux:
        system.symmetrize()
    
    # Set identically 0 amplitudes to 0 in SRs, edits for duplicate amplitudes from symmetrization
    self_conj_amps = system.find_self_conj_amps() #column indices of self conj amps within amp_pairs list
    dup_pairs = np.array(system.find_dup_pairs())
    for i,sr_mat in enumerate(sum_rules): #vanishing amplitudes are a-type for even p, s-type for odd p
        if (i % 2) == (system.p % 2):
            sr_mat[:,self_conj_amps] = 0
            if len(dup_pairs) > 0:
                sr_mat[:,dup_pairs[:,1]] = -sr_mat[:,dup_pairs[:,1]] #negated dup amplitudes are a-type for even p, etc

    sum_rules = [sr_mat.tolist() for sr_mat in sum_rules]
    return system, sum_rules, (dup_pairs+1).tolist()

def remove_dups(system,dup_pairs):
    if len(dup_pairs) > 0:
        dup_amps = np.array(dup_pairs)[:,1]-1
        for i in sorted(dup_amps, reverse=True):
            del system.amp_pairs[i]
    return system

def extract_sys(system):
    n = int(system.n)
    irreps = [map(Fraction,state) for state in system.inputs]
    return n, irreps