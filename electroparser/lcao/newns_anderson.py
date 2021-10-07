
from gpaw import GPAW
from gpaw import restart
from gpaw.utilities.dos import RestartLCAODOS, fold
import os
from gpaw.lcao.tools import get_lcao_hamiltonian
from pprint import pprint
import numpy as np
from ase.units import Hartree


class NewnsAndersonLCAO:

    def __init__(self, homedir, adsorbate_index, metal_index, spin_index=0, k_index=0):

        self.homedir = homedir
        self.adsorbate_index = adsorbate_index
        self.metal_index = metal_index
        self.spin_index = spin_index
        self.kindex = k_index

        self.H_skMM = [] # spin, kpoint, MxM matrix
        self.S_kMM = [] # kpoint, MxM matrix
        self.Ef = 0.0 # Fermi level
        self.H_MM = []
        self.S_MM = []
        self.adsorbate_bfunc_index = []
        self.metal_bfunc_index = []
        self.H = [] # subdiagonalised Hamiltonian
        self.S = [] # subdiagonalised overlap matrix
        self.eigen_ads = [] # eigenvalue of adsorbate
        self.eigen_metal = [] # eigenvalue of metal
        self.delta = [] # chemisorption function
        self.edos = [] # site projected dos energies
        self.w_ads_dos = [] # site projected dos weights

    def _get_hamiltonian_lcao(self):
        """
        Get the Hamiltonian in the form 
        that is in the basis for the LCAO orbitals
        """
        ## read in the out.gpw file
        atoms, calc = restart(os.path.join(self.homedir, 'out.gpw'))
        ## Get the Hamiltonian and overlap matrix 
        self.H_skMM, self.S_kMM = get_lcao_hamiltonian(calc)
        self.Ef = calc.get_fermi_level()

        self.H_MM = self.H_skMM[self.spin_index, self.kindex]
        self.S_MM = self.S_kMM[self.kindex]

        dos = RestartLCAODOS(calc)
        self.adsorbate_bfunc_index = dos.get_atom_indices(self.adsorbate_index)
        self.metal_bfunc_index = dos.get_atom_indices(self.metal_index)
        energies, weights_ads = dos.get_atomic_subspace_pdos(self.adsorbate_index)
        self.edos, self.w_ads_dos = fold(energies * Hartree, weights_ads, 2000, 0.1)

    def _subdiagonalise(self, indices, H_i, S_i):
        """
        Performs the subdiagonalisation of the matrix
        """
        h = self.H_MM.take(indices, 0).take(indices, 1)
        s = self.S_MM.take(indices, 0).take(indices, 1)

        # subdiadonalise by performing
        # H  = S epsilon
        # solver AX = B theen get the eigenvectors and values
        eigenval, eigenvec = np.linalg.eig(np.linalg.solve(s, h))
        # Normalise the eigenvectors taking into account the overlap matrix
        for col in eigenvec.T:
            col /= np.sqrt(np.dot(col.conj(), np.dot(s, col)))
        # Transformation matrix to convert the Hamiltonian
        t_matrix = np.identity(H_i.shape[0], dtype=complex)
        for i in range(len(indices)):
            for j in range(len(indices)):
                t_matrix[indices[i],indices[j]] = eigenvec[i,j]
        
        # Unitary transform to get the rearranged Hamiltonian and overlap
        H_r = np.dot(np.transpose(np.conj(t_matrix)), np.dot(H_i,t_matrix))
        S_r = np.dot(np.transpose(np.conj(t_matrix)), np.dot(S_i,t_matrix))

        return H_r, S_r, eigenval

    def get_subdiagonalised_hamiltonian(self):
        self._get_hamiltonian_lcao()
        H_r, S_r, eigenval_mol = self._subdiagonalise(self.adsorbate_bfunc_index, self.H_MM, self.S_MM)
        H_r, S_r, eigenval_metal = self._subdiagonalise(self.metal_bfunc_index, H_r, S_r )

        self.H = H_r
        self.S = S_r
        self.eigen_ads = eigenval_mol
        self.eigen_metal = eigenval_metal
    
    def get_chemisorption_function(self):
        """
        For each adsorbate state, get the chemistorption function 
        Delta = Sigma_{k} = | e_k * s_{ak} - v_{ak} | ^2
        """
        delta = np.zeros([len(self.adsorbate_bfunc_index), len(self.metal_bfunc_index)])
        for i, a in enumerate(self.adsorbate_bfunc_index):
            for j, k in enumerate(self.metal_bfunc_index):
                eps_k = self.eigen_metal[j]
                s_ak = self.S[a,k]
                v_ak = self.H[a,k]
                delta[i,j] = np.pi * np.abs(eps_k * s_ak - v_ak)**2 
        
        self.delta = delta














    
