import numpy as np
import Elements
from esp_fitting import *
from scipy.optimize import minimize

class Density:
    def __init__(self, 
        exponent=1.0,
        center=np.array([0.0, 0.0, 0.0]),
        coeff=1.0,
        atom_id=0,
        center_id = 0,
        element=1,
        type='s'):

        self.exp = exponent
        self.center = center
        self.coeff = coeff
        self.atom_id = atom_id
        self.center_id = center_id
        self.element = element
        self.type = type

class DensityFitting():
    def __init__(self, coords_in, nuclei_in, esp_file, charge=0, lone_pairs=False, lone_pair_k=0.005, resp_a=0.0005, resp_b=0.10) -> None:
        self._n_dens = 1
        self._nh_dens = 1
        self._n_atoms = 0
        self._numE = np.sum(nuclei_in) - charge
        self._n_atoms = len(nuclei_in)
        self._nuclei_in = np.copy(nuclei_in)
        self._coords_in = np.copy(coords_in)*1.889725989
        self._atom_map = []
        self._use_lone_pairs = lone_pairs
        self._lone_pair_k = lone_pair_k
        self._resp_a = resp_a
        self._resp_b = resp_b

        #   used to communicate results
        self._current_esp_rrms = None
        self._current_func_eval = None
        self._current_step = 0
        self._current_deriv = None

        #   esp data file
        self._esp_file = esp_file

    
    def create_lone_pairs(self, index, coords, nuclei, lp_dist=0.40):

        new_coords = []

        #   dermine conenctivity for atoms that will have lone pairs
        nuc_i = nuclei[index]
        i = index
        if nuc_i not in [7, 8]: return new_coords
        bonds = []
        for j, nuc_j in enumerate(nuclei):
            if i == j: continue
            dist = np.linalg.norm(coords[i] - coords[j])
            if nuc_j == 1 and dist < 1.2*1.88973:
                bonds.append(j)
            elif nuc_j > 1 and dist < 1.7*1.88973:
                bonds.append(j)
                
        #   nitrogen bisector
        if nuc_i == 7 and len(bonds) == 2:
            r1, r2, r3 = coords[i], coords[bonds[0]], coords[bonds[1]]
            x = r1 - (r2 + r3)*0.5
            x = x/np.linalg.norm(x)

            new_coords.append(r1 + x*lp_dist)

        return new_coords

    def assign_densities(self, coords, nuclei, add_lone_pairs=False, lp_dist=0.40, p_base_only=False):
        print(" Assigning Atomic Densities")
        num_dens = 0
        num_p_gauss = 0
        atoms = [Elements.int2name(x) for x in nuclei]
        density_list = []

        center_id = -1
        for i, atm in enumerate(atoms):
            center_id += 1
            #   add hydrogens
            if atm == 'H':
                for j in range(self._nh_dens):
                    density_list.append(
                        Density(exponent=1.0, center=coords[i], atom_id=i, 
                        element=nuclei[i], type='s', center_id=center_id))
                    num_dens += 1
                
            #   add heavy atoms
            else:
                lone_pairs = []
                if add_lone_pairs:
                    lone_pairs = self.create_lone_pairs(i, coords, nuclei, lp_dist=lp_dist)

                #   centered on nuclei
                for n in range(max(self._n_dens - 1*(len(lone_pairs) > 0), 1)):
                    density_list.append(
                        Density(exponent=1.0, center=coords[i], atom_id=i, 
                        element=nuclei[i], type='s', center_id=center_id))
                    num_dens += 1
                    #   add p-like densities
                    if (p_base_only and coords[i][2] == 0.0):
                        for tp in ['px', 'py', 'pz']:
                            g = Density(exponent=1.0, center=coords[i], atom_id=i, 
                            element=nuclei[i], type=tp, center_id=center_id)
                            density_list.append(g)
                        num_p_gauss += 1

                #   off center lone pairs
                for n in range(len(lone_pairs)):
                    center_id += 1
                    density_list.append(
                        Density(exponent=1.0, center=lone_pairs[n], atom_id=i, 
                        element=0, type='s', center_id=center_id))
                    num_dens += 1

        print("\t Total number of Densities:  ", num_dens + num_p_gauss*3)

        if p_base_only:
            print("\t Number of S-like Gaussians: ", num_dens)
            print("\t Number of P-like Gaussians: ", num_p_gauss*3)
            print("\t Requested Base only P's:    ", str(p_base_only))
            print("\t\t Note: Atoms must have z=0.0 \n\t\t to count as a base atom")
            print("\t \t Note: P's count as 3 gaussians")

        print(' -----------------------------------------------------------------------')
        print("                     ***  Density Centers ***")
        print(" {:>6s}  {:>6s}  {:>6s}  {:>3s}  {:>4s}  {:>10}  {:>10}  {:>10}".format('N', 'AtomID', 'CentID', 'Elm', 'Type', 'X (au)', 'Z (au)', 'Z (au)'))
        print(' -----------------------------------------------------------------------')
        for n, g in enumerate(density_list):
            print(" {:6d}  {:6d}  {:6d}  {:3s}  {:4s}  {:10.4f}  {:10.4f}  {:10.4f}".format(n+1, g.atom_id, g.center_id, Elements.int2name(g.element), g.type, g.center[0], g.center[1], g.center[2]))
        print(' -----------------------------------------------------------------------')
        return density_list

    def get_SLSQP_constrains(self, all_nuclei, guess, center_ids):
        '''
            Common constraints for SLSQP fitting, including hydrogen 
            exponentents, core populations, and valance population
            equalization. At minimum, it returns one constraint for the
            total number of electrons.

            Parameters
            ----------
            all_nuclei : list[int] or np.array(int)
                array of all nuclei atomic numbers.
            guess : list(float) or np.array(float)
                initial guess to used in SLSQP fitting. This used to determine
                which exponents are "core" exponents.
            center_idx: list[int] or np.array(int)
                list of integer ID's for each atom type

            Returns
            -------
            list of dictionaries:
                each dictionary has the keys 'type', 'fun', 'jac', and 'args'

        '''
        constraints = []
        dim = int(len(guess)/2)
        guess_exp = guess[dim:]
        
        #   extract core and hydrogen types
        core_idx_list = []
        hydrogen_list = []
        for id in set(center_ids):
            exp_list = [guess_exp[n] for n in range(dim) if center_ids[n] == id]
            exp_idx = [n for n in range(dim) if center_ids[n] == id]
            max_idx = exp_idx[np.argmax(exp_list)]
            if all_nuclei[max_idx] == 1:
                hydrogen_list.append(max_idx)
            elif len(exp_list) > 1:
                core_idx_list.append(max_idx)

        #   find unique types of constraints
        signatures = []
        for n in range(dim):
            if n in hydrogen_list:
                signatures.append((all_nuclei[n], 'hydro'))
            elif n in core_idx_list:
                signatures.append((all_nuclei[n], 'core'))
            else:
                signatures.append((all_nuclei[n], 'val'))

        #   add constraints for each unique type
        unique_sigs = set(signatures)
        for uniq in unique_sigs:
            if uniq[1] in ['hydro', 'core']:
            #if uniq[1] in []:
                base_idx = None
                for n, sig in enumerate(signatures):
                    if sig == uniq:
                        if base_idx is None:
                            base_idx = n
                        else:
                            jac_vec = np.zeros(dim*2)
                            jac_vec[base_idx + dim] = 1
                            jac_vec[n + dim] = -1
                            constraints.append({
                                'type': 'eq',
                                'fun': lambda x, jac_vec: np.dot(jac_vec, x),
                                'jac': lambda x, jac_vec: jac_vec,
                                'args': [jac_vec.copy()]
                            })

        #   fix core populations
        for n in range(dim):
            if n in core_idx_list:
                jac_vec = np.zeros(dim*2)
                jac_vec[n] = 1
                constraints.append({
                    'type': 'eq',
                    'fun': lambda x, jac_vec: np.dot(jac_vec, x) - 2.0,
                    'jac': lambda x, jac_vec: jac_vec,
                    'args': [jac_vec.copy()]
                })

        
            #   total electron constraint
            jac_vec = np.zeros(dim*2)
            jac_vec[0:dim] = 1
            constraints.append({
                'type': 'eq',
                'fun': lambda x, jac_vec: np.dot(jac_vec, x) - self._numE,
                'jac': lambda x, jac_vec: jac_vec,
                'args': [jac_vec.copy()]
            })

            print(" Applying {:d} constraints".format(len(constraints)))
            return constraints

    def ESP_Min(self, x, esp_norms, esp_fit, nuclei, calc_d):
        '''
            Actual Function to be minimized
        '''
        dim = int(len(x)/2)
        coeff = np.array(x[0:dim])
        logExp = x[dim:]
        exp = np.exp(logExp)

        esp_res = calc_chelp_coeff(esp_norms, esp_fit.QM_esp_elec, self._numE, exponents=exp, coeff_deriv=calc_d, exp_deriv=calc_d, coeff=coeff)
        func_eval = esp_res['rms'] / esp_fit.sum_pot_sq
        
        #   lone pair restraints:
        for n in [idx for idx, nuc in enumerate(nuclei) if nuc == 0]:
            func_eval += 0.5*self._lone_pair_k*(coeff[n]**2)

        #   hyperbolic restraint
        for n in [idx for idx, nuc in enumerate(nuclei) if nuc == 0]:
                func_eval += 0.5*self._lone_pair_k*(coeff[n]**2)

        #   charge RESP restraints
        q = [nuclei[idx[0]] - np.sum(coeff[idx]) for idx in self._atom_map]
        sqrt_rest = np.zeros_like(q)
        for n, atom_idx in enumerate(self._atom_map):
            if nuclei[atom_idx[0]] > 1:
                sqrt_rest[n] = np.sqrt(q[n]**2 + self._resp_b**2)
                func_eval += self._resp_a*(sqrt_rest[n] - self._resp_b)

        #   derivative terms
        deriv = None
        if calc_d:
            #   derivative of ESP fitting function
            deriv = np.append(esp_res['coeff_deriv'], esp_res['exp_deriv'])/ esp_fit.sum_pot_sq

            #   derivative of lone pair restraints
            for n in [idx for idx, nuc in enumerate(nuclei) if nuc == 0]:
                deriv[n] += self._lone_pair_k*coeff[n]

            #   derivative of charge RESP restraints
            for n, atom_idx in enumerate(self._atom_map):
                if nuclei[atom_idx[0]] > 1:
                    for idx in atom_idx:
                        deriv[idx] -= q[n]*self._resp_a/sqrt_rest[n]

            deriv[dim:] *= exp

        #   updates for callback
        self._current_func_eval = func_eval
        self._current_esp_rrms = esp_res['rms'] / esp_fit.sum_pot_sq
        self._current_deriv = deriv

        return (func_eval, deriv)

    def run_fitting(self):
        density_list = self.assign_densities(self._coords_in, self._nuclei_in, self._use_lone_pairs)
        dim = len(density_list)
        centers = np.array([g.center for g in density_list])
        pol_types = np.array([g.type for g in density_list])
        all_nuclei = np.array([g.element for g in density_list])
        atom_ids = np.array([g.atom_id for g in density_list])
        center_ids = np.array([g.center_id for g in density_list])
        esp_fit = ESPFit()

        esp_fit.import_esp(self._esp_file, centers, self._coords_in, self._nuclei_in, types=pol_types)

        print("\n Performing ChElP fitting BEFORE optimization: ")
        init_chelp_fit = chelp_coeff(esp_fit.points, \
            esp_fit, centers, self._numE, types = pol_types, nuclei=all_nuclei)
        esp_norms_center = init_chelp_fit['esp_norms']

        #   re-map esp_norms to the dimensions of the number of densities
        esp_norms_density = np.empty((esp_norms_center.shape[0], dim))
        for i, g in enumerate(density_list):
            esp_norms_density[:, i] = esp_norms_center[:, g.center_id]

        #   construct mapping from atoms to gaussian number
        self._atom_map = [[] for n in range(len(self._nuclei_in))]
        for n, gauss in enumerate(density_list):
            self._atom_map[gauss.atom_id].append(n)

        #   initial guess
        init_guess_coeff = np.zeros(dim)
        init_guess_exp = np.zeros(dim)
        num_centers = {}
        for id in set(center_ids):
            num_centers[id] = len(np.where(center_ids == id)[0])

        for n, dens in enumerate(density_list):
            if dens.type[0] == 'p':
                init_guess_coeff[n] = 0.0
                init_guess_exp[n] = np.array([8.0])
            else:
                val = init_chelp_fit['s_pops'][dens.center_id]
                ng = num_centers[dens.center_id]
                init_guess_coeff[n] = val/ng
                init_guess_exp[n] = np.linspace(0.7, 2.8, ng)[np.mod(n, ng)]
        init_guess = np.ndarray.flatten(np.array([init_guess_coeff, init_guess_exp]))

        #   bounds
        bounds = []
        for n in range(len(init_guess_coeff)):
            bounds.append((0, None))
        for n in range(len(init_guess_exp)):
            bounds.append((-1.7, 9.5))

        #   constraints
        constraints = self.get_SLSQP_constrains(all_nuclei, init_guess, center_ids)

        #   required arguments for fitting function
        args = (esp_norms_density, esp_fit, all_nuclei, True)
      
        #   run minimization routine
        self.ESP_Min(init_guess, *args)
        self.ESP_Min_callback(init_guess, override=True)
        res = minimize(self.ESP_Min, init_guess, 
        args=args, 
        callback=self.ESP_Min_callback,
        method="SLSQP",
        jac=True,
        bounds=bounds,
        constraints = constraints,
        options={'disp':True, 'maxiter':100})
        self.ESP_Min_callback(res.x, override=True)


        coeff = res.x[0:dim]
        exp_list = np.exp(res.x[dim:])
        for n, c in enumerate(coeff):
            density_list[n].coeff = c
            density_list[n].exp = exp_list[n]


        print("\n ####  Final fitting results  #### \n")
        print(" ESP fitting relative RMS:     {:10.5f} %".format(sqrt(self._current_esp_rrms)*100))
        print(" Value of fitting function:    {:10.5f}".format(self._current_func_eval))
        print("\n Performing ChElP fitting AFTER optimization: ")

        final_chelp_fit = chelp_coeff(esp_fit.points, \
                esp_fit, centers, self._numE, types = pol_types, nuclei=all_nuclei, exponents=exp_list, coeff=coeff)
        
        esp_res = calc_chelp_coeff(esp_norms_density, esp_fit.QM_esp_elec, self._numE, coeff=coeff, exponents=exp_list)
        #if not args.usep:
        #    esp_points.print_esp_to_points(esp_fit.points, centers, final_chelp_fit['coeff'], 'esp_diff.pdb', exponents=exp_list, vdw_ratios=final_chelp_fit['vdw_ratios'], max_pts=15000, esp_fit=esp_fit)

        print(" {:>5s}  {:>8s}  {:>8s}  {:>8s}".format("Atom", "E-chg", "N-chg", "T-chg"))
        for idx, x in enumerate(self._atom_map):
            eChg = np.sum(coeff[x])
            nChg = self._nuclei_in[idx]
            chg = nChg - eChg
            print(" {:5d}  {:8.3f}  {:8.3f}  {:8.3f}".format(idx + 1, eChg, nChg, chg))
        print()
        s_pols_idx  = np.where(pol_types == 's')[0]
        px_pols_idx = np.where(pol_types == 'px')[0]
        py_pols_idx = np.where(pol_types == 'py')[0]
        pz_pols_idx = np.where(pol_types == 'pz')[0]

        eDip = -np.sum(coeff[s_pols_idx][:, None]*centers[s_pols_idx], axis=0)
        eDip[0] += np.sum(coeff[px_pols_idx])
        eDip[1] += np.sum(coeff[py_pols_idx])
        eDip[2] += np.sum(coeff[pz_pols_idx])
        nDip = np.sum(self._nuclei_in[:, None]*self._coords_in, axis=0)

        #eDip = -np.sum(np.array([g.center*g.coeff for g in density_list]), axis=0)
        print("Total Dipole (Dyne): ", (nDip + eDip)/0.393430307)
        print("Fitted Elec Dipole:  ", eDip/0.393430307)

        print("\n\n Parameters obtained from fitting procedure")
        print(" ----------------------------------------------------------------")
        print(" {:>3s}    Atom  Type    {:>10s}  {:>10s}".format("n", "Coeff", "Exponent"))
        print(" ----------------------------------------------------------------")
        for n, g in enumerate(density_list):
            print(" {:3d}  {:3d}   {:2s}  {:3s}  {:10.4f}  {:10.4f}".format(n, density_list[n].atom_id + 1, Elements.int2name( g.element), pol_types[n], coeff[n], exp_list[n]))
        print(" ----------------------------------------------------------------")


    def ESP_Min_callback(self, x, override=False):
        dim = int(len(x)/2)
        max_d = np.max(np.abs(self._current_deriv))
        if np.mod(self._current_step, 1) == 0 or override:
            print(" Step: {:4d}; func_eval: {:12.5f}; sqrt(RRMS): {:12.5f}; max_d: {:12.5f}; "\
                .format(self._current_step, self._current_func_eval, np.sqrt(self._current_esp_rrms), max_d))

        self._current_step += 1