import numpy as np
from simtk.openmm.app import element
import Elements
from esp_fitting import *
from scipy.optimize import minimize
from lagrangians import LMLagrangian, CostLagrangian

#cost_wt = 10**np.loadtxt('tmp')

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
    def __init__(self, coords_in, nuclei_in, esp_file, charge=0, lone_pairs=False, lone_pair_k=0.005, resp_a=0.0005, resp_b=0.10, intra_constraints=None, lone_pair_dist=0.40, fitting_method='slsqp') -> None:
        self._ANG_TO_BOHR = 1.8897259885789

        self._n_dens = 1
        self._nh_dens = 1
        self._n_atoms = 0
        self._numE = np.sum(nuclei_in) - charge
        self._n_atoms = len(nuclei_in)
        self._nuclei_in = np.copy(nuclei_in)
        self._coords_in = np.copy(coords_in)*self._ANG_TO_BOHR
        self._atom_map = []
        self._use_lone_pairs = lone_pairs
        self._lone_pair_k = lone_pair_k
        self._lone_pair_dist = lone_pair_dist*self._ANG_TO_BOHR
        self._resp_a = resp_a
        self._resp_b = resp_b
        self._intra_constraints = intra_constraints
        self._bonded_to = []
        self.fitting_method = str.lower(fitting_method)

        #   used to communicate results
        self._current_esp_rrms = None
        self._current_func_eval = None
        self._current_step = 0
        self._current_deriv = None
        self._old_x = None
        self._old_f = None
        self._current_x_diff = None
        self._best_guess = None
        self._best_func = 1E50

        #   esp data file
        self._esp_file = esp_file

        self.create_bonds()

    def create_bonds(self):
        self._bonded_to = [[] for i in range(len(self._coords_in))]
        for i, coord_i in enumerate(self._coords_in):
            for j, coord_j in enumerate(self._coords_in):
                if i == j: continue
                if self._nuclei_in[i] < 0 or self._nuclei_in[j] < 0: continue
                dist = np.linalg.norm(coord_i - coord_j)
                if 1 in [self._nuclei_in[i], self._nuclei_in[j]]:
                    dist_to_accept = 1.20*self._ANG_TO_BOHR
                else:
                    dist_to_accept = 1.60*self._ANG_TO_BOHR
                if dist <= dist_to_accept:
                    self._bonded_to[i].append(j)
    
    def create_lone_pairs(self, index, coords, nuclei, lp_dist=0.80):

        new_coords = []

        #   dermine conenctivity for atoms that will have lone pairs
        nuc_i = nuclei[index]
        i = index
        if nuc_i not in [7, 8]: return new_coords
        bonds = self._bonded_to[index]
                
        #   nitrogen bisector
        if nuc_i == 7 and len(bonds) == 2:
            r1, r2, r3 = coords[i], coords[bonds[0]], coords[bonds[1]]
            x = r1 - (r2 + r3)*0.5
            x = x/np.linalg.norm(x)
            new_coords.append(r1 + x*lp_dist)

        #   oxygen bisector
        if nuc_i == 8 and len(bonds) == 1:
            r1, r2 = coords[i], coords[bonds[0]]
            x = r2 - r1
            x = x/np.linalg.norm(x)
            heavy_2 = [x for x in self._bonded_to[bonds[0]] if self._nuclei_in[x] > 1 and x != index]
            print("HEAVY_2: ", heavy_2)
            if len(heavy_2) != 2:
                print(" WARNING: can not add lone-pair to O%d" % index)
                return new_coords
            for r3 in (coords[heavy_2[0]], coords[heavy_2[1]]):
                y = r3 - r1
                z = np.cross(x, y)
                z = z/np.linalg.norm(z)
                y = np.cross(z, x)

                p = r1 + (np.cos(110*np.pi/180)*x + np.sin(110*np.pi/180)*y)*lp_dist
                new_coords.append(p)

        return new_coords

    def assign_densities(self, coords, nuclei, add_lone_pairs=False, p_base_only=False):
        print(" Assigning Atomic Densities")
        print(" Adding lone pairs with distance {:.3f} ang.".format(self._lone_pair_dist/self._ANG_TO_BOHR))
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
                    lone_pairs = self.create_lone_pairs(i, coords, nuclei, lp_dist=self._lone_pair_dist)

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
        print(" {:>6s}  {:>6s}  {:>6s}  {:>3s}  {:>4s}  {:>10}  {:>10}  {:>10}".format('N', 'AtomID', 'CentID', 'Elm', 'Type', 'X (Ang)', 'Y (Ang)', 'Z (Ang)'))
        print(' -----------------------------------------------------------------------')
        #print(len(density_list))
        #print("")
        for n, g in enumerate(density_list):
            #g.center /= self._ANG_TO_BOHR
            #print(Elements.int2name(g.element), g.center[0], g.center[1], g.center[2])
            print(" {:6d}  {:6d}  {:6d}  {:3s}  {:4s}  {:10.4f}  {:10.4f}  {:10.4f}".format(n+1, g.atom_id, g.center_id, Elements.int2name(g.element), g.type, g.center[0]/self._ANG_TO_BOHR, g.center[1]/self._ANG_TO_BOHR, g.center[2]/self._ANG_TO_BOHR))
        print(' -----------------------------------------------------------------------')
        return density_list

    def get_SLSQP_constrains(self, densities, guess):
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
        print(" -----------------------------------------------------------------")
        print("              Adding cosntraints to the system")
        print(" -----------------------------------------------------------------")


        constraints = []
        dim = int(len(guess)/2)
        guess_exp = guess[dim:]
        all_nuclei = [g.element for g in densities]
        atom_ids = np.array([g.atom_id for g in densities])
        center_ids = np.array([g.center_id for g in densities])
        
        #   extract core and hydrogen types
        core_idx_list = []
        hydrogen_list = []
        #for id in set(center_ids):
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

        #   add exponent constraints for each unique type
        unique_sigs = set(signatures)
        for uniq in unique_sigs:
            if uniq[1] in ['core']:
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
                            constraints.append(self.get_linear_constraint(jac_vec, 0.0))

        #   fix core populations
        for n in range(dim):
            if n in core_idx_list:
                jac_vec = np.zeros(dim*2)
                jac_vec[n] = 1
                constraints.append(self.get_linear_constraint(jac_vec, 2.0))

        #   total electron constraint
        print(" Applying 1 total charge constraint to the entire system")
        jac_vec = np.zeros(dim*2)
        jac_vec[0:dim] = 1
        constraints.append(self.get_linear_constraint(jac_vec, self._numE))

        #   lone pair constraints
        for i, sites in enumerate(self._atom_map):
            lp_idx = [n for n in sites if all_nuclei[n] == 0]
            if len(lp_idx) > 1:
                n_lp_constr = 0
                for n in lp_idx[1:]:
                    n_lp_constr += 1
                    for shift in []:
                        jac_vec = np.zeros(dim*2)
                        jac_vec[lp_idx[0] + shift] = 1
                        jac_vec[n + shift] = -1
                        constraints.append(self.get_linear_constraint(jac_vec, 0.0))
                print(" Applying %d lone-pair population constraints with host site O%d" % (n_lp_constr, (i+1)))
                print(" Applying %d lone-pair exponent constraints with host site O%d" % (n_lp_constr, (i+1)))

        #   equivilant hydrogens from methyl and amine groups
        for i, bonds in enumerate(self._bonded_to):
            elm_i = self._nuclei_in[i]
            hydro_atom_idx = [x for x in bonds if Elements.int2name(self._nuclei_in[x]) == 'H']
            n_hydro = len(hydro_atom_idx)
            if (elm_i == 6 and n_hydro == 3) or (elm_i == 7 and n_hydro == 2):
                for n in range(n_hydro):
                    #   first hydrogen is the reference atom
                    if n == 0: continue
                    jac_vec = np.zeros(dim*2)
                    #   density population constraints
                    for site_idx in self._atom_map[hydro_atom_idx[0]]:
                        jac_vec[site_idx] = -1
                    for site_idx in self._atom_map[hydro_atom_idx[n]]:
                        jac_vec[site_idx] = 1
                    constraints.append(self.get_linear_constraint(jac_vec, 0.0))
                    #   density exponents constraints
                    for site_idx in self._atom_map[hydro_atom_idx[0]]:
                        jac_vec[dim + site_idx] = -1
                    for site_idx in self._atom_map[hydro_atom_idx[n]]:
                        jac_vec[dim + site_idx] = 1
                    constraints.append(self.get_linear_constraint(jac_vec, 0.0))

                if elm_i == 6:
                    atom_type, group_name = "C", "methyl"
                elif elm_i == 7:
                    atom_type, group_name = "N", "amine"
                print(" Identified %s group:" % group_name)
                print(" \t{:s}{:<2d} with site index {:2d}".format(atom_type, i + 1, self._atom_map[i][0]))
                for n in hydro_atom_idx:
                    print(" \tH{:<2d} with site index {:2d}".format(n + 1, self._atom_map[n][0]))
                print(" Applying %d intra-molecular charge constraints" % (n_hydro - 1))
                print(" Applying %d intra-molecular exponent constraints\n" % (n_hydro - 1))



        #   intramoleuclar constraints
        print(" Applying %d user defined intra-molecular charge constraints" % len(self._intra_constraints))
        for constr in self._intra_constraints:
            total_nuclei = 0.0
            jac_vec = np.zeros(dim*2)
            #   convert from atom constraints to total density on each site
            for constr_idx, value in constr.constr_dict.items():
                center_idx = [n for n in range(len(atom_ids)) if atom_ids[n] == constr_idx]
                for idx in center_idx:
                    jac_vec[idx] = value
                total_nuclei += all_nuclei[center_idx[0]]
            constraints.append(self.get_linear_constraint(jac_vec, (total_nuclei - constr.charge)))
        
        
        print(" \n\n Total number of cosntraints: {:d}".format(len(constraints)))
        print(" -----------------------------------------------------------------\n")
        return constraints

    def get_linear_constraint(self, jac_vec, constr_val):
        return {
                    'type': 'eq',
                    'fun': lambda x, jac_vec: np.dot(jac_vec, x) - constr_val,
                    'jac': lambda x, jac_vec: jac_vec,
                    'args': [jac_vec.copy()]
                }

    def ESP_Min(self, x, esp_norms, esp_fit, nuclei, calc_d):
        '''
            Actual Function to be minimized
        '''
        dim = int(len(x)/2)
        coeff = np.array(x[0:dim])
        #logExp = x[dim:]
        #exp = np.exp(logExp)
        exp = x[dim:]

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

        #   hessian matrix cost function
        # corr_wt = 10**np.loadtxt('tmp')*0
        # G_inv_coeff = np.linalg.inv(esp_res['G_mat']) @ coeff
        # cost = 0.5 * coeff @ G_inv_coeff
        # #print("COST: ", cost)
        # func_eval += corr_wt * cost 

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

            #   hessian matrix cost function
            # deriv[:dim] += G_inv_coeff * corr_wt
            deriv *= 100

        #   scale function for better search performance with scipy
        func_eval *= 100
        

        #   updates for callback
        self._current_func_eval = func_eval
        self._current_esp_rrms = esp_res['rms'] / esp_fit.sum_pot_sq
        self._current_deriv = deriv

        if calc_d:
            #print("IN ESP: ", deriv)
            return (func_eval, deriv)
        else:
            return func_eval

    def _test_numerical_derivative(self, x0, args):
        eps = 1E-5
        fun, deriv = self.ESP_Min(x0, *args)
        for n in range(len(x0)):
            xp = np.copy(x0)
            xp[n] += eps
            xm = np.copy(x0)
            xm[n] -= eps

            fp = self.ESP_Min(xp, *args[:-1], calc_d=False)
            fm = self.ESP_Min(xm, *args[:-1], calc_d=False)

            num_deriv = (fp - fm)/(2*eps)
            num_deriv = (fp + fm - 2*fun)/(eps*eps)
            print(" {:10.5f}  {:10.5f}  {:10.5f}".format(x0[n], deriv[n], num_deriv))
        exit()

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
            esp_fit, centers, self._numE, all_nuclei, types = pol_types)
        esp_norms_center = init_chelp_fit['esp_norms']
        #exit()

        #   re-map esp_norms to the dimensions of the number of densities
        esp_norms_density = np.empty((esp_norms_center.shape[0], dim))
        for i, g in enumerate(density_list):
            esp_norms_density[:, i] = esp_norms_center[:, g.center_id]

        #   construct mapping from atoms to density number
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
                init_guess_exp[n] = Elements.getExponentByAtomicNumber(dens.element)
        init_guess = np.ndarray.flatten(np.array([init_guess_coeff, init_guess_exp]))

        #   bounds
        bounds = []
        for n in range(len(init_guess_coeff)):
            bounds.append((0, None))
        for n in range(len(init_guess_exp)):
            bounds.append((1.50, 15.0))

        #   constraints
        constraints = self.get_SLSQP_constrains(density_list, init_guess)

        #   required arguments for fitting function
        args = (esp_norms_density, esp_fit, all_nuclei, True)
        

        #   run minimization routine
        self.ESP_Min(init_guess, *args)
        if self.fitting_method == 'slsqp':
            print("\n Starting SciPy SLSQP minimization routine ")
            res = minimize(
                self.ESP_Min, 
                init_guess, 
                method="SLSQP",
                options={'disp':False, 'maxiter':300, 'ftol': 1e-10},
                bounds=bounds,
                constraints = constraints,
                args=args, 
                callback=self.ESP_Min_callback,
                jac=True,
                )
        elif self.fitting_method == 'bfgs':
            print("\n Starting SciPy L-BFGS-B minimization routine ")
            lagrangian = CostLagrangian(self.ESP_Min, init_guess, args, constraints, 10)
            res = minimize(
                lagrangian.min_func,
                lagrangian.get_guess(),
                method="l-bfgs-b",
                options={'disp':True, 'maxiter':500, 'gtol': 1E-2, 'ftol': 1E-14},
                args=args, 
                callback=self.ESP_Min_callback,
                jac=True,
                )
        else: raise ValueError(" Fitting function must be 'SLSQP' or 'BFGS'")

        #   overwrite with best guess found through the entire searhc process
        res.x = self._best_guess[0:2*dim]  
        #self._test_numerical_derivative(res.x, args)

        print(" ---------------------------------------------------------------------")
        print("        #### " + res.message + " #### ")
        print(" Using best guess found:")
        self.ESP_Min_callback(res.x, override=True)
        print(" ---------------------------------------------------------------------")

        coeff = res.x[0:dim]
        exp_list = np.exp(res.x[dim:])
        exp_list = res.x[dim:]
        for n, c in enumerate(coeff):
            density_list[n].coeff = c
            density_list[n].exp = exp_list[n]

        print(" ESP fitting relative RMS:     {:10.5f} %".format(sqrt(self._current_esp_rrms)*100))
        print(" Value of fitting function:    {:10.5f}".format(self._current_func_eval))
        print("\n Performing ChElP fitting AFTER optimization: ")

        final_chelp_fit = chelp_coeff(esp_fit.points, \
                esp_fit, centers, self._numE, all_nuclei, types=pol_types, exponents=exp_list, coeff=coeff, print_results=False)
        
        calc_chelp_coeff(esp_norms_density, esp_fit.QM_esp_elec, self._numE, coeff=coeff, exponents=exp_list)
        #if not args.usep:
        #    esp_points.print_esp_to_points(esp_fit.points, centers, final_chelp_fit['coeff'], 'esp_diff.pdb', exponents=exp_list, vdw_ratios=final_chelp_fit['vdw_ratios'], max_pts=15000, esp_fit=esp_fit)

        s_pols_idx  = np.where(pol_types == 's')[0]
        eDip = -np.sum(coeff[s_pols_idx][:, None]*centers[s_pols_idx], axis=0)
        if len(np.where(pol_types == 'px')[0]) != 0:
            px_pols_idx = np.where(pol_types == 'px')[0]
            py_pols_idx = np.where(pol_types == 'py')[0]
            pz_pols_idx = np.where(pol_types == 'pz')[0]
            eDip[0] += np.sum(coeff[px_pols_idx])
            eDip[1] += np.sum(coeff[py_pols_idx])
            eDip[2] += np.sum(coeff[pz_pols_idx])
        nDip = np.sum(self._nuclei_in[:, None]*self._coords_in, axis=0)

        print("Total Dipole (Dyne): ", (nDip + eDip)/0.393430307)
        print("Fitted Elec Dipole:  ", eDip/0.393430307)

        print("\n\n         Parameters obtained from fitting procedure")
        print(" --------------------------------------------------------------------")
        print("   n    Atom  Type    Atom-Charge  Site-Charge  Elec-Coeff  Exponent")
        print(" ")
        print(" --------------------------------------------------------------------")
        for n, g in enumerate(density_list):
            elm = Elements.int2name( g.element)
            if elm != 'Lp' and pol_types[n] == 's':
                atom_coeff = [coeff[i] for i in range(len(coeff)) if atom_ids[i] == atom_ids[n]]
                atm_chg = all_nuclei[n] - np.sum(atom_coeff)
                atm_chg_str = '{:10.4f}'.format(atm_chg)
            else:
                atm_chg_str = '{:10s}'.format('')
            print(" {:3d}  {:3d}   {:2s}  {:3s}  {:10s}  {:10.4f}  {:10.4f}  {:10.4f}".format(n + 1, density_list[n].atom_id + 1, elm, pol_types[n], atm_chg_str, all_nuclei[n] - coeff[n], coeff[n], exp_list[n]))
        print(" --------------------------------------------------------------------")


    def ESP_Min_callback(self, x, override=False):
        if self._current_step == 0:
            print(" ---------------------------------------------------------------------")
            print("   Cycle   Func-Eval      ESP-RRMS     Coeff-Diff    Exp-Diff")
            print(" ---------------------------------------------------------------------")
            self._old_x = np.zeros_like(x)
            self._old_f = 0.0

        if np.all(x == self._old_x):
            diff_x = np.copy(self._current_x_diff)
            diff_f = self._current_func_eval
        else:
            diff_x = np.abs(x - self._old_x)
            diff_f = self._current_func_eval - self._old_f

        if np.mod(self._current_step, 5) == 0 or override:
            dim = int(len(x)/2)
            max_diff_coeff = np.max(diff_x[:dim])
            max_diff_exp = np.max(diff_x[dim:2*dim])
            if not override:
                print(" {:5d}  {:14.10f}  {:10.5f}  {:12.2e} {:12.2e}"\
                    .format(self._current_step, self._current_func_eval, np.sqrt(self._current_esp_rrms), max_diff_coeff, max_diff_exp))
            else:
                print(" {:5d}  {:14.10f}  {:10.5f}"\
                    .format(self._current_step, self._current_func_eval, np.sqrt(self._current_esp_rrms)))



        self._current_x_diff = diff_x
        self._current_step += 1
        self._old_x = np.copy(x)
        self._old_f = self._current_func_eval

        #   keep track of best guess
        if self._current_func_eval < self._best_func and self._current_step > 5:
            self._best_func = self._current_func_eval
            self._best_guess = np.copy(x)


    def ESP_Min_callback_old(self, x, override=False):
        dim = int(len(x)/2)
        max_d = np.max(np.abs(self._current_deriv))
        if np.mod(self._current_step, 5) == 0 or override:
            print(" Step: {:4d}; func_eval: {:12.5f}; sqrt(RRMS): {:12.5f}; max_d: {:12.5f}; "\
                .format(self._current_step, self._current_func_eval, np.sqrt(self._current_esp_rrms), max_d))

        self._current_step += 1