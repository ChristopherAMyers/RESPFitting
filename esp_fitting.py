import numpy as np
from numpy.lib.arraysetops import unique
from numpy.linalg import norm
from math import ceil, pi, floor, sin, cos, sqrt
import sys
from os import devnull
from numpy.linalg.linalg import eigvals

from simtk.openmm.app import element

import Elements

class ESPFit:
    def __init__(self):
        self.G_mat = np.array(None)
        self.e_vec = np.array(None)
        self.sum_pot_sq = 1.0
        self.sum_pot_sq_elec = 1.0
        self.points = np.array(None)
        self.QM_esp = np.array(None)
        self.QM_esp_elec = np.array(None)
        self.nuclei = np.array(None)

    def import_esp(self, fileLoc, coords, nuc_coords, nuclei, types = []):
        print(" Importing ESP data")
        data = np.loadtxt(fileLoc, skiprows=4)
        self.points = data[:, 0:3]
        self.QM_esp = data[:, 3]
        print(" Done Importing ESP Data\n Forming G-Matrix elements")
        self._form_components(coords, nuc_coords, nuclei, types)
        print(" Done")

    def _form_components(self, coords, nuc_coords, nuclei, types):
        diffs = (self.points[:, None] - coords)
        norms = norm(diffs, axis=-1)

        inverse_norms = 1 / norms
        for n, pole in enumerate(types):
            if pole == 'px':
                inverse_norms[:, n] = (self.points[:, 0] - coords[n, 0])/norms[:, n]**3
            elif pole == 'py':
                inverse_norms[:, n] = (self.points[:, 1] - coords[n, 1])/norms[:, n]**3
            elif pole == 'pz':
                inverse_norms[:, n] = (self.points[:, 2] - coords[n, 2])/norms[:, n]**3

        self.G_mat = inverse_norms.T @ inverse_norms
        self.sum_pot_sq = np.sum(self.QM_esp**2)

        diffs_nuc = (self.points[:, None] - nuc_coords)
        norms_nuc = norm(diffs_nuc, axis=-1)
        nuc_pot = np.sum(nuclei / norms_nuc, axis=-1)
        self.QM_esp_elec = self.QM_esp - nuc_pot
        self.sum_pot_sq_elec = np.sum(self.QM_esp_elec**2)
        self.e_vec = inverse_norms.T @ self.QM_esp_elec


def _format_83(f):
    """Format a single float into a string of width 8, with ideally 3 decimal
    places of precision. If the number is a little too large, we can
    gracefully degrade the precision by lopping off some of the decimal
    places. If it's much too large, we throw a ValueError
        Note: Taken form simtk.openmm.app.pdbfile
    """
    if -999.999 < f < 9999.999:
        return '%8.3f' % f
    if -9999999 < f < 99999999:
        return ('%8.3f' % f)[:8]
    raise ValueError('coordinate "%s" could not be represented '
                     'in a width-8 field' % f)

def print_esp_to_points(points, coords, coeff, file_out, exponents=None, vdw_ratios=None, max_pts=None, esp_fit=None):
    ''' Print ESP evaluated at xyz points to PDB file
    '''

    diffs = points[:, None] - coords
    norms = norm(diffs, axis=-1)
    unit_potential = 1.0/norms

    #   use slater densities if exponents are provided, else use point charges
    if exponents is not None:
        alpha = np.array(exponents)
        a_R = alpha*norms
        exp_ar = np.exp(-a_R)
        unit_potential -= exp_ar*(unit_potential + alpha*0.5)
    
    esp = -np.sum(np.array(coeff) * unit_potential, axis=1)
    
    #   re-normalize so the numbers are between -99.99 and 99.99
    if esp_fit is not None:
        #   subtract out QM esp if provided and use percentage
        
        esp = 1*(esp - esp_fit.QM_esp_elec)/np.abs(esp_fit.QM_esp)
        esp = 1*(esp - esp_fit.QM_esp_elec)
        esp[esp >= 100 ] = 99.99
        esp[esp < -100 ] = -99.99
    else:
        esp = 99.99*esp/np.max(np.abs(esp))
    
    if vdw_ratios is not None:
        vdw_ratios = np.array(vdw_ratios, dtype=float)
        uniq_ratio = sorted(set(np.round(vdw_ratios, 4)))
        if len(uniq_ratio) <= 20:
            res_IDs = ['{:3d}'.format(np.argmin(np.abs(uniq_ratio - ratio))) for ratio in vdw_ratios]
            res_names = ['R{:<2d}'.format(np.argmin(np.abs(uniq_ratio - ratio))) for ratio in vdw_ratios]
    else:
        res_IDs = ['1']*len(points)
        res_names = ['R1']*len(points)

    #   limit the number of points to write to file
    if max_pts is not None:
        n_pts = len(points)
        stride = ceil(n_pts / max_pts)
        points = points[::stride]
        res_IDs = res_IDs[::stride]
        res_names = res_names[::stride]
        esp = esp[::stride]

    with open(file_out, 'w') as file:
        
        for n, point in enumerate(points):
            pt_ang = point/1.88973
            resId = res_IDs[n]
            resName = res_names[n]
            chainName = 'A'
            atomName = 'He'
            symbol = 'He'
            recordName = 'ATOM  '
            resIC = " "
            line = "%s%5d %-4s %3s %s%4s%1s   %s%s%s  1.00%6.2f          %2s  \n" % (
                        recordName, n%100000, atomName, resName, chainName, resId, resIC, _format_83(pt_ang[0]),
                        _format_83(pt_ang[1]), _format_83(pt_ang[2]), esp[n], symbol)
            file.write(line)


def calc_chelp_coeff(norms, QM_esp_elec, n_elec, coeff=None, types=None, exponents=None, coeff_deriv=False, exp_deriv=False):
    unit_potential = 1.0/norms

    if exponents is not None:
        alpha = np.array(exponents)
        a_R = alpha*norms
        exp_ar = np.exp(-a_R)
        unit_potential -= exp_ar*(unit_potential + alpha*0.5)

    G_mat = unit_potential.T @ unit_potential
    e_vec = unit_potential.T @ QM_esp_elec

    #   don't solve for coeff if already provided
    if coeff is None:
        Ginv = np.linalg.inv(G_mat)
        Ginv_e = Ginv @ e_vec
        if types is not None:
            c = np.zeros(len(e_vec)) + (types == 's')
        else:
            c = np.ones(len(e_vec))
        sum_G = c @ Ginv @ c

        lamb = (n_elec + c @ Ginv_e)/sum_G
        coeff = -Ginv_e + lamb * Ginv @ c

    res = {'rms': None, 'coeff_deriv': None, 'exp_deriv': None}
    res['G_mat'] = G_mat
    res['rms'] = coeff @ G_mat @ coeff + 2*e_vec @ coeff + np.sum(QM_esp_elec**2)
    res['coeff'] = coeff

    #   derivative of RMS w.r.t. electron coefficients
    if coeff_deriv:
        res['coeff_deriv'] = 2*G_mat @ coeff + 2*e_vec
        
    #   derivative of RMS w.r.t. electron density exponents
    if exp_deriv:
        if exponents is not None:
            #d_pot_2 = 0.5 - 0.5*norms*unit_potential + 0.25*a_R*exp_ar
            d_pot = 1.0 - norms*unit_potential - 0.5*exp_ar
            deriv = np.zeros(len(exponents))
            for a, pop in enumerate(coeff):
                G_prime = unit_potential.T @ d_pot[:, a]
                deriv[a] = 2*pop*(coeff @ G_prime + QM_esp_elec.T @ d_pot[:, a])
        res['exp_deriv'] = deriv

    #print("COEFF: ", coeff)
    #print("EXP: ", exponents)
    return res

def chelp_coeff(points, esp_fit, coords, n_elec, nuclei, types = [], exponents=None, coeff=None, print_results=True):
    '''
    Calculate ChElPG equvilant for electronic coefficients
    Note: Assumes that QM_esp_elec is the ELECTRON potential
    '''

    #   if no types are given, assume that they are all s-type
    if len(types) == 0:
        types = ['s']*len(coords)

    #   loop through centers and eliminate duplicates
    uniq_coords = []
    uniq_types = []
    uniq_coord_idx = []
    uniq_to_all = []
    for i, coord1 in enumerate(coords):
        add_to_list = True
        for j, coord2 in enumerate(uniq_coords):

            if np.array_equal(coord1, coord2) and (types[i] == uniq_types[j]):
                add_to_list = False
                uniq_to_all[j].append(i)
                break

        if add_to_list:
            #   check again so that we keep a list of unique atom centers
            if types[i] == 's':
                uniq_coord_idx.append(i)
            uniq_coords.append(coord1)
            uniq_types.append(types[i])
            uniq_to_all.append([i])

    uniq_coords = np.array(uniq_coords)
    uniq_types = np.array(uniq_types)
    uniq_nuclei = nuclei[uniq_coord_idx]
    print(" Using {:d} ESP points".format(len(points)))
    print(" Provided {:d} centers".format(len(coords)))
    print(" Found {:d} unique atomic center types".format(len(uniq_coords)))

    #   If exponents are not provided (that is, if the potential is from point charges),
    #   then the CHELP fitting will result in a singular G matrix. Therefore, we
    #   use only the unique centers for the fitting if no exponents are provided.
    if exponents is None:
        print(" Consolidating similar nuclei into unique centers")
        method_coords = uniq_coords
        method_types = uniq_types
        exponents = np.array([Elements.getExponentByAtomicNumber(x) for x in uniq_nuclei])
    else:
        method_coords = coords
        method_types = types

    diffs = (points[:, None] - method_coords)
    norms = norm(diffs, axis=-1)

    #   actual CHELP calculation
    all_res = calc_chelp_coeff(norms, esp_fit.QM_esp_elec, n_elec=n_elec, types=method_types, exponents=exponents, coeff=coeff)
    coeff = all_res['coeff']

    #   calcualte vdW ratios possible used
    # uniq_norms = norm(points[:, None] - uniq_coords, axis=-1)
    # ratios = uniq_norms/1.889725989
    # ratios = np.round(ratios, 3)
    # for n, nuc in enumerate(uniq_nuclei):
    #     if nuc == 0: continue
    #     ratios[:, n] /= Elements.getRadiiByAtomicNumber(nuc)
    # min_ratios = np.min(ratios, axis=1)


    #   calcualte vdW ratios possible used
    ratio_nuclei = [nuc for n, nuc in enumerate(uniq_nuclei) if nuc > 0]
    ratio_coords = [uniq_coords[n] for n, nuc in enumerate(uniq_nuclei) if nuc > 0]
    ratio_norms = norm(points[:, None] - ratio_coords, axis=-1)
    ratios = ratio_norms/1.889725989
    ratios = np.round(ratios, 3)
    for n, nuc in enumerate(ratio_nuclei):
        ratios[:, n] /= Elements.getRadiiByAtomicNumber(nuc)
    min_ratios = np.min(ratios, axis=1)


    #   find unique ratios and remove trailing digits from rounding
    min_ratios = np.array(['{:.4f}'.format(x) for x in min_ratios])
    unique_ratios_str = sorted(set(min_ratios))
    print( " Identified {:d} unique vdW radii ratios".format(len(unique_ratios_str)))

    if len(unique_ratios_str) <= 20:
        for ratio in unique_ratios_str:
            sub_where = (min_ratios == ratio) 
            sub_norms = norms[sub_where]

            #   norms will be inverted in calc_chelp_coeff, so we pre multiple by correct first for p-densities
            for n, pole in enumerate(method_types):
                if pole == 'px':
                    sub_norms[:, n] = 1.0/( (points[sub_where, 0] - method_coords[n, 0])/sub_norms[:, n]**3 )
                elif pole == 'py':
                    sub_norms[:, n] = 1.0/( (points[sub_where, 1] - method_coords[n, 1])/sub_norms[:, n]**3 )
                elif pole == 'pz':
                    sub_norms[:, n] = 1.0/( (points[sub_where, 2] - method_coords[n, 2])/sub_norms[:, n]**3 )

            res = calc_chelp_coeff(sub_norms, esp_fit.QM_esp_elec[sub_where], coeff=coeff, n_elec=n_elec, types=method_types, exponents=exponents)
            rrms = sqrt(res['rms']/np.sum(esp_fit.QM_esp[sub_where]**2))
            #print("     vdW ratio {:>8s} ({:7d} pts): sqrt(rms) = {:10.6f}".format(ratio, len(sub_norms), sqrt(res['rms'])))
            print("     vdW ratio {:>8s} ({:7d} pts): sqrt(rrms) = {:10.6f}".format(ratio, len(sub_norms), rrms))
        
        print("     Overall sqrt(rrms) = {:10.6f}".format(sqrt(all_res['rms']/esp_fit.sum_pot_sq)))

    rmsd = None
    if esp_fit.sum_pot_sq is not None and esp_fit.sum_pot_sq_elec is not None:
        rmsd = sqrt(all_res['rms']) / sqrt(esp_fit.sum_pot_sq)
        rmsd = np.round(rmsd, 5)

    #   if nuclei are provided, then print results
    if print_results: out = sys.stdout
    else: out = open(devnull, 'w')
    print(" Populations:", file=out)
    uniq_nuclei = nuclei[uniq_coord_idx]

    s_types_idx = np.where(uniq_types == 's')[0]
    s_pops = []
    p_pops = []
    s_dipole = np.array([0., 0., 0.])
    for i, s_idx in enumerate(s_types_idx):
        center = uniq_coords[s_idx]
        nuc = uniq_nuclei[i]
        elm = Elements.int2name(nuc)
        if len(coeff) == len(uniq_coords):
            center_pop = coeff[i]
        else:
            center_pop = np.sum(coeff[x] for x in uniq_to_all[i])
        if uniq_types[i][0] == 'p':
            p_pops.append(center_pop)
            print(" {:3d}  {:3s}  {:3s}              {:12.6f} ".format(i + 1, elm, uniq_types[i], center_pop), file=out)
        else:
            s_pops.append(center_pop)
            s_dipole += center * (nuc - center_pop)/0.3934303    # units are in Debye
            print(" {:3d}  {:3s}  {:3s}  {:12.6f} ".format(i + 1, elm, uniq_types[i], nuc - center_pop), file=out)
    
    print(" Sum of atomic charges = {:12.6f}".format(np.sum(uniq_nuclei - s_pops)), file=out)
    print(" s-dipole moment: {:8.3f} (X {:8.3f} Y {:8.3f} Z{:8.3f}) Debye\n".format(norm(s_dipole), *tuple(s_dipole)), file=out)


    res = {'s_pops': np.array(s_pops), 'p_pops': np.array(p_pops), 'esp_norms': norms, 'vdw_ratios': min_ratios}
    res.update(all_res)
    return res
    #return (s_pops, p_pops, norms, ratios)


def rel_rms(esp_fit, coeff):
    chg = coeff
    rms = chg @ esp_fit.G_mat @ chg + 2 * esp_fit.e_vec @ chg + esp_fit.sum_pot_sq_elec
    return rms / esp_fit.sum_pot_sq

def deriv_rel_rms(esp_fit, coeff):
    chg = coeff
    term1 = esp_fit.G_mat @ chg
    deriv = 2 * term1 + 2*esp_fit.e_vec
    rms = chg @ term1 + 2 * esp_fit.e_vec @ chg + esp_fit.sum_pot_sq_elec
    
    rms /= esp_fit.sum_pot_sq
    deriv /= esp_fit.sum_pot_sq
    return rms, deriv    

#   in angstroms
Bondi_radii = {
        1: 1.06,
        6: 1.53,
        7: 1.46,
        8: 1.42,
            #the following are our own radii
            #these are not in the paper.
        15: 1.86,
        16: 1.80,
}

def get_Bondi_radii(nuclei):
    return np.array([Bondi_radii.get(nuc, 2.0) for nuc in nuclei])

def neighbor_list(atoms, atom_radii, scale):
    distMat =  norm(atoms[:, None] - atoms, axis=-1)
    rrMat = (atom_radii[:, None] + atom_radii)*scale
    keepMat = distMat <= rrMat


    nbList = []
    for m in range(len(atoms)):
        tmp = []
        for n in range(len(atoms)):
            if keepMat[m, n] and n != m:
                tmp.append(n)
        nbList.append(tmp)
    return nbList 


