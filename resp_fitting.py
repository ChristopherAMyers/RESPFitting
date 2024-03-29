#!/usr/bin/env python3
import numpy as np
from openmm.app import PDBFile
from openmm.unit import *
import os
from os import chdir, makedirs, environ
import os.path as path
from distutils.util import strtobool
from shutil import copyfile, move, which, rmtree
import argparse
import warnings
import subprocess
from multiprocessing import cpu_count
import glob
import sys

from options import Options
from esp_points import ESPPointGenerator
from density_fitting import DensityFitting
from molecule import Molecule

#   set available CPU resources
if 'SLURM_NTASKS' in environ.keys():
    GLOB_n_procs = int(environ['SLURM_NTASKS'])
else:
    #   if not running a slurm job, use number of cores
    GLOB_n_procs = cpu_count()


def parse_options(opt_file):
    opts = Options()
    rem_opts = {}
    if not opt_file:
        #   no file provided means use the default options
        return opts, rem_opts
    rem_lines_in = []

    with open(opt_file, 'r') as file:
        for line in file.readlines():
            if "$" not in line:
                rem_lines_in.append(line.replace('=', ' '))

    for line in rem_lines_in:
        line = line.lower()
        sp_comment = line.split('!')
        if sp_comment[0] != '!':
            sp = line.split()
            if len(sp) == 0:
                continue
            option = sp[0].lower()
            if option == 'name':
                opts.time_step = sp[1]
            elif option == 'basis':
                opts.basis = sp[1]
            elif option == 'method':
                opts.method = sp[1]
            elif option == 'charge':
                opts.charge = int(sp[1])
            elif option == 'optimize':
                opts.optimize_first = strtobool(sp[1])
            elif option == 'density_fitting':
                opts.density_fitting = strtobool(sp[1])
            elif option == 'lone_pairs':
                opts.lone_pairs = strtobool(sp[1])
            else:
                #   all other lines must be Q-Chem rem options
                rem_opts[sp[0]] = str(sp[1])

    return opts, rem_opts

def create_qchem_job_str(options, coords, elements, total_chg = 0, spin_mult=1):

    file_str = ""
    # file_str = "$rem\n"
    # for key in rem_opts:
    #     file_str += "  {:20s}  {:20s} \n".format(key, str(rem_opts[key]))
    # file_str += "$end\n"
    for section, qc_options in options.items():
        if section in ['resp', 'intra_constraints']: continue
        file_str += '$%s\n' % section
        if section == 'rem':
            for key, value in qc_options.items():
                file_str += '   {:25s}  {:25s}\n'.format(key, value)
        else:
            for line in qc_options:
                file_str += '   {:s}\n'.format(line)
        file_str += '$end\n\n'

    file_str += "$molecule\n"
    if isinstance(coords, str):
        file_str += coords + "\n"
    else:
        file_str += "{:d} {:d} \n".format(total_chg, spin_mult)
        for n in range(len(elements)):
            file_str += '  {:2s}  {:15.8f}  {:15.8f}  {:15.8f}\n'.format(elements[n], *coords[n])
    file_str += "$end\n"
    return file_str

def write_esp_point_file(out_file_loc, esp_points):
     #   write ESP points to file file
    with open(out_file_loc, 'w') as file:
        for point in esp_points:
            file.write('{:15.8f}  {:15.8f}  {:15.8f} \n'.format(*point))

def create_qchem_job(opts, scratch, name, esp_points, coords, elements, outfile=sys.stdout):

    qc_envirn = environ.get('QC', None)
    if not qc_envirn:
        raise ValueError("QC environment varible not defined")

    #   create job directories and file names
    job_dir = path.join(scratch, name)
    makedirs(job_dir, exist_ok=True)
    input_file = path.join(job_dir, 'qchem.in')
    output_file = path.join(job_dir, 'qchem.out')
    save_flag = ""

    if opts.optimize_first:
        pass
    else:
        rem_opts = opts.input_sections.copy()
        rem_opts.pop('resp', None)
        # rem_opts = other_rem_opts.copy()
        # rem_opts['method'] = opts.method
        # rem_opts['basis'] = opts.basis
        # rem_opts['jobtype'] = 'sp'
        rem_opts['rem']['esp_grid'] = str(len(esp_points))
        rem_opts['rem']['sym_ignore'] = 'true'
        ipt_lines = create_qchem_job_str(rem_opts, coords, elements, total_chg=opts.charge)
        with open(input_file, 'w') as file:
            file.write(ipt_lines)

    #   write ESP points to file file
    write_esp_point_file(path.join(job_dir, 'ESPGrid'), esp_points)
    
    #   Change to directory containing input/output file.
    #   This is because Q-Chem looks for ESPGrid file there
    this_dir = path.abspath(path.curdir)
    chdir(path.dirname(input_file))

    #   submit Q-Chem job
    qchem_bin = path.join(qc_envirn, 'bin/qchem')
    args = [qchem_bin]
    if save_flag != '':
        args += [save_flag]
    args += ['-nt', str(GLOB_n_procs), input_file, output_file]

    print(" Starting Q-Chem job. Output will be printed to:", file=outfile)
    print("\t " + path.abspath(output_file), file=outfile)
    outfile.flush()
    try:
        pass
        output = subprocess.check_output(args)
        print(' Q-Chem Console Output:', file=outfile)
        print_subprocess_output(output, outfile)                                                   
    except subprocess.CalledProcessError as error:
        print(' Error code:', error.returncode, '. Output:', error.output.decode("utf-8"), file=outfile)

    #   check Q-Chem output
    success = False
    charges = None
    energy = None
    read_charges_line = 1E50
    with open(output_file, encoding="utf-8") as file:
        for n, line in enumerate(file.readlines()):
            if "Thank you" in line:
                success = True
            if " ChElPG" in line:
                read_charges_line = n + 4
                charges = []
            if n >= read_charges_line:
                if len(charges) == len(coords):
                    read_charges_line == 1E50
                else:
                    charges.append(float(line.split()[2]))
            if "Total energy" in line:
                energy = float(line.split()[8])

    if success:
        print(" Q-Chem job has finished successfully", file=outfile)
    else:
        print(" ERROR: Q-Chem job did finish successfully!", file=outfile)
        print("        Check Q-Chem output file at", file=outfile)
        print("\t " + path.abspath(output_file), file=outfile)
        exit()

    #   change back to correct directory
    chdir(this_dir)

    esp_file_path = path.join(job_dir, 'plot.esp')
    if path.isfile(esp_file_path):
        return (esp_file_path, energy, charges)
    else:
        raise FileNotFoundError('plot.esp does not exist.\n \
            Please check Q-Chem output for possible errors')

def crop_respin(file_loc):
    keep_lines = False
    lines = []
    with open(file_loc) as file:
        for line in file.readlines():
            if keep_lines:
                lines.append(line)
                if len(line.split()) == 0:
                    break
            elif "&end" in line:
                keep_lines = True
    return lines

def print_subprocess_output(output_info, outfile=sys.stdout):
    if type(output_info) in [str, bytes]:
        raise ValueError(' Output info must be a string or bytes')
    else:
        info = output_info.decode('-utf-8')
        for line in info.split('\n'):
            print(' \t' + line, file=outfile)

def create_amber_files(pdb, ihfree=1, ioutopt=1, irstrnt=1, outfile=sys.stdout, init_charges=[], net_charge=0):
    n_frames = pdb.getNumFrames()
    resp_lines_all = [[], []]

    for n in range(n_frames):
        print_section('Generating individual resp inputs for fragment', outfile=outfile)
        makedirs('frame_{:d}'.format(n), exist_ok=True)
        pdb_file = 'frame_{0:d}/frame_{0:d}.pdb'.format(n)
        ac_file = 'frame_{0:d}/frame_{0:d}.ac'.format(n)
        with open(pdb_file, 'w') as file:
            pdb.writeModel(pdb.topology, pdb.getPositions(frame=n), file=file)
        
        #   convert from pdb to ac format
        cmd = 'antechamber -i ' +  pdb_file + ' -fi pdb -o ' + ac_file + ' -fo ac -dr no -nc ' + str(int(net_charge))
        print(" Calling Antechamber: ", file=outfile)
        print('\t ' + cmd +  '\n', file=outfile)
        try:
            output = subprocess.check_output(cmd.split())
            print(" ANTECHAMBER OUTPUT: ", file=outfile)
            print_subprocess_output(output)
        except subprocess.CalledProcessError as error:
            print(" ANTECHAMBER FAILED\nError Code: ", error.returncode, '\nOutput: ', error.output.decode('utf-8'), file=outfile)

        #   clean up antechamber files
        rm_file_list = glob.glob('ANTECHAMBER_AC.AC*')
        rm_file_list += glob.glob('ANTECHAMBER_BOND_TYPE.AC*')
        rm_file_list += glob.glob('ATOMTYPE.INF*')
        for file in rm_file_list:
            os.remove(file)


        #   generate resp files from respgen for 2-stage fitting
        for i in (1, 2):
            cmd = 'respgen -i frame_{0:d}/frame_{0:d}.ac -o frame_{0:d}/frame_{0:d}.respin{1:d} -f resp{1:d}'.format(n, i)
            print(" Calling respgen: ", file=outfile)
            print('\t ' + cmd +  '\n', file=outfile)
            try:
                output = subprocess.check_output(cmd.split())
                print(" RESPGEN OUTPUT: ", output.decode('utf-8'), file=outfile)
            except subprocess.CalledProcessError as error:
                print(" RESPGEN FAILED\nError Code: ", error.returncode, '\nOutput: ', error.output.decode('utf-8'), file=outfile)

            resp_lines = crop_respin('frame_{0:d}/frame_{0:d}.respin{1:d}'.format(n, i))
            resp_lines_all[i-1] += resp_lines

            if irstrnt == 2:
                break

        #   combine all resp files into one fitting file
        for i in (1, 2):
            #   Note: SPACE before each line is required! Might be some weird fortran input thing
            with open('resp{:d}.in'.format(i), 'w') as resp_file:
                if i == 1 and len(init_charges) != 0:
                    iqopt = 2
                else:
                    iqopt = i
                resp_file.write(' RESP Stage {:d} input\n'.format(i))
                resp_file.write(' &cntrl\n')
                resp_file.write('  ioutopt={:d}, iqopt={:d}, nmol={:d}, ihfree={:d}, irstrnt={:d}, qwt={:.5f}\n' \
                    .format(ioutopt, iqopt, n_frames, ihfree, irstrnt, 0.0005*i))
                resp_file.write(' &end\n')
                
                for line in resp_lines_all[i-1]:
                    resp_file.write(line)

                #   write intramolecular charge constraints
                #resp_file.write('                Intra and/or inter-molecular charge restraints for atom or group of atoms\n')

                #   write intermolecular constraints
                resp_file.write('                Intermolecular charge equivalencing for similar conformations\n')
                for n in range(pdb.topology._numAtoms):
                    resp_file.write('{:5d}\n'.format(n_frames))
                    for j in range(n_frames):
                        resp_file.write('{:5d}{:5d}'.format(j+1, n+1))
                        if (j+1) % 8 == 0 and (j+1) != n_frames:
                            resp_file.write('\n')
                    resp_file.write('\n')
                
                resp_file.write('\n')

            if irstrnt == 2:
                break

        #   create file of initial charges to use in fitting
        if len(init_charges) != 0:
            with open('q_init', 'w') as file:
                for n, q in enumerate(np.array(init_charges).flatten()):
                    file.write('{:10.6f}'.format(q))
                    if n % 8 == 7:
                        file.write('\n')
                file.write('\n')

        print(" ------------------------------------------------------\n\n", file=outfile)


def create_esp_data_file(qc_esp_files, coords_in_angs, data_file='esp.dat', outfile=sys.stdout):

    print(" Formatting ESP data file for Amber RESP program", file=outfile)
    print(" \tThere are {:d} files to format".format(len(qc_esp_files)), file=outfile)

    if isinstance(qc_esp_files, str):
        qc_esp_files = [qc_esp_files]
        coords = [coords_in_angs]
    
    with open(data_file, 'w') as file:
        for n, esp_file in enumerate(qc_esp_files):
            esp_data = np.loadtxt(esp_file, skiprows=4)
            n_points = len(esp_data)
            n_atoms = len(coords_in_angs[n])
            file.write('{:5d}{:5d}\n'.format(n_atoms, n_points))
            for coord in coords_in_angs[n]:
                unit = 1.8897259886
                file.write('{:16s}{:16.6e}{:16.6e}{:16.6e}\n'.format('', coord[0]*unit, coord[1]*unit, coord[2]*unit))
            for i in range(n_points):
                x, y, z, esp = esp_data[i]
                file.write('{:16.6e}{:16.6e}{:16.6e}{:16.6e}\n'.format(esp, x, y, z))

def call_resp(outfile=sys.stdout, qinit_file=None, n_stages=2):
    print_section("Starting RESP Fitting Program", outfile=outfile)

    stages = [int(i+1) for i in range(n_stages)]
    for i in stages:
        cmd = 'resp -O -i resp{0:d}.in -o resp{0:d}.out -e esp.dat -t q{0:d}.out '.format(i)
        if i == 1 and qinit_file:
            cmd += '-q {:s}'.format(qinit_file)
        if i == 2:
            cmd += '-q q1.out'

        print(" Calling Amber RESP Program:", file=outfile)
        print("\t" + cmd, file=outfile)
        try:
            output = subprocess.check_output(cmd.split())
            print(" RESP OUTPUT: ", output.decode('utf-8'), file=outfile)
        except subprocess.CalledProcessError as error:
            print(" RESP FAILED\nError Code: ", error.returncode, '\nOutput: ', error.output.decode('utf-8'), file=outfile)

    print(" ------------------------------------------------------\n\n", file=outfile)

def is_float(number):
    try:
        float(number)
        return True
    except ValueError:
        return False

def load_charges(charge_file):
    if path.isfile(charge_file):
        with open(charge_file, 'r') as file:
            read_first_line = False
            charges = []
            format_type = 'qchem'
            for line in file.readlines():
                sp = line.split()
                if sp[0] != '#':
                    #   use the first line to determine the format of the provided charge file
                    if not read_first_line:
                        if len(sp) == 1 and is_float(sp[0]):
                            format_type = 'column'
                        elif len(sp) == 3 and is_float(sp[2]):
                            format_type = 'qchem'
                        elif np.prod([is_float(x) for x in sp]) == 1:
                            format_type = 'resp'
                        else:
                            raise ValueError(' Provided charge format cannot be determined')
                        read_first_line = True
                        print(" Provided charge format assumed to be of type '{:s}' ".format(format_type))

                    if format_type == 'qchem':
                        charges.append(float(sp[2]))
                    elif format_type == 'column':
                        charges.append(float(sp[0]))
                    else:
                        charges += [float(x) for x in sp]
        print(" Imported a total of {:d} charges".format(len(charges)))
        return np.array(charges)
    else:
        raise FileNotFoundError("Charge file not found. Exiting program.")

def print_section(message, length=54, outfile=sys.stdout):
    border = " " + "-"*length
    print("\n", file=outfile)
    print(border, file=outfile)
    print(" %s" % message, file=outfile)
    print(border, file=outfile)
    print("\n", file=outfile)

def check_amber_tools(out_file=sys.stdout):
    print("\n Checking for AmberTools executables", file=out_file)

    exe_locs = {'antechamber': None, 'resp': None, 'respgen': None}
    for exe in exe_locs:
        exe_locs[exe] = which(exe)
    
    found_all = True
    for exe, exe_loc in exe_locs.items():
        if exe_loc is None:
            found_all = False
        print("     * {:s} executable located at {:s}".format(exe, exe_loc), file=out_file)

    if not found_all:
        raise SystemError(" Could not find all required AmberTools executables")


    # antechamber_loc = which('antechamber')
    # resp_loc = which('resp')

    # if antechamber_loc is None:
    #     raise SystemError('Cannot find executable for "antechamber"')
    # else:
    #     print("     * antechamber executable located at: {:s}".format(antechamber_loc), file=out_file)

    # if resp_loc is None:
    #     raise SystemError('Cannot find executable for "resp"')
    # else:
    #     print("     * resp executable located at: {:s}".format(resp_loc), file=out_file)
    print("", file=out_file)

def main(file_args, opts, mol_file_loc, outfile):

    scratch_dir = path.abspath(path.curdir)

    #   load pdb file
    print(" Loading Molecule file", file=outfile)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mol = Molecule(mol_file_loc)
    elms = [x.element.symbol for x in mol.topology.atoms()]
    point_gen = ESPPointGenerator()
    
    qc_esp_files = []
    qc_charges   = []
    qc_energies  = []

     #   no fitting is requested, only point generation
    if file_args.pts:
        
        pts_file_split = path.splitext(file_args.pts)
        if pts_file_split[0] == '':
            pts_file_split = '.xyz'
        for n in range(mol.getNumFrames()):
            pos = mol.getPositions(asNumpy=True, frame=n)/angstroms
            esp_points = point_gen.gen_MK_points(elms, pos, outfile=outfile, intervals=opts.vdw_ratios, density=opts.mk_density)*angstrom.conversion_factor_to(bohr)

            #   number each point file with the frame number
            out_file = pts_file_split[0] + f'_{n}' + pts_file_split[1]
            point_gen.print_xyz(out_file)

        return

    if file_args.esp:
        #   load in Q-Chem generated ESP data
        print(" Loading ESP File. Q-Chem will NOT be run", file=outfile)
        print(" NOTE: Program can only use one ESP file at a time.", file=outfile)
        print("       Multiple conformation fitting is not yet implimented.", file=outfile)
        mol._positions = [mol._positions[0]]
        if path.isfile(file_args.esp):
            qc_esp_files = [file_args.esp]
        else:
            raise FileNotFoundError("ESP provided file does not exist")
    else:
        #   create ESP files with Q-Chem
        print("\n\n ------------------------------------------------------", file=outfile)
        print(" Starting ESP generation with Q-Chem", file=outfile)
        print(" ------------------------------------------------------", file=outfile)
        for n in range(mol.getNumFrames()):
            pos = mol.getPositions(asNumpy=True, frame=n)/angstroms
            esp_points = point_gen.gen_MK_points(elms, pos, outfile=outfile, intervals=opts.vdw_ratios, density=opts.mk_density)*angstrom.conversion_factor_to(bohr)
            (esp_file, energy, charges) = create_qchem_job(opts, scratch_dir, 'frame_{:d}'.format(n), esp_points, pos, elms, outfile=outfile)
            esp_file = path.join(scratch_dir, 'frame_{:d}/plot.esp'.format(n))
            qc_esp_files.append(esp_file)
            qc_charges.append(charges)
            qc_energies.append(energy)
        print("\n ------------------------------------------------------", file=outfile)
        print(" All Q-Chem jobs are finished", file=outfile)
        print(" ------------------------------------------------------\n\n", file=outfile)

    #   perform density ESP fitting
    if opts.density_fitting:
        coords = mol.getPositions(True)/angstrom
        nuclei = [atom.element.atomic_number for atom in mol.topology.atoms()]

        print_section("Starting ESP Density Fitting", outfile=outfile)
        #   overwrite std.out; density fitting does not use output file, temp fix.
        sys.stdout = outfile
        density_fitter = DensityFitting(coords, nuclei, qc_esp_files[0], 
                charge=opts.charge, lone_pairs=opts.lone_pairs, 
                intra_constraints=opts.input_sections['intra_constraints'], 
                lone_pair_dist=opts.lone_pairs_dist,
                lone_pair_k=opts.lone_pairs_k, fitting_method=opts.fitting_method)
        density_fitter.run_fitting()
        # try:
        #     density_fitter = DensityFitting(coords, nuclei, qc_esp_files[0], charge=opts.charge, lone_pairs=opts.lone_pairs, intra_constraints=opts.input_sections['intra_constraints'])
        #     density_fitter.run_fitting()
        # except: 
        #     print(" DENSITY FITTING HAS FAILED: ", sys.exc_info()[0], file=outfile)
        sys.stdout = sys.__stdout__
        print_section("Done with ESP Density Fitting", outfile=outfile)

    #   check if PDB file is provided. If not, turn off option
    pdb = None
    if opts.amber_fitting:
        check_amber_tools(outfile)
        if isinstance(mol, PDBFile):
            pdb = mol
        elif file_args.pdb:
            print(" Loading PDB file for Amber charge fitting")
            print(" Positions from -mol will be used instead")
            pdb = PDBFile(file_args.pdb)
            pdb._positions = mol._positions
        else:
            mol_file_name = path.splitext(mol_file_loc)[0]
            out_pdb_file = path.splitext(mol_file_name + '.pdb')
            PDBFile.writeFile(mol.getTopology(), mol.getPositions(), open(out_pdb_file, 'w'))
            pdb = PDBFile(out_pdb_file)
            pdb._positions = mol._positions

    #   perform RESP fitting with AmberTools
    if pdb is not None:
        coords_all = [pdb.getPositions(frame=i, asNumpy=True)/angstrom for i in range(pdb.getNumFrames())]
        create_esp_data_file(qc_esp_files, coords_all, outfile=outfile)
        if file_args.chg:
            #   load in Q-Chem pre-generated charges for fitting analysis
            print(" Loading charge file for analysis only. Q-Chem will NOT be run")
            charges = load_charges(file_args.chg)
            charges = np.reshape(charges, (pdb.getNumFrames(), pdb.topology.getNumAtoms()))
            info = create_amber_files(pdb, outfile=outfile, net_charge=opts.charge, irstrnt=2, init_charges=charges)
            call_resp(outfile=outfile, n_stages=1, qinit_file='q_init')

        else:
            create_amber_files(pdb, outfile=outfile, net_charge=opts.charge)
            call_resp(outfile=outfile)
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ipt', required=False)
    parser.add_argument('-pdb', required=False)
    parser.add_argument('-out', required=False)
    parser.add_argument('-esp', required=False)
    parser.add_argument('-chg', required=False)
    parser.add_argument('-mol', required=True)
    parser.add_argument('-pts', required=False)
    args = parser.parse_args()

    outfile=sys.stdout
    if args.out:
        outfile = open(args.out, 'w')

    #   change file args to absolute paths
    for arg, file_loc in args._get_kwargs():
        if file_loc is not None:
            setattr(args, arg, path.abspath(file_loc))

    print(" STARTING RESP PYTHON PROGRAM", file=outfile)
    #   grep job and Q-Chem options
    opts = Options(args.ipt)
    original_dir = path.abspath(path.curdir)
    print(" Creating scratch directory:", file=outfile)
    mol_file = args.mol
    if opts.name == "":
        opts.name = path.splitext(path.basename(mol_file))[0]
    
    #   check to make sure that final directory doesn't already exist
    final_dir_name = path.abspath(path.join(path.curdir, opts.name + "_resp"))
    if os.path.isdir(final_dir_name):
        print('\n\n ERROR: output directory {:s} already exists. \n Change or set the "name" option in the input file to something other than "{:s}"\n\n'.format(final_dir_name, opts.name))
        raise RuntimeError("Output directory already exists")
    
    #   set up scratch directory
    scratch_dir = path.abspath(path.join(path.curdir, "__" + opts.name + "__"))
    makedirs(scratch_dir, exist_ok=True)
    print('\t' + scratch_dir, file=outfile)

    #   copy mol file to scratch and chagne directory
    print(" Changing directory to scratch", file=outfile)
    mol_file_name = path.basename(args.mol)
    copyfile(mol_file, path.join(scratch_dir, mol_file_name))
    chdir(scratch_dir)

    #   MAIN PROGRAM
    main(args, opts, mol_file_name, outfile)

    #   clean up scratch directory
    chdir(original_dir)
    if args.pts:
        rmtree(scratch_dir)
    else:
        final_dir_name = path.abspath(path.join(path.curdir, opts.name + "_resp"))
        os.rename(scratch_dir, final_dir_name)