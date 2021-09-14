from argparse import FileType
from simtk.openmm.app import PDBFile, Topology, element
from simtk.openmm import Vec3
from simtk.unit import nanometers, angstroms, is_quantity, norm, Quantity, dot
import numpy as np
from os.path import *


class XYZFile():
    def __init__(self, file_loc) -> None:
        ''' Similar to the PDBFile object from OpenMM, but for XYZ
        files instead. This class also creates a topology object and
        is suited for multiple frames.
        
        Parameters
        ----------
        file : string
            the name of the file to load
        '''
        
        self._positions = []
        self._atoms = []

        #   load in atomic coordinates
        with open(file_loc, 'r') as file:
            n_atoms = 0
            coords = []
            for n, line in enumerate(file.readlines()):
                #   establish number of atoms from first line only
                if n == 0:
                    n_atoms = int(line.split()[0])

                if n % (n_atoms + 2) == 0:
                    new_n_atoms = int(line.split()[0])
                    #   all frames must have the same no. atoms
                    if new_n_atoms != n_atoms:
                        raise ValueError("All frames in xyz file must have the same no. of atoms")
                    if n != 0:
                        self._positions.append(coords*angstroms)
                    coords = []
                elif n % (n_atoms + 2) >= 2:
                    sp = line.split()
                    if n < (n_atoms + 2):
                        self._atoms.append(sp[0])
                    pos = [float(x) for x in sp[1:]]
                    coords.append(Vec3(pos[0], pos[1], pos[2]))

        self._positions.append(coords*angstroms)
        self.positions = self._positions[0]
        self._numpyPositions = None

        #   define a topology object
        top = Topology()
        self.topology = top

        chain = top.addChain()
        res = top.addResidue('mol', chain)
        for n, atom in enumerate(self._atoms):
            elm = element.get_by_symbol(atom)
            top.addAtom(str(n + 1) + atom, elm, res)

    def getTopology(self):
        """Get the Topology of the model."""
        return self.topology

    def getNumFrames(self):
        """Get the number of frames stored in the file."""
        return len(self._positions)

    def getPositions(self, asNumpy=False, frame=0):
        """Get the atomic positions.

        Parameters
        ----------
        asNumpy : boolean=False
            if true, the values are returned as a numpy array instead of a list
            of Vec3s
        frame : int=0
            the index of the frame for which to get positions
        """
        if asNumpy:
            if self._numpyPositions is None:
                self._numpyPositions = [None]*len(self._positions)
            if self._numpyPositions[frame] is None:
                self._numpyPositions[frame] = Quantity(np.array(self._positions[frame].value_in_unit(nanometers)), nanometers)
            return self._numpyPositions[frame]
        return self._positions[frame]


class Molecule():
    ''' Molecule importer for either a PDB or XYZ files.
        This class retures either a PDBFile or XYZFile molecule object, depending
        on the file extension provided.Each have the same methods, less the
        ability to write PDB files for the XYZ object.

        Parameters
        ----------
        file_loc : string
            The location of the file to load
        
        Returns
        -------
        mol: PDBFile or XYZFile object
            Depends in the file extension used for file_loc'''

    def __new__(self, file_loc):
        extension = splitext(file_loc)[-1]
        if extension == '.pdb':
            return PDBFile(file_loc)
        elif extension == '.xyz':
            return XYZFile(file_loc)
        else:
            raise FileType('Only .pdb or .xyz molecules files are allowed')

    def getPDB(self):
        if isinstance(self._mol, PDBFile):
            return self._mol
        else:
            raise TypeError("Molecule() was not defined as a PDBFile")




