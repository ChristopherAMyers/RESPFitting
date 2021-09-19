class IntraConstraint:
    def __init__(self, charge=0.0, mol_idx=0) -> None:
        self.charge = charge
        self.mol_idx = mol_idx
        self.constr_dict = {}
        
    def add_constraint(self, coeff, atom_idx_list):
        try:
            atom_list = list(atom_idx_list)
            for atom in atom_list:
                self.constr_dict[int(atom)] = coeff
        except TypeError:
            raise TypeError("add_constraint: atom_idx_list must be list-like, but was passed %s" % str(type(atom)))