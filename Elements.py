from simtk.openmm.app.element import Element
__name2intDict = {
        "H": 1,
        "C": 6,
        "N": 7,
        "O": 8,
        "P": 15,
        "Zn": 30,
        "Se": 34,
        "Lp": 0
    }
__int2nameDict = {
        1: "H",
        6: "C",
        7: "N",
        8: "O",
        15: "P",
        30: "Zn",
        34: "Se",
        0: "Lp"
    }
def name2int(name):
    if name == 'Lp':
        return 0
    return Element.getBySymbol(name).atomic_number
    #return __name2intDict[name]

def int2name(intIn):
    if intIn == 0:
        return 'Lp'
    return Element.getByAtomicNumber(intIn).symbol
    #return __int2nameDict[intIn]

def getRadiiByAtomicNumber(atomic_number, default=2.0):
    return vdw_radii.get(int2name(atomic_number).upper(), default)

def getRadiiBySymbol(symbol, default=2.0):
    return vdw_radii.get(symbol.upper(), default)

vdw_radii = {'H': 1.20, 'HE': 1.20,
         'LI': 1.37, 'BE': 1.45, 'B': 1.45, 'C': 1.50,
         'N': 1.50, 'O': 1.40, 'F': 1.35, 'NE': 1.30,
         'NA': 1.57, 'MG': 1.36, 'AL': 1.24, 'SI': 1.17,
         'P': 1.80, 'S': 1.75, 'CL': 1.70,
         'FE': 1.60} #  made up  