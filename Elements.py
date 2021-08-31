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
        return __name2intDict[name]

def int2name(intIn):
    return __int2nameDict[intIn]