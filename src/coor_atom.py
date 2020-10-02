"""
"""
import pandas as pd

def coord(pdb_file):
    """
    """
    with open(pdb_file, "r") as pdb_f:
        lst = []
        for line in pdb_f :
            dic = {}
            if line[0:4] == "ATOM":
                dic["atom"] = str(line[77:79].strip())
                dic["res"] = str(line[17:21].strip())
                dic["res nbr"] = int(line[22:26].strip())
                dic["x"] = float(line[30:38].strip())
                dic["y"] = float(line[38:46].strip())
                dic["z"] = float(line[46:54].strip())
                lst.append(dic)
        atoms_df = pd.DataFrame(lst)
    return atoms_df


if __name__ == "__main__":
    import coor_atom
    print(help(coor_atom))