"""
"""

import mayavi
from mayavi import mlab
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix

VDW_RADIUS = {'H':1.20, 'C':1.7, 'N':1.55, 'O':1.52, 'CL':1.75, 'F':1.47, 'P':1.80,
'S':1.80, 'CU':1.40, 'HE':1.40, 'LI':1.82, 'BE':1.53, 'B':1.92, 'NE':1.54, 'NA':2.27,
'MG':1.73, 'AL':1.84, 'SI':2.10, 'AR':1.88, 'K':2.75, 'CA':2.31, 'SC':2.11, 'NI':1.63,
'ZN':1.39, 'GA':1.87, 'GE':2.11, 'AS':1.85, 'SE':1.90, 'BR':1.85, 'KR':2.02, 'RB':3.03,
'SR':2.49, 'PD':1.63, 'AG':1.72, 'CD':1.58, 'IN':1.93, 'SN':2.17, 'SB':2.06, 'TE':2.06,
'I':1.98, 'XE':2.16, 'CS':3.43, 'BA':2.68, 'PT':1.75, 'AU':1.66, 'HG':1.55, 'TL':1.96,
'PB':2.02, 'BI':2.07, 'PO':1.97, 'AT':2.02, 'RN':2.20, 'FR':3.48, 'RA':2.83, 'U':1.86}

COL_ATM = {'H':(0.9,0.9,0.9), 'C':(0.1,0.1,0.1), 'N':(0,0.5,0.9), 'O':(0.9,0,0), 'S':(0.9,0.9,0), 
'CL':(0.5,0.5,0), 'BR':(0.6,0.1,0.1), 'I':(0.6,0,0.7)}


def atom_spheres(atoms_df):
    """
    """
    sph_lst=[]
    for row in atoms_df.iterrows():

        radius = VDW_RADIUS[row[1][0]]

        if row[1][0] in COL_ATM :
            col = COL_ATM[row[1][0]]
        else :
            col = (0.5,0.5,0.5)

        x_coor = row[1][3]
        y_coor = row[1][4]
        z_coor = row[1][5]

        [phi,theta] = np.mgrid[0:2*np.pi:12j,0:np.pi:12j]

        x = radius * np.cos(phi) * np.sin(theta) + x_coor
        y = radius * np.sin(phi) * np.sin(theta) + y_coor
        z = radius * np.cos(theta) + z_coor

        sph_lst.append((x, y, z, col, row[1][0] + row[1][1])) #, row[1][0] + row[1][1]

    return sph_lst


def visual_spacefilling(sph_lst):
    """
    """
    for sph in sph_lst:
        mlab.mesh(sph[0], sph[1], sph[2], color = sph[3], name = sph[4])



def xyz_ca_cb(atoms_df, res_name, i):
    """
    """
    dic = {}
    dic['x'] = []
    dic['y'] = []
    dic['z'] = []
    
    cc = 3
    for c in ['x', 'y', 'z'] :
        for a in atoms_df[(atoms_df["res nbr"]==atoms_df.iloc[i+1,2])&(atoms_df["ap"]==res_name)&(atoms_df["chain"]==atoms_df.iloc[i+1,7])][c]:
            x = a
        dic[c].append(x)
        dic[c].append(atoms_df.iloc[i+1,cc])
        cc += 1

    return dic


def visual_ballstick(atoms_df):
    """
    """
    for row in atoms_df.iterrows():
        if row[1][0] in COL_ATM :
            col = COL_ATM[row[1][0]]
        else :
            col = (0.5,0.5,0.5)
        mlab.points3d(row[1][3], row[1][4], row[1][5], color = col, scale_factor = 0.9, scale_mode = 'scalar')

    ## wip
    for i in range(len(atoms_df)-2):
        if atoms_df.iloc[i,2] == atoms_df.iloc[i+1,2]:
            if atoms_df.iloc[i+1,6] != "CB":
                mlab.plot3d(atoms_df.iloc[i:i+2,3], atoms_df.iloc[i:i+2,4], atoms_df.iloc[i:i+2,5], color = (0.5,0.5,0.5), tube_radius = 0.2)
            else :
                di = xyz_ca_cb(atoms_df, "CA", i)
                mlab.plot3d(di['x'], di['y'], di['z'], color = (0.5,0.5,0.5), tube_radius = 0.2)
        elif atoms_df.iloc[i,7] != atoms_df.iloc[i+1,7]:
            i += 1 
        else :
            di = xyz_ca_cb(atoms_df, "O", i)
            mlab.plot3d(di['x'], di['y'], di['z'], color = (0.5,0.5,0.5), tube_radius = 0.2)


def res_names(atoms_df):
    """
    """
    for row in atoms_df.iterrows():
        mlab.text3d(row[1][3], row[1][4], row[1][5], text = row[1][1] + " " + str(row[1][2]))



def visual_lines(atoms_df):
    """
    """
    for i in range(len(atoms_df)-2):
        if atoms_df.iloc[i,2] == atoms_df.iloc[i+1,2]:
            if atoms_df.iloc[i+1,6] != "CB":
                mlab.plot3d(atoms_df.iloc[i:i+2,3], atoms_df.iloc[i:i+2,4], atoms_df.iloc[i:i+2,5], color = (0,1,0), tube_radius = 0.2)
            else :
                di = xyz_ca_cb(atoms_df, "CA", i)
                mlab.plot3d(di['x'], di['y'], di['z'], color = (0,1,0), tube_radius = 0.2)
        elif atoms_df.iloc[i,7] != atoms_df.iloc[i+1,7]:
            i += 1 
        else :
            di = xyz_ca_cb(atoms_df, "O", i)
            mlab.plot3d(di['x'], di['y'], di['z'], color = (0,1,0), tube_radius = 0.2)


def visual_ribbon(atoms_df):
    """
    """
    df = atoms_df[atoms_df['ap'] == 'CA']
    col = 0.2
    for i in range(len(df.chain.unique())):
        c = df.chain.unique()[i]
        mlab.plot3d(df[df["chain"]==c]["x"], df[df["chain"]==c]['y'], df[df["chain"]==c]['z'], color = (col,col+0.2,col+0.3), tube_radius = 0.2)
        col += 0.3
    #mlab.plot3d(df['x'], df['y'], df['z'], color = (0,0.6,0.8), tube_radius = 0.1)


#https://docs.enthought.com/mayavi/mayavi/auto/example_pick_on_surface.html
def picker_callback(picker):
    """ Picker callback: this get called when on pick events.
    """
    #print(picker.picked_positions)
    print("*********************")
    print("mapper : ", picker.mapper_position)
    print("coord  : ", df.iloc[picker.point_id,3], df.iloc[picker.point_id,4], df.iloc[picker.point_id,5])
    #print(picker.point_id, atoms_df.iloc[picker.point_id,0])
    #print(picker.MapperPosition, atoms_df.iloc[picker_id,3], atoms_df.iloc[picker_id,4], atoms_df.iloc[picker_id,5])

#fig = mlab.figure(1)
#fig.on_mouse_pick(picker_callback)


def atom_dist_matrix(atoms_df):
    """
    """
    return pd.DataFrame(distance_matrix(atoms_df.iloc[:,3:6],atoms_df.iloc[:,3:6]), index=atoms_df.iloc[:,3:6].index, columns=atoms_df.iloc[:,3:6].index)


def threshold_dict(mtx, atoms_df):
    """
    """
    r = []
    i = []
    TRS_MAX = 4 # 4A

    for index, row in mtx.iterrows():
        r.append(row)
        i.append(index)

    dico_trsh = {}
    for a in range(len(i)):
        lst = []
        for b in range(len(r[a])-1):
            if (r[a][b] < TRS_MAX) & (atoms_df.iloc[a,2] != atoms_df.iloc[b,2]) & (atoms_df.iloc[a,6] != "CA") & (atoms_df.iloc[b,6] != "CA"):
                lst.append(b)
        dico_trsh[a] = lst

    return dico_trsh


def visual_dist(dico_trsh, atoms_df):
	"""
	"""
	for key in dico_trsh:
		for i in dico_trsh[key]:
			if i :
				mlab.plot3d(atoms_df.iloc[[key,i],3], atoms_df.iloc[[key,i],4], atoms_df.iloc[[key,i],5], color = (1,1,0), tube_radius = 0.05, opacity = 0.5)#, representation = 'wireframe')


if __name__ == "__main__":
    import firsttry
    print(help(firsttry))

#import coor_atom as ca
#import firsttry as ft
#df=ca.coord("../data/1bzv.pdb")
#ft.visual_ballstick(df)
#s=ft.atom_spheres(df)
#ft.visual_spacefilling(s)