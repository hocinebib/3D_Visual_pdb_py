import mayavi
from mayavi import mlab
import numpy as np

VDW_RADIUS = {'H':1.20, 'C':1.7, 'N':1.55, 'O':1.52, 'CL':1.75, 'F':1.47, 'P':1.80,
'S':1.80, 'CU':1.40, 'HE':1.40, 'LI':1.82, 'BE':1.53, 'B':1.92, 'NE':1.54, 'NA':2.27,
'MG':1.73, 'AL':1.84, 'SI':2.10, 'AR':1.88, 'K':2.75, 'CA':2.31, 'SC':2.11, 'NI':1.63,
'ZN':1.39, 'GA':1.87, 'GE':2.11, 'AS':1.85, 'SE':1.90, 'BR':1.85, 'KR':2.02, 'RB':3.03,
'SR':2.49, 'PD':1.63, 'AG':1.72, 'CD':1.58, 'IN':1.93, 'SN':2.17, 'SB':2.06, 'TE':2.06,
'I':1.98, 'XE':2.16, 'CS':3.43, 'BA':2.68, 'PT':1.75, 'AU':1.66, 'HG':1.55, 'TL':1.96,
'PB':2.02, 'BI':2.07, 'PO':1.97, 'AT':2.02, 'RN':2.20, 'FR':3.48, 'RA':2.83, 'U':1.86}

COL_ATM = {'H':(1,1,1), 'C':(0,0,0), 'N':(0,0.5,1), 'O':(1,0,0), 'S':(1,1,0), 'CL':(0.5,0.5,0), 
'BR':(0.6,0.1,0.1), 'I':(0.6,0,0.7)}

R = VDW_RADIUS[row[1][0]]
a = row[1][3]
b = row[1][4]
c = row[1][5]

[phi,theta] = np.mgrid[0:2*np.pi:12j,0:np.pi:12j]

x = R*np.cos(phi)*np.sin(theta) + a
y = R*np.sin(phi)*np.sin(theta) + b
z = R*np.cos(theta) + c

x1 = x + 1
y1 = y + 1
z1 = z + 1

mlab.mesh(x, y, z, color = COL_ATM[row[1][0]])
mlab.mesh(x1,y1,z1)