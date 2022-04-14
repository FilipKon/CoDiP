import logging
import numpy as np
import pandas as pd
import Bio
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pylab import figure
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import FastICA
from sklearn.cluster import KMeans
import plotly.graph_objs as go
import sys
import plotly.express as px

GlyR = ['7MLU', '7L31', '7KUY', '5BKG', '5BKF', '6PLS', '6PLT', '6PLU', '6PLV', '6PLR', '6PLO', '6PLQ', '6PLP',
        '6PLZ', '6PLW', '6PLY', '6PLX', '6PM0', '6PM1', '6PM2', '6PM3', '6PM4', '6PM5', '6PM6', '4X5T', '5VDH',
        '6PXD', '6UBT', '6UD3', '6UBS', '6VM0', '6VM2', '6VM3', '5TIO', '5TIN', '5VDI', '5VDH', '5CFB', '3JAE',
        '3JAD', '3JAF', '7MLY']
nachr = ['4BOT', '4BON', '4BOO', '4BOI', '4BOR', '4AQ5', '4AQ9', '2BG9', '7KOO', '7KOQ', '7KOX', '6UR8', '6CNK',
         '6CNJ', '5KXI', '7EKT', '7EKI', '7EKP', '6USF', '6PV7', '6PV8', '6UWZ', '7QL6', '7QL5', '7QKO']

antagonists_bound = ['6HIS', '6NP0', '6Y1Z', '6W1Y', '6W1J', '6W1M', '7KOO', '6UWZ[D-E]', '7L31',
                     '7KUY', '5CFB', '3JAD', '7PC0[E-D]', '7PC0[B-A]', '6UWZ[A-B]', '6HUK[B-A]', '6X3S[C-D]',
                     '6X3S[A-B]']
agonist_bound = ['6HIO', '6HIQ', '6HIN', '6DG8', '6DG7', '6BE1', '4PIR', '7QL5[D-C]', '7QL5[A-E]', '7QL6[A-E]',
                 '7EKP', '7KOQ', '7QL6[D-C]', '6USF[D-E]', '6PV7[D-E]', '6PV8[A-B]', '6CNK[D-E]',
                 '6CNJ[D-E]', '5KXI[D-E]', '5KXI[A-B]', '7MLY', '5BKG[E-A]', '5BKG[A-B]', '5BKG[B-C]', '5BKF[C-D]',
                 '6PLS', '5BKF[D-E]', '5BKF[E-A]', '5BKF[A-B]', '5BKF[B-C]', '6Y5A',
                 '6PLT', '6PLU', '6PLV', '6PLR', '6PLO', '6PLQ', '6PLP', '6PLZ', '6PLW', '6PLY', '6PLX', '6UBT',
                 '3JAE', '6PM0', '6PM1', '6PM2', '6PM3', '6PM4', '6PM5', '6PM6', '6USF[A-B]', '6PV7[A-B]',
                 '6X3Z[A-B]', '4COF', '6A96[B-A]', '6DW0[E-A]', '6DW0[A-B]', '6DW0[B-C]', '6PV8[D-E]', '7KOX',
                 '6X3Z[C-D]', '6UR8[D-E]', '6UR8[A-B]', '6CNK[A-B]', '7PBD[B-A]', '7PBD[E-D]', '6CNJ[A-B]']

PAM = ['5TIO', '5TIN', '5VDI', '5VDH', '3JAF']
channel_blocker = ['6UD3']
chimeras = ['4X5T', '5OSA', '5OSB', '5OSC', '5O8F', '6CDU', '6D1S', '5OJM']


GABAAR = ['6I53', '6A96', '6CDU', '6D1S', '5OSA', '5OSB', '5OSC', '4COF', '6HUP', '6HUG', '6HUO', '6HUJ', '6HUK',
          '6QFA', '6X3X', '6DW0', '6DW1', '6X40', '6X3S', '6X3V', '6X3T', '5O8F', '5OJM', '6D6U', '6D6T', '6X3Z',
          '6X3U', '6X3X', '6X3W', '7PC0', '7PBZ', '7PBD']

all = ['6HIO', '6HIQ', '6HIN', '6HIS', '6NP0', '6Y59', '6Y1Z', '6DG8', '6DG7', '6BE1', '6W1Y', '6W1J', '6W1M',
       '4PIR', '4BOT', '4BON', '4BOO', '4BOI', '4BOR', '4AQ5', '4AQ9', '2BG9', '6PLS', '6PLT', '6PLU', '6PLV',
       '6PLR', '6PLO', '6PLQ', '6PLP', '6PLZ', '6PLW', '6PLY', '6PLX', '6PM0', '6PM1', '6PM2', '6PM3', '6PM4',
       '6PM5', '6PM6', '4X5T', '6I53', '6A96', '6CDU', '6D1S', '5OSA', '5OSB', '5OSC', '4COF', '6HUP', '6HUG',
       '6HUO', '6HUJ', '6HUK', '6PV7', '6PV8', '6UWZ', '6QFA', '6X3X', '7KOO', '5VDH', '6Y5B', '6Y5A',
       '7KOO', '7KOQ', '7KOX', '6UR8', '6CNK', '6CNJ', '5KXI', '6PXD', '6UBT', '6UD3', '6UBS', '6VM0', '6VM2',
       '6VM3', '5TIO', '5TIN', '5VDI', '5VDH', '5CFB', '3JAE', '3JAD', '3JAF', '6DW0', '6DW1', '6X40', '6X3S',
       '6X3V', '6X3T', '5O8F', '5OJM', '6D6U', '6D6T', '6X3Z', '6X3U', '6X3X', '6X3W', '5BKG', '7KUY', '7L31', '7MLY',
       '7MLU', '7MLU', '7EKP']


clockwise = ['6HIO', '6HIQ', '6HIN', '6HIS', '6NP0', '6Y59', '6Y1Z', '6DG8', '6DG7', '6BE1', '6W1Y', '6W1J', '6W1M',
             '4PIR', '4BOT', '4BON', '4BOO', '4BOI', '4BOR', '4AQ5', '4AQ9', '2BG9', '6PLS', '6PLT', '6PLU', '6PLV',
             '6PLR', '6PLO', '6PLQ', '6PLP', '6PLZ', '6PLW', '6PLY', '6PLX', '6PM0', '6PM1', '6PM2', '6PM3', '6PM4',
             '6PM5', '6PM6', '4X5T', '6I53', '6A96', '6CDU', '6D1S', '5OSA', '5OSB', '5OSC', '4COF', '6HUP', '6HUG',
             '6HUO', '6HUJ', '6HUK', '6QFA', '7PC0', '7PBZ', '7PBD', '7QKO', '7QL5', '7QL6']
counter = ['6UWZ', '6X3X', '7KOO', '5VDH', '6Y5B', '6Y5A', '7KOO', '7KOQ', '7KOX', '6UR8', '6CNK',
           '6CNJ', '5KXI', '6PXD', '6UBT', '6UD3', '6UBS', '6VM0', '6VM2', '6VM3', '5TIO', '5TIN', '5VDI', '5VDH',
           '5CFB', '3JAE', '3JAD', '3JAF', '6DW0', '6X40', '6X3S', '6X3V', '6X3T', '5O8F', '5OJM', '6D6U',
           '6D6T', '6X3Z', '6X3U', '6X3X', '6X3W', '5BKF', '7EKT', '7EKI', '7EKP', '7MLU', '7MLY', '7L31', '7KUY',
           '5BKG', '6PV7', '6PV8', '6USF', '6PLR']

ligand_bound_bs3 = ['6X3W[B-C]', '6X3X[A-B]', '6HUP[B-A]', '6X3V[C-D]', '3JAF[A-B]', '3JAF[B-C]', '3JAF[C-D]',
                    '6X3W[E-A]', '6X3X[E-A]', '6X3X[C-D]', '6HUP[E-D]', '6X3V[A-B]',
                    '7EKT[E-A]', '7EKT[C-D]', '7EKT[A-B]', '7EKT[B-C]', '7EKT[D-E]', '5VDH[A-B]', '5VDH[B-C]',
                    '5VDH[C-D]', '5VDH[D-E]', '5VDH[E-A]', '5VDI[E-A]', '5VDI[A-B]', '5VDI[B-C]', '5VDI[C-D]',
                    '5VDI[D-E]', '3JAF[D-E]', '3JAF[E-A]', '6VM3[E-A]', '6VM3[A-B]', '6VM3[B-C]', '6VM3[C-D]',
                    '6VM3[D-E]', '6VM0[A-B]', '6VM0[B-C]', '6VM0[C-D]', '6VM0[D-E]', '6VM0[E-A]',
                    '6VM2[A-B]', '6VM2[B-C]', '6VM2[C-D]', '6VM2[D-E]', '6VM2[E-A]', '6X3T[C-D]', '6X3T[A-B]']

ligand_notbound_bs3 = ['6HUP[D-C]', '6HUP[C-B]', '6HUP[A-E]', '6X3W[C-D]', '6X3W[D-E]', '6X3W[A-B]', '6X3X[B-C]',
                       '6X3X[D-E]', '6X3V[B-C]', '6X3V[D-E]', '6X3V[E-A]', '6X3T[B-C]', '6X3T[D-E]', '6X3T[E-A]']

ligand_bound_bs4 = ['5OSB', '6CDU', '5O8F']

ligand_bound_bs5 = ['5OSC']

ht3r = ['6HIO', '6HIQ', '6HIN', '6HIS', '6NP0', '6Y59', '6Y1Z', '6DG8', '6DG7', '6BE1', '6W1Y', '6W1J', '6W1M',
        '4PIR', '6Y5A', '6Y5B']

PDBids = ['4BOT', '4BON', '4BOO', '4BOI', '4BOR', '4AQ5', '4AQ9', '2BG9', '6PLS', '6PLT', '6PLU', '6PLV',
          '6PLR', '6PLO', '6PLQ', '6PLP', '6PLZ', '6PLW', '6PLY', '6PLX', '6PM0', '6PM1', '6PM2', '6PM3', '6PM4',
          '6PM5', '6PM6', '4X5T', '6I53', '6A96', '6CDU', '6D1S', '5OSA', '5OSB', '5OSC', '4COF', '6HUP', '6HUG',
          '6HUO', '6HUJ', '6HUK', '6PV7', '6PV8', '6UWZ', '6QFA', '6X3X', '7KOO', '5VDH', '6USF',
          '7KOO', '7KOQ', '7KOX', '6UR8', '6CNK', '6CNJ', '5KXI', '6PXD', '6UBT', '6UD3', '6UBS', '6VM0', '6VM2',
          '6VM3', '5TIO', '5TIN', '5VDI', '5VDH', '5CFB', '3JAE', '3JAD', '3JAF', '6DW0', '6DW1', '6X40', '6X3S',
          '6X3V', '6X3T', '5O8F', '5OJM', '6D6U', '6D6T', '6X3Z', '6X3U', '6X3X', '6X3W', '7EKT', '7EKI', '7EKP',
          '7MLU', '7L31', '7KUY', '5BKG', '5BKF', '7MLY', '7PBD', '7PBZ', '7PC0', '7QKO', '7QL5', '7QL6',
          '6HIO', '6HIQ', '6HIN', '6HIS', '6NP0', '6Y59', '6Y1Z', '6DG8', '6DG7', '6BE1', '6W1Y', '6W1J', '6W1M',
          '4PIR', '6Y5A', '6Y5B', '6QFA']

nachrs_muscle = ['2BG9', '4AQ9', '4AQ5', '4BOI', '4BON', '4BOO', '4BOT', '4BOR', '6UWZ', '7QL6', '7QL5', '7QKO', '7SMT']
nachrs_heteros = ['5KXI', '6CNJ', '6CNK', '6UR8', '6PV8', '6PV7', '6USF']
nachrs_homos = ['7KOX', '7KOQ', '7KOO', '7EKP', '7EKT', '7EKI']
"""
        elif pdbid in counter:
            if '6DW0' in pdbid:
                # B - C - D - E - A
                chain_A = first_model["B"]
                chain_B = first_model["C"]
                chain_C = first_model["D"]
                chain_D = first_model["E"]
                chain_E = first_model["A"]
                chain_numbering = [chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr, chain_Anr]
                type = 3
            else:
                chain_A = first_model["A"]
                chain_B = first_model["B"]
                chain_C = first_model["C"]
                chain_D = first_model["D"]
                chain_E = first_model["E"]
                type = 4
                chain_numbering = [chain_Anr, chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr]
"""


if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def distance(residue1, residue2):
    ca1 = residue1["CA"]
    ca2 = residue2["CA"]
    distance = ca1 - ca2
    return distance


def angles(residue1, residue2, residue3):
    atom1 = residue1["CA"]
    atom2 = residue2["CA"]
    atom3 = residue3["CA"]
    vector1 = atom1.get_vector()
    vector2 = atom2.get_vector()
    vector3 = atom3.get_vector()
    angle = Bio.PDB.vectors.calc_angle(vector1, vector2, vector3)
    return angle


def read_pdb_xyz(pdb_name):
    """
    Reads atomic coordinates of C-alpha atoms from a .pdb file.

    Parameters
    ----------
    pdb_name : str
        Name of pdb file.

    Returns
    -------
    array of atomic coordinates
        [[x1 y1 z1]
         [x2 y2 z2]
         [.. .. ..]
         [xn yn zn]]
    """
    xyz = []
    with open(pdb_name, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM"):
                # extract x, y, z coordinates for carbon alpha atoms
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                if line[12:16].strip() == "CA":
                    xyz.append([x, y, z])
    return xyz


def create_axis_1(pdbid):
    filename = 'data//Structures//' + pdbid + '.pdb'
    #filename = 'C:\\Users\\TARIQOPLATA\\PycharmProjects\\Con_pLGIC\\data\\Structures\\' + pdbid + '.pdb'
    axis_coords_h = []
    axis_coords_v1 = []
    axis_coords_v2 = []
    xyz = read_pdb_xyz(filename)
    coord = np.array(xyz, float)
    center = np.mean(coord, 0)
    coord = coord - center
    inertia = np.dot(coord.transpose(), coord)
    e_values, e_vectors = np.linalg.eig(inertia)
    order = np.argsort(e_values)
    eval3, eval2, eval1 = e_values[order]
    axis3, axis2, axis1 = e_vectors[:, order].transpose()
    i = -2.5
    while i <= 5:
        point_h = i * 20 * axis1 + center
        axis_coords_h.append(point_h)
        # the medium vector is the second principal axis
        point_v1 = i * 20 * axis2 + center
        axis_coords_v1.append(point_v1)
        # the small vector is the third principal axis
        point_v2 = i * 20 * axis3 + center
        axis_coords_v2.append(point_v2)
        i += 0.1
    y = 0
    for item in axis_coords_h:
        #aSetPos[dum_a[1], [62.473, 64.107, 157.783]];
        t = 'aSetPos[dum_a[' + str(y) + '], ' + str(item) + '];'
        print(t)
        y += 1


def PS_site(pdbid, chain_A, chain_B, chain_C, chain_D, chain_E, PDBids):
    df_A = []
    df_B = []
    df_C = []
    df_D = []
    df_E = []
    df_A_nr = []
    df_B_nr = []
    df_C_nr = []
    df_D_nr = []
    df_E_nr = []
    if pdbid.upper() in PDBids:
        parser = PDBParser()
        #filename = '/home/ftk/PycharmProjects/Con_pLGIC_new/data/Structures/' + pdbid + '.pdb'
        filename = 'C:\\Users\\TARIQOPLATA\\PycharmProjects\\Con_pLGIC\\data\\Structures\\' + pdbid + '.pdb'
        #filename = 'C:\\Users\\filip\\PycharmProjects\\plgic\\Structures\\' + pdbid + ".pdb"
        structure = parser.get_structure(pdbid, filename)
        first_model = structure[0]
        chainA = first_model["A"]
        chainB = first_model["B"]
        chainC = first_model["C"]
        chainD = first_model["D"]
        chainE = first_model["E"]
        chains = [chainA, chainB, chainC, chainD, chainE]
        chainsnrs = [chain_A, chain_B, chain_C, chain_D, chain_E]
        i = 0
        j = 0
        x = 0
        for item in chainsnrs:
            if x == 0:
                txt = pdbid + '_chain_A'
                df_A.append(txt)
                df_A_nr.append(txt)
            elif x == 1:
                txt = pdbid + '_chain_B'
                df_B.append(txt)
                df_B_nr.append(txt)
            elif x == 2:
                txt = pdbid + '_chain_C'
                df_C.append(txt)
                df_C_nr.append(txt)
            elif x == 3:
                txt = pdbid + '_chain_D'
                df_D.append(txt)
                df_D_nr.append(txt)
            elif x == 4:
                txt = pdbid + '_chain_E'
                df_E.append(txt)
                df_E_nr.append(txt)
            i = 0
            while i < len(item):
                j = i + 1
                while j < len(item):
                    print([pdbid, item[i], item[j]])
                    residue1 = chains[x][int(item[i])]
                    residue2 = chains[x][int(item[j])]
                    d1 = distance(residue1, residue2)
                    #if j == (len(item) - 1):
                    #    residue2 = chains[x][int(item[0])]
                    #    residue3 = chains[x][int(item[1])]
                    #elif j == (len(item) - 2):
                    #    residue3 = chains[x][int(item[0])]
                    #else:
                    #    residue3 = chains[x][int(item[j+2])]
                    #a1 = angles(residue1, residue2, residue3)
                    #print([chains[x], item[i], item[j]])
                    if x == 0:
                        t = str(chains[x]) + str(residue1) + ' - ' + str(chains[x]) + str(residue2)
                        df_A_nr.append(t)
                        df_A.append(d1)
                        #df_A.append(a1)
                    elif x == 1:
                        t = str(chains[x]) + str(residue1) + ' - ' + str(chains[x]) + str(residue2)
                        df_B_nr.append(t)
                        df_B.append(d1)
                        #df_B.append(a1)
                    elif x == 2:
                        t = str(chains[x]) + str(residue1) + ' - ' + str(chains[x]) + str(residue2)
                        df_C_nr.append(t)
                        df_C.append(d1)
                        #df_C.append(a1)
                    elif x == 3:
                        t = str(chains[x]) + str(residue1) + ' - ' + str(chains[x]) + str(residue2)
                        df_D_nr.append(t)
                        df_D.append(d1)
                        #df_D.append(a1)
                    elif x == 4:
                        t = str(chains[x]) + str(residue1) + ' - ' + str(chains[x]) + str(residue2)
                        df_E_nr.append(t)
                        df_E.append(d1)
                        #df_E.append(a1)
                    j += 1
                i += 1
            x += 1
    return df_A, df_B, df_C, df_D, df_E, df_A_nr, df_B_nr, df_C_nr, df_D_nr, df_E_nr


def upperTMD_site(pdbid, chain_Anr, chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr, PDBids):
    df_A = []
    df_B = []
    df_C = []
    df_D = []
    df_E = []
    df_A_nr = []
    df_B_nr = []
    df_C_nr = []
    df_D_nr = []
    df_E_nr = []
    if pdbid.upper() in PDBids:
        parser = PDBParser()
        filename = 'C:\\Users\\TARIQOPLATA\\PycharmProjects\\Con_pLGIC\\data\\Structures\\' + pdbid + '.pdb'
        #filename = '/home/ftk/PycharmProjects/Con_pLGIC_new/data/Structures/' + pdbid + '.pdb'
        #filename = 'C:\\Users\\filip\\PycharmProjects\\plgic\\Structures\\' + pdbid + ".pdb"
        structure = parser.get_structure(pdbid, filename)
        first_model = structure[0]
        if pdbid in clockwise:
            if '6HU' in pdbid or '6I53' in pdbid or '6A96' in pdbid or '4COF' in pdbid:
                chain_A = first_model["B"]
                chain_B = first_model["A"]
                chain_C = first_model["E"]
                chain_D = first_model["D"]
                chain_E = first_model["C"]
                chain_numbering = [chain_Bnr, chain_Anr, chain_Enr, chain_Dnr, chain_Cnr]
                type = 1
            else:
                chain_A = first_model["A"]
                chain_B = first_model["E"]
                chain_C = first_model["D"]
                chain_D = first_model["C"]
                chain_E = first_model["B"]
                chain_numbering = [chain_Anr, chain_Enr, chain_Dnr, chain_Cnr, chain_Bnr]
                type = 2
        elif pdbid in counter:
            #if '6DW0' in pdbid:
                # B - C - D - E - A
            #    chain_A = first_model["B"]
            #    chain_B = first_model["C"]
            #    chain_C = first_model["D"]
            #    chain_D = first_model["E"]
            #    chain_E = first_model["A"]
            #    chain_numbering = [chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr, chain_Anr]
            #    type = 3
            #else:
            chain_A = first_model["A"]
            chain_B = first_model["B"]
            chain_C = first_model["C"]
            chain_D = first_model["D"]
            chain_E = first_model["E"]
            type = 4
            chain_numbering = [chain_Anr, chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr]
        else:
            print("Failed to find correct direction")
            print(pdbid)
        chains = [chain_A, chain_B, chain_C, chain_D, chain_E]
        x = 0
        for item in chain_numbering:
            if x == 0:
                if type == 1:
                    txt = pdbid + '[B-A]'
                elif type == 2:
                    txt = pdbid + '[A-E]'
                elif type == 3:
                    txt = pdbid + '[B-C]'
                elif type == 4:
                    txt = pdbid + '[A-B]'
                df_A.append(txt)
                df_A_nr.append(txt)
            elif x == 1:
                if type == 1:
                    txt = pdbid + '[A-E]'
                elif type == 2:
                    txt = pdbid + '[E-D]'
                elif type == 3:
                    txt = pdbid + '[C-D]'
                elif type == 4:
                    txt = pdbid + '[B-C]'
                df_B.append(txt)
                df_B_nr.append(txt)
            elif x == 2:
                if type == 1:
                    txt = pdbid + '[E-D]'
                elif type == 2:
                    txt = pdbid + '[D-C]'
                elif type == 3:
                    txt = pdbid + '[D-E]'
                elif type == 4:
                    txt = pdbid + '[C-D]'
                df_C.append(txt)
                df_C_nr.append(txt)
            elif x == 3:
                if type == 1:
                    txt = pdbid + '[D-C]'
                elif type == 2:
                    txt = pdbid + '[C-B]'
                elif type == 3:
                    txt = pdbid + '[E-A]'
                elif type == 4:
                    txt = pdbid + '[D-E]'
                df_D.append(txt)
                df_D_nr.append(txt)
            elif x == 4:
                if type == 1:
                    txt = pdbid + '[C-B]'
                elif type == 2:
                    txt = pdbid + '[B-A]'
                elif type == 3:
                    txt = pdbid + '[A-B]'
                elif type == 4:
                    txt = pdbid + '[E-A]'
                df_E.append(txt)
                df_E_nr.append(txt)
            # PRINCIPAL = 12
            # COMPLEMENTARY = 10
            if x < 4:
                # FIRST PRINCIPAL
                residue1 = chains[x][int(chain_numbering[x][7])]
                residue2 = chains[x][int(chain_numbering[x][8])]
                residue3 = chains[x][int(chain_numbering[x][11])]
                residue4 = chains[x][int(chain_numbering[x][12])]
                residue5 = chains[x][int(chain_numbering[x][14])]
                residue6 = chains[x][int(chain_numbering[x][15])]
                residue7 = chains[x][int(chain_numbering[x][16])]
                residue8 = chains[x][int(chain_numbering[x][17])]
                residue9 = chains[x][int(chain_numbering[x][18])]
                residue10 = chains[x][int(chain_numbering[x][19])]
                residue11 = chains[x][int(chain_numbering[x][20])]
                # SECOND COMPLEMENTARY
                residue12 = chains[x+1][int(chain_numbering[x+1][0])]
                residue13 = chains[x+1][int(chain_numbering[x+1][1])]
                residue14 = chains[x+1][int(chain_numbering[x+1][2])]
                residue15 = chains[x+1][int(chain_numbering[x+1][3])]
                residue16 = chains[x+1][int(chain_numbering[x+1][4])]
                residue17 = chains[x+1][int(chain_numbering[x+1][5])]
                residue18 = chains[x+1][int(chain_numbering[x+1][6])]
                residue19 = chains[x+1][int(chain_numbering[x+1][9])]
                residue20 = chains[x+1][int(chain_numbering[x+1][10])]
                residue21 = chains[x+1][int(chain_numbering[x+1][13])]
                aa_list = [residue1, residue2, residue3, residue4, residue5, residue6, residue7, residue8,
                           residue9, residue10, residue11, residue12, residue13, residue14, residue15, residue16,
                           residue17, residue18, residue19, residue20, residue21]
            else:
                # FIRST PRINCIPAL
                residue1 = chains[x][int(chain_numbering[x][7])]
                residue2 = chains[x][int(chain_numbering[x][8])]
                residue3 = chains[x][int(chain_numbering[x][11])]
                residue4 = chains[x][int(chain_numbering[x][12])]
                residue5 = chains[x][int(chain_numbering[x][14])]
                residue6 = chains[x][int(chain_numbering[x][15])]
                residue7 = chains[x][int(chain_numbering[x][16])]
                residue8 = chains[x][int(chain_numbering[x][17])]
                residue9 = chains[x][int(chain_numbering[x][18])]
                residue10 = chains[x][int(chain_numbering[x][19])]
                residue11 = chains[x][int(chain_numbering[x][20])]
                # SECOND COMPLEMENTARY
                residue12 = chains[0][int(chain_numbering[0][0])]
                residue13 = chains[0][int(chain_numbering[0][1])]
                residue14 = chains[0][int(chain_numbering[0][2])]
                residue15 = chains[0][int(chain_numbering[0][3])]
                residue16 = chains[0][int(chain_numbering[0][4])]
                residue17 = chains[0][int(chain_numbering[0][5])]
                residue18 = chains[0][int(chain_numbering[0][6])]
                residue19 = chains[0][int(chain_numbering[0][9])]
                residue20 = chains[0][int(chain_numbering[0][10])]
                residue21 = chains[0][int(chain_numbering[0][13])]
                aa_list = [residue1, residue2, residue3, residue4, residue5, residue6, residue7, residue8,
                           residue9, residue10, residue11, residue12, residue13, residue14, residue15, residue16,
                           residue17, residue18, residue19, residue20, residue21]
            i = 0
            j = 0
            nr = 0
            while i < len(aa_list):
                j = i + 1
                while j < len(aa_list):
                    d = distance(aa_list[i], aa_list[j])
                    #print([aa_list[0].resname, chains[x], chain_numbering[x][7]])
                    # ANPASSEN AN DIE RICHTIGE RESIDUES
                    if j == (len(aa_list) - 1):
                        v = angles(aa_list[i], aa_list[j], aa_list[0])
                    else:
                        v = angles(aa_list[i], aa_list[j], aa_list[j + 1])
                    #print([i, j])
                    if x == 0:
                        df_A.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_A_nr.append(t)
                    elif x == 1:
                        df_B.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_B_nr.append(t)
                    elif x == 2:
                        df_C.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x + 1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x + 1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_C_nr.append(t)
                    elif x == 3:
                        df_D.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_D_nr.append(t)
                    elif x == 4:
                        df_E.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[0]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[0]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[0]) + str(aa_list[i]) + ' - ' + str(chains[0]) + str(aa_list[j])
                        df_E_nr.append(t)
                    #print(t)
                    j += 1
                i += 1
            x += 1
    return df_A, df_B, df_C, df_D, df_E, df_A_nr, df_B_nr, df_C_nr, df_D_nr, df_E_nr


def steroid_interface(pdbid, chain_Anr, chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr, PDBids):
    df_A = []
    df_B = []
    df_C = []
    df_D = []
    df_E = []
    df_A_nr = []
    df_B_nr = []
    df_C_nr = []
    df_D_nr = []
    df_E_nr = []
    #print([len(chain_Anr), len(chain_Bnr), len(chain_Cnr), len(chain_Dnr), len(chain_Enr)])
    if len(chain_Anr) < 12 or len(chain_Bnr) < 12 or len(chain_Cnr) < 12 or len(chain_Dnr) < 12 or len(chain_Enr) < 12:
        PDBids = [None]
    if pdbid.upper() in PDBids:
        parser = PDBParser()
        #filename = '/home/ftk/PycharmProjects/Con_pLGIC/data/Structures/' + pdbid + '.pdb'
        filename = 'C:\\Users\\TARIQOPLATA\\PycharmProjects\\Con_pLGIC\\data\\Structures\\' + pdbid + '.pdb'
        #filename = 'C:\\Users\\filip\\PycharmProjects\\Con_pLGIC\\data\\Structures\\' + pdbid + ".pdb"
        structure = parser.get_structure(pdbid, filename)
        first_model = structure[0]
        if pdbid in clockwise:
            if '6HU' in pdbid or '6I53' in pdbid or '6A96' in pdbid or '4COF' in pdbid:
                chain_A = first_model["B"]
                chain_B = first_model["A"]
                chain_C = first_model["E"]
                chain_D = first_model["D"]
                chain_E = first_model["C"]
                chain_numbering = [chain_Bnr, chain_Anr, chain_Enr, chain_Dnr, chain_Cnr]
                type = 1
            else:
                chain_A = first_model["A"]
                chain_B = first_model["E"]
                chain_C = first_model["D"]
                chain_D = first_model["C"]
                chain_E = first_model["B"]
                chain_numbering = [chain_Anr, chain_Enr, chain_Dnr, chain_Cnr, chain_Bnr]
                type = 2
        elif pdbid in counter:
            if '6DW0' in pdbid:
                # B - C - D - E - A
                chain_A = first_model["B"]
                chain_B = first_model["C"]
                chain_C = first_model["D"]
                chain_D = first_model["E"]
                chain_E = first_model["A"]
                chain_numbering = [chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr, chain_Anr]
                type = 3
            else:
                chain_A = first_model["A"]
                chain_B = first_model["B"]
                chain_C = first_model["C"]
                chain_D = first_model["D"]
                chain_E = first_model["E"]
                type = 4
                chain_numbering = [chain_Anr, chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr]
        else:
            print("Failed to find correct direction")
            print(pdbid)
        chains = [chain_A, chain_B, chain_C, chain_D, chain_E]
        x = 0
        for item in chain_numbering:
            if len(item) < 10:
                continue
            if x == 0:
                if type == 1:
                    txt = pdbid + '[B-A]'
                elif type == 2:
                    txt = pdbid + '[A-E]'
                elif type == 3:
                    txt = pdbid + '[B-C]'
                elif type == 4:
                    txt = pdbid + '[A-B]'
                df_A.append(txt)
                df_A_nr.append(txt)
            elif x == 1:
                if type == 1:
                    txt = pdbid + '[A-E]'
                elif type == 2:
                    txt = pdbid + '[E-D]'
                elif type == 3:
                    txt = pdbid + '[C-D]'
                elif type == 4:
                    txt = pdbid + '[B-C]'
                df_B.append(txt)
                df_B_nr.append(txt)
            elif x == 2:
                if type == 1:
                    txt = pdbid + '[E-D]'
                elif type == 2:
                    txt = pdbid + '[D-C]'
                elif type == 3:
                    txt = pdbid + '[D-E]'
                elif type == 4:
                    txt = pdbid + '[C-D]'
                df_C.append(txt)
                df_C_nr.append(txt)
            elif x == 3:
                if type == 1:
                    txt = pdbid + '[D-C]'
                elif type == 2:
                    txt = pdbid + '[C-B]'
                elif type == 3:
                    txt = pdbid + '[E-A]'
                elif type == 4:
                    txt = pdbid + '[D-E]'
                df_D.append(txt)
                df_D_nr.append(txt)
            elif x == 4:
                if type == 1:
                    txt = pdbid + '[C-B]'
                elif type == 2:
                    txt = pdbid + '[B-A]'
                elif type == 3:
                    txt = pdbid + '[A-B]'
                elif type == 4:
                    txt = pdbid + '[E-A]'
                df_E.append(txt)
                df_E_nr.append(txt)
            if x < 4:
                # FIRST PRINCIPAL
                print([pdbid, chains[x], chain_numbering[x]])
                residue1 = chains[x][int(chain_numbering[x][4])]
                residue2 = chains[x][int(chain_numbering[x][5])]
                residue3 = chains[x][int(chain_numbering[x][6])]
                residue4 = chains[x][int(chain_numbering[x][7])]
                residue5 = chains[x][int(chain_numbering[x][8])]
                residue6 = chains[x][int(chain_numbering[x][9])]
                # SECOND COMPLEMENTARY
                residue7 = chains[x + 1][int(chain_numbering[x + 1][0])]
                residue8 = chains[x + 1][int(chain_numbering[x + 1][1])]
                residue9 = chains[x + 1][int(chain_numbering[x + 1][2])]
                residue10 = chains[x + 1][int(chain_numbering[x + 1][3])]
                residue11 = chains[x + 1][int(chain_numbering[x + 1][10])]
                residue12 = chains[x + 1][int(chain_numbering[x + 1][11])]
                aa_list = [residue1, residue2, residue3, residue4, residue5, residue6, residue7, residue8,
                           residue9, residue10, residue11, residue12]
            else:
                # FIRST PRINCIPAL
                residue1 = chains[x][int(chain_numbering[x][4])]
                residue2 = chains[x][int(chain_numbering[x][5])]
                residue3 = chains[x][int(chain_numbering[x][6])]
                residue4 = chains[x][int(chain_numbering[x][7])]
                residue5 = chains[x][int(chain_numbering[x][8])]
                residue6 = chains[x][int(chain_numbering[x][9])]
                # SECOND COMPLEMENTARY
                residue7 = chains[0][int(chain_numbering[0][0])]
                residue8 = chains[0][int(chain_numbering[0][1])]
                residue9 = chains[0][int(chain_numbering[0][2])]
                residue10 = chains[0][int(chain_numbering[0][3])]
                residue11 = chains[0][int(chain_numbering[0][10])]
                residue12 = chains[0][int(chain_numbering[0][11])]
                aa_list = [residue1, residue2, residue3, residue4, residue5, residue6, residue7, residue8,
                           residue9, residue10, residue11, residue12]
            i = 0
            j = 0
            while i < len(aa_list):
                j = i + 1
                while j < len(aa_list):
                    d = distance(aa_list[i], aa_list[j])
                    if j == (len(aa_list) - 1):
                        v = angles(aa_list[i], aa_list[j], aa_list[0])
                    else:
                        v = angles(aa_list[i], aa_list[j], aa_list[j + 1])
                    if x == 0:
                        df_A.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_A_nr.append(t)
                    elif x == 1:
                        df_B.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_B_nr.append(t)
                    elif x == 2:
                        df_C.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x + 1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x + 1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_C_nr.append(t)
                    elif x == 3:
                        df_D.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x+1]) + str(aa_list[i]) + ' - ' + str(chains[x+1]) + str(aa_list[j])
                        df_D_nr.append(t)
                    elif x == 4:
                        df_E.append(d)
                        if i < 11:
                            if j < 11:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[x]) + str(aa_list[i]) + ' - ' + str(chains[0]) + str(aa_list[j])
                        else:
                            if j < 10:
                                t = str(chains[0]) + str(aa_list[i]) + ' - ' + str(chains[x]) + str(aa_list[j])
                            else:
                                t = str(chains[0]) + str(aa_list[i]) + ' - ' + str(chains[0]) + str(aa_list[j])
                        df_E_nr.append(t)
                    j += 1
                i += 1
            x += 1
    print('---')
    print(len(df_A_nr))
    print(len(df_A))
    print('---')
    return df_A, df_B, df_C, df_D, df_E, df_A_nr, df_B_nr, df_C_nr, df_D_nr, df_E_nr


def global_confo(pdbid, chain_Anr, chain_Bnr, chain_Cnr, chain_Dnr, chain_Enr, PDBids):
    df = []
    df_nr = []
    if pdbid.upper() in PDBids:
        parser = PDBParser()
        filename = 'data/Structures/' + pdbid + '.pdb'
        #filename = 'C:\\Users\\filip\\PycharmProjects\\plgic\\Structures\\' + pdbid + ".pdb"
        #filename = 'C:\\Users\\TARIQOPLATA\\PycharmProjects\\Con_pLGIC\\data\\Structures\\' + pdbid + '.pdb'
        structure = parser.get_structure(pdbid, filename)
        first_model = structure[0]
        chain_A = first_model["A"]
        chain_B = first_model["B"]
        chain_C = first_model["C"]
        chain_D = first_model["D"]
        chain_E = first_model["E"]
        chains = [chain_A, chain_B, chain_C, chain_D, chain_E]
        chain_numbering = chain_Anr + chain_Bnr + chain_Cnr + chain_Dnr + chain_Enr
        df.append(pdbid)
        df_nr.append(pdbid)
        aa_list = []
        x = 0
        i = 0
        j = 1
        while i < len(chain_numbering):
            j = i + 1
            while j < len(chain_numbering):

                if 'A' in chain_numbering[i]:
                    nr = chain_numbering[i].split('-')
                    nr = nr[1]
                    residue1 = chains[0][int(nr)]
                    t1 = str(chains[0]) + str(residue1)
                elif 'B' in chain_numbering[i]:
                    nr = chain_numbering[i].split('-')
                    nr = nr[1]
                    residue1 = chains[1][int(nr)]
                    t1 = str(chains[1]) + str(residue1)
                elif 'C' in chain_numbering[i]:
                    nr = chain_numbering[i].split('-')
                    nr = nr[1]
                    residue1 = chains[2][int(nr)]
                    t1 = str(chains[2]) + str(residue1)
                elif 'D' in chain_numbering[i]:
                    nr = chain_numbering[i].split('-')
                    nr = nr[1]
                    residue1 = chains[3][int(nr)]
                    t1 = str(chains[3]) + str(residue1)
                elif 'E' in chain_numbering[i]:
                    nr = chain_numbering[i].split('-')
                    nr = nr[1]
                    residue1 = chains[4][int(nr)]
                    t1 = str(chains[4]) + str(residue1)
                if 'A' in chain_numbering[j]:
                    nr = chain_numbering[j].split('-')
                    nr = nr[1]
                    residue2 = chains[0][int(nr)]
                    t2 = str(chains[0]) + str(residue2)
                elif 'B' in chain_numbering[j]:
                    nr = chain_numbering[j].split('-')
                    nr = nr[1]
                    residue2 = chains[1][int(nr)]
                    t2 = str(chains[1]) + str(residue2)
                elif 'C' in chain_numbering[j]:
                    nr = chain_numbering[j].split('-')
                    nr = nr[1]
                    residue2 = chains[2][int(nr)]
                    t2 = str(chains[2]) + str(residue2)
                elif 'D' in chain_numbering[j]:
                    nr = chain_numbering[j].split('-')
                    nr = nr[1]
                    residue2 = chains[3][int(nr)]
                    t2 = str(chains[3]) + str(residue2)
                elif 'E' in chain_numbering[j]:
                    nr = chain_numbering[j].split('-')
                    nr = nr[1]
                    residue2 = chains[4][int(nr)]
                    t2 = str(chains[4]) + str(residue2)
                d = distance(residue1, residue2)
                #print(['DISTANCE', d])
                df.append(d)
                t = t1 + ' - ' + t2
                #print(t)
                df_nr.append(t)
                j += 1
            i += 1
    return df, df_nr


def read_seq(data, title, site, PDBids):
    i = 0
    sequ_number = []
    all_sequ_nrs = []
    for line in data:
        if '@<TRIPOS>MOLECULE' in line:
            if len(sequ_number) > 0:
                all_sequ_nrs.append(sequ_number)
            sequ_number = []
            i = 2
            continue
        if i == 2:
            header = line
            if '\n' in header:
                header = header[:-1]
            # print(header)
            sequ_number.append(header)
            i = 0
        if '@<TRIPOS>SUBSTRUCTURE' in line:
            i = 4
            continue
        if '# MOE 2020.09 (io_trps.svl 2020.11)' in line:
            i = 0
            continue
        if i == 4:
            if line == '\n':
                i = 0
                continue
            items = line.split(' ')
            items = list(filter(None, items))
            number = items[1]
            chain = items[5]
            result = number + '.' + chain
            if result not in sequ_number:
                sequ_number.append(result)
    chainsA = []
    chainsB = []
    chainsC = []
    chainsD = []
    chainsE = []
    df = []
    df_full = []
    labels = []
    pdbs_labels = []
    data_dict = {}
    sub_data = {}
    data_dict_nr = {}
    sub_data = {}
    for item in all_sequ_nrs:
        i = 0
        #if item[0] in ht3r:
        #    continue
        #if item[0] in GlyR:
        #    continue
        #if item[0] in nachr:
        #    continue
        #print([item[0], len(item)])
        while i < len(item):
            # print(item[i])
            if i == 0:
                pdbid = item[0]
            else:
                items = item[i].split('.')
                # print(items)
                nr = items[0]
                nr = nr[3:]
                if items[1] == 'A':
                    if site == 99:
                        txt = 'A-' + str(nr)
                        chainsA.append(txt)
                    else:
                        chainsA.append(nr)
                elif items[1] == 'B':
                    if site == 99:
                        txt = 'B-' + str(nr)
                        chainsA.append(txt)
                    else:
                        chainsB.append(nr)
                elif items[1] == 'C':
                    if site == 99:
                        txt = 'C-' + str(nr)
                        chainsA.append(txt)
                    else:
                        chainsC.append(nr)
                elif items[1] == 'D':
                    if site == 99:
                        txt = 'D-' + str(nr)
                        chainsA.append(txt)
                    else:
                        chainsD.append(nr)
                else:
                    if site == 99:
                        txt = 'E-' + str(nr)
                        chainsA.append(txt)
                    else:
                        chainsE.append(nr)
            i += 1
        if site == 1:
        ################## PS SITE
            df_X = PS_site(pdbid, chainsA, chainsB, chainsC, chainsD, chainsE, PDBids)
            df_A = df_X[0]
            df_B = df_X[1]
            df_C = df_X[2]
            df_D = df_X[3]
            df_E = df_X[4]
            df_A_nr = df_X[5]
            df_B_nr = df_X[6]
            df_C_nr = df_X[7]
            df_D_nr = df_X[8]
            df_E_nr = df_X[9]
            df_X = [df_A, df_B, df_C, df_D, df_E]
            df_Z = [df_A_nr, df_B_nr, df_C_nr, df_D_nr, df_E_nr]
            # PER SUBUNIT
            for item in df_X:
                #print(['df_X', len(item)])
                if pdbid == '6D6U':
                    df_X = []
                if pdbid == '6D6T':
                    df_X = []
                if len(item) < 66:
                    continue
                labels.append(item[0])
                df.append(item[1:])
                xx = 0
                for it in item:
                    # if xx == 0:
                    #    continue
                    sub_data[xx] = it
                    xx = xx + 1
                n = item[0]
                data_dict[n] = sub_data
                sub_data = {}
            sub_data = {}
            for item in df_Z:
                print(['df_Z', len(item)])
                if pdbid == '6D6U':
                    df_Z = []
                elif pdbid == '6D6T':
                    df_Z = []
                elif len(item) < 66:
                    continue
                else:
                    xx = 0
                    for it in item:
                        print(['IT', it])
                        sub_data[xx] = it
                        xx = xx + 1
                    n = item[0]
                    print(sub_data)
                    data_dict_nr[n] = sub_data
                    print(data_dict_nr)
                    sub_data = {}
        ################# upper TMD SITE
        elif site == 2:
            df_X = upperTMD_site(pdbid, chainsA, chainsB, chainsC, chainsD, chainsE, PDBids)
            df_A = df_X[0]
            df_B = df_X[1]
            df_C = df_X[2]
            df_D = df_X[3]
            df_E = df_X[4]
            df_A_nr = df_X[5]
            df_B_nr = df_X[6]
            df_C_nr = df_X[7]
            df_D_nr = df_X[8]
            df_E_nr = df_X[9]
            df_X = [df_A, df_B, df_C, df_D, df_E]
            df_Z = [df_A_nr, df_B_nr, df_C_nr, df_D_nr, df_E_nr]
            #print('Z: ----------------')
            #print(df_Z)
            for item in df_X:
                print([pdbid, len(item)])
                if pdbid == '6D6U':
                    df_X = []
                elif pdbid == '6D6T':
                    df_X = []
                elif len(item) < 211:
                    continue
                else:
                    labels.append(item[0])
                    df.append(item[1:])
                    xx = 0
                    for it in item:
                        #if xx == 0:
                        #    continue
                        sub_data[xx] = it
                        xx = xx + 1
                    n = item[0]
                    data_dict[n] = sub_data
                    sub_data = {}
            sub_data = {}
            for item in df_Z:
                if pdbid == '6D6U':
                    df_Z = []
                elif pdbid == '6D6T':
                    df_Z = []
                elif len(item) < 211:
                    continue
                else:
                    xx = 0
                    for it in item:
                        sub_data[xx] = it
                        xx = xx + 1
                    n = item[0]
                    data_dict_nr[n] = sub_data
                    sub_data = {}
        elif site == 3:
            ############## interface steroid site
            df_X = steroid_interface(pdbid, chainsA, chainsB, chainsC, chainsD, chainsE, PDBids)
            df_A = df_X[0]
            df_B = df_X[1]
            df_C = df_X[2]
            df_D = df_X[3]
            df_E = df_X[4]
            df_A_nr = df_X[5]
            df_B_nr = df_X[6]
            df_C_nr = df_X[7]
            df_D_nr = df_X[8]
            df_E_nr = df_X[9]
            df_X = [df_A, df_B, df_C, df_D, df_E]
            df_Z = [df_A_nr, df_B_nr, df_C_nr, df_D_nr, df_E_nr]
            # PER DIMER
            for item in df_X:
                if '6D6U' in item:
                    df_X = []
                elif '6D6T' in item:
                    df_X = []
                elif len(item) < 66:
                    continue
                else:
                    labels.append(item[0])
                    df.append(item[1:])
                    xx = 0
                    for it in item:
                        # if xx == 0:
                        #    continue
                        sub_data[xx] = it
                        xx = xx + 1
                    n = item[0]
                    data_dict[n] = sub_data
                    sub_data = {}
            sub_data = {}
            for item in df_Z:
                if pdbid == '6D6U':
                    df_Z = []
                elif pdbid == '6D6T':
                    df_Z = []
                elif len(item) < 66:
                    continue
                else:
                    xx = 0
                    for it in item:
                        sub_data[xx] = it
                        xx = xx + 1
                    n = item[0]
                    data_dict_nr[n] = sub_data
                    sub_data = {}
        elif site == 99:
            df_X = global_confo(pdbid, chainsA, chainsB, chainsC, chainsD, chainsE, PDBids)
            df_Z = df_X[1]
            df_X = df_X[0]
            print(pdbid)
            print(len(df_Z))
            print(len(df_X))
            if len(df_X) == 0:
                continue
            labels.append(df_X[0])
            df.append(df_X[1:])
            xx = 0
            for it in df_X:
                sub_data[xx] = it
                xx = xx + 1
            n = item[0]
            data_dict[n] = sub_data
            sub_data = {}
            xx = 0
            for it in df_Z:
                sub_data[xx] = it
                xx = xx + 1
            n = item[0]
            data_dict_nr[n] = sub_data
            sub_data = {}
        chainsA = []
        chainsB = []
        chainsC = []
        chainsD = []
        chainsE = []
    data = np.array([np.array(xi) for xi in df])
    labels = np.array(labels)
    # For nested
    df_ne2 = pd.DataFrame(data_dict).T
    df_nr2 = pd.DataFrame(data_dict_nr).T
    print(df_nr2)
    print(df_ne2)
    print('-----------------------------------------------------------------')
    #print('VERSION 2')
    #print(df_ne2)
    #df_ne3 = pd.DataFrame.from_dict(data_dict, orient='index')
    #print('VERSION 3')
    #print(df_ne3)
    #fig2 = pca_2d(data, labels, title, site)
    fig = pca(data, labels, title, site)
    df_ne2 = df_ne2.drop(df_ne2.columns[0], axis=1)
    df_ne2['pdbid'] = df_ne2.index
    df_ne2.index = range(len(df_ne2))
    pdbids = df_ne2['pdbid']
    df_ne2.drop(labels=['pdbid'], axis=1, inplace=True)
    df_ne2.insert(0, 'pdbid', pdbids)
    df_nr2 = df_nr2.drop(df_nr2.columns[0], axis=1)
    df_nr2['pdbid'] = df_nr2.index
    df_nr2.index = range(len(df_nr2))
    pdbids = df_nr2['pdbid']
    df_nr2.drop(labels=['pdbid'], axis=1, inplace=True)
    df_nr2.insert(0, 'pdbid', pdbids)
    print('DF !!----')
    print(df_nr2)
    print(df_ne2)
    df, df_nrs = MinMax(df_ne2, df_nr2)
    #df_ne2.to_csv('Distances_v2.csv')
    #df_nr2.to_csv('Names_v2.csv')
    np.savetxt('Data_bs3_apos_1.txt', data)
    np.savetxt('Labels_bs3_apos_1.txt', labels, delimiter=" ", fmt="%s")
    return fig, df, df_nrs


def pca_2d(data, labels, title, site):
    pca = PCA(n_components=2)
    X = pca.fit_transform(data)
    df2 = pd.DataFrame(X)
    df2['labels'] = labels.tolist()
    symbols = []
    for item in labels.tolist():
        #t = item[:4]
        if item in ligand_bound_bs3:
            symbols.append('diamond')
        elif item in ligand_notbound_bs3:
            symbols.append('circle-open')
        else:
            symbols.append('circle')
    df2['markers'] = symbols
    df2.columns = ['PC1', 'PC2', 'labels', 'markers']
    df2 = add_families(df2, site)
    #txt = title + '_PCA_Result_TMD.csv'
    df2.to_csv('BS3_all_data_new.csv', sep=';')
    Scene = dict(xaxis=dict(title='PC1'), yaxis=dict(title='PC2'))
    trace = go.Scatter(x=df2['PC1'], y=df2['PC2'], mode='markers+text', hovertext=df2['labels'],
                         marker_symbol=df2['markers'], #symbol=df2['markers'],
                         marker=dict(color=df2['Colors'], size=10, line=dict(color='black', width=10)))
    data = [trace]
    layout = go.Layout(margin=dict(l=0, r=0), scene=Scene, height=800, width=1200)
    fig = go.Figure(data=data, layout=layout)
    #fig = px.scatter_3d(df2, x='PC1', y='PC2', z='PC3', title=title, hover_name='labels', color='Colors')
    #txt = 'C:\\Users\\filip\\Desktop\\' + title + '.html'
    fig.write_html('upperTMD_newColors_newSymbols_allFamilies_2d.html')
    fig.show()
    return fig


def pca(data, labels, title, site):
    pca = PCA(n_components=3)
    X = pca.fit_transform(data)
    df2 = pd.DataFrame(X)
    df2['labels'] = labels.tolist()
    symbols = []
    for item in labels.tolist():
        t = item[:4]
        if t in ligand_bound_bs4:
            symbols.append('diamond')
        #elif item in ligand_notbound_bs3:
        #    symbols.append('diamond-cross')
        else:
            symbols.append('circle')
    df2['markers'] = symbols
    df2.columns = ['PC1', 'PC2', 'PC3', 'labels', 'markers']
    df2 = add_families(df2, site)
    #txt = title + '_PCA_Result_TMD.csv'
    df2.to_csv('BS3_PCA_data_new.csv', sep=';')
    Scene = dict(xaxis=dict(title='PC1', showgrid=True, gridwidth=1, gridcolor='black'),
                 yaxis=dict(title='PC2', showgrid=True, gridwidth=1, gridcolor='black'),
                 zaxis=dict(title='PC3', showgrid=True, gridwidth=1, gridcolor='black'))
    trace = go.Scatter3d(x=df2['PC1'], y=df2['PC2'], z=df2['PC3'], mode='markers+text', hovertext=df2['labels'],
                         marker_symbol=df2['markers'], #symbol=df2['markers'],
                         marker=dict(color=df2['Colors'], size=10, line=dict(color='black', width=10)))
    data = [trace]
    layout = go.Layout(margin=dict(l=0, r=0), scene=Scene, height=800, width=1200)
    fig = go.Figure(data=data, layout=layout)
    #fig = px.scatter_3d(df2, x='PC1', y='PC2', z='PC3', title=title, hover_name='labels', color='Colors')
    #txt = 'C:\\Users\\filip\\Desktop\\' + title + '.html'
    fig.write_html('BS3_apos_newSymbols_nachrsMarked.html')
    #fig.show()
    return fig


def color_interfaces(df, site):
    families = []
    colors = []
    a7 = ['']
    for index, row in df.iterrows():
        if site == 1:
            txt = row['labels'].split('_')
            txt = txt[0]
        elif site == 2:
            txt = row['labels']
            txt = txt[:-5]
        elif site == 3:
            txt = row['labels']
            txt = txt[:-5]
        if row['labels'] in a7:
            colors.append('rgb(128,0,0)')


def add_families(df, site):
    families = []
    colors = []
    for index, row in df.iterrows():
        if site == 1:
            txt = row['labels'].split('_')
            txt = txt[0]
        elif site == 2:
            txt = row['labels']
            txt = txt[:-5]
        elif site == 3:
            txt = row['labels']
            txt = txt[:-5]
        if txt in GABAAR:
            if txt in chimeras:
                families.append('GABAAR_chimera')
                colors.append('rgb(106, 128, 156)')
            else:
                families.append('GABAAR')
                colors.append('rgb(43, 113, 255)')
                """
                if row['labels'] in agonist_bound:
                    families.append('GABAAR_agonist')
                    colors.append('rgb(2, 49, 160)')
                elif txt in agonist_bound:
                    families.append('GABAAR_agonist')
                    colors.append('rgb(2, 49, 160)')
                elif row['labels'] in antagonists_bound:
                    families.append('GABAAR_antagonist')
                    colors.append('rgb(127, 170, 255)')
                elif txt in antagonists_bound:
                    families.append('GABAAR_antagonist')
                    colors.append('rgb(127, 170, 255')
                else:
                    families.append('GABAAR')
                    colors.append('rgb(43, 113, 255)')"""
        elif txt in GlyR:
            if txt in chimeras:
                families.append('GlyR_chimera')
                colors.append('rgb(119, 150, 142)')
            else:
                families.append('GlyR')
                colors.append('rgb(0, 204, 150)')
                """
                if row['labels'] in agonist_bound:
                    families.append('GlyR_agonist')
                    colors.append('rgb(0, 134, 150)')
                elif txt in agonist_bound:
                    families.append('GlyR_agonist')
                    colors.append('rgb(0, 134, 150)')
                elif row['labels'] in antagonists_bound:
                    families.append('GlyR_antagonist')
                    colors.append('rgb(95, 255, 255)')
                elif txt in antagonists_bound:
                    families.append('GlyR_antagonist')
                    colors.append('rgb(95, 255, 255)')
                else:
                    families.append('GlyR')
                    colors.append('rgb(0, 204, 150)')"""
        elif txt in nachr:
            if txt in nachrs_muscle:
                families.append('nAChR_muscle-type')
                colors.append('rgb(254, 132, 132)')
            elif txt in nachrs_homos:
                families.append('nAChR_alpha7-homopentamers')
                colors.append('rgb(83, 0, 0)')
            elif txt in nachrs_heteros:
                families.append('nAChR_heteropentamers')
                colors.append('rgb(254, 32, 32)')
            else:
                print('NOT FOUNDS BOYS')
                print(txt)
                families.append('nAChR')
                colors.append('rgb(239, 85, 59)')
            """
            # AGONIST / ANTAGONISTS
            if row['labels'] in agonist_bound:
                families.append('nAChR_agonist')
                colors.append('rgb(151, 21, 21)')
            elif txt in agonist_bound:
                families.append('nAChR_agonist')
                colors.append('rgb(151, 21, 21)')
            elif row['labels'] in antagonists_bound:
                families.append('nAChR_antagonist')
                colors.append('rgb(255, 172, 194)')
            elif txt in antagonists_bound:
                families.append('nAChR_antagonist')
                colors.append('rgb(255, 172, 194)')
            else:
                families.append('nAChR')
                colors.append('rgb(239, 85, 59)')"""
        elif txt in ht3r:
            families.append('5HT3R')
            colors.append('rgb(239, 255, 0)')
            """
            if txt in agonist_bound:
                families.append('5HT3R_agonist')
                colors.append('rgb(143, 154, 0)')
            elif txt in antagonists_bound:
                families.append('5HT3R_antagonist')
                colors.append('rgb(249, 255, 148)')
            else:
                families.append('5HT3R')
                colors.append('rgb(239, 255, 0)')"""
        else:
            print(['CANT FIND !!!!!!', txt])
    df['Family'] = families
    df['Colors'] = colors
    return df


def main():
    #data = open('/home/ftk/PycharmProjects/Con_pLGIC_new/data/upperTMD_newest_aas.mol2', 'r')
    data = open('C:\\Users\\TARIQOPLATA\\PycharmProjects\\Con_pLGIC\\data\\upperTMD_newest_aas.mol2', 'r')
    #data = open('/home/ftk/Conformations/Con_pLGIC/data/upperTMD.mol2', 'r')
    data = data.readlines()
    title = ''
    fig, df, df_nrs = read_seq(data, title, 2, PDBids)
    fig.write_html('BS3_apos_newSymbols_nachrsMarked.html')
    fig.show()
    #MinMax(df, df_nrs)


def MinMax(df, df_nrs):
    pdbids = df['pdbid']
    df_a = pd.Series(['Max', 'Max_Names'])
    pdbids = pdbids.append(df_a, ignore_index=True)
    df.set_index('pdbid', inplace=True)
    df = df.apply(pd.to_numeric, errors='ignore')
    maxValues = df.max(axis=0)
    maxValueIndex = df.idxmax(axis=0)
    minValues = df.min(axis=0)
    minValueIndex = df.idxmin(axis=0)
    df.index = range(len(df))
    df['pdbid'] = pdbids
    maxs = ['Max']
    mins = ['Min']
    column_names = df.columns
    print(df_nrs)
    #print(column_names)
    #print(len(column_names))
    for idx, item in enumerate(column_names):
        if 'pdb' in str(item):
            continue
        print(item)
        my_max_ind = df[item].idxmax()
        my_min_ind = df[item].idxmin()
        print(my_max_ind)
        max_pdbname = df_nrs.iloc[my_max_ind]['pdbid']
        min_pdbname = df_nrs.iloc[my_min_ind]['pdbid']
        max_name = df_nrs.iloc[my_max_ind][item]
        min_name = df_nrs.iloc[my_min_ind][item]
        max_txt = ' - ' + max_name + '- of PDB_ID[Chains] ' + max_pdbname
        min_txt = ' - ' + min_name + '- of PDB_ID[Chains] ' + min_pdbname
        maxs.append(max_txt)
        mins.append(min_txt)
    a_series = pd.Series(maxs, index=df_nrs.columns)
    df_nrs = df_nrs.append(a_series, ignore_index=True)
    a_series = pd.Series(mins, index=df_nrs.columns)
    df_nrs = df_nrs.append(a_series, ignore_index=True)
    a_series = pd.Series(maxValues, index=df.columns)
    df = df.append(a_series, ignore_index=True)
    a_series = pd.Series(maxValueIndex, index=df.columns)
    df = df.append(a_series, ignore_index=True)
    a_series = pd.Series(minValues, index=df.columns)
    df = df.append(a_series, ignore_index=True)
    a_series = pd.Series(minValueIndex, index=df.columns)
    df = df.append(a_series, ignore_index=True)
    df.at[495, 'pdbid'] = 'Max'
    df.at[496, 'pdbid'] = 'Max_Names'
    df.at[497, 'pdbid'] = 'Min'
    df.at[498, 'pdbid'] = 'Min_Names'
    df.to_csv('Apos_full_data_bs3.csv')
    df_nrs.to_csv('Apos_full_labels_bs3.csv')
    return df, df_nrs


if __name__ == '__main__':
    main()
