import sys
import os
import numpy as np
from ase.geometry import get_distances
from ase.db import connect
from ase.io import write
from ase.build import make_supercell
from ase.data import covalent_radii as cradii
from scipy.spatial import Voronoi


def get_ads_dist(atoms, ads0, ads1='H'):
    index0 = [i for i in range(len(atoms)) if atoms[i].symbol == ads0]
    index1 = [i for i in range(len(atoms)) if atoms[i].symbol == ads1]
    dist = []
    D, D_len = get_distances(atoms.positions[index0],
                             atoms.positions[index1],
                             atoms.cell, pbc=True)

    return np.max(D_len)


def get_ads_slab_dist(atoms, ads0):
    index0 = [i for i in range(len(atoms)) if atoms[i].symbol == ads0]
    index1 = range(12)
    dist = []
    D, D_len = get_distances(atoms.positions[index0],
                             atoms.positions[index1],
                             atoms.cell, pbc=True)

    return np.min(D_len)


def check_adsorbate(A, B):
    dissociated = False
    if not len(B) > 13:  # only one adsorbate
        return dissociated, A, B

    adsatoms = [atom for atom in B[12:]]
    ads0, ads1 = set(atom.symbol for atom in adsatoms)
    dist_A = get_ads_dist(A, ads0, ads1)
    dist_B = get_ads_dist(B, ads0, ads1)

    if dist_B > 1.2 * dist_A:  # dissociation
        print('DISSOCIATED')
        dissociated = True

    adsatom = adsatoms[np.argmin([atom.position[2]
                                  for atom in adsatoms])].symbol
    removeads = [atom.symbol for atom in adsatoms
                 if not atom.symbol == adsatom]
    for ra in removeads:
        del A[[atom.index for atom in A if atom.symbol == ra]]
        del B[[atom.index for atom in B if atom.symbol == ra]]

    return dissociated, A, B


def is_desorbed(B):
    desorbed = False
    D, D_len = get_distances(B.positions, cell=B.cell, pbc=True)

    indexM = np.argmin(D_len[12, :12])
    dist_S = D_len[12, indexM]

    if dist_S > (cradii[B[12].number] + cradii[B[indexM].number]) * 2:
        print('DESORBED FROM SLAB')
        print(cradii[B[12].number] + cradii[B[indexM].number], dist_S)
        desorbed = True
    return desorbed


def is_reconstructed(B):
    a1 = B.get_angle(0, 7, 8, mic=True)
    a2 = B.get_angle(1, 6, 9, mic=True)
    cond1 = 170 < a1 < 190
    cond2 = 170 < a1 < 190
    a3 = B.get_angle(0, 1, 2)
    a4 = B.get_angle(8, 9, 10)
    cond3 = a3 - 10 < a4 < a3 + 10

    if cond1 and cond2 and cond3:
        return False
    else:
        return True


def compare_slab(A, B):
    hollow_dist = np.linalg.norm(B.positions[8, :2] - B.positions[0, :2])
    a = A.positions[-5:-1]
    b = B.positions[-5:-1]

    allowed_z_movement = 0.4 * cradii[A.get_atomic_numbers()[-5:-1]]

    d_xy = np.empty(4)
    d_z = np.empty(4)

    D, D_len = get_distances(p1=a, p2=b, cell=A.cell, pbc=True)
    d_xy = np.linalg.norm(np.diagonal(D)[:2], axis=0)

    d_z = np.diagonal(D)[2:][0]

    cond1 = np.all(d_xy < 0.15 * hollow_dist)
    cond2 = np.all([d_z[i] < allowed_z_movement[i]
                    for i in range(len(a))])

    if cond1 and cond2:
        return True

    else:
        return False


def is_subsurface(A):

    pos0 = A.positions[-5:-1][:, 2]
    pos1 = A.positions[-1][2]
    metal_covalent_radii = cradii[A.get_atomic_numbers()[-5:-1]]
    ads_covalent_radii = cradii[A.get_atomic_numbers()[-1]]

    if np.any([(pos0[i] - pos1) > 0.5 * metal_covalent_radii[i]
               for i in range(len(pos0))]):
        return True
    else:
        return False


def get_hollow(A):

    from ase.build import make_supercell
    P = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]
    SC = make_supercell(A, P)

    second_layer_Z = A.positions[4][2]
    ads_pos = SC.positions[25][:2]

    second_layer = []
    for pos in SC.positions:
        if abs(pos[2] - second_layer_Z) < 0.05:
            second_layer += [pos[:2]]

    if np.any(np.linalg.norm(second_layer - ads_pos, axis=1) <
              0.25 * cradii[A.get_atomic_numbers()[0]]):
        return 'HCP'
    else:
        return 'FCC'


def get_site_dict(B, ads_pos0):
    C = B.copy()
    C = C[-5:-1]

    SC = C * (3, 3, 1)

    D = {}
    cell = C.get_cell()
    t_dist = 1

    """Top sites: just all atoms """
    t_ind = 0
    for atom in [atom for atom in SC if
                 np.linalg.norm(atom.position[:2] - ads_pos0) < t_dist]:
        k = 'top_site-' + str(t_ind)
        D[k] = {}
        D[k]['pos'] = atom.position
        D[k]['sym'] = atom.symbol
        t_ind += 1

    b_dist = 0.5 * cell[0][0]
    b_ind = 1

    """Bridge sites: bewteen pairs """
    for atomi in [atom for atom in SC if
                  np.linalg.norm(atom.position[:2] - ads_pos0) < b_dist]:
        pos = atomi.position
        for atom in [atom for atom in SC if
                     np.linalg.norm(atom.position - pos) < 1.5 * b_dist and
                     np.linalg.norm(atom.position[:2] - ads_pos0) < b_dist
                     and not (atom.position == pos).all()]:
            ads_pos = atom.position
            temp = 0.5 * (pos + ads_pos)
            k = 'bridge_site-' + str(b_ind)
            D[k] = {}
            D[k]['pos'] = temp
            D[k]['sym'] = atomi.symbol + '_' + atom.symbol
            b_ind += 1

    h_dist = np.linalg.norm(B.positions[11, :2] - B.positions[4, :2])
    h_ind = 1
    h4_ind = 1

    """Hollow sites: Voronoi vertices """
    vor = Voronoi(SC.positions[:, :2])
    vertices = vor.vertices
    close_v = [v for v in vertices if np.linalg.norm(v - ads_pos0[:2]) < 1]
    for v in close_v:
        close_v_v = [v0 for v0 in close_v if np.linalg.norm(v0 - v) < 0.5]
        if len(close_v_v) > 1:
            v_mean = np.mean(close_v_v, axis=0)
            k = '4fold_' + str(h4_ind)
            D[k] = {}
            D[k]['pos'] = v_mean

            # Delete bridge sites overlapping with 4fold
            for key in [key for key in list(D.keys()) if 'bridge' in key]:
                bridge_pos = D[key]['pos']
                if np.linalg.norm(bridge_pos[:2] - v) < 0.3:
                    del D[key]

            min_dist = sorted([np.linalg.norm(v_mean - m[:2])
                               for m in SC.positions])[4]
            symb = [atom.symbol for atom in
                    SC if np.linalg.norm(v_mean - atom.position[:2]) < min_dist]
            D[k]['sym'] = '_'.join(symb)
            h4_ind += 1
        else:
            k = 'hollow_site-' + str(h_ind)
            D[k] = {}
            D[k]['pos'] = v
            min_dist = sorted([np.linalg.norm(v - m[:2])
                               for m in SC.positions])[3]
            symb = [atom.symbol for atom in
                    SC if np.linalg.norm(v - atom.position[:2]) < min_dist]
            D[k]['sym'] = '_'.join(symb)
            h_ind += 1

    return D


def get_under_bridge(B):
    SC0 = B * (3, 3, 1)
    ads_pos = SC0.positions[64]

    B = B[np.r_[4:8]]
    SC = B * (3, 3, 1)
    cell = B.get_cell()

    ret = None
    dis = cell[0][0] * 2

    for ele in SC:
        new_dis = np.linalg.norm(ads_pos - ele.position)
        if new_dis < dis:
            dis = new_dis
            ret = ele.symbol

    return ret


def get_under_hollow(B):
    SC0 = B * (3, 3, 1)
    ads_pos = SC0.positions[64]

    B = B[np.r_[4:8]]
    SC = B * (3, 3, 1)
    cell = B.get_cell()

    ret = 'FCC'

    if np.any([np.linalg.norm(ads_pos[:2] - ele.position[:2]) < 0.5 *
               cradii[ele.number] for ele in SC]):
        ret = 'HCP'

    return ret


def get_site(B):

    C = B.copy()
    C = C * (3, 3, 1)

    ads_pos = C.positions[64]

    Dict = get_site_dict(B, ads_pos[:2])

    # View sampled sites
    # from ase.visualize import view
    # from ase.build import add_adsorbate
    # X = B.copy()
    # X = X * (3, 3, 1)
    # del X[12]
    # for pos in Dict.values():
    #     add_adsorbate(X, 'X', position=(pos['pos'][:2]), height=0.2)
    # view(X)

    f_a_s = None
    dis = B.get_cell()[0][0]
    Kind = None

    values = [np.linalg.norm(ads_pos[:2] - d['pos'][:2])
              for d in list(Dict.values())]
    if len(values) == 0:
        return 'N/A', ''
    idx = np.argmin(values)
    dis = values[idx]
    k = list(Dict.keys())[idx]
    f_a_s = k.split('_')[0]
    Kind = k

    if f_a_s == 'top':
        s_t = Dict[Kind]['sym']

    if f_a_s == 'bridge':
        s_t = Dict[Kind]['sym'] + '|' + get_under_bridge(B)

    elif f_a_s == 'hollow':
        s_t = Dict[Kind]['sym'] + '|' + get_under_hollow(B)

    elif f_a_s == '4fold':
        s_t = Dict[Kind]['sym']

    if dis > 0.5:
        f_a_s += '-tilt'
        print('Warnign: A strong site match could not be found!')
        print('  structure labeled as {}'.format(f_a_s))

    return f_a_s, s_t


def get_info(A, B):

    dissociated, A, B = check_adsorbate(A, B)

    reconstructed = not compare_slab(A, B)

    if dissociated:
        return reconstructed, 'dissociated', ''

    if is_desorbed(B):
        return reconstructed, 'desorbed', ''

    if is_subsurface(B):
        return reconstructed, 'subsurface', ''

    site, site_type = get_site(B)

    return reconstructed, site, site_type
