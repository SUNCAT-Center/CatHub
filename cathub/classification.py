import sys
import os
import numpy as np
import ase
from ase.geometry import get_distances
from ase.db import connect
from ase.io import write
from ase.visualize import view
from ase.build import add_adsorbate

from ase.data import covalent_radii as cradii
from scipy.spatial import Voronoi


class SiteClassification:
    """ Determine surface reconstruction (True/False) and adsorption
    site for chemisorption.

    A: ASE Atoms object
        Initial surface structure. Adsorbate atoms must be last indices
    B: ASE Atoms object
        Final surface structure. Adsorbate atoms must be last indices
    natoms_top_layer: int
        Number of atoms in top layer of slab
    natoms_slab: int
        Number of atoms in slab
    """

    def __init__(self, B, natoms_top_layer=4, natoms_slab=12, A=None):
        self.ntop = natoms_top_layer
        self.nslab = natoms_slab
        self.B = B

        # Check if multiatomic adsorbates dissociates on surface
        self.dissociated = self.check_dissociated()

        # Only keep the adsorbate closest to surface
        B = self.remove_extra_atoms(B)

        # Sort layers according to z-position
        layer_indices = np.argsort(B[:-1].positions[:, 2])
        self.B = B[:-1][layer_indices] + B[-1]
        if A:
            A = self.remove_extra_atoms(A)
            self.A = A
            self.A = A[:-1][layer_indices] + A[-1]

    def get_info(self):
        """Return surface reconstruction as well as primary and
        secondary adsorption site labels"""

        reconstructed = self.is_reconstructed()

        site, site_type = self.get_site()

        return reconstructed, site, site_type

    def check_dissociated(self, cutoff=1.2):
        """Check if adsorbate dissociates"""
        dissociated = False
        if not len(self.B) > self.nslab + 1:  # only one adsorbate
            return dissociated

        adsatoms = [atom for atom in self.B[self.nslab:]]
        ads0, ads1 = set(atom.symbol for atom in adsatoms)
        bond_dist = get_ads_dist(self.B, ads0, ads1)

        Cradii = [cradii[atom.number]
                  for atom in [ase.Atom(ads0), ase.Atom(ads1)]]
        bond_dist0 = sum(Cradii)

        if bond_dist > cutoff * bond_dist0:
            print('DISSOCIATED: {} Ang > 1.2 * {} Ang'
                  .format(bond_dist, bond_dist0))
            dissociated = True

            return dissociated

    def remove_extra_atoms(self, slab):
        if not len(slab) > self.nslab + 1:
            return slab

        adsatoms = [atom for atom in slab[self.nslab:]]
        adsindex = np.argmin([atom.position[2]
                              for atom in adsatoms])

        del slab[[atom.index + self.nslab for atom in adsatoms
                  if not atom.index == adsindex]]

        return slab

    def is_desorbed(self):
        desorbed = False
        D, D_len = get_distances(self.B.positions, cell=self.B.cell, pbc=True)

        indexM = np.argmin(D_len[-1, :-1])
        dist_S = D_len[-1, indexM]

        if dist_S > (cradii[self.B[-1].number] +
                     cradii[self.B[indexM].number]) * 2:
            print('DESORBED FROM SLAB')
            desorbed = True
        return desorbed

    def is_reconstructed(self, xy_cutoff=0.3, z_cutoff=0.4):
        """Compare initial and final slab configuration 
        to determine if slab reconstructs during relaxation

        xy_cutoff: Allowed xy-movement is determined from 
                   the covalent radii as:
                   xy_cutoff * np.mean(cradii)
        z_cutoff:  Allowed z-movement is determined as
                   z_cutoff * cradii_i
        """

        assert self.A, \
            'Initial slab geometry needed to classify reconstruction'

        # remove adsorbate
        A = self.A[:-1].copy()
        B = self.B[:-1].copy()

        # Order wrt x-positions
        x_indices = np.argsort(A.positions[:, 0])
        A = A[x_indices]
        B = B[x_indices]
        a = A.positions
        b = B.positions

        allowed_z_movement = z_cutoff * cradii[A.get_atomic_numbers()]
        allowed_xy_movement = \
            xy_cutoff * np.mean(cradii[A.get_atomic_numbers()])

        D, D_len = get_distances(p1=a, p2=b, cell=A.cell, pbc=True)
        d_xy = np.linalg.norm(np.diagonal(D)[:2], axis=0)
        d_z = np.diagonal(D)[2:][0]

        cond1 = np.all(d_xy < allowed_xy_movement)
        cond2 = np.all([d_z[i] < allowed_z_movement[i]
                        for i in range(len(a))])

        if cond1 and cond2:  # not reconstructed
            return False
        else:
            return True

    def is_subsurface(self):
        pos0 = self.B.positions[:-1][:, 2]
        pos1 = self.B.positions[-1][2]
        metal_covalent_radii = cradii[self.B.get_atomic_numbers()[:-1]]
        ads_covalent_radii = cradii[self.B.get_atomic_numbers()[-1]]

        if np.any([(pos0[i] - pos1) > 0.5 * metal_covalent_radii[i]
                   for i in range(len(pos0))]):
            return True
        else:
            return False

    def get_site_dict(self, ads_pos):
        """Get dictionary with high symmetry sites close to adsorbate
        position. 
        Top sites: Optained from the atomic positions of
                   the top layer of the slab. 
        Bridge sites: The average position of atomic pairs
        Hollow sites: Optained as the Voronoi vertices
        4-fold sites: Assigned when Voronoi vertives overlap
        """

        # Get top layer
        C = self.B[-self.ntop - 1:-1]
        SC = C * (3, 3, 1)

        D = {}
        cell = C.get_cell()
        top_dist = 1

        """Top sites: atomic positions """
        for i, atom in \
            enumerate([atom for atom in SC if
                       np.linalg.norm(atom.position[:2] - ads_pos)
                       < top_dist]):
            k = 'top_site-{}'.format(i)
            D[k] = {}
            D[k]['pos'] = atom.position
            D[k]['sym'] = atom.symbol

        """Bridge sites: bewteen atomic pairs """
        bridge_dist = 0.5 * cell[0][0]
        b_index = 0
        for atomi in [atom for atom in SC if
                      np.linalg.norm(atom.position[:2] - ads_pos)
                      < bridge_dist]:
            pos1 = atomi.position
            # Atom bairs should be close to each other and close to
            # adsorbate
            for atom in \
                [atom for atom in SC if np.linalg.norm(atom.position - pos1)
                 < 1.5 * bridge_dist and
                 np.linalg.norm(atom.position[:2] - ads_pos) < bridge_dist
                 and not (atom.position == pos1).all()]:
                pos2 = atom.position
                # bridge site is average position
                bridge_pos = 0.5 * (pos1 + pos2)
                k = 'bridge_site-{}'.format(b_index)
                D[k] = {}
                D[k]['pos'] = bridge_pos
                D[k]['sym'] = atomi.symbol + '_' + atom.symbol
                b_index += 1

        """Hollow sites: Voronoi vertices """
        hollow_dist = 1
        vor = Voronoi(SC.positions[:, :2])
        vertices = vor.vertices
        # map vertices close to adsorbate
        close_v = [v for v in vertices if np.linalg.norm(v - ads_pos[:2])
                   < hollow_dist]
        h_index = 1
        ff_index = 1
        for v in close_v:
            # Check if vertices overlap
            close_v_v = [v0 for v0 in close_v if np.linalg.norm(v0 - v) < 0.5]
            if len(close_v_v) > 1:
                v_mean = np.mean(close_v_v, axis=0)
                k = '4fold_{}'.format(ff_index)
                D[k] = {}
                D[k]['pos'] = v_mean

                # Delete bridge sites overlapping with 4fold
                for key in [key for key in list(D.keys()) if 'bridge' in key]:
                    bridge_pos = D[key]['pos']
                    if np.linalg.norm(bridge_pos[:2] - v) < 0.3:
                        del D[key]

                ffold_max_dist = sorted([np.linalg.norm(v_mean - m[:2])
                                         for m in SC.positions])[4]
                symb = [atom.symbol for atom in SC if
                        np.linalg.norm(v_mean - atom.position[:2])
                        < ffold_max_dist]

                D[k]['sym'] = '_'.join(symb)
                ff_index += 1
            else:  # Regular hollow site
                k = 'hollow_site-{}'.format(h_index)
                D[k] = {}
                D[k]['pos'] = v
                hollow_max_dist = sorted([np.linalg.norm(v - m[:2])
                                          for m in SC.positions])[3]
                symb = [atom.symbol for atom in SC
                        if np.linalg.norm(v - atom.position[:2])
                        < hollow_max_dist]
                D[k]['sym'] = '_'.join(symb)
                h_index += 1

        return D

    def get_subsurface_layer(self):
        return self.B[np.argsort(self.B.positions[:, 2])][-self.ntop * 2 - 1:
                                                          -self.ntop - 1]

    def get_under_bridge(self):
        """Return element closest to the adsorbate in the subsurface layer"""
        C0 = self.B[-1:] * (3, 3, 1)
        ads_pos = C0.positions[4]

        C = self.get_subsurface_layer() * (3, 3, 1)
        dis = self.B.cell[0][0] * 2

        ret = None

        for ele in C:
            new_dis = np.linalg.norm(ads_pos - ele.position)
            if new_dis < dis:
                dis = new_dis
                ret = ele.symbol

        return ret

    def get_under_hollow(self):
        """ Return HCP if an atom is present below the adsorbate in the 
        subsurface layer and FCC if not"""
        C0 = self.B[-1:] * (3, 3, 1)
        ads_pos = C0.positions[4]

        C = self.get_subsurface_layer() * (3, 3, 1)

        ret = 'FCC'
        if np.any([np.linalg.norm(ads_pos[:2] - ele.position[:2]) < 0.5 *
                   cradii[ele.number] for ele in C]):
            ret = 'HCP'

        return ret

    def get_site(self, plot_voronoi_sites=False):
        """Return primaty (top, bridge, hollow, 4fold) and
        secondary (chemical elements in close environment) site designation"""

        if self.dissociated:
            return 'dissociated', ''

        if self.is_desorbed():
            return 'desorbed', ''

        if self.is_subsurface():
            return 'subsurface', ''

        C0 = self.B[-1:] * (3, 3, 1)
        ads_pos = C0.positions[4]

        C = self.B.copy() * (3, 3, 1)

        # Use top layer and adsorbate to map sites
        Dict = self.get_site_dict(ads_pos[:2])

        primary_site = None
        dis = self.B.get_cell()[0][0]
        Kind = None

        values = [np.linalg.norm(ads_pos[:2] - d['pos'][:2])
                  for d in list(Dict.values())]
        if len(values) == 0:
            return 'N/A', ''
        idx = np.argmin(values)
        dis = values[idx]
        kind = list(Dict.keys())[idx]
        primary_site = kind.split('_')[0]

        if plot_voronoi_sites:  # View sampled sites
            X = self.B.copy()
            X = X * (3, 3, 1)
            del X[-1]
            for pos in Dict.values():
                add_adsorbate(X, 'X', position=(pos['pos'][:2]), height=0.2)
            view(X)

        if primary_site == 'top':
            site_type = Dict[kind]['sym']

        if primary_site == 'bridge':
            site_type = Dict[kind]['sym'] + '|' + self.get_under_bridge()

        elif primary_site == 'hollow':
            site_type = Dict[kind]['sym'] + '|' + self.get_under_hollow()

        elif primary_site == '4fold':
            site_type = Dict[kind]['sym']

        if dis > 0.5:
            primary_site += '-tilt'
            print('Warning: A strong site match could not be found!')
            print('  structure labeled as {}'.format(primary_site))

        return primary_site, site_type


def get_ads_dist(atoms, ads0, ads1='H'):
    index0 = [i for i in range(len(atoms)) if atoms[i].symbol == ads0]
    index1 = [i for i in range(len(atoms)) if atoms[i].symbol == ads1]
    dist = []
    D, D_len = get_distances(atoms.positions[index0],
                             atoms.positions[index1],
                             atoms.cell, pbc=True)

    return np.max(D_len)
