import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from scipy.optimize import curve_fit

sns.set_style('white')
sns.set_palette(sns.hls_palette(6, h=0.5,l=0.4,s=0.5))
font = {'size': 18}
plt.rc('font', **font)

# Etot(H) - 1/2 Etot(H2), DFT, QE, BEEF-vdW
H2_dissociation_energy = 3.3087431266500005
H2_zpe = 0.276
dG_H2_diss = H2_dissociation_energy - 0.5*H2_zpe

# # https://doi.org/10.1016/0022-2852(61)90111-4
De_H2_exp = 4.7469853


def file_to_df(file_name):
    """Read in file and return pandas data_frame.

    Parameters
    ----------
    filename : Filename including path.

    Returns
    -------
    df : pandas data frame
    """
    filename, file_extension = os.path.splitext(file_name)
    if file_extension=='.csv':
        df = pd.read_csv(file_name, sep=',', header=0).iloc[:,:]
    elif file_extension=='.tsv':
        df = pd.read_csv(file_name, sep='\t', header=0).iloc[:,:]
    else:
        print('Please provide valid csv or tsv file format with header names.')
    return df


def morse_norm(r, a):
    return (1-np.exp(-a*r))**2


def morse_diff(a,b,r):
    y1 = (1-np.exp(-a*r))**2
    y2 = (1-np.exp(-b*r))**2
    return y1-y2


class PES:
    """Set up a 1D potential energy surface (PES) for a state (hydrogen position).
    Args:
        filepath (str) : Path to data to which the morse potential will be fitted.
        plot: Whether to plot normalized data after read in.
        position: 'left' or 'right, determines curve symmetry:
                If the transferred hydrogen is above the position: 'left'.
                If the transferred hydrogen is above the position: 'right'
                Special key words: 'H2O' and 'H3O+':
                Automatically plots curves of 'H2O' and 'H3O+'.
        dHeq: Equilibrium distance between hydrogen and position in Angstroem.
        G_rel: Gibbs free energy of reactants/products incl. e*U, (H+ + e-), OH-, H2O in eV,
        with respect to the hydrogen atom in vacuum.
    Attributes:
        filepath: Path to data, containing tsv with first column distances
        and second column energy, no header.
        position: left or right.
        dHeq: see attributes dHeq.
        gH: State absolute potential well depth in eV.
        De: Depth of morse potential when gH is 0 eV.
        df: Normalized distance and energy values.
        a: Morse potential constant.
        morse_fit_error: Error of morse fit.
        U: Potential if state contains electron in eV,
        supply only for states which contain electron.
    """

    def __init__(self,
                 a=None,
                 De=None,
                 filepath=None,
                 position='left',
                 dHeq=1.,
                 G_rel=0.,
                 U=None,
                 plot=False,
                 preset_donor=None):
        self.preset_donor = preset_donor
        self.filepath = filepath
        self.position = position
        self.dHeq = dHeq
        self.df = None
        self.U = U
        if U is not None:
            self.gH = G_rel + self.U
        else:
            self.gH = G_rel

        self.a = None
        self.De = None
        self.morse_fit_error = None

        # Check if position H2O/H3O and set up curves automatically
        if self.preset_donor is not None:
            self.De = dG_H2_diss
            if self.preset_donor == 'H2O':
                self.a = 1.70815017
            elif self.preset_donor == 'H3O+':
                self.a = 2.14555189
            else:
                print('Possible preset donors are: H2O or H3O+')
            self.make_df()

        # Set up from file
        elif filepath is not None:
            self.preprocess()
            self.normalize()
            self.fit_morse()

        # Set up from parameters
        elif De is not None and a is not None:
            self.De = De
            self.a = a
            self.make_df()

        # Ensure data supply
        else:
            print('Provide Morse De and a, or filepath to data, see docstrings. ')

        # Adjust potential well depth
        self.gH -= self.De

        if plot:
            self.plot_morse()

    def plot_data(self):
        """Plots the read-in data."""
        fig, ax = plt.subplots()
        if os.path.isfile(self.filepath):
            ax.plot(self.df.distance, self.df.energy*self.De)
            ax.set_xlabel('Distance d-d$_{eq}$ (Å)')
            ax.set_ylabel('Energy E (eV)')
            ax.legend()
            ax.set_title(self.position + ' read in data')
            ax.set_ylim((0, 7))
            ax.set_xlim((0, 4))
            plt.show()
        return None

    def make_df(self):
        """Creates a pandas DataFrame from the morse-potential"""
        self.df = pd.DataFrame()
        self.df['distance'] = np.linspace(0, 10, 500) + self.dHeq
        self.df['energy'] = self.morse(self.df['distance'] - self.dHeq)
        self.normalize()
        return None

    def preprocess(self):
        """Pre-processes the data by smoothing vacuum energy fluctuations and computing De.
        """
        # Read in and preprocess data
        df = pd.read_csv(self.filepath, sep='\t', header=None)
        df.columns = ['distance', 'energy']

        # H2O ata are noisy at large d.
        if not self.position == 'left':
            for i in range(10):
                df = self._smoothen(df)

        # Get dissociation energy and set Emin = 0
        df.energy = df.energy - df.energy.min()
        self.De = df.energy.max()

        self.df = df
        return self.df

    def normalize(self):
        """Normalized the data by normalizing the total energy to energy/De (eV) (0--1 eV).
        The distance is transformed to distance relative to equilibrium distance.
        """
        df = self.df.copy(deep=True)

        # Min-Max scaling to 0-->1
        df.energy = df.energy / self.De

        # Get distance at minimum energy and set to 0.
        dmin = df[df.energy == df.energy.min()].distance.get_values()[0]
        if self.position is 'right':
            df.distance = dmin - df.distance
            dmax = df.distance.max() / 1.1
            df = df[df.distance < dmax]
        else:
            df.distance = df.distance - dmin
        self.df = df
        return self.df

    def fit_morse(self):
        """Fits a morse potential to the normalized data."""
        # Fit to get a
        popt, pcov = curve_fit(morse_norm, self.df.distance, self.df.energy)
        perr = np.sqrt(np.diag(pcov))
        self.a = popt
        self.morse_fit_error = perr

    def morse(self, r=None):
        """
        Returns absolute morse potential in eV.
        :param r: Array or list of distance points.
        :return: Morse potential values.
        """
        # Define the morse potential for this proton position.
        if r is None:
            r = self.df.distance
        if not self.position == 'left':
            return self.De * (1 - np.exp(-self.a * (self.dHeq - r)))**2 + self.gH
        else:
            return self.De * (1 - np.exp(-self.a * (r - self.dHeq)))**2 + self.gH

    def morse_norm(self, r):
        """
        Returns normalized morse potential.
        :param r: Array or list of distance points.
        :return: Normalized morse potential.
        """
        return (1 - np.exp(-self.a * r))**2

    def plot_morse(self):
        """Plots morse potential."""
        fig, ax = plt.subplots()

        if self.filepath is not None:
            if self.position == 'left':
                ax.plot(self.df.distance + self.dHeq, self.df.energy * self.De + self.gH)
            else:
                ax.plot(self.dHeq - self.df.distance, self.df.energy * self.De + self.gH)

        x = np.linspace(-10, 10, 500)
        ax.plot(x, self.morse(r=x), '--',
                label='Morse: a=%5.3f' % self.a)
        ax.axhline(color='k', zorder=0, lw=0.6)
        ax.set_xlabel('Distance d (Å)')
        ax.set_ylabel('Energy E (eV)')
        ax.legend()
        ax.set_title(self.position + ' morse fit')
        ax.set_ylim((-10, 10))
        plt.show()
        return None

    def _smoothen(self, df):
        """
        Helper funciton to smooth out vacuum energy fluctuations (DFT problem).
        :param df: A DataFrame with distance and energy values.
        :return: DataFrame with smooth data.
        """
        df_new = df.copy(deep=True)
        for idx, entry in enumerate(df.energy):
            if idx == 0:
                previos_entry = entry
            else:
                dE = entry - previos_entry
                if dE < 0:
                    df_new = df_new.drop([idx], axis=0)
                previos_entry = entry
        df_new = df_new.reset_index(drop=True)
        return df_new


class Energy:
    """Compute diabatic and adiabatic energy interceptions (activation energies).
    Args:
        pes1 and pes2 : PES objects.
    Attributes:
        left: Left state
        right: Right state (Note that this has nothing to do
            with the orientation of the fitting data as in PES,
            but simply denotes where to position the morse potential.
        x: Distance points
        r_corr: Distance points between the two states
        morse_left: Returns left morse energy values given distance values.
        morse_right: Returns right morse energy values given distance values.
        Ea_left: Forward diabatic barrier in eV.
        Ea_right: Backwards diabatic barrier in eV.
        xint: Distance between hydrogen and left state at diabatic transition state in Angstroem.
        Ea_ad_left: Forward adiabatic barrier in eV.
        Ea_ad_right: Backwards adiabatic barrier in eV.
        xint_ad: Distance between hydrogen and left state at adiabatic transition state in Angstroem.
    """

    def __init__(self, left, right):
        self.x = np.linspace(-10., 10., 5000)
        self.left = left
        self.right = right
        self.r_corr = np.linspace(self.left.dHeq, self.right.dHeq, 1000)
        self.adiabatic_correction()
        self._beta = None

    def morse_left(self, r=None):
        """
        Returns Morse potential for left proton donor.
        :param r:  Array or list of distance points.
        :return: Morse potential energy values.
        """
        if r is None:
            r = self.x
        return self.left.De * (1 - np.exp(-self.left.a * (r - self.left.dHeq))) ** 2 + self.left.gH

    def morse_right(self, r=None):
        """
        Returns Morse potential for right proton donor.
        :param r: Array or list of distance points.
        :return: Morse potential energy values.
        """
        if r is None:
            r = self.x
        return self.right.De * (1 - np.exp(-self.right.a * (self.right.dHeq - r))) ** 2 + self.right.gH

    def interception(self,
                     adiabatic=False,
                     plot=False,
                     xlim=(0.5, 4.),
                     ylim=(-0.5, 4.)):
        """
        Calculates and, if desired, plots interception between Morse potentials.
        :param adiabatic: Whether to calculate the interception with adiabatic correction.
        :param plot: Whether to plot resulting potential energy curves and intercepts.
        :param xlim: Adjust x limits for plots, supply tuple like this: xlim=(1,2).
        :param ylim: Adjust y limits for plots, supply tuple like this: ylim=(1,2).
        :return: xint is the distance of the transition state,
        yint is the y position of the transition state (refer to potential well depth to
        get the activation energy)
        """
        a = self.morse_left()
        b = self.morse_right()

        # Find intercept
        xint_list = list(self.x[np.argwhere(np.diff(np.sign(a - b))).flatten()])
        yint_list = []
        for xi in xint_list:
            yint_list.append(self.morse_left(r=xi))

        if len(xint_list) == 1 and xint_list[0] < self.left.dHeq:
            self.xint = self.left.dHeq
            self.yint = -self.left.De
        else:
            yint_list = [a[0] for a in yint_list]
            val, idx = min((val, idx) for (idx, val) in enumerate(yint_list))
            self.xint = xint_list[idx]
            self.yint = val

        self.Ea_left = self.yint + self.left.De
        self.Ea_right = self.yint + self.right.De

        if plot:
            fig, ax = plt.subplots(figsize=(8, 5))

            # Plot diabatic curves and intercept
            ax.plot(self.x, a, label='dia-left')
            ax.plot(self.x, b, label='dia-right')
            ax.plot([self.xint], [self.yint], marker='o', markersize=6.0,
                    c='grey', markeredgecolor='k', ls='',
                    label='E$^{a}_{left}$ = %5.2f eV' % self.Ea_left)
            ax.plot([self.xint], [self.yint], marker='o', markersize=6.0,
                    c='grey', markeredgecolor='k', ls='',
                    label='E$^{a}_{right}$ = %5.2f eV' % self.Ea_right)

            # Plot adiabatic curves and intercept
            if adiabatic:
                ax.plot(self.r_corr, self.adia_left, label='adia-left')
                ax.plot(self.r_corr, self.adia_right, label='adia-right')
                ax.plot([self.xint_ad], [self.yint_ad], marker='o', markersize=6.0,
                        c='w', markeredgecolor='b', ls='',
                        label='E$^{a}_{left}$ = %5.2f eV' % (self.Ea_ad_left))
                ax.plot([self.xint_ad], [self.yint_ad], marker='o', markersize=6.0,
                        c='w', markeredgecolor='b', ls='',
                        label='E$^{a}_{right}$ = %5.2f eV' % (self.Ea_ad_right))

            ax.set_xlabel('Distance to right hydrogen position (Å)')
            ax.set_ylabel('Energy (eV)')

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.legend(loc=2, bbox_to_anchor=(1.0, 1.0))
            plt.show()
        return self.xint, self.yint

    @ property
    def beta_left(self):
        """Calculates charge transfer coefficient for forward reaction."""
        if not self._beta:
            _ = self.interception()
            self._beta_left = self.left.morse_norm(self.xint-self.left.dHeq)
        return self._beta_left

    @ property
    def beta_right(self):
        """Calculates charge transfer coefficient for backward reaction."""
        if not self._beta:
            _ = self.interception()
            self._beta_right = self.right.morse_norm(self.left.dHeq-self.xint)
        return self._beta_right

    def adiabatic_correction(self):
        """
        Calculates the adiabatic curves, intercepts and activation barriers.
        :return: Forward and backward adiabatic activation energy.
        """
        # Get Gamma values between 0 and 1 for all distances of interest
        gamma_left = self.left.morse_norm(self.r_corr - self.left.dHeq)       # 0-->0.66
        gamma_right = self.right.morse_norm(self.right.dHeq - self.r_corr)    # 0.66-->0

        # Morse values
        E_stretch_left = gamma_left * self.left.De + self.left.gH
        E_stretch_right = gamma_right * self.right.De + self.right.gH

        # E_hyb
        self.adia_left = E_stretch_left + gamma_left * (1-gamma_right) * (self.right.gH - E_stretch_left)
        self.adia_right = E_stretch_right + gamma_right * (1-gamma_left) * (self.left.gH - E_stretch_right)

        yts_ad = [min(self.adia_left[i],self.adia_right[i]) for i, _ in enumerate(self.adia_right)]
        self.yint_ad, idx = max((val, idx) for (idx, val) in enumerate(yts_ad))
        self.xint_ad = self.r_corr[idx]

        self.Ea_ad_left = self.yint_ad + self.left.De
        self.Ea_ad_right = self.yint_ad + self.right.De

        return self.Ea_ad_left, self.Ea_ad_right


if __name__ == '__main__':
    print('Executed without errors.')
