#%%

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pathlib import Path
import seaborn as sns
from scipy.optimize import curve_fit

sns.set_style('white')
sns.set_palette(sns.hls_palette(6, h=0.5, l=0.4, s=0.5))
font = {'size': 18}
plt.rc('font', **font)

# Etot(H) - 1/2 Etot(H2), DFT, QE, BEEF-vdW
H2_dissociation_energy = 3.3087431266500005
H2_zpe = 0.276
dG_H2_diss = H2_dissociation_energy - 0.5*H2_zpe

# # https://doi.org/10.1016/0022-2852(61)90111-4
De_H2_exp = 4.7469853

proton_donors = {
    'H2O':  {'a': 1.70815017,
             'De': dG_H2_diss},
    'H3O+': {'a': 2.14555189,
             'De': dG_H2_diss}
    }


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


def preprocess(filepath, position):
    """Pre-processes the data by smoothing vacuum energy fluctuations and normalizing the data."""
    # Read in and preprocess data
    df = pd.read_csv(filepath, sep='\t', header=None)
    df.columns = ['distance', 'energy']

    # H2O ata are noisy at large d.
    if not position == 'left':
        for i in range(10):
            df = smoothen(df)

    # Get dissociation energy and set Emin = 0
    df.energy = df.energy - df.energy.min()
    De = df.energy.max()

    # Min-Max scaling to 0-->1
    df.energy = df.energy / De

    # Get distance at minimum energy and set to 0.
    dmin = df[df.energy == df.energy.min()].distance.get_values()[0]
    if position is 'right':
        df.distance = dmin - df.distance
        dmax = df.distance.max() / 1.1
        df = df[df.distance < dmax]
    else:
        df.distance = df.distance - dmin

    return df, De


def smoothen(df):
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


def morse_norm(r, a):
    return (1-np.exp(-a*r))**2


def morse_diff(a,b,r):
    y1 = (1-np.exp(-a*r))**2
    y2 = (1-np.exp(-b*r))**2
    return y1-y2


def fit_morse(df):
    """Fits a morse potential to the normalized dataframe."""
    popt, pcov = curve_fit(morse_norm, df.distance, df.energy)
    perr = np.sqrt(np.diag(pcov))
    a = popt
    morse_fit_error = perr
    return a, morse_fit_error


class PES:
    """Set up a 1D potential energy surface (PES) for a state (hydrogen position).

    Args:
        position: 'left' or 'right, determines curve symmetry:
                If the transferred hydrogen is above the position: 'left'.
                If the transferred hydrogen is above the position: 'right'
                Special key words: 'H2O' and 'H3O+':
                Automatically plots curves of 'H2O' and 'H3O+'.
        d_Heq: Equilibrium distance between hydrogen and position in Angstroem.
        g_rel: Gibbs free energy of reactants/products incl. e*U, (H+ + e-), OH-, H2O in eV,
        with respect to the hydrogen atom in vacuum.
        potential: Potential vs. SHE of this state.

    Attributes:
        filepath: Path to data, containing tsv with first column distances
        and second column energy, no header.
        position: left or right.
        d_Heq: see attributes d_Heq.
        g_H: State absolute potential well depth in eV.
        De: Depth of morse potential when g_H is 0 eV.
        df: Normalized distance and energy values.
        a: Morse potential constant.
        morse_fit_error: Error of morse fit.
        U: Potential if state contains electron in eV,
        supply only for states which contain electron.
    """

    def __init__(
            self,
            position='left',
            d_Heq=1.0,
            potential=0.0,
            g_rel=0.0,
            a=None,
            De_U0=None,
            filepath=None,
            proton_donor=None,
            df=None,
            morse_fit_error=None
    ):
        self.position = position
        self.d_Heq = d_Heq
        self.U = potential
        self.g_rel = g_rel

        self.a = a
        self.De_U0 = De_U0
        self.filepath = filepath
        self.proton_donor = proton_donor
        self._df = df
        self.morse_fit_error = morse_fit_error

        self.g_H = self.g_rel + self.U
        self.De = self.De_U0 - self.g_H

    @classmethod
    def init_from_file(cls,
                       filepath: str = None,
                       position: str = 'left',
                       **kw):
        if Path(filepath).exists():
            df, De_U0 = preprocess(filepath, position)
            a, morse_fit_error = fit_morse(df)
            return cls(filepath=filepath,
                       position=position,
                       De_U0=De_U0,
                       df=df,
                       a=a,
                       morse_fit_error=morse_fit_error,
                       **kw)

    @classmethod
    def init_from_parameters(cls,
                             De_U0: float,
                             a: float,
                             **kw):
        return cls(De_U0=De_U0,
                   a=a,
                   **kw)

    @classmethod
    def init_from_database(cls,
                           proton_donor: str ='H2O',
                           **kw):
        h_dict = proton_donors[proton_donor]
        return cls(a=h_dict['a'],
                   De_U0=h_dict['De'],
                   **kw)

    @property
    def df(self):
        if self._df is None:
            df = pd.DataFrame()
            df['distance'] = np.linspace(-10, 10, 5000)
            df['energy'] = self.morse(r=df['distance'])
            self._df = df
        return self._df

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
            return self.De * (1 - np.exp(-self.a * (self.d_Heq - r)))**2 - self.De
        else:
            return self.De * (1 - np.exp(-self.a * (r - self.d_Heq)))**2 - self.De

    def morse_norm(self, r):
        """
        Returns normalized morse potential.
        :param r: Array or list of distance points.
        :return: Normalized morse potential.
        """
        return (1 - np.exp(-self.a * r))**2

    def plot_morse(self,
                   xlim=(-0,5),
                   ylim=(-5,5),
                   title=None):
        """
        Plots morse potential.
        :param title: Figure title.
        :param xlim: x-axis limits.
        :param ylim: y-axis limits.
        :return:
        """
        fig, ax = plt.subplots()
        x = np.linspace(-10, 10, 500)
        ax.plot(x, self.morse(r=x), '--',
                label='Morse: a=%5.3f' % self.a)
        ax.axhline(color='k', zorder=0, lw=0.6)
        ax.set_xlabel('Distance d (Å)')
        ax.set_ylabel('Energy E (eV)')
        ax.legend()
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.show()
        return None


class Energy:
    """Compute diabatic and adiabatic energy interceptions (activation energies).
    Args:
        pes1 and pes2 : PES objects.
    Attributes:
        left: Left state
        right: right state (Note that this has nothing to do
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
        self.r_corr = np.linspace(self.left.d_Heq, self.right.d_Heq, 1000)
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
        return self.left.De * (1 - np.exp(-self.left.a * (r - self.left.d_Heq))) ** 2 - self.left.De

    def morse_right(self, r=None):
        """
        Returns Morse potential for right proton donor.
        :param r: Array or list of distance points.
        :return: Morse potential energy values.
        """
        if r is None:
            r = self.x
        return self.right.De * (1 - np.exp(-self.right.a * (self.right.d_Heq - r))) ** 2 - self.right.De

    def interception(self,
                     adiabatic=False,
                     plot=False,
                     xlim=(0.5, 4.),
                     ylim=(-7, 2.)):
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

        if len(xint_list) == 1 and xint_list[0] < self.left.d_Heq:
            self.xint = self.left.d_Heq
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
                        label='E$^{a}_{left}$ = %5.2f eV' % self.Ea_ad_left)
                ax.plot([self.xint_ad], [self.yint_ad], marker='o', markersize=6.0,
                        c='w', markeredgecolor='b', ls='',
                        label='E$^{a}_{right}$ = %5.2f eV' % self.Ea_ad_right)

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
            self._beta_left = self.left.morse_norm(self.xint-self.left.d_Heq)
        return self._beta_left

    @ property
    def beta_right(self):
        """Calculates charge transfer coefficient for backward reaction."""
        if not self._beta:
            _ = self.interception()
            self._beta_right = self.right.morse_norm(self.xint-self.right.d_Heq)
        return self._beta_right

    def adiabatic_correction(self):
        """
        Calculates the adiabatic curves, intercepts and activation barriers.
        :return: Forward and backward adiabatic activation energy.
        """
        # Get Gamma values between 0 and 1 for all distances of interest
        gamma_left = self.left.morse_norm(self.r_corr - self.left.d_Heq)       # 0-->0.66
        gamma_right = self.right.morse_norm(self.right.d_Heq - self.r_corr)    # 0.66-->0

        # Morse values
        E_stretch_left = (gamma_left - 1) * self.left.De
        E_stretch_right = (gamma_right - 1) * self.right.De

        # E_hyb
        self.adia_left = E_stretch_left - gamma_left * (1-gamma_right) * (self.right.De + E_stretch_left)
        self.adia_right = E_stretch_right - gamma_right * (1-gamma_left) * (self.left.De + E_stretch_right)

        yts_ad = [min(self.adia_left[i],self.adia_right[i]) for i, _ in enumerate(self.adia_right)]
        self.yint_ad, idx = max((val, idx) for (idx, val) in enumerate(yts_ad))
        self.xint_ad = self.r_corr[idx]

        self.Ea_ad_left = self.yint_ad + self.left.De
        self.Ea_ad_right = self.yint_ad + self.right.De

        return self.Ea_ad_left, self.Ea_ad_right


if __name__ == '__main__':
    print('Executed without errors.')
