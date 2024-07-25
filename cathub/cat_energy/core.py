'''
Module to perform free energetic analysis for CatHub ASE database
'''
import pandas as pd
from .io import make_mkm_input_files, write_columns
from .conversion import u_she_to_u_rhe

class EnergeticAnalysis:
    """
    Perform energetic analysis for CatHub ASE database.

    This class handles the computation of energetics for gases, adsorbates, and
    transition states based on the provided parameters. It integrates various
    corrections and external effects to compute the free energy changes and 
    optionally writes the corrected energy data to input files for microkinetic 
    modeling (MKM).

    Attributes:
    -----------
    db_filepath : pathlib.Path
        Path to the ASE database file containing the reaction data.
    system_parameters : dict
        Dictionary containing system parameters such as:
            'desired_surface' (str): The surface material being analyzed.
            'desired_facet' (str): The facet of the surface material being analyzed.
            'temp' (float): Temperature in Kelvin.
            'pH' (float): pH value of the system.
            'u_she' (float): Standard Hydrogen Electrode potential.
    reference_gases : list of str
        List of reference gas species.
    dummy_gases : list of str
        List of dummy gas species used for computations.
    dft_corrections_gases : dict
        Dictionary of DFT corrections for gas species.
    fake_ads : dict
        Dictionary of fake adsorbates used for testing or placeholder purposes.
    beef_dft_helmholtz_offset : dict
        Dictionary of Helmholtz offset values for BEEF-DFT corrections, keyed by species.
    external_effects : dict
        Dictionary of external effects to be considered in the calculations.
    facet_conditional : str
        String to indicate if facet-specific calculations are required.
    gas_jsondata_filepath : pathlib.Path
        Path to the JSON file containing gas species data.
    ads_jsondata_filepath : pathlib.Path
        Path to the JSON file containing adsorbate species data.
    ts_jsondata_filepath : pathlib.Path
        Path to the JSON file containing transition state species data.
    rxn_expressions_filepath : pathlib.Path
        Path to the file containing reaction expressions for transition states.
    ts_data : dict
        Dictionary containing data for transition states, structured as:
            'barrier_fits' (dict): Nested dictionaries containing backward barrier fits for each reaction.
            'ts_states' (dict): Dictionary containing transition state information, with keys:
                'wf_data' (list of float): Workfunction data for initial, transition, and final states.
                'neb_image_id_range' (list of int): Range of NEB image IDs for the transition state.
            'extrapolation' (dict): Dictionary containing extrapolation parameters:
                'perform' (bool): Whether to perform extrapolation.
                'method' (str): Method of extrapolation (e.g., 'charge').
                'potential_type' (str): Type of potential ('fixed' or 'potential-dependent').
                'potential_value' (float): Value of the potential if 'fixed' type is used.
            'phi_correction' (float): Phi correction value.
            'alk_corr' (float): Alkaline correction value.
    write_mkm_input_files : bool
        Boolean to indicate if MKM input files should be generated.
    verbose : bool, optional
        If True, print detailed information about the calculations. Default is True.
    latex : bool, optional
        If True, format the output for LaTeX. Default is True.

    Methods:
    --------
    write_energies(write_species)
        Perform computations and return energetics for specified species.
    """
    def __init__(self, db_filepath, system_parameters, reference_gases, dummy_gases,
                 dft_corrections_gases, fake_ads, beef_dft_helmholtz_offset, external_effects,
                 facet_conditional, gas_jsondata_filepath, ads_jsondata_filepath,
                 ts_jsondata_filepath, rxn_expressions_filepath, ts_data, write_mkm_input_files=True,
                 verbose=True, latex=True):
        self.db_filepath = db_filepath
        self.system_parameters = system_parameters
        self.reference_gases = reference_gases
        self.dummy_gases = dummy_gases
        self.dft_corrections_gases = dft_corrections_gases
        self.fake_ads = fake_ads
        self.beef_dft_helmholtz_offset = beef_dft_helmholtz_offset
        self.external_effects = external_effects
        self.facet_conditional = facet_conditional
        self.gas_jsondata_filepath = gas_jsondata_filepath
        self.ads_jsondata_filepath = ads_jsondata_filepath
        self.ts_jsondata_filepath = ts_jsondata_filepath
        self.rxn_expressions_filepath = rxn_expressions_filepath
        self.ts_data = ts_data
        self.write_mkm_input_files = write_mkm_input_files
        self.verbose = verbose
        self.latex = latex

        self._calculate_u_rhe()

    def _calculate_u_rhe(self):
        """ Calculate and set u_rhe in system parameters """
        u_she = self.system_parameters['u_she']
        temp = self.system_parameters['temp']
        self.system_parameters['u_rhe'] = u_she_to_u_rhe(u_she, temp, self.system_parameters['pH'])

    def _print_verbose_info(self):
        """ Print detailed information about the calculations if verbose is True """
        system_header = (f'###### {self.system_parameters["desired_surface"]}'
                         f'_{self.system_parameters["desired_facet"]}: '
                         f'SHE Potential = {self.system_parameters["u_she"]:.2f} V ######')
        print('-' * len(system_header))
        print(system_header)
        print('-' * len(system_header))
        print()
        print('Term1 = Calculated Electronic Energy + DFT Correction')
        print('Term2 = Enthalpic Temperature Correction + Entropy Contribution')
        print('Term3 = RHE-scale Dependency')
        print('Term4 = External Effect Corrections')
        print('Chemical Potential, µ = Term1 + Term2 + Term3 + Term4')
        print('Free Energy Change, ∆G = µ_species - µ_ref. For example, '
              '∆G_CH4 = µ_CH4 - (µ_CO + 3 * µ_H2(ref) - µ_H2O)')
        print('∆G at U_RHE=0.0 V = ∆G - Term3')
        print()


def write_energies(analysis, write_species):
    '''
    Perform computations and return energetics for specified species.

    This function handles the computation of energetics for gases, adsorbates, 
    and transition states based on the provided parameters. It integrates various 
    corrections and external effects to compute the free energy changes and 
    optionally writes the corrected energy data to input files for microkinetic 
    modeling (MKM).

    Parameters:
    -----------
    write_species : dict
        Dictionary indicating which species (gases, adsorbates, transition states) to process.

    Returns:
    --------
    pandas.DataFrame
        Dataframe containing the computed energetics of the requested species.
    '''
    df_out = pd.DataFrame(columns=write_columns)

    if analysis.verbose:
        analysis._print_verbose_info()

    if write_species['gases']:
        from .gas_species import GasAnalysis
        gas_analysis = GasAnalysis(analysis)
        df_out = gas_analysis.write_gas_energies(df_out)

    if write_species['adsorbates']:
        from .ads_species import AdsAnalysis
        ads_analysis = AdsAnalysis(analysis)
        df_out = ads_analysis.write_adsorbate_energies(df_out)

    if write_species['transition_states']:
        from .ts_species import TSAnalysis
        ts_analysis = TSAnalysis(analysis)
        exec(compile(open(analysis.rxn_expressions_filepath, 'rb').read(), '<string>', 'exec'))
        df_out = ts_analysis.write_ts_energies(df_out, locals()['rxn_expressions'])

    # write corrected energy data to mkm input file
    if analysis.write_mkm_input_files:
        make_mkm_input_files(analysis.db_filepath, analysis.system_parameters, df_out)

    return df_out
