'''
Module to perform free energetic analysis for CatHub ASE database
'''
import pandas as pd
from .io import make_mkm_input_files, write_columns
from .ts_species import write_ts_energies
from .ads_species import write_adsorbate_energies
from .gas_species import write_gas_energies
from .conversion import u_she_to_u_rhe

def write_energies(
        db_filepath, reference_gases, dummy_gases, dft_corrections_gases,
        fake_ads, beef_dft_helmholtz_offset, external_effects, system_parameters,
        facet_conditional, write_species, gas_jsondata_filepath,
        ads_jsondata_filepath, ts_jsondata_filepath, rxn_expressions_filepath,
        ts_data, write_mkm_input_files, verbose=True, latex=True):
    '''
    Function to delegate computation and returning energetics of requested
    species
    '''
    df_out = pd.DataFrame(columns=write_columns)
    u_she = system_parameters['u_she']
    temp = system_parameters['temp']
    u_rhe = u_she_to_u_rhe(u_she, temp, system_parameters['pH'])
    system_parameters['u_rhe'] = u_rhe

    if verbose:
        system_header = (f'###### {system_parameters["desired_surface"]}'
                         f'_{system_parameters["desired_facet"]}: '
                         f'SHE Potential = {u_she:.2f} V ######')
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

    if write_species['gases']:
        df_out = write_gas_energies(db_filepath, df_out, gas_jsondata_filepath,
                                    system_parameters, reference_gases,
                                    dummy_gases, dft_corrections_gases,
                                    beef_dft_helmholtz_offset, external_effects,
                                    verbose, latex)

    if write_species['adsorbates']:
        df_out = write_adsorbate_energies(db_filepath, df_out,
                                          ads_jsondata_filepath,
                                          system_parameters, facet_conditional,
                                          reference_gases, fake_ads,
                                          dft_corrections_gases,
                                          external_effects, verbose, latex)

    if write_species['transition_states']:
        exec(compile(open(rxn_expressions_filepath, 'rb').read(), '<string>',
                     'exec'))
        df_out = write_ts_energies(db_filepath, df_out, ts_jsondata_filepath,
                                   locals()['rxn_expressions'], ts_data,
                                   system_parameters, reference_gases,
                                   external_effects, verbose, latex)

    # write corrected energy data to mkm input file
    if write_mkm_input_files:
        make_mkm_input_files(db_filepath, system_parameters, df_out)
    return df_out
