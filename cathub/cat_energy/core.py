'''
Module to perform free energetic analysis for CatHub ASE database
'''
import pandas as pd
from .io import make_mkm_input_files, write_columns
from .ts_species import write_ts_energies
from .ads_species import write_adsorbate_energies
from .gas_species import write_gas_energies


def write_energies(
        db_filepath, reference_gases, dummy_gases, dft_corrections_gases,
        beef_dft_helmholtz_offset, field_effects, adsorbate_parameters,
        facet_conditional, write_gases, write_adsorbates,
        write_transition_states, gas_jsondata_filepath, ads_jsondata_filepath,
        ts_jsondata_filepath, ts_data, temp, pH, write_mkm_input_files,
        verbose=True, latex=True):
    '''
    Function to delegate computation and returning energetics of requested
    species
    '''
    df_out = pd.DataFrame(columns=write_columns)
    if verbose:
        system_header = (f'###### {adsorbate_parameters["desired_surface"]}'
                         f'_{adsorbate_parameters["desired_facet"]}: '
                         f'Electric Field Strength = '
                         f'{field_effects["epsilon"]:.2f} V/A ######')
        print('-' * len(system_header))
        print(system_header)
        print('-' * len(system_header))

        print()
        print('Term1 = Calculated Electronic Energy + DFT Correction')
        print('Term2 = Enthalpic Temperature Correction + Entropy Contribution')
        print('Term3 = RHE-scale Dependency')
        print('Term4 = Solvation Correction + Electric Field Correction')
        print('Chemical Potential, µ = Term1 + Term2 + Term3 + Term4')
        print('Free Energy Change, ∆G = µ_species - µ_ref. For example, '
              '∆G_CH4 = µ_CH4 - (µ_CO + 3 * µ_H2(ref) - µ_H2O)')
        print('∆G at U_RHE=0.0 V = ∆G - Term3')
        print()

    if write_gases:
        df_out = write_gas_energies(db_filepath, df_out, gas_jsondata_filepath,
                                    reference_gases, dummy_gases,
                                    dft_corrections_gases,
                                    beef_dft_helmholtz_offset, field_effects,
                                    temp, verbose, latex)

    if write_adsorbates:
        df_out = write_adsorbate_energies(db_filepath, df_out,
                                          ads_jsondata_filepath,
                                          adsorbate_parameters,
                                          facet_conditional, reference_gases,
                                          dft_corrections_gases, field_effects,
                                          temp, verbose, latex)

    if write_transition_states:
        if verbose:
            ts_phase_header = 'Transition State Free Energy Correction:'
            print(ts_phase_header)
            print('-' * len(ts_phase_header))

            print('Term1 = Backward Electronic Activation Energy '
                  '+ DFT Correction')
            print('Term2 = Enthalpic Temperature Correction '
                   '+ Entropy Contribution')
            print('Term3 = RHE-scale Dependency')
            if ts_data['extrapolation']:
                print('Term4 = Solvation Correction + Electric Field Correction'
                      ' + Alkaline Correction + Charge Extrapolation Correction'
                      ' + Final Adsorbate Energy')
            else:
                print('Term4 = Solvation Correction + Electric Field Correction'
                      ' + Alkaline Correction + Final Adsorbate Energy')
            print('Free Energy Change, ∆G = Term1 + Term2 + Term3 + Term4')
            print('∆G at U_RHE=0.0 V = ∆G - Term3')
            print()
        # exec(compile(open(rxn_expressions_filepath, 'rb').read(), '<string>',
        #              'exec'))
        df_out = write_ts_energies(db_filepath, df_out, ts_jsondata_filepath,
                                   locals()['rxn_expressions'], ts_data,
                                   adsorbate_parameters, field_effects,
                                   temp, pH, verbose, latex)

    # write corrected energy data to mkm input file
    if write_mkm_input_files:
        make_mkm_input_files(db_filepath, adsorbate_parameters, field_effects,
                             df_out)
    return df_out
