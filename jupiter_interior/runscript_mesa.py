################################################
# MESA Run Script for make_planet suite        #
# mugalino                                     #
################################################

import math
import numpy as np
import os
import shutil
import sys
from utilities import *
import constants
import warnings

# Initializing dictionary of relevant parameters
params = {}

# fiducial parameters
params['mp'] = 1.67 * constants.Mjup_cgs          # planet mass in mjup
params['rp'] = 2.00 * constants.Rjup_cgs           # inital planet radius in rjup
params['mcore'] = 10.0 * constants.Mearth_cgs      # planet core mass in mearth
params['rhocore'] = 10.0                           # core density in g/cc

# mass of planet before core put in
params['mp_wo_core'] = params['mp'] - params['mcore'] 

params['z'] = 0.02                                # metallicity of both planet and star
params['y'] = 0.24                                # helium fraction of both planet and (initial) star
params['maxage_1'] = 50.                          # ending age
params['maxage_2'] = 2.e3                         # ending age
params['maxage_3'] = 1.e10                        # ending age

# parameters having to do with irradiation
ms = 0.903                                         # star mass in msun
rs = 0.910                                         # star radius in rsun
Teff_star = 5617                                  # stellar Teff
orb_sep = 0.02335                                   # orbital separation in AU (to set day-side flux)
params['irrad_col'] = 300.0                       # column depth for depositing stellar radiation as heat

# flux hitting planet's dayside
params['flux_dayside'] = constants.sigSB_cgs * Teff_star**4 *   \
                        ( rs * constants.Rsolar_cgs / orb_sep / constants.au_cgs )**2 
params['Teq'] = (params['flux_dayside'] / 4.0 / constants.sigSB_cgs)**0.25

params['heating_gamma'] = 1.e-4                   # lowercase gamma, Gamma / Lstar
params['heating_Gamma'] = params['heating_gamma'] * constants.sigSB_cgs * Teff_star**4 *  \
                        ( 4 * np.pi * (rs * constants.Rsolar_cgs) ** 2.) * (params['rp'] / (2.0 * orb_sep * constants.au_cgs))**2.0

# pressure depth, Pdep
pdep = 100. # pressure in bars
params['pdep'] = pdep * constants.bar_cgs # cgs units for pressure

def main(params_new = None, verbose = True, run_sims = True, depth_heating=False):
    # updating some parameters of interest
    if (params_new != None) and depth_heating:
        params.update({'pdep' : params_new[0]})
        params.update({'heating_gamma': params_new[1]})
        params.update({'z' : params_new[2]})
    elif (params_new != None):
        params.update({'heating_gamma': params_new[0]})
        params.update({'z' : params_new[1]})
    else:
        None

    # flags to skip steps
    do_create_planet = True
    do_put_in_core = False
    do_evolve_planet = True

    # Choose between uniform heating or depth-dependent heating
    if depth_heating:
        do_use_other_energy = '.true.'
        params.update({'heating_Gamma':0.e0})
        heating_mode = 'depthdependent_{}bar_pdep'.format(params['pdep']/constants.bar_cgs)
    else:
        do_use_other_energy = '.false.'
        heating_mode = 'uniform'

    if verbose:
        print_parameters(params, do_use_other_energy)

    createmodel = "planet_create_{:.2f}_MJ_{:.2f}_ME_{:.2f}_RJ_{:.5f}_gamma_{}_Zmetal_{}.mod".format(params['mp'] / constants.Mjup_cgs, \
                                                                            params['mcore'] / constants.Mearth_cgs, params['rp'] / constants.Rjup_cgs, \
                                                                            params['heating_gamma'], params['z'], heating_mode)
    coremodel = "planet_core_{:.2f}_MJ_{:.2f}_ME_{:.2f}_RJ_{:.5f}_gamma_{}_Zmetal_{}.mod".format(params['mp'] / constants.Mjup_cgs, \
                                                                            params['mcore'] / constants.Mearth_cgs, params['rp'] / constants.Rjup_cgs, \
                                                                            params['heating_gamma'], params['z'], heating_mode)
    evolvemodel = "planet_evolve_{:.2f}_MJ_{:.2f}_ME_{:.2f}_RJ_{:.5f}_gamma_{}_Zmetal_{}.mod".format(params['mp'] / constants.Mjup_cgs, \
                                                                            params['mcore'] / constants.Mearth_cgs, params['rp'] / constants.Rjup_cgs, \
                                                                            params['heating_gamma'], params['z'], heating_mode)
    if do_create_planet:
        inlist1 = "inlist_create_{:.2f}_MJ_{:.2f}_ME_{:.2f}_RJ_{}_Zmetal_{}".format(params['mp'] / constants.Mjup_cgs, params['mcore'] / constants.Mearth_cgs, \
                                                                            params['rp'] / constants.Rjup_cgs, params['z'], heating_mode)
        photdir = inlist1[7:] + "_photos"
        logdir = inlist1[7:] + "_logs"

        if not os.path.exists(photdir): os.mkdir(photdir)
        if not os.path.exists(logdir): os.mkdir(logdir)

        if run_sims: run_time = create_planet(params, inlist1, createmodel, run_sims, logdir, photdir)

    if not os.path.exists(createmodel) and run_sims:
        print("[do_create_planet] MESA failed creating the planet. Check your parameters.")

    if do_put_in_core and (params['mcore'] > 0.0):
            inlist2 = "inlist_core_{:.2f}_MJ_{:.2f}_ME_{:.2f}_RJ_{}_Zmetal_{}".format(params['mp'] / constants.Mjup_cgs, params['mcore'] / constants.Mearth_cgs, \
                                                                            params['rp'] / constants.Rjup_cgs, params['z'], heating_mode)

            photdir = inlist2[7:] + "_photos"
            logdir = inlist2[7:] + "_logs"

            if not os.path.exists(photdir): os.mkdir(photdir)
            if not os.path.exists(logdir): os.mkdir(logdir)

            if run_sims: run_time = put_core_in_planet(params, inlist2, createmodel, coremodel, run_sims, logdir, photdir)

    else:
        shutil.copyfile(createmodel, coremodel)
        warnings.warn('Code did not generate an actual core model because Mcore =< 0.0.') 

    if not os.path.exists(coremodel) and run_sims:
        print("[do_put_in_core] MESA failed creating core model. Check your parameters.")

    if do_evolve_planet:
        inlist3 = "inlist_evolve_{:.2f}_MJ_{:.2f}_ME_{:.2f}_RJ_{:.5f}_gamma_{}_Zmetal_{}".format(params['mp'] / constants.Mjup_cgs, \
                                                                            params['mcore'] / constants.Mearth_cgs, params['rp'] / constants.Rjup_cgs, \
                                                                            params['heating_gamma'], params['z'], heating_mode)
        photdir = inlist3[7:] + "_photos"
        logdir = inlist3[7:] + "_logs"

        if not os.path.exists(photdir): os.mkdir(photdir)
        if not os.path.exists(logdir): os.mkdir(logdir)

        if run_sims: run_time = evolve_planet(params, inlist3, coremodel, evolvemodel, run_sims, logdir, photdir, do_use_other_energy)

    if not os.path.exists(evolvemodel) and run_sims:
        print("[do_evolve_planet] MESA failed evolving core model. Check your parameters.")


# main()
# Uncomment for different planet parameters
depth_heating = False

if depth_heating:
    for gamma in [0.01, 0.001, 0.1]:
        for pdepth in [10., 100., 1000.]:
            for metallicity in [0.01, 0.015, 0.02]:
                model_params = [pdepth * constants.bar_cgs, gamma, metallicity]
                main(params_new=model_params, depth_heating=depth_heating)
else:
    for gamma in [0.01, 0.001, 0.1]:
            for metallicity in [0.01, 0.015, 0.02]:
                model_params = [gamma, metallicity]
                main(params_new=model_params)
