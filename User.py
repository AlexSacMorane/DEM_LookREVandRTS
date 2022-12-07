# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file where the user can change the different parameters for the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

def All_parameters():
    '''this function is called in main.py to have all the parameters needed in the simulation'''

    #---------------------------------------------------------------------------
    #Geometric parameters

    N_grain_disk = 300 #number of grains
    R_mean = 350 #µm radius to compute the grain distribution. Then recomputed
    L_R = [1.2*R_mean, 1.1*R_mean, 0.9*R_mean, 0.8*R_mean] #from larger to smaller
    L_percentage_R = [1/6, 1/3, 1/3, 1/6] #distribution of the different radius
    #Recompute the mean radius
    R_mean = 0
    for i in range(len(L_R)):
        R_mean = R_mean + L_R[i]*L_percentage_R[i]

    #write dict
    dict_geometry = {
    'N_grain_disk' : N_grain_disk,
    'R_mean' : R_mean,
    'L_R' : L_R,
    'L_percentage_R' : L_percentage_R,
    }

    #---------------------------------------------------------------------------
    #Material parameters

    Y = 70*(10**9)*(10**6)*(10**(-12)) #Young Modulus µN/µm2
    nu = 0.3 #Poisson's ratio
    rho = 2500*10**(-6*3) #density kg/µm3
    mu_friction_gg = 0.5 #grain-grain
    mu_friction_gw = 0 #grain-wall
    coeff_restitution = 0.2 #1 is perfect elastic

    #write dict
    dict_material = {
    'Y' : Y,
    'nu' : nu,
    'rho' : rho,
    'mu_friction_gg' : mu_friction_gg,
    'mu_friction_gw' : mu_friction_gw,
    'coeff_restitution' : coeff_restitution
    }

    #---------------------------------------------------------------------------

    #Box définition
    x_box_min = 0 #µm
    x_box_max = 2*R_mean*math.sqrt(N_grain_disk/0.6) #µm 0.6 from Santamarina, 2014 to avoid boundaries effect
    y_box_min = 0 #µm

    #write dict
    dict_sample = {
    'x_box_min' : x_box_min,
    'x_box_max' : x_box_max,
    'y_box_min' : y_box_min
    }

    #---------------------------------------------------------------------------
    #External sollicitations
    
    #Confinement
    Vertical_Confinement_Linear_Force = Y*2*R_mean/1000 #µN/µm used to compute the Vertical_Confinement_Force
    Vertical_Confinement_Force = Vertical_Confinement_Linear_Force*(x_box_max-x_box_min) #µN
    gravity = 0 #µm/s2
    
    #Dissolution
    frac_dissolved = 0.15 #Percentage of grain dissolved
    frac_Rmean0 = 0.05
    DR_dissolution = frac_Rmean0*R_mean #Reduction of the grain radius at eact iteration

    #write dict
    dict_sollicitations = {
    'Vertical_Confinement_Force' : Vertical_Confinement_Force,
    'gravity' : gravity,
    'frac_dissolved' : frac_dissolved,
    'DR_dissolution' : DR_dissolution
    }

    #---------------------------------------------------------------------------
    #Algorithm parameters

    #DEM parameters
    dt_DEM_crit = math.pi*min(L_R)/(0.16*nu+0.88)*math.sqrt(rho*(2+2*nu)/Y) #s critical time step from O'Sullivan 2011
    dt_DEM = dt_DEM_crit/8 #s time step during DEM simulation
    Spring_type = 'Ponctual' #Kind of contact
    factor_neighborhood = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods = 200 #the frequency of the update of the neighborhood of the grains and the walls
    
    #Stop criteria of the DEM
    i_DEM_stop = 3000 #maximum iteration for one DEM simulation
    Ecin_ratio = 0.0002
    n_window_stop = 50
    dk0_stop = 0.03
    dy_box_max_stop = 0.5

    #PF-DEM
    n_t_PFDEM = 20 #number of cycle PF-DEM

    #Debugging
    i_print_plot = 100 #frenquency of the print and plot (if Debug_DEM) in DEM step
    Debug = True #plot configuration before and after DEM simulation
    Debug_DEM = False #plot configuration inside DEM
    SaveData = True #save simulation
    main_folder_name = 'Data_RTS' #where data are saved
    template_simulation_name = 'f_'+str(int(1000*frac_Rmean0))+'_Run_' #template of the simulation name

    #write dict
    dict_algorithm = {
    'dt_DEM_crit' : dt_DEM_crit,
    'dt_DEM' : dt_DEM,
    'i_update_neighborhoods': i_update_neighborhoods,
    'i_DEM_stop' : i_DEM_stop,
    'Ecin_ratio' : Ecin_ratio,
    'n_window_stop' : n_window_stop,
    'dk0_stop' : dk0_stop,
    'dy_box_max_stop' : dy_box_max_stop,
    'n_t_PFDEM' : n_t_PFDEM,
    'Spring_type' : Spring_type,
    'Debug' : Debug,
    'Debug_DEM' : Debug_DEM,
    'SaveData' : SaveData,
    'main_folder_name' : main_folder_name,
    'template_simulation_name' : template_simulation_name,
    'i_print_plot' : i_print_plot,
    'factor_neighborhood' : factor_neighborhood
    }

    #---------------------------------------------------------------------------
    #Initial condition parameters

    #Generation of grains
    N_test_max = 5000 # maximum number of tries to generate a grain without overlap
    factor_ymax_box = 2.5 #margin to generate grains
    n_generation = 2 #number of grains generation
    #/!\ Work only for 2 /!\
    
    #DEM
    dt_DEM_IC = dt_DEM_crit/5 #s time step during IC
    factor_neighborhood_IC = 1.5 #margin to detect a grain into a neighborhood
    i_update_neighborhoods_gen = 5 #the frequency of the update of the neighborhood of the grains and the walls during IC generations
    i_update_neighborhoods_com = 100 #the frequency of the update of the neighborhood of the grains and the walls during IC combination
    #Stop DEM
    i_DEM_stop_IC = 3000 #stop criteria for DEM during IC
    Ecin_ratio_IC = 0.0005
    
    #Plot
    Debug_DEM_IC = False #plot configuration inside DEM during IC
    i_print_plot_IC = 100 #frequency of the print and plot (if Debug_DEM_IC) for IC
    
    #write dict
    dict_ic = {
    'n_generation' : n_generation,
    'i_update_neighborhoods_gen': i_update_neighborhoods_gen,
    'i_update_neighborhoods_com': i_update_neighborhoods_com,
    'factor_ymax_box' : factor_ymax_box,
    'i_DEM_stop_IC' : i_DEM_stop_IC,
    'Debug_DEM' : Debug_DEM_IC,
    'dt_DEM_IC' : dt_DEM_IC,
    'Ecin_ratio_IC' : Ecin_ratio_IC,
    'i_print_plot_IC' : i_print_plot_IC,
    'factor_neighborhood_IC' : factor_neighborhood_IC,
    'N_test_max' : N_test_max
    }

    #---------------------------------------------------------------------------

    return dict_algorithm, dict_geometry, dict_ic, dict_material, dict_sample, dict_sollicitations

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm):
    '''Criteria to stop simulation (PF and DEM)'''
    Criteria_Verified = False
    if dict_algorithm['i_PF'] >= dict_algorithm['n_t_PFDEM']:
        Criteria_Verified = True
    return Criteria_Verified
