# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between two grains
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math
import random
import Grain

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G1, G2, dict_material):
    '''
    defining the contact grain-grain
    ...
    '''
    #Id of the contact
    self.id = ID
    self.g1 = G1
    self.g2 = G2
    #Material properties
    self.mu = dict_material['mu_friction_gg']
    self.coeff_restitution = dict_material['coeff_restitution']
    Y_eq = 1/((1-G1.nu*G1.nu)/G1.y+(1-G2.nu*G2.nu)/G2.y)
    R_eq = 1/(1/G1.radius+1/G2.radius)
    k = 4/3*Y_eq*math.sqrt(R_eq) #unlinear spring
    self.k = k
    gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
    mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
    eta = 2 * gamma * math.sqrt(mass_eq*k)
    self.eta = eta
    #Initialisation of the tangetial behavior
    self.tangential_old_statut = False
    self.overlap_tangential = 0
    self.ft = 0

#-------------------------------------------------------------------------------

  def init_contact(self,L_g):
    '''
    initialize the contact with updating the grains,
                                putting at 0 the tangential reaction
                                saying the boolean at False (new contact grain-grain)
    '''
    self.g1 = L_g[self.g1.id]
    self.g2 = L_g[self.g2.id]
    self.ft = 0
    self.tangential_old_statut = False
    self.overlap_tangential = 0

#-------------------------------------------------------------------------------

  def normal(self):
    '''compute the normal reaction of the contact

    ...'''
    #unlinear stiffness
    n12 = (self.g2.center-self.g1.center)/np.linalg.norm(self.g2.center-self.g1.center)
    self.n12 = n12
    overlap = self.g1.radius + self.g2.radius - np.linalg.norm(self.g1.center - self.g2.center)
    self.overlap = overlap
    F12_n = self.k*overlap**(3/2)
    F12 = F12_n*n12
    self.F12_n = F12_n
    self.g1.update_f(-F12[0],-F12[1])
    self.g2.update_f( F12[0], F12[1])
    #damping
    F12_damp = -np.dot(self.g2.v - self.g1.v,n12)*self.eta*n12
    self.F12_damp_n = np.linalg.norm(F12_damp)
    self.g1.update_f(-F12_damp[0],-F12_damp[1])
    self.g2.update_f( F12_damp[0], F12_damp[1])

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
    '''compute the tangential reaction of the contact
    ...
    '''
    if self.mu > 0 :
        #unlinear stiffness
        G_eq = 1/((1-self.g1.nu)/self.g1.g+(1-self.g2.nu)/self.g2.g)
        R_eq = 1/(1/self.g1.radius+1/self.g2.radius)
        kt0 = 8 * G_eq *math.sqrt(R_eq*abs(self.overlap))
        kt = kt0*math.sqrt(max(1-2/3*kt0*abs(self.overlap_tangential)/self.mu/abs(self.F12_n),0))
        self.kt = kt

        t12 = np.array([-self.n12[1], self.n12[0]])
        self.t12 = t12
        if self.tangential_old_statut:
            #if a reaction has been already computed
            #need to project the tangential reaction on the new tangential plane
            self.ft = self.ft*np.dot(self.t12_old,self.t12)
        else:
            self.tangential_old_statut = True
        Delta_Us = np.dot(self.g1.v-self.g2.v,self.t12) * dt_DEM
        self.overlap_tangential = self.overlap_tangential + Delta_Us
        self.ft = self.ft - self.kt*Delta_Us
        self.t12_old = self.t12
        if abs(self.ft) > abs(self.mu*self.F12_n) or kt == 0: #Coulomb criteria
            self.ft = self.mu * abs(self.F12_n) * np.sign(self.ft)
        F12 = -self.ft*t12
        self.g1.update_f(-F12[0],-F12[1])
        self.g2.update_f( F12[0], F12[1])
        #damping
        gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
        mass_eq = self.g1.mass*self.g2.mass/(self.g1.mass+self.g2.mass)
        eta = 2 * gamma * math.sqrt(mass_eq*kt)
        F12_damp = -np.dot(self.g2.v - self.g1.v,t12)*eta/2*t12
        self.f12_damp_t = np.linalg.norm(F12_damp)
        self.g1.update_f(-F12_damp[0],-F12_damp[1])
        self.g2.update_f( F12_damp[0], F12_damp[1])
    else :
        self.kt = 0
        self.t12 = np.array([-self.n12[1], self.n12[0]])

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Update_Neighborhoods(dict_algorithm,dict_sample):
    '''
    determine a neighborhoods for each grain. This function is called every x time step
    grain contact is determined by Grains_Polyhedral_contact_Neighborhoods

    notice that if there is a potential contact between grain_i and grain_j
    grain_i is not in the neighborhood of grain_j
    whereas grain_j is in the neighborhood of grain_i
    with i_grain < j_grain
    '''
    for i_grain in range(len(dict_sample['L_g'])-1) :
        neighborhood = []
        for j_grain in range(i_grain+1,len(dict_sample['L_g'])):
            if np.linalg.norm(dict_sample['L_g'][i_grain].center-dict_sample['L_g'][j_grain].center) < dict_algorithm['factor_neighborhood']*(dict_sample['L_g'][i_grain].radius+dict_sample['L_g'][j_grain].radius):
                neighborhood.append(dict_sample['L_g'][j_grain])
        dict_sample['L_g'][i_grain].neighbourood = neighborhood

#-------------------------------------------------------------------------------

def Grains_Disks_contact_Neighborhoods(dict_material,dict_sample):
    '''
    detect contact between a grain and grains from its neighborhood
    the neighborhood is updated with Update_Neighborhoods_f()
    '''
    for i_grain in range(len(dict_sample['L_g'])-1) :
        for neighbor in dict_sample['L_g'][i_grain].neighbourood:
            j_grain = neighbor.id
            if Intersection(dict_sample['L_g'][i_grain],dict_sample['L_g'][j_grain]):
                if (i_grain,j_grain) not in dict_sample['L_ij_contact']:  #contact not detected previously
                   #creation of contact
                   dict_sample['L_ij_contact'].append((i_grain,j_grain))
                   dict_sample['L_contact'].append(Contact(dict_sample['id_contact'], dict_sample['L_g'][i_grain], dict_sample['L_g'][j_grain], dict_material))
                   dict_sample['id_contact'] = dict_sample['id_contact'] + 1

            else :
                if (i_grain,j_grain) in dict_sample['L_ij_contact'] : #contact detected previously is not anymore
                       dict_sample['L_contact'].pop(dict_sample['L_ij_contact'].index((i_grain,j_grain)))
                       dict_sample['L_ij_contact'].remove((i_grain,j_grain))

#-------------------------------------------------------------------------------

def Intersection(g1,g2):
    '''verify if there is an overlap between two grains or not'''
    d_12 = np.linalg.norm(g1.center - g2.center)
    return d_12 < g1.radius + g2.radius
