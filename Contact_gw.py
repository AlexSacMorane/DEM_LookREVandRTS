# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the contact between a wall and a grain
"""

#-------------------------------------------------------------------------------
#Librairies
#-------------------------------------------------------------------------------

import numpy as np
import math

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Contact_gw:

#-------------------------------------------------------------------------------

  def __init__(self, ID, G, dict_material, Nature, Limit):
    '''defining the contact grain-wall
    ...
    '''
    #Id of the contact
    self.id = ID
    self.nature = Nature
    self.g = G
    self.limit = Limit
    #Material properties
    self.mu = dict_material['mu_friction_gw']
    self.coeff_restitution = dict_material['coeff_restitution']
    factor = 5 #factor just to increase the stiffness
    k = factor*4/3*G.y/(1-G.nu*G.nu)*math.sqrt(G.radius) #unlinear spring
    self.k = k
    self.kt = 0
    gamma = -math.log(self.coeff_restitution)/math.sqrt(math.pi**2+math.log(self.coeff_restitution)**2)
    mass_eq = self.g.mass
    eta = 2 * gamma * math.sqrt(mass_eq*k)
    self.eta = eta
    #Initialisation of the tangetial behavior
    self.tangential_old_statut = False
    self.overlap_tangential = 0
    self.ft = 0

#-------------------------------------------------------------------------------

  def init_contact_gw(self,L_g):
    '''
    initialize the contact with updating the grain,
                                putting at 0 the tangential reaction
                                saying the boolean at False (new contact grain-wall)
    '''
    self.g = L_g[self.g.id]
    self.ft = 0
    self.tangential_old_statut = False

#-------------------------------------------------------------------------------

  def update_overlap(self,new_overlap):
    '''
    update the overlap of a contact already created.
    '''
    self.overlap = new_overlap

#-------------------------------------------------------------------------------

  def normal(self):
    '''compute the normal reaction of the contact

    one "elif" by wall'''
    if self.nature == 'gwy_min':
        #unlinear stiffness
        nwg = np.array([0,1])
        self.nwg = nwg
        overlap = self.limit-(self.g.center[1]-self.g.radius)
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1])
        #damping
        Fwg_damp = -np.dot(self.g.v,nwg)*self.eta*nwg
        self.Fwg_damp_n = np.linalg.norm(Fwg_damp)
        self.g.update_f(Fwg_damp[0],Fwg_damp[1])

    elif self.nature == 'gwy_max':
        #unlinear stiffness
        nwg = np.array([0,-1])
        self.nwg = nwg
        overlap = (self.g.center[1]+self.g.radius)-self.limit
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = np.linalg.norm(Fwg_n)
        self.g.update_f(Fwg[0],Fwg[1])
        #damping
        self.Fwg_damp_n = 0

    elif self.nature == 'gwx_min':
        #unlinear stiffness
        nwg = np.array([1,0])
        self.nwg = nwg
        overlap = self.limit-(self.g.center[0]-self.g.radius)
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1])
        #damping
        Fwg_damp = -np.dot(self.g.v,nwg)*self.eta*nwg
        self.Fwg_damp_n = np.linalg.norm(Fwg_damp)
        self.g.update_f(Fwg_damp[0],Fwg_damp[1])

    elif self.nature == 'gwx_max':
        #linear stiffness
        nwg = np.array([-1,0])
        self.nwg = nwg
        overlap = self.g.center[0]+self.g.radius-self.limit
        self.overlap = overlap
        Fwg_n = self.k*overlap**(3/2)
        Fwg = Fwg_n*nwg
        self.Fwg_n = Fwg_n
        self.g.update_f(Fwg[0],Fwg[1])
        #damping
        Fwg_damp = -np.dot(self.g.v,nwg)*self.eta*nwg
        self.Fwg_damp_n = np.linalg.norm(Fwg_damp)
        self.g.update_f(Fwg_damp[0],Fwg_damp[1])

#-------------------------------------------------------------------------------

  def tangential(self, dt_DEM):
    '''compute the tangential reaction of the contact

    one "elif" by wall'''
    if self.nature == 'gwy_min':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([-1, 0])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.update_f(Fwg[0],Fwg[1])
        else :
            twg = np.array([-1, 0])
            self.twg = twg

    elif self.nature == 'gwy_max':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([1, 0])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.update_f(Fwg[0],Fwg[1])
        else :
            twg = np.array([1, 0])
            self.twg = twg

    elif self.nature == 'gwx_min':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([0, 1])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.update_f(Fwg[0],Fwg[1])
        else :
            twg = np.array([0, 1])
            self.twg = twg

    elif self.nature == 'gwx_max':
        if self.mu > 0 :
            #unlinear stiffness
            twg = np.array([0, -1])
            self.twg = twg
            Delta_Us = np.dot(self.g.v,self.twg) * dt_DEM
            self.overlap_tangential = self.overlap_tangential + Delta_Us
            self.ft = self.ft - self.kt*Delta_Us
            if abs(self.ft) > abs(self.mu*self.Fwg_n) : #Coulomb criteria
                self.ft = self.mu * abs(self.Fwg_n) * np.sign(self.ft)
            Fwg = self.ft*twg
            self.g.update_f(Fwg[0],Fwg[1])
        else :
            twg = np.array([0, 1])
            self.twg = twg

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def Update_wall_Neighborhoods(dict_algorithm, dict_sample):
    '''
    determine a neighborhoods for wall. This function is called every x time step

    grain_wall contact is determined by Grains_Polyhedral_Wall_contact_Neighborhood
    '''
    #factor determines the size of the neighborhood window
    wall_neighborhood = []
    for grain in dict_sample['L_g']:

        p_x_min = grain.center[0]-grain.radius
        p_x_max = grain.center[0]+grain.radius
        p_y_min = grain.center[1]-grain.radius
        p_y_max = grain.center[1]+grain.radius

        #grain-wall x_min
        if abs(p_x_min-dict_sample['x_box_min']) < dict_algorithm['factor_neighborhood']*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall x_max
        if abs(p_x_max-dict_sample['x_box_max']) < dict_algorithm['factor_neighborhood']*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall y_min
        if abs(p_y_min-dict_sample['y_box_min']) < dict_algorithm['factor_neighborhood']*grain.radius :
            wall_neighborhood.append(grain)
        #grain-wall y_max
        if abs(p_y_max-dict_sample['y_box_max']) < dict_algorithm['factor_neighborhood']*grain.radius :
            wall_neighborhood.append(grain)

    #add element in dict
    dict_sample['wall_neighborhood'] = wall_neighborhood

#-------------------------------------------------------------------------------

def Grains_Disk_Wall_contact_Neighborhood(dict_material,dict_sample):
  '''
  detect contact grain in the neighborhood of the wall and  the wall

  the neighborhood is updated with Update_wall_Neighborhoods()
  we realize iterations on the grain list and compare with the coordinate of the different walls
  '''
  for grain in dict_sample['wall_neighborhood']:

      # contact grain-wall x_min
      if grain.center[0] < dict_sample['x_box_min'] + grain.radius and (grain.id,-1) not in dict_sample['L_ij_contact_gw']:
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwx_min', dict_sample['x_box_min']))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
          dict_sample['L_ij_contact_gw'].append((grain.id,-1))
      elif grain.center[0] > dict_sample['x_box_min'] + grain.radius and (grain.id,-1) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-1))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
      # contact grain-wall x_max
      if grain.center[0] > dict_sample['x_box_max'] - grain.radius and (grain.id,-2) not in dict_sample['L_ij_contact_gw']:
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwx_max', dict_sample['x_box_max']))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
          dict_sample['L_ij_contact_gw'].append((grain.id,-2))
      elif grain.center[0] < dict_sample['x_box_max'] - grain.radius and (grain.id,-2) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-2))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
      # contact grain-wall y_min
      if grain.center[1] < dict_sample['y_box_min'] + grain.radius and (grain.id,-3) not in dict_sample['L_ij_contact_gw']:
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwy_min', dict_sample['y_box_min']))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
          dict_sample['L_ij_contact_gw'].append((grain.id,-3))
      elif grain.center[1] > dict_sample['y_box_min'] + grain.radius and (grain.id,-3) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-3))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
      # contact grain-wall y_max
      if grain.center[1] > dict_sample['y_box_max'] - grain.radius and (grain.id,-4) not in dict_sample['L_ij_contact_gw']:
          dict_sample['L_contact_gw'].append(Contact_gw(dict_sample['id_contact_gw'], grain, dict_material, 'gwy_max', dict_sample['y_box_max']))
          dict_sample['id_contact_gw'] = dict_sample['id_contact_gw'] + 1
          dict_sample['L_ij_contact_gw'].append((grain.id,-4))
      elif grain.center[1] < dict_sample['y_box_max'] - grain.radius and (grain.id,-4) in dict_sample['L_ij_contact_gw']:
          i_contact = dict_sample['L_ij_contact_gw'].index((grain.id,-4))
          dict_sample['L_contact_gw'].pop(i_contact)
          dict_sample['L_ij_contact_gw'].pop(i_contact)
