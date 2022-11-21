# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to define a new class.
The new class is about the grains
"""

#-------------------------------------------------------------------------------
#Libs
#-------------------------------------------------------------------------------

import numpy as np
import math
import random

#-------------------------------------------------------------------------------
#Class
#-------------------------------------------------------------------------------

class Grain:

#-------------------------------------------------------------------------------

  def __init__(self, ID, Center, Radius, dict_material):
    '''defining the grain

    each grain is described ...
    '''
    #Id of the grain
    self.id = ID
    #Material property
    self.dissolvable = False
    self.y = dict_material['Y']
    self.nu = dict_material['nu']
    self.g = dict_material['Y']/2/(1+dict_material['nu']) #shear modulus
    self.rho = dict_material['rho']
    #other characteristic
    self.radius = Radius
    self.surface = math.pi*self.radius**2
    self.rho_surf = 4/3*self.rho*self.radius
    self.mass = self.rho_surf*self.surface
    self.inertia = self.mass*self.radius**2
    #kinematic
    self.v = np.array([0,0])
    self.a = np.array([0,0])
    self.theta = 0 #no rotation assumed
    self.w = 0 #dtheta/dt
    #position
    self.center = Center
    self.plot_preparation()

#-------------------------------------------------------------------------------

  def update_geometry_kinetic(self, V, A, DT):
    '''update the acceleration and the velocity of a grain
    update geometrical parameters as border and center nodes'''
    #translation
    self.v = V
    self.a = A
    for i in range(len(self.l_border)):
        self.l_border[i] = self.l_border[i] + self.v*DT
        self.l_border_x[i] = self.l_border_x[i] + self.v[0]*DT
        self.l_border_y[i] = self.l_border_y[i] + self.v[1]*DT
    self.center = self.center + self.v*DT

#-------------------------------------------------------------------------------

  def init_f_control(self,dict_sollicitations):
      '''initialize the force applied to the grain
      a gravity of g is applied'''
      #-------------------------------------------------------------------------
      #load data needed
      g = dict_sollicitations['gravity']
      #-------------------------------------------------------------------------

      self.fx = 0
      self.fy = -g*self.mass
      self.f = np.array([self.fx,self.fy])
      self.mz = 0

#-------------------------------------------------------------------------------

  def update_f(self, Fx, Fy):
    '''add a force (an array [Fx,Fy]) to the grain'''
    self.fx = self.fx + Fx
    self.fy = self.fy + Fy
    self.f = np.array([self.fx,self.fy])

#-------------------------------------------------------------------------------

  def Update_characteristic(self):
      '''to write...'''
      self.surface = math.pi*self.radius**2
      self.rho_surf = 4/3*self.rho*self.radius
      self.mass = self.rho_surf*self.surface
      self.inertia = self.mass*self.radius**2

#-------------------------------------------------------------------------------

  def plot_preparation(self):
    '''prepare the plot of a grain'''
    N_p_border = 180 #number of nodes
    L_border_x = []
    L_border_y = []
    L_border = []
    for i in range(N_p_border):
       theta = 2*math.pi*i/N_p_border
       L_border_x.append(self.center[0]+self.radius*math.cos(theta))
       L_border_y.append(self.center[1]+self.radius*math.sin(theta))
       L_border.append(np.array([self.center[0]+self.radius*math.cos(theta),self.center[1]+self.radius*math.sin(theta)]))
    L_border_x.append(L_border_x[0])
    L_border_y.append(L_border_y[0])
    L_border.append(L_border[0])
    self.l_border_x = L_border_x
    self.l_border_y = L_border_y
    self.l_border = L_border

#-------------------------------------------------------------------------------
#Function
#-------------------------------------------------------------------------------
