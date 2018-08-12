#This python script contains the construct for a component
import numpy as np
from sympy import *
#from sympy.vector import CoordSysCartesian
import quaternion

class component(object):
    def __init__(self, compname, mass=0.,shape='box',dims={"l":0.,"w":0.,"h":0.},\
        PelectricalIn=[0.,0.,0.,0.,0.],MaxTemp=[1.,1.,1.,1.,1.],MinTemp=[0.,0.,0.,0.,0.],\
        specific_heat_capacity=0.,emissivity=0.,absorptivity=0.):#SCBody, 

        self.compname = 'box1'#args.get('compname')#compname
        
        self.componentFrame = np.asarray([0.,0.,0.,])#self.componentFrame = CoordSysCartesian(self.compname)
        #C = self.componentFrame
        self.mass = mass
        self.shape = shape
        self.dims = dims
        self.PelectricalIn = PelectricalIn
        self.MaxTemp = MaxTemp
        self.MinTemp = MinTemp
        self.specific_heat_capacity = specific_heat_capacity
        self.emissivity = emissivity
        self.absorptivity = absorptivity

        #What shape is the component
        if self.shape == 'box':
            self.l = dims['l']
            self.w = dims['w']
            self.h = dims['h']

            #self.r_component_Body = SCBody.origin.locate_new('Body',0.*SCBody.i + 0.*SCBody.j + 0.*SCBody.k)#np.asarray([0.,0.,0.])
            self.vertices = list()
            self.vertices.append(np.asarray([0.,    0.,    0.]))# C.origin.locate_new('A',0.*C.i     + 0.*C.j        + 0.*C.k))
            self.vertices.append(np.asarray([self.l,0.,    0.]))#C.origin.locate_new('B',self.l*C.i + 0.*C.j        + 0.*C.k))
            self.vertices.append(np.asarray([0.,    self.w,0.]))#C.origin.locate_new('C',0.*C.i     + self.w*C.j    + 0.*C.k))
            self.vertices.append(np.asarray([0.,    0.,    self.h]))#C.origin.locate_new('D',0.*C.i     + 0.*C.j        + self.h*C.k))
            self.vertices.append(np.asarray([self.l,self.w,0.]))#C.origin.locate_new('E',self.l*C.i + self.w*C.j    + 0.*C.k))
            self.vertices.append(np.asarray([self.l,0.,    self.h]))#C.origin.locate_new('F',self.l*C.i + 0.*C.j        + self.h*C.k))
            self.vertices.append(np.asarray([0.,    self.w,self.h]))#C.origin.locate_new('G',0.*C.i     + self.w*C.j    + self.h*C.k))
            self.vertices.append(np.asarray([self.l,self.w,self.h]))#C.origin.locate_new('H',self.l*C.i + self.w*C.j    + self.h*C.k))
        elif self.shape == 'cylinder':
            self.h = dims['h']
            self.r = dims['r']
            if dims.haskey('theta'):#Angular size of cylinder, i.e. 2pi is a half cylinder
                self.theta = dims['theta']
            else:
                self.theta = 2*np.pi

            self.vertices = list()
            self.vertices.append(np.asarray([0.,    0.,    0.]))
            self.vertices.append(np.asarray([0.,    0.,    self.h]))
        elif self.shape == 'cone':
            self.h = dims['h']
            self.r = dims['r']
            if dims.haskey('theta'):#Angular size of cylinder, i.e. 2pi is a half cylinder
                self.theta = dims['theta']
            else:
                self.theta = 2*np.pi

            self.vertices = list()
            self.vertices.append(np.asarray([0.,    0.,    0.]))
            self.vertices.append(np.asarray([0.,    0.,    self.h]))

        for vertex in self.vertices:
            vertex = vertex + self.r_component_Body

        self.q_component_Body = np.quaternion(0.,0.,0.,0.)

