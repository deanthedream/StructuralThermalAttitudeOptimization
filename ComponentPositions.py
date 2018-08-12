#This python script plots 

from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import json
import numpy as np


class ComponentPositions(object):
    def __init__(self):
        #all distances in units of meters
        #Fairing Dimensions
        self.FairingR = 5.
        self.FairingH = 12.

        self.InertialFrameAxes = {'x':[1.,0.,0.],'y':[0.,1.,0.],'z':[0.,0.,1.]} #the frame from which everything is defined
        self.InertialFrameOrigin = np.asarray([0.,0.,0.])
        self.r_SCbody_I = np.asarray([0.,0.,0.])#the position of the spacecraft body frame 
            #w.r.t the inertial frame #just some arbitrary definition for now
        self.q_SCbody_I = self.eulerAnglesToQuaternion(0.,0.,0.)#spacecraft body frame orientation
            #w.r.t inertial frame #just some arbitrary instantiation

        #TODO SPACECRAFT PARTS NEED TO BE SEPARATE OBJECTS THAT ARE INSTANTIATED


        #### Read in the spacecraft parts ####################
        fname = 'parts.json'
        with open(fname, 'rb') as f:#load STLfile
            self.comps = json.load(f)#f.readlines()
        ######################################################

        #### iterate over each part ##########################
        for part in self.comps.keys():
            #check if part DOES NOT HAVE body_vector assigned to it
            if not self.comps[part].has_key('part_position'):
                #assign body vector
                self.comps[part]['part_position'] = self.randomPartPositionVector()

            #check if part has part_body_frame
            if not self.comps[part].has_key('part_body_frame'):
                self.comps[part]['part_body_frame'] = self.randomPartBodyFrame()

            #check if part overlaps with any other part

        ######################################################

    def randomPartPositionVector(self):
        """
        return:
            part_position - vector describing the part location of a component
                has form [x,y,z] 
        """
        r = np.random.uniform(low=0.,high=self.FairingR)
        theta = np.random.uniform(low=0.,high=2.*np.pi)
        h = np.random.uniform(low=0.,high=self.FairingH)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        body_vector = np.array([x,y,h]) + self.r_SCbody_I#this should produce a spacecraft body_vector position that updates if the spacecraft body updates
        return body_vector

    def eulerAnglesToQuaternion(self, pitch, roll, yaw):
        """
        Args:
            roll - 
            pitch - 
            yaw - 
        return:
            quaternion - 
        """
        cy = np.cos(yaw * 0.5);
        sy = np.sin(yaw * 0.5);
        cr = np.cos(roll * 0.5);
        sr = np.sin(roll * 0.5);
        cp = np.cos(pitch * 0.5);
        sp = np.sin(pitch * 0.5);

        q = np.asarray([0,0,0,0])
        q[0] = cy * cr * cp + sy * sr * sp #w
        q[1] = cy * sr * cp - sy * cr * sp #x
        q[2] = cy * cr * sp + sy * sr * cp #y
        q[3] = sy * cr * cp - cy * sr * sp #z
        return q

    def quaterionToEulerAngles(self, q):
        """
        Args:
            q - 
        return:
            roll, pitch, yaw
        """
        roll = np.atan2(2*(q[0]*q[1] + q[2]*q[3]),1-2*(q[1]**2 + q[2]**2))#phi
        pitch = np.arcsin(2*(q[0]*q[1] - q[3]*q[1]))#theta
        yaw = np.atan2(2*(q[0]*q[3] + q[1]*q[2]),1-2*(q[2]**2 + q[3]**2))#psi
        EulerAngles = [pitch, roll, yaw]
        return EulerAngles

    def randomPartBodyFrame(self):
        """
        return:
            part_body_frame - the body frame of the part

        """
        #https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
        #Generate Euler Angles
        yaw = np.random.uniform(low=0,high=2*np.pi)
        roll = np.random.uniform(low=0,high=2*np.pi)
        pitch = np.random.uniform(low=0,high=2*np.pi)
        q = self.eulerAnglesToQuaternion(pitch,roll,yaw)
        part_body_frame = q + self.q_SCbody_I
        return part_body_frame


