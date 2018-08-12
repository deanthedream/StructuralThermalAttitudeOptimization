#run ComponentPositions

from component import component
from ComponentPositions import ComponentPositions as CP
import numpy as np

from sympy import *
import quaternion

SCBody = CoordSysCartesian('SCBody')
SCBody.origin
comp1 = component(SCBody,'box1')#=SCBody)#compname='box1',
comp1.r_component_Body += 1.*SCBody.i + 1.*SCBody.j + 1.*SCBody.k
#comp1.r_component_Body += [1.,1.,1.]#np.asarray([1.,1.,1.])
comp2 = component(SCBody,'box2')
comp2.r_component_Body += comp1.r_component_Body

print comp1.vertices
print comp2.vertices
comp1.r_component_Body += 1.*SCBody.i + 1.*SCBody.j + 1.*SCBody.k
print comp1.vertices
print comp2.vertices

# cp = CP()

# print cp.comps['box1']['part_position']
# cp.r_SCbody_I += np.asarray([1,1,1])
# print cp.comps['box1']['part_position']

