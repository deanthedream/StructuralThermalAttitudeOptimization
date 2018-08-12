#This python script reads in the cassini STL and plots it in matplotlib

from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.optimize import fmin
import math, random

class spacecraftModelPoints(object):
    def __init__(self,fname):
        with open(fname, 'rb') as f:#load STLfile
            lines = f.readlines()

        self.CADstructure = list()

        iSolid = -1#index for which solid is being iterated over
        firstIteration = True
        for line in lines:
            if 'solid ' in line:
                if firstIteration == False:
                    self.CADstructure.append({'iSolid':iSolid, 'facets':facetList})#'ifacet':ifacet, 'facetNormal':facetNormal, 'vertices':vertices})
                else:
                    firstIteration = False
                iSolid += 1
                ifacet = -1#reset facet count to 0
            if 'facet normal' in line:
                if ifacet == -1:
                    facetList = list()
                else:#If not first facet
                    facetList.append({'ifacet':ifacet, 'facetNormal':facetNormal, 'vertices':vertices})
                ifacet += 1
                facetNormal = list()
                facetNormal.append([float(item) for item in line.rstrip().split(' ') if isFloat(item)])
                vertices = list()
            if 'vertex' in line:
                vertices.append([float(item) for item in line.rstrip().split(' ') if isFloat(item)])
        self.CADstructure.append({'iSolid':iSolid, 'facets':facetList})

        print 'Number of Components: ' + str(len(self.CADstructure))
        print 'Number of Facets: ' + str(len(self.CADstructure[0]['facets']))

        #Extract all Vertices
        #ASSUMES ONLY 1 SOLID
        self.allVertices = list()
        for i in np.arange(len(self.CADstructure[0]['facets'])):
            if i%50 == 0:#DOWNSELECT NUMBER OF POINTS
                tmpVertices = self.CADstructure[0]['facets'][i]['vertices']
                for K in tmpVertices:
                    self.allVertices.append(np.asarray(K))
        #We now have a stack of all vertices in one list

        #Extract all X, Y, Z
        self.X = list()
        self.Y = list()
        self.Z = list()
        for i in self.allVertices:
            self.X.append(i[0])
            self.Y.append(i[1])
            self.Z.append(i[2])

        #Center and Normalize all X,Y,Z
        meanX = np.mean(self.X)
        meanY = np.mean(self.Y)
        meanZ = np.mean(self.Z)
        self.X -= meanX
        self.Y -= meanY
        self.Z -= meanZ
        rangeX = np.max(self.X) - np.min(self.X)
        rangeY = np.max(self.Y) - np.min(self.Y)
        rangeZ = np.max(self.Z) - np.min(self.Z)
        maxRange = np.max([rangeX,rangeY,rangeZ])
        self.X = self.X/maxRange
        self.Y = self.Y/maxRange
        self.Z = self.Z/maxRange

        #REDEFINE ALL VERTICES to have normalized points
        self.allVertices = list()
        for (tx, ty, tz) in zip(self.X, self.Y, self.Z):
            arr = np.asarray([tx,ty,tz])
            self.allVertices.append(arr)

        #### Pick 100 random points
        n = 100
        rndpts = random.sample(range(1, len(self.allVertices)), n)
        distances = list()
        meanDist = list()
        stdDist = list()
        for ind in np.arange(len(rndpts)):
            rndptInd = rndpts[ind]
            #Calculate distance between point and all other points
            refVect = np.asarray(self.allVertices[rndptInd])
            otherInds = np.arange(len(self.allVertices))
            otherInds = delete(otherInds,rndptInd)
            distPT = list()
            for otherInd in otherInds:
                tmp = np.linalg.norm(np.asarray(self.allVertices[otherInd])-refVect)
                distPT.append(tmp)
            meanDist.append(np.mean(distPT))
            stdDist.append(np.std(distPT))
            distances.append(distPT)
            print 'Done Rnd Pt %d of %d' %(ind, n)
        indOfMinMean = np.argmin(meanDist)
            #hist(distances)
            #show(block=False)
        print 'Done Picking Reference Distance'
        ################################################

        ####
        Rk = 0.5*meanDist[indOfMinMean]#+1.*stdDist[indOfMinMean]#arbitrarily try Rk as meanDist
        pInds = np.arange(len(self.allVertices)).tolist()# An Array of All Indices
        #Pick Random pInd From pInds
        PVInds = [random.sample(pInds,1)[0]]#Add first point to PVInds, these are the indices of points to keep
        sampleInds = [ind for ind in pInds if ind not in PVInds]#All points in pInds and not in PVInds
        #Start Loop
        while len(pInds) > 0:#len(PVInds):#While at least 1 pInds is not in PVInds
            #Randomly pick one ind not in Picked Inds but still in pInds for testing
            testingInd = random.sample(sampleInds,1)[0]
            for ind in pInds:#Iterate over all Inds
                # if ind == testingInd:#Skip if current ind is the testing ind
                #     continue
                # else:
                dist = np.linalg.norm(self.allVertices[ind]-self.allVertices[testingInd])
                if dist < Rk:#If distance between points is less than reference distance
                    pInds.remove(ind)#delete point from pInds
                    #sampleInds.remove(ind)#delete point from sampleInds
                    print "dist: %f, Rk: %f" %(dist,Rk)
                    break
            else:
                #No Break Occured in For Loop So Add to PVInds
                PVInds.append(testingInd)
                sampleInds.remove(testingInd)#delete point from sampleInds
                print 'Appending to PVInds'
            print "len(pInds):%d, len(PVInds):%d" %(len(pInds),len(PVInds))
        #Setup All Points
        self.pts = list()
        for ptInd in PVInds:
            self.pts.append(self.allVertices[ptInd])

        #Re-Extract all X, Y, Z
        self.X = list()
        self.Y = list()
        self.Z = list()
        for i in self.pts:
            self.X.append(i[0])
            self.Y.append(i[1])
            self.Z.append(i[2])

        self.PVInds = PVInds
        self.pInds = pInds
        ####

        # #Setup All Points
        # self.pts = list()
        # self.pts.append(np.asarray([self.X[0],self.Y[0],self.Z[0]]))
        # cnt = 0
        # for (tx, ty, tz) in zip(self.X, self.Y, self.Z):
        #     arr = np.asarray([tx,ty,tz])
        #     distPT = list()
        #     for ind in np.arange(len(self.pts)):
        #         distPT.append(np.linalg.norm(np.asarray(self.pts[ind])-arr))
        #     TMPmeanDist = np.mean(distPT)
        #     if TMPmeanDist > meanDist[indOfMinMean]+1.*stdDist[indOfMinMean]:#The mean calcualted for pt is less than other pts
        #         self.pts.append(arr)
        #     cnt += 1
        #     print "Pt Num: %d of %d, numpts: %d, TMPmeanDist: %f, meanDist: %f, stdDist: %f" %(cnt,len(self.X),len(self.pts),TMPmeanDist,meanDist[indOfMinMean],stdDist[indOfMinMean])
        # #self.pts = np.asarray(self.pts)      


        #Calculate Means and Ranges of Points
        self.meanX = np.mean(self.X)
        self.meanY = np.mean(self.Y)
        self.meanZ = np.mean(self.Z)
        self.rangeX = np.max(self.X) - np.min(self.X)
        self.rangeY = np.max(self.Y) - np.min(self.Y)
        self.rangeZ = np.max(self.Z) - np.min(self.Z)
        self.maxRange = np.max([self.rangeX,self.rangeY,self.rangeZ])

def isFloat(string):
    try:
        float(string)
        return True
    except:
        return False

if __name__ == '__main__':
    close('all')
    fname = 'cassiniSTLFROM3DBUILDER3ascii.STL'
    scMDLpts = spacecraftModelPoints(fname)
    
    fig = plt.figure()
    ax= fig.gca(projection='3d')
    #ax = fig.add_subplot(111, projection='3d')


    ax.scatter(scMDLpts.X,scMDLpts.Y,scMDLpts.Z, c='k', marker='o', s=1)
    show(block=False)

    #### Minimize Central Axis ######################
    def erf(w0,scMDLpts):
        pos = np.asarray(w0[3:6])
        bhat = np.asarray(w0[0:3])/np.linalg.norm(np.asarray(w0[0:3]))
        error = 0.
        #for i in np.arange(len(scMDLpts.pts)):
        #    error += np.linalg.norm(scMDLpts.pts[i]-dot(scMDLpts.pts[i],bhat)*bhat + pos)**2.
        tmpPTS = scMDLpts.pts-pos
        error = np.sum(np.linalg.norm(tmpPTS - outer(np.einsum('ij,j',tmpPTS,bhat),bhat),axis=1)**2.)
        #np.einsum('ij,j',tmp,tmp[0])    
        return error

    def fibonacci_sphere(samples=1,randomize=True):
        rnd = 1.
        if randomize:
            rnd = random.random() * samples

        points = []
        offset = 2./samples
        increment = math.pi * (3. - math.sqrt(5.));

        for i in range(samples):
            y = ((i * offset) - 1) + (offset / 2);
            r = math.sqrt(1 - pow(y,2))

            phi = ((i + rnd) % samples) * increment

            x = math.cos(phi) * r
            z = math.sin(phi) * r

            points.append([x,y,z])

        return points

    points = fibonacci_sphere(samples=200)
    initERF = list()
    w0 = [np.sqrt(2.)/2.,0.,np.sqrt(2.)/2.,scMDLpts.meanX,scMDLpts.meanY,scMDLpts.meanZ]
    for i in np.arange(len(points)):
        w0[0:3] = [points[i][0],points[i][1],points[i][2]]
        f = erf(w0,scMDLpts)
        initERF.append(f)
    indMin = np.argmin(initERF)

    w0[0:3] = [points[indMin][0],points[indMin][1],points[indMin][2]]
    #w0 = [np.sqrt(2.)/2.,0.,np.sqrt(2.)/2.,scMDLpts.meanX,scMDLpts.meanY,scMDLpts.meanZ]
    
    #pos0 = [meanX,meanY,meanZ]
    #WORKS#out = fmin(erf,w0,args=(scMDLpts,),xtol=1e-3,ftol=1e-3,full_output=True)
    from scipy.optimize import minimize
    # from scipy.optimize import LinearConstraint
    # linear_constraint = LinearConstraint(\
    #     [[1.,1.,1.,0.,0.,0.],\
    #     [1.,1.,1.,0.,0.,0.],\
    #     [1.,1.,1.,0.,0.,0.],\
    #     [0.,0.,0.,1.,0.,0.],\
    #     [0.,0.,0.,0.,1.,0.],\
    #     [0.,0.,0.,0.,0.,1.]],\
    #     [1., 1., 1.,-np.inf, -np.inf, -np.inf],\
    #     [1., 1., 1., np.inf, np.inf, np.inf])
    

    # const1 = {'type':'ineq','fun': lambda x: np.linalg.norm(np.asarray(x[0],x[1],x[2]))-1.}
    # const2 = {'type':'ineq','fun': lambda x: -np.linalg.norm(np.asarray(x[0],x[1],x[2]))+1.}
    # out = minimize(erf,w0,args=(scMDLpts,),options={'gtol': 1e-3, 'disp': True},constraints=[const1,const2])#constraints=[{'type':'eq','fun':linear_constraint}])
    from scipy.optimize import NonlinearConstraint
    def myConstraint(x):
        out = np.linalg.norm(np.asarray([x[0],x[1],x[2]]))-1.
        print "%f, %f, %f, out:%f" %(x[0],x[1],x[2],out)
        return out#np.asarray([out,out,out,out,out,out]) 
    lb = 0.
    ub = 0.
    NLCon1 = NonlinearConstraint(myConstraint,lb,ub)
    con1 = {'type':'eq','fun': myConstraint}
    con2 = {'type':'ineq','fun': lambda x: x[0]+1}
    con3 = {'type':'ineq','fun': lambda x: -x[0]+1}
    con4 = {'type':'ineq','fun': lambda x: x[1]+1}
    con5 = {'type':'ineq','fun': lambda x: -x[1]+1}
    con6 = {'type':'ineq','fun': lambda x: x[2]+1}
    con7 = {'type':'ineq','fun': lambda x: -x[2]+1}
    const1 = {'type':'ineq','fun': lambda x: np.linalg.norm(np.asarray(x[0],x[1],x[2]))-1.}
    const2 = {'type':'ineq','fun': lambda x: -np.linalg.norm(np.asarray(x[0],x[1],x[2]))+1.}#constraints=(const1,const2)
    #A = np.asarray([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0]])
    #lb = 
    #LinearConstraint()
    #const1 = {{'type':'eq',}
    out = minimize(erf,w0,args=(scMDLpts,),options={'tol': 1e-3, 'disp': True},constraints=(con1,con2,con3,con4,con5,con6,con7))
    #constraints=[{'type':'eq','fun':linear_constraint}])



    # from scipy.optimize import Bounds
    # bounds = Bounds([0, -0.5], [1.0, 2.0])

    print out
    #################################################


    #### Plot
    bhat = np.asarray([out['x'][0],out['x'][1],out['x'][2]])
    #scMDLpts.maxRange = np.max([scMDLpts.rangeX,scMDLpts.rangeY,scMDLpts.rangeZ])
    pos0 = np.asarray([out['x'][3],out['x'][4],out['x'][5]])
    ax.scatter(pos0[0],pos0[1],pos0[2],color='red')
    ax.plot([-0.5*scMDLpts.maxRange*bhat[0]+pos0[0],0.5*scMDLpts.maxRange*bhat[0]+pos0[0]],[-0.5*scMDLpts.maxRange*bhat[1]+pos0[1],0.5*scMDLpts.maxRange*bhat[1]+pos0[1]],[-0.5*scMDLpts.maxRange*bhat[2]+pos0[2],0.5*scMDLpts.maxRange*bhat[2]+pos0[2]],color='blue')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    show(block=False)
