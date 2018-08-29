#runFitCirclesToPoints2.py
import numpy as np
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn import datasets
import heapq
import math

def udot(vect1,vect2):
    #Calculates unit vectors before doing dot product
    vect1 = vect1/np.linalg.norm(vect1)
    vect2 = vect2/np.linalg.norm(vect2)
    dp = np.dot(vect1,vect2)
    if dp > 1.0:
        dp = 1.0
    elif dp < -1.0:
        dp = -1.0
    return dp

close('all')
#### Plot Filtered Model ###################
fig = plt.figure()
ax= fig.gca(projection='3d')
#ax = fig.add_subplot(111, projection='3d')
ax.scatter(scMDLpts.X,scMDLpts.Y,scMDLpts.Z, c='k', marker='o', s=1)
#############################################

#### Select specific Point For Testing Circle Fit##############
xInds = np.where(np.asarray(scMDLpts.X) <-0.05)[0]
zInds = np.where(np.asarray(scMDLpts.Z) >0.10)[0]
selectInds = np.intersect1d(xInds,zInds)

# print xInds
# print zInds
# print selectInds

ax.scatter(np.asarray(scMDLpts.X)[selectInds],np.asarray(scMDLpts.Y)[selectInds],np.asarray(scMDLpts.Z)[selectInds], c='r', marker='x', s=5)
show(block=False)
###############################################################


#### FIND Closest 6 or less points to each point ##############
try:
    cPS[0]
    print 'Using previously defined cPS'
except NameError:
    ptInds = np.arange(len(scMDLpts.pts))#create array of points
    cPS = {}#Closest Points Struct
    for i in ptInds:#Instantiate struct of closest points for each ind
        cPS[i] = [i]

    k=8 #number of points to find#10 is horrible, 7 is ok
    for i in ptInds:
        if len(cPS[i]) >= k:#Only looking for k closest points.... might want to rethink this Maybe calculate how many based on filtering distance and the minimum spatial distance between points on the unit sphere
            continue
        tmp = map(lambda x: np.linalg.norm(x - scMDLpts.pts[i]), scMDLpts.pts)#Calculate distance between points and all other points
        #numPts = k - len(cPS[i]) + 1#num of points to find
        TMP1 = np.argsort(tmp)
        ptsToAdd = np.delete(TMP1,np.argwhere(TMP1==i))#np.asarray(TMPptsToAdd)#heapq.nsmallest(numPts,np.nditer(tmp),key=tmp.__getitem__)#heapq.nsmallest(numPts, tmp)
        if len(cPS[i]) == 1:
            if cPS[i][0] == i:#If cPS has only 1 item and that item is itself
                numPts = k - len(cPS[i]) - 1 + 1#num of points to find
                cPS[i] = ptsToAdd[0:numPts].tolist()#Reassign points
            else:#Otherwise append all of them
                numPts = k - len(cPS[i]) + 1#num of points to find
                cPS[i] += ptsToAdd[0:numPts].tolist()#data format might be messy
        else:#Otherwise append all of them
            numPts = k - len(cPS[i]) + 1#num of points to find
            cPS[i] += ptsToAdd[0:numPts].tolist()#data format might be messy

        # numPts = k - len(cPS[i]) + 1#num of points to find
        # for indG in ptsToAdd[1:numPts]:#Add current point to other points
        #     if cPS[indG] == [indG]:#If cPS has only 1 item and that item is itself
        #         cPS[indG] = [i]
        #     else:
        #         cPS[indG].append(i)
    print 'Done Generating cPS'
#################################################################

#### Calc Vectors to Nearest Points #############################
try:
    cPSvectors[0]
    print 'Using previously defined cPSvectors'
except:
    cPSvectors = {}#Closest Points Struct
    for i in ptInds:#Instantiate struct of closest points for each ind
        cPSvectors[i] = [i]

    for i in ptInds:
        pt = scMDLpts.pts[i]#x,y,z of current points
        ptVectors = list()#will contain vectors to app adjoining points
        for j in cPS[i]:#Inds of all adjoining points
            tmp = scMDLpts.pts[j] - pt
            if math.isnan(tmp[0]):
                print saltyburrito
            if math.isnan((tmp/np.linalg.norm(tmp))[0]):
                print saltyburrito
            ptVectors.append(tmp/np.linalg.norm(tmp))
        cPSvectors[i] = ptVectors
    print 'Done Calculating cPSvectors'
#################################################################


#### Create Matrix of Dot Products ##########################################
connList = list()
THTOL = 1.*np.pi/180.
for iInd in ptInds:
    for tmpjjInd1 in np.arange(len(cPS[iInd])):#iterate over connecting points
        jInd1 = cPS[iInd][tmpjjInd1]
        if iInd == jInd1:
            continue
        vect1 = cPSvectors[iInd][tmpjjInd1]
        for tmpjjInd2 in np.arange(len(cPS[jInd1])):#iterate over connecting points
            jInd2 = cPS[jInd1][tmpjjInd2]
            if jInd2 == jInd1 or jInd2 == iInd:
                continue
            vect2 = cPSvectors[jInd1][tmpjjInd2]
            th1 = np.pi - np.arccos(udot(vect1,vect2))
            phi1 = np.arccos(udot(vect1,vect2))
            for tmpjjInd3 in np.arange(len(cPS[jInd2])):#iterate over connecting points
                jInd3 = cPS[jInd2][tmpjjInd3]
                if jInd3 == jInd2 or jInd3 == iInd or jInd3 == jInd2:
                    continue
                vect3 = cPSvectors[jInd2][tmpjjInd3]
                th2 = np.pi - np.arccos(udot(vect2,vect3))
                phi2 = np.arccos(udot(vect2,vect3))
                if np.abs(th1-th2) > THTOL or th2 < 0.75*np.pi:
                    continue
                else:
                    for tmpjjInd4 in np.arange(len(cPS[jInd3])):#iterate over connecting points
                        jInd4 = cPS[jInd3][tmpjjInd4]
                        if jInd4 == jInd2 or jInd4 == iInd or jInd4 == jInd2 or jInd4 == jInd3:
                            continue
                        vect4 = cPSvectors[jInd3][tmpjjInd4]
                        th3 = np.pi - np.arccos(udot(vect3,vect4))
                        phi3 = np.arccos(udot(vect3,vect4))
                        if np.abs(th2-th3) > THTOL or th3 < 0.75*np.pi:
                            continue
                        else:
                            connList.append([iInd,jInd1,jInd2,jInd3,jInd4])
    if iInd%100 == 0:
        print(str(float(iInd)/len(ptInds)))
print(len(connList))
#### Plot on Figure
rotColor = ['cyan','red','blue','green','orange','purple','yellow']
cInd = 0
for item in connList:
    if len(item)>1:
        ind0 = np.arange(len(item)-1)
        ind1 = np.arange(len(item)-1)+1
        for i in np.arange(len(item)-1):
            ax.plot([scMDLpts.X[item[ind0[i]]],scMDLpts.X[item[ind1[i]]]],\
                [scMDLpts.Y[item[ind0[i]]],scMDLpts.Y[item[ind1[i]]]],
                [scMDLpts.Z[item[ind0[i]]],scMDLpts.Z[item[ind1[i]]]],color=rotColor[cInd])
    cInd = (cInd + 1)%len(rotColor)

show(block=False)
print 'Done Plotting'


