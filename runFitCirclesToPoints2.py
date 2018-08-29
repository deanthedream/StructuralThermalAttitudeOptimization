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
cPSconn = {}
cPSdpmat = {}
cPSdpmatL = {}#Contains angle between vectiIndjInd and vectL
cPSdpmatR = {}#Contains angle between vectjIndiInd and vectR
cPSdpmatInd = {}
for i in ptInds:#Instantiate struct of closest points for each ind
    cPSconn[i] = [i]
    cPSdpmat[i] = [i]#Contains angle between vectL and vectR
    cPSdpmatL[i] = [i]#Contains angle between vectiIndjInd and vectL
    cPSdpmatR[i] = [i]#Contains angle between vectjIndiInd and vectR
    cPSdpmatInd[i] = [i]

for i in ptInds:#Iterate over all Points i
    iInd = i
    mat2 = list()
    mat2L = list()
    mat2R = list()
    mat2Ind = list()
    for jInd in cPS[iInd]:#iterate over connecting points
        if (jInd == iInd) or not (iInd in cPS[jInd]) or not (jInd in cPS[iInd]):#Dont count the connecting vectors
            mat2.append(None)
            mat2L.append(None)
            mat2R.append(None)
            mat2Ind.append(None)
            continue
        vectiIndjInd = scMDLpts.pts[jInd] - scMDLpts.pts[iInd]#cPSvectors[iInd][cPS[iInd].index(jInd)]#vector from iInd to jInd
        vectjIndiInd = scMDLpts.pts[iInd] - scMDLpts.pts[jInd]#cPSvectors[jInd][cPS[jInd].index(iInd)]#vector from jInd to iInd
        mat = list()
        matL = list()
        matR = list()
        matInds = list()
        for tmpiiInd in np.arange(len(cPS[iInd])):#vectors connecting to iInd
            iiInd = cPS[iInd][tmpiiInd]
            vectL = cPSvectors[iInd][tmpiiInd]#vectors of iInd that are not iInd->jInd
            vectRrow = list()
            vectjIndiIndDvectR = list()
            #vectiIndjIndDvectL = list()
            vectRrowinds = list()
            for tmpjjInd in np.arange(len(cPS[jInd])):#Vectors connecting to jInd
                jjInd = cPS[jInd][tmpjjInd]
                vectR = cPSvectors[jInd][tmpjjInd]#vectors of jInd that are not jInd->iInd
                vectRrow.append(np.arccos(udot(vectL,vectR)))#angles between vectL and vectR in radians
                if tmpiiInd == 0:
                    vectjIndiIndDvectR.append(np.arccos(udot(vectjIndiInd,vectR)))
                #vectiIndjIndDvectL.append(np.arccos(udot(vectiIndjInd,vectL)))
                if math.isnan(np.arccos(udot(vectL,vectR))):
                    print saltyburrito
                vectRrowinds.append([iInd, jInd, tmpiiInd, iiInd, tmpjjInd, jjInd])
            mat.append(vectRrow)
            matL.append(np.arccos(udot(vectiIndjInd,vectL)))#has length len(cPS[iInd])
            if tmpiiInd == 0:
                matR = vectjIndiIndDvectR#Each matR has length len(cPS[jInd])
            matInds.append(vectRrowinds)
        mat2.append(mat)
        mat2L.append(matL)#has length
        mat2R.append(matR)
        mat2Ind.append(matInds)
    cPSdpmatL[iInd] = mat2L#Contains angle between vectiIndjInd and vectL
    cPSdpmatR[iInd] = mat2R#Contains angle between vectjIndiInd and vectR
    cPSdpmat[iInd] = mat2
    cPSdpmatInd[iInd] = mat2Ind
print('Done creating Matrix of Dot Products')
########################################################################

#cPSdpmatL has len(ptInds)
#cPSdpmatL[iInd] has len(cPS[iInd])
#cPSdpmatL[iInd][iInd] has len(cPS[iInd])

#cPSdpmatR has len(ptInds)
#cPSdpmatR[iInd] has len(cPS[iInd])
#cPSdpmatR[iInd][iInd] has len(cPS[iInd])

#cPSdpmat has len(ptInds)
#cPSdpmat[iInd] has len(cPS[iInd])
#cPSdpmat[iInd][iInd] has len(cPS[iInd])
#cPSdpmat[iInd][iInd][jInd] has len(cPS[jInd])#I think


connList = list()
for iInd in ptInds:
    #cPSdpmat[iInd]
    #cPSdpmatL[iInd]
    #cPSdpmatR[iInd] 
    for jInd in cPS[iInd]:
        if (jInd == iInd) or not (iInd in cPS[jInd]) or not (jInd in cPS[iInd]):#Dont count the connecting vectors
            continue
        #do stuff
        for tmpiiInd in np.arange(len(cPS[iInd])):#vectors connecting to iInd
            iiInd = cPS[iInd][tmpiiInd]
            for tmpjjInd in np.arange(len(cPS[jInd])):#Vectors connecting to jInd
                jjInd = cPS[jInd][tmpjjInd]

                psiD = cPSdpmat[iInd][cPS[iInd].index(jInd)][tmpiiInd][tmpjjInd] #angle between vectL and vectR
                th1 = cPSdpmatL[iInd][cPS[iInd].index(jInd)][tmpiiInd]
                th2 = cPSdpmatR[iInd][cPS[iInd].index(jInd)][tmpjjInd]
                tol = 4.*np.pi/180.#5 deg tolerance
                psiC1 = 2*np.pi - 2.*th1
                psiC2 = 2*np.pi - 2.*th2
                psiC3 = 2*np.pi - th1 - th2
                if np.abs(psiD-psiC1) < tol and np.abs(psiD-psiC2) < tol and np.abs(psiD-psiC3) < tol:#Any 2 conditions produce the same angle
                    connList.append([iInd,jInd,cPS[jInd][tmpjjInd]])
                # if (np.abs(psiD-psiC1) < tol and np.abs(psiD-psiC2) < tol) or\
                #     (np.abs(psiD-psiC1) < tol and np.abs(psiD-psiC3) < tol) or\
                #     (np.abs(psiD-psiC2) < tol and np.abs(psiD-psiC3) < tol):#Any 2 conditions produce the same angle
                #     connList.append([iInd,jInd,cPS[jInd][tmpjjInd]])
                    #([iInd,jInd,cPS[iInd][tmpiiInd],cPS[jInd][tmpjjInd]])
print('Done generating connList')


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