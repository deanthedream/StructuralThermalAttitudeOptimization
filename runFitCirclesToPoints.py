#runFitCirclesToPoints.py
import numpy as np
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.cluster import SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn import datasets
import heapq
import math

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

    k=6 #number of points to find#10 is horrible, 7 is ok
    for i in ptInds:
        if len(cPS[i]) >= k:#Only looking for k closest points.... might want to rethink this Maybe calculate how many based on filtering distance and the minimum spatial distance between points on the unit sphere
            continue
        tmp = map(lambda x: np.linalg.norm(x - scMDLpts.pts[i]), scMDLpts.pts)#Calculate distance between points and all other points
        #numPts = k - len(cPS[i]) + 1#num of points to find
        ptsToAdd = np.argsort(tmp)#heapq.nsmallest(numPts,np.nditer(tmp),key=tmp.__getitem__)#heapq.nsmallest(numPts, tmp)
        if len(cPS[i]) == 1:
            if cPS[i][0] == i:#If cPS has only 1 item and that item is itself
                numPts = k - len(cPS[i]) - 1 + 1#num of points to find
                cPS[i] = ptsToAdd[1:numPts].tolist()#Reassign points
            else:#Otherwise append all of them
                numPts = k - len(cPS[i]) + 1#num of points to find
                cPS[i] += ptsToAdd[1:numPts].tolist()#data format might be messy
        else:#Otherwise append all of them
            numPts = k - len(cPS[i]) + 1#num of points to find
            cPS[i] += ptsToAdd[1:numPts].tolist()#data format might be messy

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
            ptVectors.append(tmp/np.linalg.norm(tmp))
        cPSvectors[i] = ptVectors
    print 'Done Calculating cPSvectors'
#################################################################

#### Calc Vector Dot Products
def calcRowDots(j, jInd, cPS, iInd, vectij, cPSvectors):
    """Calcs point angles from vectij to all other vectors
    Args:
        cPS - dict of all points and all connecting points
        iInd - point ind of originating vectij
        vectij - vector from point iInd to j
        cPSvectos - dict of all vectors from all points to all connecting points
    return:
        angleRow - dot products between originating vector vectij and outgoing vectors vectij2
        indRow
    """
    angleRow = list()
    indRow = list()
    for j2 in np.arange(len(cPS[iInd])):
        j2Ind = cPS[iInd][j2]
        if j == j2:#skip if its referencing ij
            angleRow.append(None)
        else:
            vectij2 = cPSvectors[iInd][j2]
            angleRow.append(dot(vectij,vectij2))
        indRow.append([jInd,j2Ind])
    return angleRow, indRow

def calcMatrixDots(cPS,cPSvectors,iInd):
    """Create a matrix of all primary and secondary connections
    """
    angleMatrix = list()
    indMatrix = list()
    for j in np.arange(len(cPS[iInd])):#Iterate over all points j
        jInd = cPS[iInd][j]#ind of Point j
        vectij = cPSvectors[iInd][j]#Grab vector connecting ith point to a jth point
        tmpAngleRow, tmpIndRow = calcRowDots(j, jInd, cPS, iInd, vectij, cPSvectors)
        angleMatrix.append(tmpAngleRow)
        indMatrix.append(tmpIndRow)
    return angleMatrix, indMatrix


cPSdp = {}#Closest Points dot product
cPSindMat = {}
for i in ptInds:#Instantiate struct of closest points for each ind
    cPSdp[i] = [i]
    cPSindMat[i] = [i]

for i in ptInds:#Iterate over all Points i
    iInd = i#ind of point1
    angleMatrix, indMatrix = calcMatrixDots(cPS,cPSvectors,iInd)
    cPSdp[iInd] = angleMatrix
    cPSindMat[iInd] = indMatrix
    #print saltyburrito
print 'Done calcMatrixDots'

cPSconn = {}
cPSconn2 = {}
for i in ptInds:#Instantiate struct of closest points for each ind
    cPSconn[i] = [i]
    cPSconn2[i] = [i]

#Here we go through each angle combination of a point, for each angle combination, we propagate outwards and create a data structure of the best fitting parameters
for i in ptInds:#Iterate over all Points i
    iInd = i#ind of point1
    for j in np.arange(len(cPSdp[iInd])):#Iterate over each dot product combo rows
        for k in np.arange(len(cPSdp[iInd][j])):#iterate over each dot product combos cols
            if j == k:#Don't want to include identical dot products
                continue#Technically value will be none
            else:#Explore outwards 1
                currDPval = cPSdp[iInd][j][k]#save current dp value
                if math.isnan(currDPval):
                    continue
                currDPind = cPSindMat[iInd][j][k]#will have left vector and right vector
                minErrorLR = None
                minErrorjj = None
                minErrorkk = None
                err = 100.
                for lrInd in currDPind:#Explore left and right indices
                    #We are passing in an iInd to lrInd vector
                    #We want to know of there is a complementary vector with same dot product for lrInd
                    minErrorjj_lrInd = None
                    minErrorkk_lrInd = None
                    err_lrInd = 100.#arbitrarily large value
                    for jj in np.arange(len(cPSdp[lrInd])):#Iterate over each dot product combo rows
                        for kk in np.arange(len(cPSdp[lrInd][jj])):#iterate over each dot product combos cols
                            if jj == kk:#Don't want to include identical dot products
                                continue#Technically value will be none
                            else:#evaluate
                                # ttInd = cPS[iInd].index(lrInd)#getting ind of for finding vector
                                # ttVect = cPSvectors[iInd][ttInd]
                                # cPS[lrInd]
                                if np.min([np.abs(cPSdp[lrInd][jj][kk]),np.abs(currDPval)]) == 0.:
                                    tmperr = np.abs(cPSdp[lrInd][jj][kk] - currDPval)/np.min([np.abs(cPSdp[lrInd][jj][kk]),np.abs(currDPval)])#difference in dp
                                elif np.min([np.abs(cPSdp[lrInd][jj][kk]),np.abs(currDPval)]) == None:
                                    tmperr = 100.
                                else:
                                    tmperr = np.abs(cPSdp[lrInd][jj][kk] - currDPval)#difference in dp
                                # print saltyburrito
                                if tmperr < err_lrInd:#Update because we found the most in common angle
                                    minErrorjj_lrInd = jj
                                    minErrorkk_lrInd = kk
                                    err_lrInd = tmperr
                    if err_lrInd < err:
                        minErrorLR = lrInd
                        minErrorjj = minErrorjj_lrInd
                        minErrorkk = minErrorkk_lrInd
                        err = err_lrInd
                # if minErrorLR == None:
                #     print saltyburrito
                if not minErrorLR == None:
                    tmpList = [minErrorLR, minErrorjj, minErrorkk, err, cPSdp[minErrorLR][minErrorjj][minErrorkk]]
                    if cPSconn[iInd] == [[iInd]]:
                        cPSconn[iInd] = [tmpList]
                    else:
                        cPSconn[iInd].append(tmpList)
    print(str(float(i)/float(len(ptInds))))
print 'Done finding Most Similar Angles between two Points'
                
#### Distill data down
for i in ptInds:#Iterate over all Points i
    iInd = i#ind of point1
    minErrorInd = np.argmin([err[3] for err in cPSconn[iInd]])
    cPSconn2[iInd] = cPSconn[iInd][minErrorInd]


tmpInds1 = ptInds.tolist()
connList = list()#contains a list of connections between nodes
ind1 = 0
while len(tmpInds1) > 1:#Create another strain
    tmpList = list()
    ind = tmpInds1[0]
    tmpList.append(ind)
    tmpInds1.remove(ind)
    tmpInds2 = ptInds.tolist()
    tmpInds2.remove(ind)
    while cPSconn2[ind][0] in tmpInds2:#Should iterate until there are no points left in strain
        ind = cPSconn2[ind][0]
        tmpList.append(ind)
        tmpInds2.remove(ind)
    connList.append(tmpList)
print 'Done Creating connection List'

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


#### Pull out Max length connList #################
tmpMaxInd = np.argmax([len(sl) for sl in connList])
connInds = connList[tmpMaxInd]
ax.scatter(np.asarray(scMDLpts.X)[connInds],np.asarray(scMDLpts.Y)[connInds],np.asarray(scMDLpts.Z)[connInds], c='k', marker='o', s=50)
show(block=False)
###################################################

    #In work for same direction but must propagate twice
    # for j in np.arange(len(cPS[iInd])):#Iterate over all points j
    #     jInd = cPS[iInd][j]#ind of Point j
    #     vectij = cPSvectors[iInd][j]#Grab vector connecting ith point to a jth point
        
    #     dprowk = list()
    #     tmpkInd = list()
    #     for k in np.arange(len(cPS[jInd])):#iterate over all points k
    #         kInd = cPS[jInd][k]#ind of Point k
    #         if kInd == iInd:#Skip the ind that is already being considered
    #             continue
    #         vectjk = cPSvectors[jInd][k]
    #         dprowk.append(abs(dot(vectij,vectjk)))#add dot product to row of dot products
    #         tmpkInd.append(kInd)
    #     ksameInd = tmpkInd[np.argmax(dprowk)]#find ind of vector in closest to the same direction as previous
    #     ksameVal = np.max(dprowk)
    #     #WOULD DO SAME FOR ORTHOGONAL HERE BUT MIGHT BE MORE CHALLENGING since 2 vectors could be orthogonal and should be mutually orthongonal
    #     print saltyburrito  

    # dp = cPSdp[i]
#print saltyburrito
#################################################################


# #### KMeans Clustering ########################################
# n=8
# estimators=KMeans(n_clusters=n)
# estimators.fit(np.asarray(scMDLpts.allVertices))
# labels = estimators.labels_
# fig5 = figure(num=1000+n,figsize=(4, 3))
# ax5 = gca(projection='3d')#Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
# ax5.scatter(np.asarray(scMDLpts.allVertices)[:, 0], np.asarray(scMDLpts.allVertices)[:, 1], np.asarray(scMDLpts.allVertices)[:, 2],\
#     c=labels.astype(np.float), edgecolor='k')
# title('KMeans estimator num: ' + str(n))
# show(block=False)
# ###############################################################


# Performs Ok but not great
# #### KMeans Clustering ########################################
# estimators = [0,0]#Produces a pretty good fit to some features [KMeans(n_clusters=8)]
# for n in [2,3,4,5,6,7,8,9,10,11,12]:
#     #estimators[n].fit(scMDLpts.pts)
#     estimators.append(KMeans(n_clusters=n))
#     estimators[n].fit(scMDLpts.pts)
#     labels = estimators[n].labels_
#     fig5 = figure(num=1000+n,figsize=(4, 3))
#     ax5 = gca(projection='3d')#Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
#     ax5.scatter(np.asarray(scMDLpts.pts)[:, 0], np.asarray(scMDLpts.pts)[:, 1], np.asarray(scMDLpts.pts)[:, 2],\
#                    c=labels.astype(np.float), edgecolor='k')
#     title('KMeans estimator num: ' + str(n))
#     show(block=False)
#     savefig('KMeansEstimatorNum'+ str(n) + '.png')
# ###############################################################

# #### Spectral Clustering ########################################
# estimators = [0,0]#Produces a pretty good fit to some features [KMeans(n_clusters=8)]
# for n in [2,3,4,5,6,7,8,9,10,11,12]:
#     #estimators[n].fit(scMDLpts.pts)
#     estimators.append(SpectralClustering(n_clusters=n))
#     estimators[n].fit(scMDLpts.pts)
#     labels = estimators[n].labels_
#     fig5 = figure(num=2000+n,figsize=(4, 3))
#     ax5 = gca(projection='3d')#Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
#     ax5.scatter(np.asarray(scMDLpts.pts)[:, 0], np.asarray(scMDLpts.pts)[:, 1], np.asarray(scMDLpts.pts)[:, 2],\
#                    c=labels.astype(np.float), edgecolor='k')
#     title('spectral estimator num: ' + str(n))
#     show(block=False)
# ###############################################################

# #### Gaussian Mixture ########################################
# estimators = [0,0]#Produces a pretty good fit to some features [KMeans(n_clusters=8)]
# for n in [2,3,4,5,6,7,8,9,10,11,12]:
#     #estimators[n].fit(scMDLpts.pts)
#     estimators.append(GaussianMixture(n_components=n))
#     estimators[n].fit(scMDLpts.pts)
#     labels = estimators[n].labels_
#     fig5 = figure(num=2000+n,figsize=(4, 3))
#     ax5 = gca(projection='3d')#Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
#     ax5.scatter(np.asarray(scMDLpts.pts)[:, 0], np.asarray(scMDLpts.pts)[:, 1], np.asarray(scMDLpts.pts)[:, 2],\
#                    c=labels.astype(np.float), edgecolor='k')
#     title('spectral estimator num: ' + str(n))
#     show(block=False)
# ###############################################################


# #### Equation For Circle Center Vector ########################
# #https://www.physicsforums.com/threads/equation-of-a-circle-through-3-points-in-3d-space.173847/
# testpts = np.random.choice(selectInds,3,replace=False)#pick 3 random points from selectInds
# p1 = np.asarray([np.asarray(scMDLpts.X)[testpts[0]],np.asarray(scMDLpts.Y)[testpts[0]],np.asarray(scMDLpts.Z)[testpts[0]]])
# p2 = np.asarray([np.asarray(scMDLpts.X)[testpts[1]],np.asarray(scMDLpts.Y)[testpts[1]],np.asarray(scMDLpts.Z)[testpts[1]]])
# p3 = np.asarray([np.asarray(scMDLpts.X)[testpts[2]],np.asarray(scMDLpts.Y)[testpts[2]],np.asarray(scMDLpts.Z)[testpts[2]]])


# # #Solve For Circle Center
# # AA = np.asarray([p2-p1,p3-p1,p3-p2])
# # bb = np.asarray([dot(p1,p1)-dot(p2,p2),dot(p1,p1)-dot(p3,p3),dot(p2,p2)-dot(p3,p3)]).T
# # cc = 0.5*np.matmul(np.linalg.inv(AA),bb)

# #Get vector representing circle plane
# r_p2_p1 = p2-p1
# r_p3_p1 = p3-p1
# planeNormal = np.cross(r_p2_p1,r_p3_p1)/np.linalg.norm(np.cross(r_p2_p1,r_p3_p1))#Vector Normal to surface or plane passing through all 3 points

# #Find subset of points within 

# #Plot Circle
# ax.scatter(cc[0],cc[1],cc[2],color='green',s=10,marker='o')#plots circle center
# num = 100
# theta = np.linspace(0,2*np.pi,num=num)
# e1 = (p1-cc)/np.linalg.norm(p1-cc)
# e3 = np.cross(e1,p2-cc)/np.linalg.norm(np.cross(e1,p2-cc))#unit vector perpendicular to circle
# e2 = np.cross(e3,e1)
# Rc = np.linalg.norm(p1-cc)
# circlePts = np.transpose(np.tile(cc,(num,1))) + Rc*outer(e1,np.cos(theta)) + Rc*outer(e2,np.sin(theta))
# ax.plot(circlePts[0],circlePts[1],circlePts[2],color='red')
# show(block=False)
# ###############################################################




#-------------------------------------------------------------------------------
# FIT CIRCLE 2D
# - Find center [xc, yc] and radius r of circle fitting to set of 2D points
# - Optionally specify weights for points
#
# - Implicit circle function:
#   (x-xc)^2 + (y-yc)^2 = r^2
#   (2*xc)*x + (2*yc)*y + (r^2-xc^2-yc^2) = x^2+y^2
#   c[0]*x + c[1]*y + c[2] = x^2+y^2
#
# - Solution by method of least squares:
#   A*c = b, c' = argmin(||A*c - b||^2)
#   A = [x y 1], b = [x^2+y^2]
#-------------------------------------------------------------------------------
def fit_circle_2d(x, y, w=[]):
    
    A = array([x, y, ones(len(x))]).T
    b = x**2 + y**2
    
    # Modify A,b for weighted least squares
    if len(w) == len(x):
        W = diag(w)
        A = dot(W,A)
        b = dot(W,b)
    
    # Solve by method of least squares
    c = linalg.lstsq(A,b,rcond=None)[0]
    
    # Get circle parameters from solution c
    xc = c[0]/2
    yc = c[1]/2
    r = sqrt(c[2] + xc**2 + yc**2)
    return xc, yc, r

#-------------------------------------------------------------------------------
# RODRIGUES ROTATION
# - Rotate given points based on a starting and ending vector
# - Axis k and angle of rotation theta given by vectors n0,n1
#   P_rot = P*cos(theta) + (k x P)*sin(theta) + k*<k,P>*(1-cos(theta))
#-------------------------------------------------------------------------------
def rodrigues_rot(P, n0, n1):
    
    # If P is only 1d array (coords of single point), fix it to be matrix
    if P.ndim == 1:
        P = P[newaxis,:]
    
    # Get vector of rotation k and angle theta
    n0 = n0/linalg.norm(n0)
    n1 = n1/linalg.norm(n1)
    k = cross(n0,n1)
    k = k/linalg.norm(k)
    theta = arccos(dot(n0,n1))
    
    # Compute rotated points
    P_rot = zeros((len(P),3))
    for i in range(len(P)):
        P_rot[i] = P[i]*cos(theta) + cross(k,P[i])*sin(theta) + k*dot(k,P[i])*(1-cos(theta))

    return P_rot

#-------------------------------------------------------------------------------
# ANGLE BETWEEN
# - Get angle between vectors u,v with sign based on plane with unit normal n
#-------------------------------------------------------------------------------
def angle_between(u, v, n=None):
    if n is None:
        return arctan2(linalg.norm(cross(u,v)), dot(u,v))
    else:
        return arctan2(dot(n,cross(u,v)), dot(u,v))

    
#-------------------------------------------------------------------------------
# - Make axes of 3D plot to have equal scales
# - This is a workaround to Matplotlib's set_aspect('equal') and axis('equal')
#   which were not working for 3D
#-------------------------------------------------------------------------------
def set_axes_equal_3d(ax):
    limits = array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    spans = abs(limits[:,0] - limits[:,1])
    centers = mean(limits, axis=1)
    radius = 0.5 * max(spans)
    ax.set_xlim3d([centers[0]-radius, centers[0]+radius])
    ax.set_ylim3d([centers[1]-radius, centers[1]+radius])
    ax.set_zlim3d([centers[2]-radius, centers[2]+radius])

#### Of This Point Subset, fit circle ######

############################################