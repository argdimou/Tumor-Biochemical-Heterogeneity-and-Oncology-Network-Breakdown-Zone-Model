# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 20:49:50 2022

@author: argdi
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 12:44:58 2022

@author: argdi
"""

import networkx as nx
import numpy as np
import sys
import copy


def Pk(kList): #k and Pk values 
    d = {}
    kAll = []
    PkAll = []
    for i in range(len(kList)):    
        if kList[i] not in d:
            d[kList[i]] = 1
        else:        
            d[kList[i]] += 1
    Pktemp = np.array(list(d.values()))
    PkAll.append(Pktemp/sum(Pktemp))
    kAll.append(np.array(list(d.keys())))
    return [PkAll, kAll]

if __name__=='__main__':
    N = 10000
    nei0 = 20
    p0 = nei0/(N-1)
    np.random.seed(int(sys.argv[1]))
    graphInitial = nx.fast_gnp_random_graph(N, p0, seed=int(sys.argv[1]))
    r0List = np.arange(0, 1.1, 0.1)
    sizeMax = np.zeros(( len(r0List), N) )
    sizeMax[:, 0] = 1.
    remNodes = []
    clustersNum = []
    maxCluster = 0

    ktemp = np.array(graphInitial.degree())
            
    kList = ktemp[:,1]
    nodes = ktemp[:,0]
    kListSort = np.sort(kList)
    kLimInd = int(len(kList)/2)

    kLim = kListSort[kLimInd]
    totalNodes = np.array(graphInitial.nodes())

    if (int(sys.argv[1]) < 501): # each zone will contain the kLim once for half of the realizations
        zone1Ind = np.where(kList<=kLim)[0]
        zone1Nodes = nodes[zone1Ind]
        zone0Ind = np.where(kList>kLim)[0]
        zone0Nodes = nodes[zone0Ind]
    else:
        zone1Ind = np.where(kList<kLim)[0]
        zone1Nodes = nodes[zone1Ind]
        zone0Ind = np.where(kList>=kLim)[0]
        zone0Nodes = nodes[zone0Ind]

    
    cntsAll = []
    for indR0, r0 in enumerate(r0List):
        numRem = 0
        graph = graphInitial.copy()
        
        kappa = 100
        zoneP = np.ones(2)    
        
        zoneP[0] = zoneP[1]*(1-r0)
        
        while kappa>2:
            while True:
                randomNode = np.random.choice(np.array(graph.nodes())) 
                tpl = np.where(zone0Ind==randomNode)
                if len(tpl[0] > 0):
                    zoneNum=0
                else:
                    tpl = np.where(zone1Ind==randomNode)
                    if len(tpl[0] > 0):
                        zoneNum=1
                    else:
                        sys.exit("error")
                        
                        
                pRemoval = zoneP[zoneNum]
                if np.random.rand() < pRemoval:
                    numRem +=1 
                    break
            graph.remove_node(randomNode)
            ktemp = list(graph.degree())
            
            kList = np.array([i[1] for i in ktemp])
            kM = np.mean(kList)
            k2List = kList**2
            k2M = np.mean(k2List)
            kappa = k2M/kM
            clustSizes = [len(c) for c in nx.connected_components(graph)]

            relSize = max(clustSizes)/len(graph.nodes())
            sizeMax[indR0, numRem-1] = relSize
            
        clustersNum.append(len(clustSizes))    
        remNodes.append(N - len(graph.nodes()))
        

        cnts = [0]*max(clustSizes)
        for i in range(1, max(clustSizes) + 1):
            cnts[i-1] = clustSizes.count(i)
        cntsAll.append(cnts)
        if max(clustSizes) > maxCluster:
            maxCluster = max(clustSizes)


    for i in range(len(cntsAll)):
        if len(cntsAll[i]) < maxCluster:
            cntsAll[i].extend([0]*(maxCluster - len(cntsAll[i])))
        


    numFold = int(sys.argv[1])
    saveArray = np.array([clustersNum, remNodes]) 
    np.savetxt('ErdClusts_RemNodes' + str(numFold) + '.txt', saveArray)


    saveArray = np.array(cntsAll, dtype=object,) 
    np.savetxt('ErdosDistr' + str(numFold) + '.txt', saveArray, delimiter='  ',  comments='', fmt='%s')

    np.savetxt('maxSize' + str(numFold) + '.txt', sizeMax)
