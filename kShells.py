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

if __name__ == '__main__':

    N = 10000
    nei0 = 20
    p0 = nei0/(N-1)
    np.random.seed(int(sys.argv[1])) #give the seed of the realization
    graph0 = nx.fast_gnp_random_graph(N, p0, seed=int(sys.argv[1]))
    kappaValues = []
    k=1
    graph = copy.deepcopy(graph0)
    clustersNum = []
    maxCluster = 0
    relSizeAll = []
    shells = []

    while True: #k-shell decomposition
        shellsTemp = []
        while True:
            kTemp = np.array(graph.degree())
            if len(kTemp)<1:
                break
            nodeList = kTemp[:,0]
            kList = kTemp[:,1]
            remNodesInd = np.where(kList<=k)[0]
            remNodes = nodeList[remNodesInd]
            shellsTemp.extend(remNodes)
            if len(remNodes)>0:
                graph.remove_nodes_from(remNodes)
                
            else:
                shells.append(shellsTemp)
                k += 1
                shellsTemp = []
        break

    shells.append(shellsTemp)    
    zone0Nodes = shells[-1]
    zone1Nodes = []
    for i in range(len(shells)-1):
        zone1Nodes.extend(shells[i])
        
    zone0Nodes, zone1Nodes = np.array(zone0Nodes), np.array(zone1Nodes)
    r0List = np.arange(0, 1.1, 0.1)
    sizeMax = np.zeros(( len(r0List), N) )

    remNodes = []
    cntsAll = []

    ktemp = list(graph0.degree())            
    kList = np.array([i[1] for i in ktemp])
    clustSizes = [len(c) for c in nx.connected_components(graph0)]
    relSize = max(clustSizes)/len(graph0.nodes())
    relSizeAll.append(relSize)

    for indR0, r0 in enumerate(r0List):
        numRem = 0
        graph = graph0.copy()    
        kappa = 100
        zoneP = np.ones(2)        
            
        zoneP[0] = zoneP[1]*(1-r0)
        
        kappaValues = []
        
        while kappa>2: #removing process
            while True:
                randomNode = np.random.choice(np.array(graph.nodes())) 
                tpl = np.where(zone0Nodes==randomNode)  #internal zone 
                if len(tpl[0] > 0):
                    zoneNum=0
                else:
                    tpl = np.where(zone1Nodes==randomNode) #external zone
                    if len(tpl[0] > 0):
                        zoneNum=1
                    else:
                        sys.exit("error")                        
                        
                pRemoval = zoneP[zoneNum]
                if np.random.rand() < pRemoval:
                    numRem += 1
                    break
            graph.remove_node(randomNode)
            ktemp = list(graph.degree())            
            kList = np.array([i[1] for i in ktemp])        
            kM = np.mean(kList)
            k2List = kList**2
            k2M = np.mean(k2List)
            kappa = k2M/kM
            kappaValues.append(kappa)
            clustSizes = [len(c) for c in nx.connected_components(graph)]
            relSize = max(clustSizes)/len(graph.nodes())         
            sizeMax[indR0, numRem-1] = relSize # relative size of the largest cluster 
            
        clustersNum.append(len(clustSizes))  #number of clusters 
        remNodes.append(N - len(graph.nodes())) #removed nodes
        
        cnts = list(np.zeros(max(clustSizes))) 
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
    #
    saveArray = np.array(cntsAll, dtype=object,) 
    np.savetxt('ErdosDistr' + str(numFold) + '.txt', saveArray, delimiter='  ',  comments='', fmt='%s')


    np.savetxt('maxSize' + str(numFold) + '.txt', sizeMax)
