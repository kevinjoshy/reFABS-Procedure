import pandas as pd 
import numpy as np
from numpy import random, array
import time 
import os # Dont think this is used 
from multiprocessing import Pool

from ontologyPackage.ontologySTATanalysis import *
from ontologyPackage.ontologySTATanalysis import go2gene

pd.set_option('display.max_colwidth', -1) # Values in columns won't be shortned 
pd.set_option('chained_assignment',None) # Disabling chained assignments 

geneColID = ["chromosome","source","type","start","end","score","strand","phase","gene_symbol","gene_ensID","length","entrezid"]
geneAnnotationDF = pd.read_csv('entrez_id/geneAnnotationsDF_Selected_entrezID.csv', sep=',', comment='#', low_memory=False, header=0, names=geneColID)
chromosomeColID = ['chromosome','source','type','start','end','score','strand','phase']
chromosomesDF = pd.read_csv('chromosomesDF.csv', sep=',', comment='#', low_memory=False, header=0, names=chromosomeColID)


import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests as padjust


#Dropping unneccesary stuff and resetting index from 0
geneAnnotationDF = geneAnnotationDF[geneAnnotationDF.gene_symbol != 'CTD-2207O23.3']
geneAnnotationDF = geneAnnotationDF.sort_values(by=['chromosome', 'start'])
chromosomesDF = chromosomesDF.sort_values(by=['chromosome'])
chromosomesDF = chromosomesDF.reset_index()
geneAnnotationDF = geneAnnotationDF.reset_index()
chromosomesDF = chromosomesDF.drop(columns=['score', 'strand', 'phase', 'index'])
geneAnnotationDF = geneAnnotationDF.drop(columns=['score', 'phase', 'index'])


# Making a list of DataFrames to be used in addwindow tss so that if - gene the start is its end and its end is its start 
u = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','MT','X','Y'] # Chromosome ids 
geneDFL = []
for c in u:
    geneDF = geneAnnotationDF.loc[geneAnnotationDF['chromosome'] == c].copy()
    for i, row in geneDF.iterrows():
        if row['strand'] == '-':
            start = row['end']
            end = row['start']
            geneDF.at[i, 'start'] = start
            geneDF.at[i, 'end'] = end
    geneDF = geneDF.sort_values(by=['start'])
    geneDF = geneDF.reset_index()
    geneDF = geneDF.drop(columns=['index'])
    geneDFL.append(geneDF)

# Making a list of DataFrames to be used in future functions in var geneDFL
u = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','MT','X','Y']
geneDFLBODY = []
for cc in u:
    geneDFbod = geneAnnotationDF.loc[geneAnnotationDF['chromosome'] == cc].copy()
    geneDFbod = geneDFbod.sort_values(by=['start'])
    geneDFbod = geneDFbod.reset_index()
    geneDFbod = geneDFbod.drop(columns=['index'])
    geneDFLBODY.append(geneDFbod)


# probability that a randomly chosen gene is in a certain chromosome. Ordered from 1 - MT,X,Y | As shown chromosome 1 has highest prob of selection
probabilityL = [0.08069467597786323, 0.07850260775517634, 0.06427388269957329, 0.06165457287524136, 0.05884231002379223, 0.05536363751707386, 0.05164908592141782, 0.0470440371698401, 0.044858119022639725, 0.043367989830914035, 0.043785860460049814, 0.04319875643982657, 0.037069107410936775, 0.034696265410969616, 0.033058580434273406, 0.029281523923645837, 0.026986378244674356, 0.026051531764496164, 0.018999829086634144, 0.02088839915494138, 0.015140187376094325, 0.016471877735809576, 5.370214538878023e-06, 0.05057780526426391, 0.017537607961181506]


# Generates random sites
dl = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'] # 23, 24, 25
chrD = {'1': [1, 248956422],'2': [1, 242193529],'3': [1, 198295559],'4': [1, 190214555],'5': [1, 181538259],'6': [1, 170805979],'7': [1, 159345973],'8': [1, 145138636],'9': [1, 138394717],'10': [1, 133797422],'11': [1, 135086622],'12': [1, 133275309],'13': [1, 114364328],'14': [1, 107043718],'15': [1, 101991189],'16': [1, 90338345],'17': [1, 83257441],'18': [1, 80373285],'19': [1, 58617616],'20': [1, 64444167],'21': [1, 46709983],'22': [1, 50818468],'23': [1, 16569],'24': [1, 156040895],'25': [2781480, 56887902]}
def nRandSites(data):
    nSite = data[0]
    rand = data[1] 
    sitesbyC =[[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

    random.seed(rand)
        
    for site in range(nSite):
        chrN = random.choice(dl, p=probabilityL)
        randnum = round(random.uniform(chrD[chrN][0], chrD[chrN][1]))
        sitesbyC[int(chrN) - 1].append(randnum)
    for l in sitesbyC:
        l.sort()
    return array([array(l) for l in sitesbyC])

def nRandSitesSim(nSite, nSim):
    totalsites = [nSite] * nSim
    chunks = []
    for i in totalsites:
        chunks.append([i])
    for l in chunks:
        l.append(random.randint(100000))
    
    pool = Pool(5)
    result = pool.map(nRandSites, chunks)
    pool.close()
    return(result) 

##### Adds midpoint in gene and upper and lower bound values to geneDFLtss based on size of window or nearby genes. Adding around geneTSS.
def addWindowTSS(window):
    chrD = {'1': [1, 248956422],'2': [1, 242193529],'3': [1, 198295559],'4': [1, 190214555],'5': [1, 181538259],'6': [1, 170805979],'7': [1, 159345973],'8': [1, 145138636],'9': [1, 138394717],'10': [1, 133797422],'11': [1, 135086622],'12': [1, 133275309],'13': [1, 114364328],'14': [1, 107043718],'15': [1, 101991189],'16': [1, 90338345],'17': [1, 83257441],'18': [1, 80373285],'19': [1, 58617616],'20': [1, 64444167],'21': [1, 46709983],'22': [1, 50818468],'23': [1, 16569],'24': [1, 156040895],'25': [2781480, 56887902]} 
    geneWindowDFL = []
    for i in range(len(geneDFL)):
        df = geneDFL[i].copy()
        geneWindowDFL.append(df)
    outList = []
    chrIDNUM = 0
    for dfG in geneWindowDFL: # ordered by chromosome 1- MT,X,Y
        chrIDNUM += 1 
        dfIndexMax = len(dfG) - 1   
        for i, row in dfG.iterrows():
            start = row['start']
            
            dfG.at[i, 'midBody'] = (((row['length'] / 2)) + start) # Middle of gene body 
            if i == 0 or i == dfIndexMax: # Case for start and end genes
                if i == 0:
                    dfG.at[i, 'lowB'] = max(start - window, chrD[str(chrIDNUM)][0]) # chrD[chrIDNUM][0] = start val of the chromosome
                    dfG.at[i, 'upperB'] = min(start + window, ((dfG.iat[i+1, 3] - start) / 2) + start)
                else: # index max case 
                    dfG.at[i, 'lowB'] = max(start - window, start - ((start - dfG.iat[i-1, 3]) / 2)) # dfG.iat[i-1, 3] = start of gene below
                    dfG.at[i, 'upperB'] = min(start + window, chrD[str(chrIDNUM)][1]) # chrD[chrIDNUM][1] = end val of chromosome 
            else: # i > 0 or i < len(geneDF) - 1 
                dfG.at[i, 'lowB'] = max(start - window, start - ((start - dfG.iat[i-1, 3]) / 2)) # dfG.iat[i-1, 3] = start of gene below
                dfG.at[i, 'upperB'] = min(start + window, ((dfG.iat[i+1, 3] - start) / 2) + start) # Start of the gene above position [i,3]
        dfG['midBody'] = dfG['midBody'].astype(int)
        dfG['lowB'] = dfG['lowB'].astype(int)
        dfG['upperB'] = dfG['upperB'].astype(int)
#         dfG['domainLEN'] = dfG['upperB'] - dfG['lowB']
        LBandUB = []
        LBandUB.append(dfG.lowB.values)
        LBandUB.append(dfG.upperB.values)
        outList.append(LBandUB)
    
    return(outList)


##### Adds midpoint in gene and upper and lower bound values to geneDFLtss based on size of window or nearby genes. Adding around geneTSS.
def addWindowTSS20(window):
    chrD = {'1': [1, 248956422],'2': [1, 242193529],'3': [1, 198295559],'4': [1, 190214555],'5': [1, 181538259],'6': [1, 170805979],'7': [1, 159345973],'8': [1, 145138636],'9': [1, 138394717],'10': [1, 133797422],'11': [1, 135086622],'12': [1, 133275309],'13': [1, 114364328],'14': [1, 107043718],'15': [1, 101991189],'16': [1, 90338345],'17': [1, 83257441],'18': [1, 80373285],'19': [1, 58617616],'20': [1, 64444167],'21': [1, 46709983],'22': [1, 50818468],'23': [1, 16569],'24': [1, 156040895],'25': [2781480, 56887902]} 
    geneWindowDFL = []
    for i in range(len(geneDFL)):
        df = geneDFL[i].copy()
        geneWindowDFL.append(df)
    chrIDNUM = 0
    for dfG in geneWindowDFL: # ordered by chromosome 1- MT,X,Y
        chrIDNUM += 1 
        dfIndexMax = len(dfG) - 1   
        for i, row in dfG.iterrows():
            start = row['start']
            
            dfG.at[i, 'midBody'] = (((row['length'] / 2)) + start) # Middle of gene body 
            if i == 0 or i == dfIndexMax: # Case for start and end genes
                if i == 0:
                    dfG.at[i, 'lowB'] = max(start - window, chrD[str(chrIDNUM)][0]) # chrD[chrIDNUM][0] = start val of the chromosome
                    dfG.at[i, 'upperB'] = min(start + window, ((dfG.iat[i+1, 3] - start) / 2) + start)
                else: # index max case 
                    dfG.at[i, 'lowB'] = max(start - window, start - ((start - dfG.iat[i-1, 3]) / 2)) # dfG.iat[i-1, 3] = start of gene below
                    dfG.at[i, 'upperB'] = min(start + window, chrD[str(chrIDNUM)][1]) # chrD[chrIDNUM][1] = end val of chromosome 
            else: # i > 0 or i < len(geneDF) - 1 
                dfG.at[i, 'lowB'] = max(start - window, start - ((start - dfG.iat[i-1, 3]) / 2)) # dfG.iat[i-1, 3] = start of gene below
                dfG.at[i, 'upperB'] = min(start + window, ((dfG.iat[i+1, 3] - start) / 2) + start) # Start of the gene above position [i,3]
        dfG['midBody'] = dfG['midBody'].astype(int)
        dfG['lowB'] = dfG['lowB'].astype(int)
        dfG['upperB'] = dfG['upperB'].astype(int)
        dfG['domainLEN'] = dfG['upperB'] - dfG['lowB']
    
    return(geneWindowDFL)


# Adds midpoint in gene and upper and lower bound values to geneDFL based on size of window or nearby genes. Adding around geneBODY.
def addWindowBODY(window):
    chrD = {'1': [1, 248956422],'2': [1, 242193529],'3': [1, 198295559],'4': [1, 190214555],'5': [1, 181538259],'6': [1, 170805979],'7': [1, 159345973],'8': [1, 145138636],'9': [1, 138394717],'10': [1, 133797422],'11': [1, 135086622],'12': [1, 133275309],'13': [1, 114364328],'14': [1, 107043718],'15': [1, 101991189],'16': [1, 90338345],'17': [1, 83257441],'18': [1, 80373285],'19': [1, 58617616],'20': [1, 64444167],'21': [1, 46709983],'22': [1, 50818468],'23': [1, 16569],'24': [1, 156040895],'25': [2781480, 56887902]} 
    chrIDNUM = 0
    geneWindowDFL = []
    for i in range(len(geneDFLBODY)):
        df = geneDFLBODY[i].copy()
        geneWindowDFL.append(df)
    
    outList = []
    for dfG in geneWindowDFL: # ordered by chromosome 1- MT,X,Y
        chrIDNUM += 1 
        dfIndexMax = len(dfG) - 1         
        for i, row in dfG.iterrows():
            start = row['start']
            end = row['end']
            dfG.at[i, 'midBody'] = (((row['length'] / 2)) + row['start']) # Middle of gene body 
            
            if i == 0 or i == dfIndexMax: # Case for start and end sites
                if i == 0:
                    dfG.at[i, 'lowB'] = max(start - window, chrD[str(chrIDNUM)][0]) # end of the gene below position
                    dfG.at[i, 'upperB'] = min(end + window, ((dfG.iat[i+1, 3] - end) / 2) + end)
                else: # index max case 
                    dfG.at[i, 'lowB'] = max(start - window, start - (abs(start - dfG.iat[i-1, 4]) / 2)) 
                    dfG.at[i, 'upperB'] = min(row['end'] + window, chrD[str(chrIDNUM)][1])
            else: # i > 0 or i < len(geneDF) - 1 
                dfG.at[i, 'lowB'] = max(start - window, start - (abs(start - dfG.iat[i-1, 4]) / 2)) # end of the gene below position [i,4]
                dfG.at[i, 'upperB'] = min(end + window, ((dfG.iat[i+1, 3] - end) / 2) + end) # Start of the gene above position [i,3]
        dfG['midBody'] = dfG['midBody'].astype(int)
        dfG['lowB'] = dfG['lowB'].astype(int)
        dfG['upperB'] = dfG['upperB'].astype(int)
#         dfG['domainLEN'] = dfG['upperB'] - dfG['lowB']
        LBandUB = []
        LBandUB.append(dfG.lowB.values)
        LBandUB.append(dfG.upperB.values)
        outList.append(LBandUB)
    
    return(outList)


# Adds midpoint in gene and upper and lower bound values to geneDFL based on size of window or nearby genes. Adding around geneBODY.
def addWindowBODY20(window):
    chrD = {'1': [1, 248956422],'2': [1, 242193529],'3': [1, 198295559],'4': [1, 190214555],'5': [1, 181538259],'6': [1, 170805979],'7': [1, 159345973],'8': [1, 145138636],'9': [1, 138394717],'10': [1, 133797422],'11': [1, 135086622],'12': [1, 133275309],'13': [1, 114364328],'14': [1, 107043718],'15': [1, 101991189],'16': [1, 90338345],'17': [1, 83257441],'18': [1, 80373285],'19': [1, 58617616],'20': [1, 64444167],'21': [1, 46709983],'22': [1, 50818468],'23': [1, 16569],'24': [1, 156040895],'25': [2781480, 56887902]} 
    chrIDNUM = 0
    geneWindowDFL = []
    for i in range(len(geneDFLBODY)):
        df = geneDFLBODY[i].copy()
        geneWindowDFL.append(df)
    
    for dfG in geneWindowDFL: # ordered by chromosome 1- MT,X,Y
        chrIDNUM += 1 
        dfIndexMax = len(dfG) - 1         
        for i, row in dfG.iterrows():
            start = row['start']
            end = row['end']
            
            dfG.at[i, 'midBody'] = (((row['length'] / 2)) + row['start']) # Middle of gene body 
            
            if i == 0 or i == dfIndexMax: # Case for start and end sites
                if i == 0:
                    dfG.at[i, 'lowB'] = max(start - window, chrD[str(chrIDNUM)][0]) # end of the gene below position
                    dfG.at[i, 'upperB'] = min(end + window, ((dfG.iat[i+1, 3] - end) / 2) + end)
                else: # index max case 
                    dfG.at[i, 'lowB'] = max(start - window, start - (abs(start - dfG.iat[i-1, 4]) / 2)) 
                    dfG.at[i, 'upperB'] = min(row['end'] + window, chrD[str(chrIDNUM)][1])
            else: # i > 0 or i < len(geneDF) - 1 
                dfG.at[i, 'lowB'] = max(start - window, start - (abs(start - dfG.iat[i-1, 4]) / 2)) # end of the gene below position [i,4]
                dfG.at[i, 'upperB'] = min(end + window, ((dfG.iat[i+1, 3] - end) / 2) + end) # Start of the gene above position [i,3]
        dfG['midBody'] = dfG['midBody'].astype(int)
        dfG['lowB'] = dfG['lowB'].astype(int)
        dfG['upperB'] = dfG['upperB'].astype(int)
        dfG['domainLEN'] = dfG['upperB'] - dfG['lowB']
    
    return(geneWindowDFL)


def findSMLdomain(window, method):
    if method == 'TSS':
        windowLens = addWindowTSS20(window)
    elif method == 'BODY':
        windowLens = addWindowBODY20(window)
    else:
        return 'Error in method'
    allVALS = []
    for df in windowLens:
        df['domainLEN'] = df['upperB'] - df['lowB']
        values = list(df.domainLEN.values)
        allVALS += values
    
    return sorted(allVALS)



# Creating list of all the entrezids 
entrezIDLtss = []
for geneDF in geneDFL:
    entrezIDLtss.append(geneDF.entrezid.values)


entrezIDLbody = []
for geneDF in geneDFLBODY:
    entrezIDLbody.append(geneDF.entrezid.values)


# CLOSEST Genes potentionally regulated by sites in window with respect to the genes body. k = 1 
def geneReadSites(bsL, geneWindow, method='TSS'): # geneWindow = windows lists
    #     geneWindow = addWindowBODY(window)
    chromosomeI = 0 # 0 == chromosome group 1          
    entrezidOutL = [] # Output list of entrez ids 
    if method == 'TSS':
        entrezIDL = entrezIDLtss
    if method == 'BODY':
        entrezIDL = entrezIDLbody
    if method != 'TSS'and method != 'BODY':
        return 'method must be TSS or BODY'
    
    for bounds in geneWindow:
        lowB = bounds[0]
        upperB = bounds[1]
        geneIDS = entrezIDL[chromosomeI]
        sitesL = bsL[chromosomeI]
        
        for site in sitesL:
            for i in range(len(lowB)):
                if site < lowB[i]:
                    lowB = lowB[i:] # Getting rid of the gene windows because it has been mapped. Genes at a lower positions are removed casuse they are small
                    upperB = upperB[i:]
                    geneIDS = geneIDS[i:]
                    break
                if lowB[i] <= site and upperB[i] > site:
                    entrezidOutL.append(geneIDS[i])
#                     lowB = lowB[i+1:] # Getting rid of the gene windows because it has been mapped. Genes at a lower positions are removed casuse they are small
#                     upperB = upperB[i+1:]
#                     geneIDS = geneIDS[i+1:]
#                     break # Done with the current site so break 
        chromosomeI += 1   
    return entrezidOutL


# Given a list of lists with go term and its odds-ratio,p-val return a list of just GO terms 
def getGOfromAnalysis(goAnalysis):
    goTermL = []
    for key in goAnalysis:
        goTermL.append(key)
    return goTermL



# Given original analysis and simulated analysis. Adds count to go term in original list if the random go term has a lower p-val and higher odds ratio 
def compareGOAnalysis(origAnalysis, simAnalysis):
    for k in simAnalysis:
        if (simAnalysis[k][1] < origAnalysis[k][1]): #  (simAnalysis[k][0] > origAnalysis[k][0]) and 
            origAnalysis[k][2] += 1
    return origAnalysis     



# converts analysis to list format and divides by nSim to get refabs p value
def convertAnalysistoFormat(analysis, nSim):
    analist = []
    for k in analysis.keys():
        analysis[k][2] = (analysis[k][2] + 1) / (nSim + 1) # + 1 to get rid of p vals of 0 for correction
        analist.append([k, analysis[k]])
    return analist



# Compares original analysis with simulated analysis for nSim
def simulation(GOlist, geneL, nSim, nSites, origAnalysis, method): # GOList is the new go list that we analyze through since we don't want whole GO list, geneL is windowDF
#     sitesL = nRandSitesSim(nSites, nSim)
    sitesL = nRandSitesSim(nSites, nSim)
    outputAnalysis = origAnalysis # We will keep updating this dictionary and return when all sims are done 
    
    for sim in range(nSim):
        mapped = geneReadSites(sitesL[sim], geneL, method)
        randAnalysis = conductAnalysisFAST(mapped, OntologyL=GOlist) # Updating OntologyL to make the calculations quicker 
        outputAnalysis = compareGOAnalysis(outputAnalysis, randAnalysis)
        
    return sorted(convertAnalysistoFormat(outputAnalysis, nSim), key = lambda x: x[1][1])


# First run that sets up for the simulation
def firstRun(userSitesL, window, method='TSS', simN=100): # userSites must be ordered by chromosome number (List of lists) # Default simN = 10000
    if method == 'TSS':
        geneL = addWindowTSS(window)
    elif method == 'BODY': 
        geneL = addWindowBODY(window)
    else:
        return 'Only TSS or BODY method allowed!'
    
    mappedGenes = geneReadSites(userSitesL, geneL, method)
    numMapped = len(mappedGenes)
    numsites = 0
    for l in userSitesL:
        numsites += len(l)
    
    goAnalysis = conductAnalysisFAST(mappedGenes)
    newGOlist = getGOfromAnalysis(goAnalysis)
    
    # estimate num sites to sample: numSitesSamp
    nBSL = nRandSitesSim(numsites,5)
    nPrime2 = 0 # Total num genes mapped to figure out nPrime value eventually for estimation
    for BSL in nBSL:
        nPrime2 += len(geneReadSites(BSL, geneL, method))
    nPrime = int(nPrime2 / 5)
    numSitesSamp = int((numMapped * numsites) / nPrime) # Num mapped genes * num sites inputted divided by nPrime 
    
    return simulation(newGOlist, geneL, simN, numSitesSamp, goAnalysis, method) # GeneL is list of genes with window bounds 


def adjP(dataSET):
    GOdataSet = dataSET
    goGrp = [[],[],[],[],[],[],[],[],[],[]]
    numAssGrp = [[1], [2], [3], [4], [5, 6], [7, 8], [9, 10, 11, 12], [13, 14, 15, 16, 17, 18], [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33], [34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 160, 161, 163, 164, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 182, 183, 185, 187, 188, 189, 190, 191, 193, 194, 196, 198, 200, 201, 202, 203, 204, 205, 206, 208, 210, 211, 212, 213, 214, 216, 217, 218, 219, 221, 225, 226, 229, 232, 234, 235, 236, 238, 239, 240, 242, 244, 249, 250, 251, 252, 254, 262, 263, 264, 265, 268, 272, 273, 278, 280, 282, 285, 290, 295, 297, 301, 307, 308, 310, 311, 312, 313, 315, 324, 329, 330, 334, 343, 344, 347, 354, 360, 362, 364, 366, 368, 371, 374, 375, 385, 387, 395, 397, 405, 409, 411, 413, 420, 424, 427, 431, 439, 441, 457, 458, 467, 474, 481, 483, 491, 496, 502, 504, 505, 524, 529, 537, 540, 547, 548, 552, 565, 578, 605, 610, 611, 667, 670, 689, 703, 706, 812, 816, 859, 927, 960, 973, 980, 1059, 1132, 1134, 1150, 1227, 1304, 1381, 1388, 1456, 1473, 1546, 1778, 1911, 1987, 2163, 2292, 3168, 3637, 4438, 4482, 5026, 5627, 9691]]
    numGO = [5184, 2895, 1743, 1241, 1551, 1018, 1219, 1028, 1020, 1275] # Number of GO terms in each group above.
    
    for go in GOdataSet:
        numAss = len(go2gene[go[0]])
        for i in range(len(numAssGrp)):
            if numAss in numAssGrp[i]:
                goGrp[i].append(go)
    
    pos = 0
    for GOdata in goGrp:
        fishersP = []
#         refabsP = []
        for l in GOdata:
            fishersP.append(l[1][1])
#             refabsP.append(l[1][2])
         
        if len(fishersP) < numGO[pos]:
            numAppend = numGO[pos] - len(fishersP)
            fishersP += [1] * numAppend
#             refabsP += [1] * numAppend
        pos += 1         
        
        reject, fishersPadj, alphacSidak, alphacBonf = padjust(fishersP, method='fdr_bh', is_sorted=False)
#         reject2, refabsPadj, alphacSidak2, alphacBonf2 = padjust(refabsP, method='fdr_bh', is_sorted=False)
        correctedfishers = []
#         correctedrefabs = []
        for i in range(len(GOdata)):
            correctedfishers.append(float(fishersPadj[i]))
#             correctedrefabs.append(float(refabsPadj[i]))
        
        for i in range(len(GOdata)):
            GOdata[i][1].append(correctedfishers[i])
#             GOdata[i][1].append(correctedrefabs[i])
    groupedGO = []
    for i in goGrp:
        for j in i:
            groupedGO.append(j)
    groupedGO = sorted(groupedGO, key = lambda x: x[1][1])
    refabsP = []
    for go in groupedGO:
        refabsP.append(go[1][2])
    reject2, refabsPadj, alphacSidak2, alphacBonf2 = padjust(refabsP, method='fdr_bh', is_sorted=False)
    
    for i in range(len(groupedGO)):
        groupedGO[i][1].append(refabsPadj[i])
    
    return groupedGO # go: [odds, fishersP, refabsP, correctedfishers, correctedrefabs]
t = time.time()
data1 = firstRun(nRandSites([10000, random.randint(100000)]), 100000, method='BODY', simN=2000)
analysed = adjP(data1)
e = time.time()
for l in analysed:
    print(l)
print('------')
print(len(analysed))
print(e-t)