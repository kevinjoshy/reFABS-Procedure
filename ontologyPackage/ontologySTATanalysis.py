import time 
import scipy.stats as stats
import numpy as np 
from multiprocessing import Pool
import json
from operator import itemgetter


# This is the old way to get data. New is stored in a txt file. gene2go and go2gene.
# from goatools.base import download_ncbi_associations # Only needed to download: 'gene2go'
# from goatools.anno.genetogo_reader import Gene2GoReader
# fin_gene2go = download_ncbi_associations() # Downloading associations file: 'gene2go'
# fin_gene2go = 'gene2go'
# objanno = Gene2GoReader(fin_gene2go, taxids=[9606]) # Reading human associasions 
# ns2assoc = objanno.get_ns2assc() # Sorted by NS 
# assoc = objanno.get_id2gos_nss() # Not sorted by NS 


with open('go2gene.txt', 'r') as file:
     go2gene = json.loads(file.read())
        
with open('genes2Ontology.txt', 'r') as file:
     genes2go = json.loads(file.read())
genes2go = {int(k):(v) for k,v in genes2go.items()} # Converts key back to int 



genel = [51150,
 6339,
 3104,
 54206,
 5195,
 57035,
 23585,
 378807,
 9064,
 199870,
 4146,
 101929464,
 64769,
 64064,
 128209,
 55929,
 5052,
 55624,
 127435,
 101927034,
 3716,
 6121,
 494115,
 23285,
 55599,
 2703,
 653513,
 56957,
 51177,
 23126,
 100191040,
 284486,
 6281,
 6274,
 127579,
 55974,
 84283,
 164118,
 81875,
 477,
 844,
 2214,
 127933,
 60676,
 730102,
 26092,
 23179,
 7175,
 22874,
 79098,
 101930114,
 29920,
 56997,
 100287814,
 81466,
 401994,
 151354,
 343930,
 1788,
 11117,
 6683,
 6546,
 805,
 644093,
 3344,
 55704,
 57223,
 5966,
 3099,
 25885,
 101928371,
 7850,
 6549,
 400997,
 26525,
 54520,
 8685,
 80097,
 10752,
 375323,
 27258,
 8292,
 9975,
 27303,
 7048,
 114884,
 64689,
 199223,
 10289,
 7123,
 83598,
 646424,
 11334,
 1795,
 57406,
 131177,
 55079,
 389136,
 257144,
 100874246,
 285282,
 165631,
 8971,
 7200,
 255403,
 5158,
 7884,
 118,
 6002,
 94031,
 1816,
 9079,
 55300,
 100506444,
 79644,
 55236,
 9407,
 1446,
 401138,
 132660,
 113510,
 345275,
 93627,
 64850,
 2169,
 4085,
 166378,
 170690,
 1007,
 4883,
 114899,
 253260,
 5019,
 493869,
 54505,
 643155,
 2297,
 411,
 51752,
 101927023,
 10455,
 285782,
 10048,
 6310,
 7172,
 221656,
 3018,
 10385,
 696,
 170954,
 135656,
 221527,
 10665,
 6890,
 6046,
 2289,
 1026,
 10695,
 101926934,
 167,
 83741,
 56479,
 100507381,
 441161,
 23033,
 5238,
 81491,
 56975,
 340260,
 541472,
 10643,
 9805,
 223075,
 273,
 136647,
 58498,
 107,
 641977,
 100130849,
 285905,
 83451,
 7462,
 100505767,
 55971,
 23660,
 2055,
 1667,
 504180,
 286097,
 137492,
 9108,
 64760,
 23221,
 4017,
 64641,
 157574,
 1846,
 6867,
 7336,
 96764,
 3444,
 9373,
 10850,
 26206,
 881,
 5079,
 9925,
 169693,
 116224,
 9414,
 5125,
 389762,
 389763,
 257019,
 23560,
 55526,
 101928834,
 1326,
 79741,
 79290,
 240,
 83849,
 219623,
 5552,
 6865,
 10367,
 283008,
 6050,
 143503,
 340980,
 53840,
 120796,
 387751,
 338322,
 10819,
 79733,
 2188,
 8539,
 55761,
 60529,
 4038,
 5702,
 79841,
 440040,
 390155,
 81029,
 50813,
 1822,
 8076,
 64581,
 3821,
 79370,
 341350,
 84070,
 613227,
 6778,
 3458,
 55508,
 3747,
 196477,
 1634,
 5250,
 429,
 29902,
 160760,
 6910,
 79794,
 283455,
 8843,
 80324,
 23141,
 7637,
 221143,
 9818,
 90627,
 160857,
 7178,
 55213,
 81624,
 100874196,
 122742,
 328,
 390443,
 51222,
 4323,
 10278,
 10548,
 283629,
 84312,
 390539,
 283768,
 9824,
 6263,
 89978,
 5888,
 643338,
 26015,
 8773,
 197131,
 1854,
 89941,
 388199,
 9894,
 5310,
 899,
 527,
 102724927,
 7023,
 8651,
 730013,
 51704,
 54988,
 112479,
 1039,
 124446,
 6299,
 23322,
 51716,
 2806,
 8883,
 10725,
 23035,
 5713,
 497190,
 22879,
 57687,
 101927839,
 162514,
 1973,
 3000,
 6844,
 84314,
 5376,
 55090,
 125170,
 10750,
 3965,
 146861,
 246176,
 6871,
 100132476,
 728318,
 83899,
 100507608,
 2118,
 51629,
 2896,
 2535,
 79170,
 5245,
 5164,
 84687,
 3131,
 56155,
 146771,
 6329,
 3384,
 2186,
 353174,
 400627,
 124565,
 6182,
 116729,
 2837,
 11031,
 143471,
 57536,
 4200,
 51320,
 23239,
 5055,
 1991,
 404665,
 729359,
 8192,
 79230,
 5300,
 93145,
 53637,
 10498,
 90580,
 388507,
 84292,
 8677,
 10755,
 30817,
 80726,
 163227,
 79156,
 9745,
 58510,
 84063,
 148268,
 126433,
 9149,
 284323,
 126526,
 90273,
 5452,
 101928063,
 24139,
 1175,
 84215,
 6822,
 106144526,
 342918,
 3817,
 7728,
 162967,
 163033,
 23619,
 374928,
 54807,
 55321,
 7053,
 57506,
 9770,
 101929225,
 57325,
 440757,
 1471,
 92086,
 140706,
 51654,
 128866,
 84557,
 101926987,
 51526,
 55861,
 391253,
 57580,
 140731,
 8771,
 64092,
 29980,
 7267,
 100874006,
 10785,
 1292,
 51807,
 728226,
 128989,
 440822,
 157,
 8563,
 550631,
 253143,
 51493,
 4627,
 5816,
 23466,
 9145,
 644186,
 170062,
 8573,
 55634,
 79917,
 158835,
 56159,
 4841,
 27160,
 101928259,
 105373251,
 389874,
 51132,
 9452,
 643486,
 9075,
 11043,
 1288,
 2182,
 1641,
 63932,
 9016,
 55796,
 100131434,
 2564,
 8260,
 140032,
 441543]


# Returns all GO terms associated with geneL input. Has no duplicates 
def getOntologyID(geneIDL):
    GOList = []
    for gene in geneIDL:
        if gene in genes2go.keys():
            terms = genes2go[gene]
            GOList += terms 
    return list(set(GOList)) 


# Finds common between two lists and returns its length 0 if none 
def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if (a_set & b_set): 
        return len(a_set & b_set) 
    else: 
        return 0 


# Returns genes that have go term 
# This isn't needed anymore! slow 
def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item in listOfItems:
        if valueToFind in list(item[1]):
            listOfKeys.append(item[0])
    return listOfKeys


# Conducts GO analysis to find odds ratio and p-value for each go term. 
def conductAnalysis(geneL, OntologyL = None):
    geneNum = 21352 - 58 # Num of genes intotal in DF - the duplicates
    mappedGeneNum = len(geneL)
    if OntologyL == None:
        OntologyL = getOntologyID(geneL)
    goTermAnalysis = {}
    
    for go in OntologyL:
        genesAssocGO = go2gene[go] 
        A = common_member(genesAssocGO, geneL) # Num of common members between two lists  
        B = mappedGeneNum - A  # Num of genes in mapped list but not in go gene list 
        C = len(genesAssocGO) - A # Num genes in go gene list but not in mapped list 
        D = geneNum - (A+B+C) # Num total genes in genome - (A+B+C) = D 
        oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
        out = [oddsratio, pvalue, 0]
        goTermAnalysis.update({go:out})
    
#     goTermAnalysis = {k: v for k, v in goTermAnalysis.items() if v[1] <= 1} # less than 0.05
    return goTermAnalysis


def analysis(inputs):
    goTerms = inputs[0]
    mappedGeneNum = inputs[1]
    geneL=inputs[2]
    geneNum = 21294 
    goTermAnalysis = {}
    
    for go in goTerms:
        genesAssocGO = go2gene[go] 
        A = common_member(genesAssocGO, geneL) # Num of common members between two lists  
        B = mappedGeneNum - A  # Num of genes in mapped list but not in go gene list 
        C = len(genesAssocGO) - A # Num genes in go gene list but not in mapped list 
        D = geneNum - (A+B+C) # Num total genes in genome - (A+B+C) = D 
        oddsratio, pvalue = stats.fisher_exact([[A, B], [C, D]])
        out = [oddsratio, pvalue, 0]
        goTermAnalysis.update({go:out})
    return goTermAnalysis


# Conducts GO analysis to find odds ratio and p-value for each go term. 
def conductAnalysisFAST(geneL, OntologyL = None):
    geneNum = 21352 - 58 # Num of genes intotal in DF - the duplicates
    mappedGeneNum = len(geneL)
    if OntologyL == None:
        OntologyL = getOntologyID(geneL)
    goTermAnalysis = {}
    
    x = int(len(OntologyL) / 5)
    x1 = OntologyL[:x]
    x2 = OntologyL[x:x*2]
    x3 = OntologyL[x*2:x*3]
    x4 = OntologyL[x*3:x*4]
    x5 = OntologyL[x*4:]
    pool = Pool(5)
    result = pool.map(analysis, [[x1,mappedGeneNum,geneL],[x2,mappedGeneNum,geneL],[x3,mappedGeneNum,geneL],[x4,mappedGeneNum,geneL],[x5,mappedGeneNum,geneL]])
    pool.close()
    
    for res in result:
        goTermAnalysis.update(res)
    
#     goTermAnalysis = {k: v for k, v in goTermAnalysis.items() if v[1] <= 1} # less than 0.05
    return goTermAnalysis

