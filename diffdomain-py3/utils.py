#-*- coding : utf-8-*-
# coding:unicode_escape

from collections import defaultdict
import numpy as np
import hicstraw as straw
import TracyWidom as tw
import gzip
import pandas as pd
import os
import sys
import pandas as pd
from docopt import docopt
from multiprocessing import Pool
import warnings
warnings.filterwarnings('ignore')


def finftypes(filepath):
    """
    filepath can be stdin, normal txt, and gz file
    """
    if(filepath == 'stdin'):
        fin = sys.stdin
    elif(filepath[-3:] == '.gz'):
        fin = gzip.open(filepath)
    else:
        fin = open(filepath,'r')
    return fin


def loadtads(inpath, sep, chrnum, reso, min_nbin): # SyntaxError:non-default argument follows dafault argument
    tadb = []
    fin = finftypes(inpath)
    next(fin)
    for line in fin:
        #col = line.rstrip().split(sep)
        col = line.rstrip().split(sep)
        tadb.append(col[:3])
    fin.close()

    tadb = pd.DataFrame(tadb)
    tadb.iloc[:,1:3] = tadb.iloc[:,1:3].astype(int)

#     format the chromosomal names
#     if 'chr' in tadb.iloc[0,0]:
#         tadb[0] = tadb[0].str.replace('chr', '')

    if chrnum != 'ALL':
        ind = tadb[0] == chrnum
        tadb = tadb[ind] #tadb = [_ for _ in tadb if _[0]==chrnum]

    # select domain with at least min_nbin
    nbins = np.ceil(tadb[2] /reso) - np.floor(tadb[1] / reso) 
    tadb = tadb[nbins>=min_nbin] #huadm

    return tadb

def makewindow(x, reso):
    
    k = int(np.ceil(x / reso))
    
    wins = [_*reso for _ in range(k+1) ]

    return wins

def makewindow2(start, end, reso):
    if start > end :
        print('The end of TAD should be bigger than its start')
        return None
    k1 = int(np.floor(start / reso))
    k2 = int(np.ceil(end / reso))

    wins = [ _*reso for _ in range(k1, k2)]
    return wins


def compute_nbins(start, end, reso):
    domwin = makewindow2(start, end, reso)
    # create index for domwin
    nb = len(domwin)
    return nb




def contact_matrix_from_hic(chrn, start, end, reso, fhic, hicnorm):
    # handle some regions cannot be loaded by straw
    # find the bins for a domain
    domwin = makewindow2(start, end, reso)
    # create index for domwin
    nb = len(domwin)
    domwin_dict = {domwin[_]:_ for _ in range(nb)}
    # create empty matrix
    mat = np.empty(shape=(nb, nb))
    mat[:] = np.nan
    # find the edgelist from .hic
    region = "{0}:{1}:{2}".format(chrn, domwin[0], domwin[-1])
    if fhic[-4:] == '.hic':
        try:
            el = straw.straw('observed',hicnorm, fhic, region, region, 'BP', reso)
            for i in range(len(el)):
                bin0 = el[i].binX
                bin1 = el[i].binY
                k=domwin_dict[bin0]
                l=domwin_dict[bin1]
                if k == l:
                    mat[k, l] = el[i].counts
                else:
                    mat[k, l] = el[i].counts
                    mat[l, k] = el[i].counts

            return mat
        except IOError:
            print("Sorry,{%s} does't exist."% fhic )
            mat = None
            return mat

    ## This module is just written to test the .h5 sample data (chr1_50M_GM12878.h5, chr1_50M_K562.h5)
    elif fhic[-3:] == '.h5':
        import h5py
        hf = h5py.File(fhic,'r')
        regions = "{0}:{1}:{2}".format(chrn, start, end)
        el = hf.get(regions)
        if el is None:
            print(' {regions} not found in {fhic}')
            return
        el = np.array(el) # transform to numpy format

        for i in range(el.shape[0]):
            bin0 = el[0][i]
            bin1 = el[1][i]
            k=domwin_dict[bin0]
            l=domwin_dict[bin1]
            if k == l:
                mat[k, l] = el[2][i]
            else:
                mat[k, l] = el[2][i]
                mat[l, k] = el[2][i]

        return mat

    # load a sparse matrix with three columns
    else:
        data = pd.read_table(fhic,sep='\t')
        el =[[],[],[]]
        el[0] = data[data.columns[0]].tolist()
        el[1] = data[data.columns[1]].tolist()
        el[2] = data[data.columns[2]].tolist()
        
        for i in range(len(el[2])):
            bin0 = el[0][i]
            bin1 = el[1][i]
            k=domwin_dict[bin0]
            l=domwin_dict[bin1]
            if k == l:
                mat[k, l] = el[2][i]
            else:
                mat[k, l] = el[2][i]
                mat[l, k] = el[2][i]

        return mat


def extractKdiagonalCsrMatrix(spsCsrMat):
    nbin = spsCsrMat.shape[0]
    nonZeroIndex = spsCsrMat.nonzero()
    contactByDistance = defaultdict(list)
    for _ in range(len(nonZeroIndex[0])):
        rowIndex = nonZeroIndex[0][_]
        colIndex = nonZeroIndex[1][_]
        if colIndex >= rowIndex:
            dist = colIndex - rowIndex
            contactByDistance[dist].append(spsCsrMat[rowIndex, colIndex])

    numNonzero = []
    median = []
    for k, v in contactByDistance.items():
        numNonzero.append(len(v))
        median.append(np.nanmedian(v))


    #print numNonzero[:10], numNonzero[-10:]
    #print median[:10], median[-10:]

    return contactByDistance

def normDiffbyMeanSD(D):
    # log transformation
    D = np.log(D)
    contactByDistanceDiff = extractKdiagonalCsrMatrix(D)
    # imputation of nan and inf
    # part0: get the median and maximum for each off-diagonal
    a, b = defaultdict(list), defaultdict(list)
    m, sd = defaultdict(float), defaultdict(float)
    for k, val in contactByDistanceDiff.items():
        indnan = np.isnan(val)
        indinf = np.isinf(val)
        if np.any(indnan) and np.any(indinf):
            val = np.array(val)
            ind = np.logical_or(indnan, indinf)
            val1 = val[np.logical_not(ind)]
            a[k] = np.median(val1)
            b[k] = np.max(val1)
            m[k] = np.mean(val1)
            sd[k] = np.std(val1)
        elif np.any(indnan) and not np.any(indinf):
            val = np.array(val)
            val1 = val[np.logical_not(indnan)]
            a[k] = np.median(val1)
            m[k] = np.mean(val1)
            sd[k] = np.std(val1)
        elif not np.any(indnan) and np.any(indinf):
            val = np.array(val)
            val1 = val[np.logical_not(indinf)]
            b[k] = np.max(val1)
            m[k] = np.mean(val1)
            sd[k] = np.std(val1)
        else:
            m[k] = np.mean(val)
            sd[k] = np.std(val)

    #part1: impute nan and inf
    indnan = np.isnan(D)
    indnan.astype(int)
    indr, indc = np.nonzero(indnan)
    for _ in range(len(indr)):
        k = abs(indr[_] - indc[_])
        D[indr[_], indc[_]] = a[k]

    indinf = np.isinf(D)
    indr, indc = np.nonzero(indinf)
    for _ in range(len(indr)):
        k = abs(indr[_] - indc[_])
        D[indr[_], indc[_]] = b[k]

    # subtract mean and dividing by sd
    for k, v in sd.items():
        if np.isnan(v) or v == 0:
            sd[k] = 1

    nd = D.shape[0]
    for i in range(nd):
        for j in range(nd):
            k = abs(i-j)
            D[i, j] = (D[i,j] - m[k]) / sd[k]


    return D

def twtest_formula(Mat):
    nd = Mat.shape[0]
    # normalize it to a Wigner matrix
    Mat = Mat / np.sqrt(nd)
    try:
        w, v = np.linalg.eigh(Mat)
        #U, s, V = np.linalg.svd(Mat, full_matrices=True)
        #w = np.abs(w)
        #w.sort()
        #lambdan = w[-1]
        lambdan = np.amax(np.abs(w))
        # normalize by sqrt(n)
        lambdan = np.power(nd, 2.0/3) * (lambdan -2)

        twdist = tw.TracyWidom(beta=1)
        p = 1 - twdist.cdf(lambdan)
    except:
        nd = np.nan
        lambdan = np.nan
        p = np.nan

    return nd, lambdan, p

def visualization(chrn, start, end, reso, hicnorm, fhic0, fhic1):

        mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm)
        mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm)
        return mat0,mat1


def comp2domins_by_twtest(chrn, start, end, reso, hicnorm, fhic0, fhic1, min_nbin, f):
        mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm)
        mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm)

        if not mat0 is None and not mat1 is None:
            # rlmove rows that have more than half np.nan
            nbins = compute_nbins(start, end, reso)
            ind0 = np.sum(np.isnan(mat0), axis=0) < nbins * (1-float(f))
            ind1 = np.sum(np.isnan(mat1), axis=0) < nbins * (1-float(f))
            ind = ind0 & ind1
            if ind.sum() >= min_nbin:
                indarray = np.array(ind)
                mat0rmna = mat0[indarray, :]
                mat0rmna = mat0rmna[:, indarray]
                mat1rmna = mat1[indarray, :]
                mat1rmna = mat1rmna[:, indarray]

                # compute the differnece matrix
                Diffmat = mat0rmna / mat1rmna
                #print Diffmat
                ind1 = np.sum(np.isnan(Diffmat), axis=0)
                #print ind1

                Diffmatnorm = normDiffbyMeanSD(D=Diffmat)
                #print Diffmatnorm

                result = twtest_formula(Diffmatnorm)

                domname = '%s:%s-%s' % (chrn, start, end)
                result = [chrn, start, end, domname, result[1], result[2], result[0]]
            else:
                domname = '%s:%s-%s' % (chrn, start, end)
                print('The length of this TAD is too small at this resolution to be calculated !')
                result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
        else:
            domname = '%s:%s-%s' % (chrn, start, end)
            result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
            print('The matrix is too spase at this resolution to be calculated !')

        return result



