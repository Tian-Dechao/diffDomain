from collections import defaultdict
import numpy as np
import straw
import TracyWidom as tw
import gzip
import pandas as pd

def finftypes(filepath):
    """
    filepath can be stdin, normal txt, and gz file
    """
    if(filepath == 'stdin'):
        fin = sys.stdin
    elif(filepath[-3:] == '.gz'):
        fin = gzip.open(filepath)
    else:
        fin = open(filepath)

    return fin

def loadtads(inpath, sep, chrnum, min_nbin, reso):
    tadb = []
    fin = finftypes(inpath)
    next(fin)
    for line in fin:
        #col = line.rstrip().split(sep)
        col = line.decode().rstrip().split('\t')
        tadb.append(col[:3])

    fin.close()

    tadb = pd.DataFrame(tadb)
    tadb[[1,2]] = tadb[[1,2]].astype(int)
    # format the chromosomal names
    if 'chr' in tadb.iloc[0,0]:
        tadb[0] = tadb[0].str.replace('chr', '')

    if chrnum != 'ALL':
        ind = tadb[0] == chrnum
        tadb = tadb[ind] #tadb = [_ for _ in tadb if _[0]==chrnum]

    # select domain with at least min_nbin
    nbins = np.ceil(tadb[2] /reso) - np.floor(tadb[1] / reso)
    tadb = tadb[nbins>=min_nbin]

    return tadb

def makewindow(x, reso):

    k = x / reso
    l = x % reso
    wins = [_*reso for _ in range(k) ]
    if(l>0):
        wins.append(k*reso)

    return wins

def makewindow2(start, end, reso):
    k1 = start / reso
    k2 = end / reso
    l2 = end % reso
    if l2 > 0:
        k2 += 1

    wins = [ _*reso for _ in range(k1, k2)]

    return wins

def compute_nbins(start, end, reso):
    domwin = makewindow2(start, end, reso)
    # create index for domwin
    nb = len(domwin)

    return nb


def contact_matrix_from_hic(chrn, start, end, reso, fhic, hicnorm='KR'):
    # handle some regions cannot be loaded by straw
    try:
        # find the bins for a domain
        domwin = makewindow2(start, end, reso)
        # create index for domwin
        nb = len(domwin)
        domwin_dict = {domwin[_]:_ for _ in range(nb) }
        # create empty matrix
        mat = np.empty(shape=(nb, nb))
        mat[:] = np.nan
        # find the edgelist from .hic
        region = "{0}:{1}:{2}".format(chrn, domwin[0], domwin[-1])
        el = straw.straw(hicnorm, fhic, region, region, 'BP', reso)

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
    except:
        mat = None

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
    for k, v in contactByDistance.iteritems():
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
    for k, val in contactByDistanceDiff.iteritems():
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
    for k, v in sd.iteritems():
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
#     print('wigner D:',Mat)
    Mat_ = pd.DataFrame(Mat)
#     Mat_.to_csv('WignerD_',sep='\t',header=None,index=None)
    try:
        w, v = np.linalg.eigh(Mat)
        # print('w:',np.abs(w))
        
        #U, s, V = np.linalg.svd(Mat, full_matrices=True)
        #w = np.abs(w)
        #w.sort()
        #lambdan = w[-1]
        lambdan = np.amax(np.abs(w))
#         print('\nlambdan:',lambdan)
        # normalize by sqrt(n)
        lambdan = np.power(nd, 2.0/3) * (lambdan -2)
#         print('theta:',lambdan)
        twdist = tw.TracyWidom(beta=1)
        p = 1 - twdist.cdf(lambdan)
    except:
        nd = np.nan
        lambdan = np.nan
        p = np.nan

    return nd, lambdan, p


def comp2domins_by_twtest(chrn, start, end, reso, hicnorm, fhic0, fhic1):
        mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm=hicnorm)
        # gu print the process to generate graphics
#         print('\nmat0:',mat0)
        mat0_ = pd.DataFrame(mat0)
#         mat0_.to_csv('mat0_',sep='\t',header=None,index=None)
        mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm=hicnorm)
#         print('\nmat1:',mat1)
        mat1_ = pd.DataFrame(mat1)
#         mat1_.to_csv('mat1_',sep='\t',header=None,index=None)
        if not mat0 is None and not mat1 is None:
            # remove rows that have more than half np.nan
            nbins = compute_nbins(start, end, reso)
            ind0 = np.sum(np.isnan(mat0), axis=0) < nbins * 0.5
            ind1 = np.sum(np.isnan(mat1), axis=0) < nbins * 0.5
            ind = ind0 & ind1
            if ind.sum() >=10:
                indarray = np.array(ind)
                mat0rmna = mat0[indarray, :]
                mat0rmna = mat0rmna[:, indarray]
                mat1rmna = mat1[indarray, :]
                mat1rmna = mat1rmna[:, indarray]

                # compute the differnece matrix
                Diffmat = mat0rmna / mat1rmna
#                 print('\nDiffmat:',Diffmat)
                Diffmat_ = pd.DataFrame(Diffmat)
                Diffmat_.to_csv('Diffmat_',sep='\t',header=None,index=None)
                
                #print Diffmat
                ind1 = np.sum(np.isnan(Diffmat), axis=0)
                #print ind1

                Diffmatnorm = normDiffbyMeanSD(D=Diffmat)
                #print Diffmatnorm
#                 print('\nDiffmatnorm:',Diffmatnorm)
                Diffmatnorm_ = pd.DataFrame(Diffmatnorm)
#                 Diffmatnorm_.to_csv('Diffmatnorm_',sep='\t',header=None,index=None)

                result = twtest_formula(Diffmatnorm)

                domname = 'chr%s:%s-%s' % (chrn, start, end)
                result = [chrn, start, end, domname, result[1], result[2], result[0]]
            else:
                domname = 'chr%s:%s-%s' % (chrn, start, end)
                result = [chrn, start, end, domname, np.nan, np.nan, np.nan]

        else:
            domname = 'chr%s:%s-%s' % (chrn, start, end)
            result = [chrn, start, end, domname, np.nan, np.nan, np.nan]

        return result

