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
import random
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
    if start >= end :
        print('The end of TAD should be bigger than its start')
        return None
    # k1 = int(np.floor(start / reso))
    # k2 = int(np.ceil(end / reso))
    k1 = start // reso

    k2 = end // reso 

    wins = [ _*reso for _ in range(k1, k2+1)]
    return wins


def compute_nbins(start, end, reso):
    domwin = makewindow2(start, end, reso)
    # create index for domwin
    nb = len(domwin)
    return nb


# def load_cool(fhic,reso):
#     import subprocess as sp
#     # translate the cool file into h5
#     name = fhic.replace('.cool','').replace('.mcool','')

#     cmd1 = f'hicConvertFormat -m {fhic} --inputFormat cool \
#         --outputFormat ginteractions -o {name} --reso {reso}'
#     cp1 = sp.run(cmd1,shell=True,capture_output=True,encoding='utf-8')

#     if cp1.returncode == 0:
#         # sp.run(f'rm {name_h5}',shell=True,capture_output=False,encoding='utf-8')
#         return name+'.tsv'
#     else:
#         raise(Exception(cp1.stderr))


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
    region2 = "{0}:{1}-{2}".format(chrn, domwin[0], domwin[-1])
    print(region)

    def load_h5(fhic,chrn, region):
            import h5py
            hf = h5py.File(fhic,'r')
            el = hf.get(region)
            if el is None:
                print(f'{region} is not found in {fhic}')
                return None
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
    if fhic.endswith('.hic'):
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
    # elif fhic[-3:] == '.h5' or fhic.split('.')[-1] == 'hdf5':
        
    #     return load_h5(fhic,region)
    
    
    elif fhic.endswith('cool') : # .cool or .mcool
        import cooler
        import subprocess as sp
        # Normalization

        suffix = f'_{int(reso/1000)}k_{hicnorm}.cool'
        if fhic.endswith('.cool'):
            hicfile = fhic
            hic_norm = fhic.replace('.cool',suffix)
        elif fhic.endswith('.mcool'):
            hicfile = f'{fhic}::resolutions/{reso}'
            hic_norm = fhic.replace('.mcool',suffix)
                
        try:
            c = cooler.Cooler(f'{hic_norm}')

        except FileNotFoundError:
            # normalization
            # KR or ICE
            cmd = f'hicConvertFormat -m {hicfile} --inputFormat cool \
                --outputFormat cool --correction_name {hicnorm} \
                --reso {reso} -o {hic_norm}'
            cp1 = sp.run(cmd,shell=True,capture_output=True,encoding='utf-8')
            
            if cp1.returncode == 0: # True
                c = cooler.Cooler(f'{hic_norm}')
            else: # False
                raise(Exception(cp1.stderr))

        mat = c.matrix(balance=False).fetch(region2)
        return mat
        
    # load a sparse matrix with three columns
    else:
        print('Trying to sparse the files as three columns by "\t" ')
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
            # remove rows that have more than half np.nan
            nbins = compute_nbins(start, end, reso)
            ind0 = np.sum(np.isnan(mat0), axis=0) < nbins * (1-float(f))
            ind1 = np.sum(np.isnan(mat1), axis=0) < nbins * (1-float(f))
            ind = ind0 & ind1
            
            # change made to suit SVs - replace the nan with 1 in mat1 and mat0 not changed.
            #ind0 = np.sum(np.isnan(mat1), axis=0) == 0
            #ind1 = np.sum(np.isnan(mat1), axis=1) == 0
            #indc = np.where(ind0 == False)
            #indr = np.where(ind1 == False)
            
            # change mat1 to suit SVs  
            #mat1_changed = mat1[indr,indc] = 1
            #X = [random.expovariate(5) for i in range(nbins)]
            #for i in range(len(indr[0])):
            #    x = random.choice(X)
            #    mat1[indr[0][i],indc[0][i]] = 1+x
            #mat1_changed = mat1
           
            if ind.sum() >= min_nbin:
            #if ind.sum() > 0:
                indarray = np.array(ind)
                mat0rmna = mat0[indarray, :]
                mat0rmna = mat0rmna[:, indarray]
                mat1rmna = mat1[indarray, :]
                mat1rmna = mat1rmna[:, indarray]
        
                # compute the differnece matrix
                Diffmat = mat0rmna / mat1rmna
                #Diffmat = mat0 / mat1
                #Diffmat = mat0 / mat1_changed

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
                print('The matrix is too spase or too small at this resolution to be calculated !')
                result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
        else:
            domname = '%s:%s-%s' % (chrn, start, end)
            result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
            print('The matrixes failed to be load!')

        return result




def comp2domins_by_twtest_changed(chrn, start, end, reso, hicnorm, fhic0, fhic1, min_nbin, f):
        mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm)
        mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm)
        #print(mat0)
        #print(mat1)
        #print(sum(np.isnan(mat1)))
        #print(np.isnan(mat1).sum())
        #print(len(mat1)*len(mat1))

        if not mat0 is None and not mat1 is None:
            # rlmove rows that have more than half np.nan
            nbins = compute_nbins(start, end, reso)
            ind0 = np.sum(np.isnan(mat0), axis=0) < nbins * (1-float(f))
            ind1 = np.sum(np.isnan(mat1), axis=0) < nbins * (1-float(f))
            ind = ind0 & ind1
            
            row,col = [],[]
            for i in range(nbins):
                for j in range(nbins):
                    if np.isnan(mat1[i][j]):
                        row.append(i)
                        col.append(j)

            #X = [random.expovariate(5) for i in range(nbins)]

            for i in range(len(row)):
                #x = random.choice(X)
                mat1[row[i],col[i]] = 1

            if ind.sum() >= min_nbin:
                #indarray = np.array(ind)
                #mat0rmna = mat0[indarray, :]
                #mat0rmna = mat0rmna[:, indarray]
                #mat1rmna = mat1[indarray, :]
                #mat1rmna = mat1rmna[:, indarray]

                # compute the differnece matrix
                #Diffmat = mat0rmna / mat1rmna
                Diffmat = mat0 / mat1
                
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

## ADD Functions to get svs data and make adjustment to diffdomain

def SVs_from_bedfile(fsvs):
    SV_data = pd.read_table(fsvs,header = None)
    SV_data = SV_data.rename(columns={0:'chrom',1:'breakpoint1',2:'breakpoint2',3:'type'})
    SV_data['chrom'] = [str(x) for x in SV_data['chrom']]
    return SV_data

def find_large_missing_blocks(matrix, threshold):
    # threshold is a ratio,if nan is larger than threshold of the whole matrix, then find the missing_blocks.
    rows, cols = len(matrix), len(matrix[0])
    missing_blocks = []
    visited = set()  # To keep track of visited positions
    
    total_elements = rows * cols
    threshold_nan_count = int(threshold * total_elements)
    
    def is_visited(i, j, block_rows, block_cols):
        return any((i + k, j + l) in visited for k in range(block_rows) for l in range(block_cols))
    
    for i in range(rows):
        for j in range(cols):
            if (i, j) in visited:  # Skip if position has been visited
                continue
            
            # Check for missing value at position (i, j)
            if np.isnan(matrix[i][j]):
                block_rows = 0
                block_cols = 0
                
                # Check how large the missing block is in rows
                while i + block_rows < rows and all(np.isnan(matrix[i + k][j]) for k in range(block_rows + 1)):
                    block_rows += 1
                
                # Check how large the missing block is in columns
                while j + block_cols < cols and all(np.isnan(matrix[i][j + l]) for l in range(block_cols + 1)):
                    block_cols += 1
                    
                block_nan_count = block_rows * block_cols
                
                if block_nan_count >= threshold_nan_count:
                    if is_visited(i, j, block_rows, block_cols):
                        continue
                    
                    missing_blocks.append(((i, j), block_rows, block_cols))
                    
                    # Mark positions as visited
                    for k in range(block_rows):
                        for l in range(block_cols):
                            visited.add((i + k, j + l))
    
    # Sort the missing blocks by area (number of NaN values) in descending order
    missing_blocks.sort(key=lambda block: block[1] * block[2], reverse=True)
    
    return missing_blocks

def make_symmetric(matrix):
    rows, cols = matrix.shape
    for i in range(rows):
        for j in range(i+1, cols):
            matrix[j, i] = matrix[i, j]
    
    return matrix

def fill_missing_blocks(matrix, missing_blocks, fill_value = 1):
    filled_matrix = np.copy(matrix)
    
    for block_info in missing_blocks:
        (i, j), block_rows, block_cols = block_info
        
        for k in range(block_rows):
            for l in range(block_cols):
                filled_matrix[i + k][j + l] = fill_value
        filled_matrix = make_symmetric(filled_matrix)
    
    return filled_matrix

def missing_block_with_SVs(blockstart,blockend,SVsdataset,chrn):
    df = SVsdataset[SVsdataset['type'] == '+-']
    df = df[df['chrom'] == chrn]
    related_num = 0
    for i in range(df.shape[0]):
        df = df[((df['breakpoint1'] >= blockstart) & (df['breakpoint1'] <= blockend))
                |((df['breakpoint2'] >= blockstart) & (df['breakpoint2'] <= blockend))
                |((df['breakpoint1'] <= blockstart) & (df['breakpoint2'] >= blockend))
                ]
    if df.shape[0] >0:
        return True,
                
    return False


def is_large_missing_block(missing_blocks, total_elements, threshold = 0.8):
    total_missing_elements = sum(block_rows * block_cols for (_, _), block_rows, block_cols in missing_blocks)
    if (total_missing_elements / total_elements) > threshold:
        return True
    else:
        return False


def comp2domins_by_twtest_withSVs(chrn, start, end, reso, hicnorm, fhic0, fhic1,min_nbin,f,min_ratio = 0.2,SVs_file_condition2 = 'nosvsfile',fsvs_condition2 = None):
    # SVs_file_condition2 and fsvs_condition2 are all under condition2 
    mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm)
    mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm)
    domname = '%s:%s-%s' % (chrn, start, end)
    filled_matrix = mat1

    if not mat0 is None and not mat1 is None:
        if SVs_file_condition2 == 'svsfile':
            SVs = SVs_from_bedfile(fsvs_condition2)
            # if the large loss in the mat1 resulted from deletion type Svs in conditon2(fhic1)
            if missing_block_with_SVs(start,end,SVs,chrn):
                # imputation1 SV related
                missing_blocks = find_large_missing_blocks(matrix = mat1, threshold =min_ratio)
                if missing_blocks:
                    filled_matrix = fill_missing_blocks(matrix = mat1, missing_blocks = missing_blocks)
        else:
            missing_blocks = find_large_missing_blocks(matrix = mat1, threshold = min_ratio)
            if missing_blocks:
                filled_matrix = fill_missing_blocks(matrix = mat1, missing_blocks = missing_blocks)
                
        #mat1_changed = filled_matrix 
        nbins = compute_nbins(start, end, reso)
        if is_large_missing_block(missing_blocks = missing_blocks,total_elements=nbins*nbins):
            result = [chrn, start, end, domname, np.nan, 2.22e-16, np.nan]
        else:
            # if rows that have more than half np.nan
            ind0 = np.sum(np.isnan(mat0), axis=0) < nbins * (1-float(f))
            ind1 = np.sum(np.isnan(filled_matrix), axis=0) < nbins * (1-float(f))
            ind = ind0 & ind1       
            if ind.sum() >= min_nbin:
                Diffmat = mat0 / filled_matrix
                #print(Diffmat)
                # imputation2
                Diffmatnorm = normDiffbyMeanSD(D=Diffmat)

                result = twtest_formula(Diffmatnorm)

                result = [chrn, start, end, domname, result[1], result[2], result[0]]
            else:
                print('The length of this TAD is too small at this resolution to be calculated !')
                result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
    else:
        result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
        print('The matrix is too spase at this resolution to be calculated !')

    return result

