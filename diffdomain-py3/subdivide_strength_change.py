import numpy as np
import pandas as pd
import os
import hicstraw as straw
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='python diffdomain-py3/subdivide_strength_change.py <-f> <-h1> <-h2> <-t1> <-t2> [options]')
parser.add_argument('-f','--file', type=str, default = None,help="Reorganized TADs identified by classification.py")
parser.add_argument('-h1','--hic1',type=str, default=None,help='Condition 1 hic file')
parser.add_argument('-h2','--hic2',type=str,default=None,help='Condition 2 hic file')
parser.add_argument('-t1','--tad1',type=str, default=None,help='condition1 TADs List')
parser.add_argument('-t2','--tad2',type=str, default=None,help='condition2 TADs List')
parser.add_argument('-o','--out',type=str,default=None,help='the output path')
parser.add_argument('-r','--reso',type=int,default=100000,help='resolution')
parser.add_argument('--sep',type=str,default='\t',help="The separator for file,tad1, and tad2")
args = parser.parse_args()

# load the files
data_types = pd.read_table(args.file,sep=args.sep)
cols = data_types.columns
data_types.rename(columns={cols[0]:'chr',cols[1]:'start',cols[2]:'end',cols[3]:'range',cols[4]:'type',cols[5]:'origin',cols[6]:'subtype',cols[7]:'significant'},inplace=True,errors='raise')
t1 = pd.read_table(args.tad1,sep=args.sep)
cols = t1.columns
t1.rename(columns={cols[0]:'chr',cols[1]:'start',cols[2]:'end'},inplace=True,errors='raise')
t2 = pd.read_table(args.tad2,sep=args.sep)
cols = t2.columns
t2.rename(columns={cols[0]:'chr',cols[1]:'start',cols[2]:'end'},inplace=True,errors='raise')

if t1['chr'].str.contains('chr').all() and t2['chr'].str.contains('chr').all() and not data_types['chr'].str.contains('chr').any():
    data_types['chr'] = 'chr' + data_types['chr'].astype(str)

hicpath1 = args.hic1
hicpath2 = args.hic2

def transCool(h,reso = args.reso):
    
    if h.endswith('cool') : # .cool or .mcool
        import cooler
        import subprocess as sp

        # Normalization
        hicnorm = 'KR'
        suffix = f'_{int(reso/1000)}k_{hicnorm}.cool'
        if h.endswith('.cool'):
            hicfile = h
            hic_norm = h.replace('.cool',suffix)
        elif h.endswith('.mcool'):
            hicfile = f'{h}::resolutions/{reso}'
            hic_norm = h.replace('.mcool',suffix)
            
        try:
            c = cooler.Cooler(f'{hic_norm}')

        except FileNotFoundError:
            # normalization
            # KR or ICE
            cmd = f'hicConvertFormat -m {hicfile} --inputFormat cool \
                --outputFormat cool --correction_name KR \
                --reso {reso} -o {hic_norm}'
            cp1 = sp.run(cmd,shell=True,capture_output=True,encoding='utf-8')

            if cp1.returncode == 0: # True
                print('Successful KR normalized.cool file:',hic_norm)
                return hic_norm
            else: # False
                raise(Exception(cp1.stderr))
        else:
            return hic_norm
        
    elif h.endswith('.hic'):
        return h
    else:
        print(r'Trying to sparse the files as three columns by "\t" ')
        print('Expect three columns :start[int] end[int] counts[float]')
        return h

def getSumTADs(h1, h2, t1, t2, reso=args.reso):

    def get_sum(row, h):

        region = f"{row['chr']}:{row['start']}:{row['end']}"
        region2 = "{0}:{1}-{2}".format(row['chr'], row['start'], row['end'])
        
        if h.endswith('.hic'):
            try:
                values = straw.straw('observed', 'KR', h, region, region, 'BP', reso)
                # print(region,np.sum([i.counts for i in values]))
                return np.sum([i.counts for i in values])
            except Exception as e:
                print(f'Error processing {region}: {e}')
                return 0  
            
        elif h.endswith('cool') : # .cool or .mcool

            import cooler

            try:
                c = cooler.Cooler(f'{h}')
                mat = c.matrix(balance=False).fetch(region2)
                return np.sum(mat)
            except Exception as e:
                print('Error:', e)
    
        # load a sparse matrix with three columns
        else:
            df = pd.read_table(h,sep=args.sep)
            start =int(row['start'])
            end =int(row['end'])
            td = df.loc[(df[df.columns[0]] >= start) & (df[df.columns[1]] <= end), df.columns[2]]

            return np.sum(td)

    s1 = t1.apply(lambda row: get_sum(row, h1), axis=1).sum()
    print('The sum of the KR-normalized Hi-C contact frequencies for all condition 1 TADs has been calculated.')

    s2 = t2.apply(lambda row: get_sum(row, h2), axis=1).sum()
    print('The sum of the KR-normalized Hi-C contact frequencies for all condition 2 TADs has been calculated.')

    return s1 / s2 if s2 != 0 else 0  

def subdivide_types(data_types, h1, h2, scaling_factor, reso=args.reso):

    types = data_types.copy()
    types.insert(types.shape[1]-1, 'subdivide_strength_change', np.nan)
    
    def getc(h,loc,region,reso=args.reso):

        if h.endswith('.hic'):

            c = straw.straw('observed','KR', h, loc, loc, 'BP', reso)
            con_hic_value = np.nanmedian([i.counts for i in c])
            return con_hic_value
        
        elif h.endswith('cool'):

            import cooler
                
            try:
                c = cooler.Cooler(f'{h}')
                mat = c.matrix(balance=False).fetch(region)
                return np.nanmedian(mat)
            except Exception as e:
                print('Error:', e)
        
        else:
            df = pd.read_table(h,sep=args.sep)
            index =  loc.split(':', 2)
            start =int(index[1])
            end =int(index[2])
            td = df.loc[(df[df.columns[0]] >= start) & (df[df.columns[1]] <= end), df.columns[2]]
            return np.nanmedian(td)
        
    for i in range(types.shape[0]):

        if types['subtype'][i] == 'strength' and types['significant'][i] == 1:

            loc1 = f"{types['chr'][i]}:{types['start'][i]}:{types['end'][i]}"
            loc2 = f"{types['chr'][i+1]}:{types['start'][i+1]}:{types['end'][i+1]}"

            region1 = "{0}:{1}-{2}".format(types['chr'][i], types['start'][i], types['end'][i])
            region2 = "{0}:{1}-{2}".format(types['chr'][i+1], types['start'][i+1], types['end'][i+1])
            
            # Compute Hi-C values using straw
            con1_hic_value = getc(h1,loc1,region1,reso)
            con2_hic_value = getc(h2,loc2,region2,reso)

            # Determine if strength up or down
            if con1_hic_value < con2_hic_value * scaling_factor:
                types['subdivide_strength_change'][i] = 'strength-change up'
            elif con1_hic_value > con2_hic_value * scaling_factor:
                types['subdivide_strength_change'][i] = 'strength-change down'

    types.fillna('nan', inplace=True)
    print("The median was successfully calculated")

    return types


try:
    hicpath1 = transCool(hicpath1,reso = args.reso)
    hicpath2 = transCool(hicpath2,reso = args.reso)
    scaling_factor = getSumTADs(hicpath1, hicpath2, t1, t2, reso=args.reso)
    subdivide_result = subdivide_types(data_types, hicpath1, hicpath2, scaling_factor=scaling_factor, reso=args.reso)
    if args.out != None: 

        if os.path.isdir(args.out): # directory path
            os.chdir(args.out) 
            subdivide_result.to_csv(f'subdivide_strength_change.txt', sep='\t', index=False)
            print(f'Results saved to: {os.path.join(os.getcwd(), "subdivide_strength_change.txt")}')
        else: # file path
            subdivide_result.to_csv(f'{args.out}', sep='\t', index=False)
            print(f'Results saved to: {args.out}')

    else: # default
        directory = os.path.dirname(args.file) + '/'
        subdivide_result.to_csv(f'{directory}subdivide_strength_change.txt', sep='\t', index=False)
        print(f'Results saved to: {directory}subdivide_strength_change.txt')

    print("Completed processing")
except Exception as e:
    print('Error:', e)