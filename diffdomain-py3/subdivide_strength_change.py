import numpy as np
import pandas as pd
import os
import hicstraw as straw
import argparse
import warnings
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='python diffdomain-py3/subdivide_strength_change.py <-f> <-h1> <-h2> <-t1> <-t2> [options]')
parser.add_argument('-f','--file', type=str, default = None,help="Reorganized TADs identified by classification.py")
parser.add_argument('-h1','--hic1',type=str, default=None,help='Condition 1.hic file')
parser.add_argument('-h2','--hic2',type=str,default=None,help='Condition 2.hic file')
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
hicpath1 = args.hic1
hicpath2 = args.hic2

def getSumTADs(h1, h2, t1, t2, reso=args.reso):

    def get_sum(row, h):
        region = f"{row['chr']}:{row['start']}:{row['end']}"
        try:
            values = straw.straw('observed', 'KR', h, region, region, 'BP', reso)
            return np.sum([i.counts for i in values])
        except Exception as e:
            print(f'Error processing {region}: {e}')
            return 0  

    s1 = t1.apply(lambda row: get_sum(row, h1), axis=1).sum()
    s2 = t2.apply(lambda row: get_sum(row, h2), axis=1).sum()

    print('the sum of the KR-normalized Hi-C contact frequencies across all condition 1 TADs:',s1)
    print('the sum of the KR-normalized Hi-C contact frequencies across all condition 2 TADs:',s2)
    return s1 / s2 if s2 != 0 else 0  

def subdivide_types(data_types, h1, h2, scaling_factor, reso=args.reso):

    types = data_types.copy()
    types.insert(types.shape[1]-1, 'subdivide_strength_change', np.nan)
    
    for i in range(types.shape[0]):

        if types['subtype'][i] == 'strength' and types['significant'][i] == 1:

            loc1 = f"{types['chr'][i]}:{types['start'][i]}:{types['end'][i]}"
            loc2 = f"{types['chr'][i+1]}:{types['start'][i+1]}:{types['end'][i+1]}"

            # Compute Hi-C values using straw
            c1 = straw.straw('observed','KR', h1, loc1, loc1, 'BP', reso)
            c2 = straw.straw('observed','KR', h2, loc2, loc2, 'BP', reso)
            con1_hic_value = np.nanmedian([i.counts for i in c1])
            con2_hic_value = np.nanmedian([i.counts for i in c2])

            # Determine if strength up or down
            if con1_hic_value < con2_hic_value * scaling_factor:
                types['subdivide_strength_change'][i] = 'strength-change up'
            elif con1_hic_value > con2_hic_value * scaling_factor:
                types['subdivide_strength_change'][i] = 'strength-change down'
    types.fillna('nan', inplace=True)
    print("The median was successfully calculated")

    return types


try:
    scaling_factor = getSumTADs(hicpath1, hicpath2, t1, t2, reso=args.reso)
    subdivide_result = subdivide_types(data_types, hicpath1, hicpath2, scaling_factor=scaling_factor, reso=args.reso)
    if args.out != None: 
        subdivide_result.to_csv(f'{args.out}_subdivide_strength_change.bed', sep='\t', index=False)
        print(f'Results saved to: {args.out}_subdivide_strength_change.bed')
    else:
        directory = os.path.dirname(args.file) + '/'
        subdivide_result.to_csv(f'{directory}subdivide_strength_change.bed', sep='\t', index=False)
        print(f'Results saved to: {directory}subdivide_strength_change.bed')
    print("Completed processing")
except Exception as e:
    print('Error:', e)