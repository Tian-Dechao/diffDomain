# all domains
# merge/split common boundary x = max(3bin,0.1 TAD Length)
# region < agrs.remote
# less complex
# zoom
# to filter the strength first
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import os
import warnings
warnings.filterwarnings('ignore')

# the arguments from command line
parser = argparse.ArgumentParser(description='python scriptname <-d> <-t> [options]')
parser.add_argument('-d','--diff', type=str, default = None,help="path/ the text of diffdoamin's outcome")
parser.add_argument('-t','--tad',type=str, default=None,help='path/ the other tadlist')
# parser.add_argument('-m','--method',type=str,default='BH',help='to choose the method that judge differential tads in "--diff"')
parser.add_argument('-o','--out',type=str,default=None,help='the output path')
parser.add_argument('-l','--limit',type=int,default=30000,help='the range(length of bases) to judge the common boundary')
parser.add_argument('-k','--kpercent',type=int,default=10,help='the common boundareis are within max(l*bin,k% TAD length)')
parser.add_argument('-r','--remote',type=int,default=1000000,help='the limitation of the biggest region')
parser.add_argument('-s1','--skip1',type=int,default=0,help='to skip the sth rows in "--diff" file (like 1,2,3')
parser.add_argument('-s2','--skip2',type=int,default=0,help='to skip the sth rows in the other tadlist file')
parser.add_argument('--sep1',type=str,default='\t',help="the seperator of the diffdoamin's outcome (like ',')")
parser.add_argument('--sep2',type=str,default='\t',help="the seperator of the other tadlist")
args = parser.parse_args()

# load the files

data = pd.read_table(args.diff,skiprows=range(args.skip1),sep=args.sep1)
tad = pd.read_table(args.tad,skiprows=range(args.skip2),sep=args.sep2)
# print(data)
#preprocessing
cols = data.columns
data.rename(columns={cols[0]:'chr',cols[1]:'start',cols[2]:'end'},inplace=True,errors='raise')
tad = tad.iloc[:,0:3]
tad.columns = ['chr','start','end']

data['chr'] = [str(i).replace('chr','').strip() for i in data['chr']]
tad['chr'] = [str(i).replace('chr','').strip() for i in tad['chr']]

data_diff = data.loc[data['adj_pvalue']<0.05,['chr','start','end']]
data_diff['significant'] = 1 
data_diff.reset_index(inplace=True,drop=True)

tad.sort_values(by=['chr','start','end'],inplace=True)
tad.reset_index(inplace=True,drop = True)
tad['range'] = [*map(lambda a,b:(a,b) , tad.start,tad.end)]

# preparation
chrs = list(map(str,[*range(1,24)]))+['X']
colnames = ['chr','start','end','range','type','origin','subtype','significant']
tad_ = data_main = loss = single = merge = split = multi = pd.DataFrame(columns=colnames)

tad_ = pd.concat([tad_,tad],axis=0)
tad = tad_
data_main = pd.concat([data_main,data.iloc[:,0:3]],axis=0)
data_main['significant'] = 0
data_main = pd.concat([data_main,data_diff],axis=0)
data_main.drop_duplicates(subset=['chr','start','end'],keep='last',inplace=True)
data_main['range'] = [*map(lambda a,b:(a,b) , data_main.start,data_main.end)]
data_main['origin'] = 'condition1'
data_main.sort_values(by=['chr','start','end'],inplace=True)
data_main.reset_index(inplace=True,drop=True)


def identical(boundary1,boundary2):
    # to judge the "common boundary"
    if int(boundary1) <= int(boundary2)+limit and int(boundary1) >= int(boundary2)-limit:
        return True 
    else: 
        return False

def cross(main,vise): 
    # main is the protagnist tad
    # to find the tads related to main in vise
    note=pd.DataFrame(columns=colnames)
    for i in range(vise.shape[0]):
        if (int(main[2])-limit > int(vise.iloc[i,1]) and int(main[1])+limit < int(vise.iloc[i,2]) ):
            note=pd.concat([note,pd.DataFrame(vise.iloc[i,:].values.reshape(1,-1),columns=colnames)],axis=0)
    return note

def n_of_region(outcome):
    # to count the number of regions in the dataframe 
    n_region = 0
    if len(outcome) != 0 :
        n_region = 1
        for i in range(2,len(outcome)):
            if outcome['origin'].values[i]=='conition1' and outcome['origin'].values[i-1]=='condition2':
                n_region = n_region+1
    return n_region

def n_diffdomain(outcome):
    n_diff = outcome.loc[outcome['origin']=='conition1',:].shape[0]
    return n_diff


# the 4th virsion+ bin
try:
    for c in tqdm(data_main.chr.unique()):
        temp = data_main.loc[data_main['chr']==c,:].copy()
        # print(c,temp.shape[0])
        tadlist = tad.loc[tad['chr']==c,:].copy()
        tadlist['origin'] = 'condition2'
        temp.reset_index(inplace=True,drop=True)
        tadlist.reset_index(inplace=True,drop=True)

        # filter the strength-change diffdomains and other non-significantly differentail tads with common boudaries in vise tadlist
        tad_index = []
        cross_index = []
        for i in range(temp.shape[0]): 
            # the i th TADs in the result of DiffDomain
            # to filter the TADs with common boundaries in different conditions

            # initialize the variables 
            note_tad = note_cross = pd.DataFrame(columns=colnames)
            # set the "limit" for judging the common boundaries
            limit = max(args.limit,args.kpercent*0.01*(temp['end'][i]-temp['start'][i]))
            note_tad = pd.concat([note_tad,pd.DataFrame(temp.loc[i,:].values.reshape(1,-1),columns=colnames)],axis=0)

            for k in range(tadlist.shape[0]):
                if (identical(temp.loc[i,'start'],tadlist.loc[k,'start'])) and (identical(temp.loc[i,'end'],tadlist.loc[k,'end'])) :
                    note_cross = pd.concat([note_cross,pd.DataFrame(tadlist.loc[k,:].values.reshape(1,-1),columns=colnames)],
                                           axis=0,ignore_index = True)
                    cross_index.append(k)
                    tad_index.append(i)

            n_cross = note_cross.shape[0]

            if n_cross !=0 : 
                # in case that there are TADs in temp having common boundaries but not in tadlist 

                for j in range(i+1,temp.shape[0]):
                    # to find the TADs (in the result of DiffDomain) located on the same boundaries with the i th TADs
                    if (identical(temp.loc[i,'start'],temp.loc[j,'start'])) and (identical(temp.loc[i,'end'],temp.loc[j,'end'])):
                        note_tad = pd.concat([note_tad,pd.DataFrame(temp.loc[j,:].values.reshape(1,-1),columns=colnames)],
                                             axis=0,ignore_index = True)
                        tad_index.append(i)
                        tad_index.append(j)

            note_tad.drop_duplicates(subset=['chr','start','end'],inplace=True)
            note_cross.drop_duplicates(subset=['chr','range','end'],inplace=True)
            n_tad = note_tad.shape[0]


            if n_tad ==1 and n_cross ==1 :
                note_tad['type'] = 'single'
                note_tad['subtype'] = 'strength'
                single = pd.concat([single,note_tad,note_cross],axis=0,
                                  ignore_index = True)

            elif n_tad == 1 and n_cross >=2 :
                note_tad['type'] = 'split'
                split = pd.concat([split,note_tad,note_cross],axis=0,
                                 ignore_index = True)

            elif n_tad >= 2 and n_cross ==1 :
                note_tad['type'] = 'merge'
                merge = pd.concat([merge,note_tad,note_cross],axis=0,
                                 ignore_index = True)

            elif n_tad >= 2 and n_cross >= 2 :
                if n_tad == n_cross :
                    note_tad['type'] = 'single'
                    note_tad['subtype'] = 'strength'
                    single = pd.concat([single,note_tad,note_cross],axis=0,
                                      ignore_index = True)
                else:
                    note_tad['type'] = 'complex'
                    multi = pd.concat([multi,note_tad,note_cross],axis=0,
                                     ignore_index = True)

        temp.drop(tad_index,inplace=True)
        temp.reset_index(drop=True,inplace=True)
        tadlist.drop(cross_index,inplace = True)
        tadlist.reset_index(drop=True,inplace=True)
    #         temp_sig = temp.loc[temp['significant']==1,:].copy()
    #         temp_sig.reset_index(drop = True,inplace=True)

        for i in range(temp.shape[0]):
            # to adjust the longest distance between "common boundaries"
            # to find the related TADs without common boundaries in different conditions

            limit = max(args.limit,(temp['end'][i]-temp['start'][i])*args.kpercent*0.01)
            note_cross = note_tad = pd.DataFrame(columns=colnames)

            # to find the tads in tadlist related to the significantly differential tad
            note_cross = pd.concat([note_cross,cross(temp.iloc[i,:].values,tadlist)],axis=0,
                                  ignore_index = True)

            note_tad = pd.concat([note_tad,pd.DataFrame(temp.iloc[i,:].values.reshape(1,-1),columns=colnames)],
                                 axis=0,ignore_index = True)

            n_cross = note_cross.shape[0]

            if n_cross == 0:
            # the significantly differential tad grew out of nothing
                note_tad['type'] = 'loss'
                loss = pd.concat([loss,note_tad],axis=0
                                 ,ignore_index = True)

            elif n_cross >=1:
                flag = 1
                while flag == 1 and note_tad['end'].max()-note_tad['start'].min() <= args.remote:
                    for p in range(note_cross.shape[0]):
                        # to find TADs related to the tads found in vise tadlist
                        note_tad = pd.concat([note_tad,cross(note_cross.iloc[p,:].values,temp)],axis=0,ignore_index = True)
                    for q in range(note_tad.shape[0]):
                        # to find TADs in the tadlist related to the TADs in the result of DiffDomain
                        note_cross = pd.concat([note_cross,cross(note_tad.iloc[q,:].values,tadlist)],axis=0,ignore_index = True)

                    first_tad =  note_tad.loc[note_tad.start == min(note_tad.start),:]
                    last_tad = note_tad.loc[note_tad.end == max(note_tad.end),:]
                    first_cross = note_cross.loc[note_cross.start == min(note_cross.start),:]
                    last_cross = note_cross.loc[note_cross.end == max(note_cross.end),:] 

                    thres1 = pd.concat([cross(first_tad.values[0],tadlist),cross(last_tad.values[0],tadlist)],axis=0)
                    thres2 = pd.concat([cross(first_cross.values[0],temp),cross(last_cross.values[0],temp)],axis=0)

                    if (thres1['range'].isin(note_cross['range'])).all() and thres2['range'].isin(note_tad['range']).all():
                        flag = 2


                note_tad.drop_duplicates(subset=['chr','start','end'],inplace=True)
                note_cross.drop_duplicates(subset=['chr','start','end'],inplace=True)
                note_tad.reset_index(inplace=True,drop=True)
                note_cross.reset_index(inplace=True,drop=True)
                n_tad = note_tad.shape[0]
                n_cross  = note_cross.shape[0]

                if n_tad == 1 and n_cross == 1: 
                    note_tad['type'] = 'single'
                    note_tad['subtype'] = 'zoom'
                    single = pd.concat([single,note_tad,note_cross],axis=0,ignore_index = True)

                elif n_tad == 1 and n_cross >= 2:
                    note_tad['type'] = 'split'
                    split = pd.concat([split,note_tad,note_cross],axis=0,ignore_index = True)

                elif n_tad >= 2 and n_cross ==1:
                    note_tad['type'] = 'merge'
                    merge = pd.concat([merge,note_tad,note_cross],axis=0,ignore_index = True)

                elif n_tad >=2 and n_cross >=2:
                    note_tad['type'] = 'complex'
                    multi = pd.concat([multi,note_tad,note_cross],axis=0,ignore_index = True)


except Exception as e:
    print(e)
    print('Interrupted!')
else:
    result = pd.DataFrame(columns=colnames)
    # there will be duplicates because of travelsal
    result = pd.concat([result,loss,single,merge,split,multi],axis=0,ignore_index=True)
    result.drop_duplicates(subset=['chr','start','end','origin','significant'],
                           inplace=True)

    if args.out != None:
        os.chdir(args.out)    
    filename = os.path.split(args.diff)[1]
    result.fillna('nan',inplace=True)
    result.to_csv('{}_types.txt'.format(filename),sep='\t',index=False)
    print(filename)

    # whether all of the significantly defferential tads have been selected
    len_sig = result.loc[result['significant']==1,:].shape[0]
    if len_sig == data_diff.shape[0]:
        print(f'Complete.There are {data_diff.shape[0]} reorganized TADs.')
    else:
        print('Error!')
        print(f'there are {data_diff.shape[0]} reorganized TADs,but the length of result is {len_sig}.')
        
        


    
