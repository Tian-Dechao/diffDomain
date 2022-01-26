
def diffdomain_one(chr,start,end,fhic0,fhic1,
                   reso=10000,hicnorm='KR',min_nbin=10,f=0.5):
    """
    Test if one domain is significantly rewried in another biological condition
    Usage:
    diffdomain_one(<chr>,<start>,<end>,<hic0>,<hic1>,[options])
    
    --chr chromosome number
    --start the start position of the domain tested by this func
    --end the end position of the domain
    --fhic0 the path of the first hic file
    --fhic1 the path of the second hic file
    
    Opts:
    --reso resolution for hicfile  [default: 100000]
    --hicnorm hic matrix normalization method [default: 'KR']
    --min_nbin effective number of bin  [default: 10]
    --f parameters for filtering the null values of the matrix[0~1)  [default: 0.5]

"""

    result = comp2domins_by_twtest(chrn=chr, start=int(start),
                                   end=int(end), reso=int(reso),
                                   hicnorm=hicnorm,
                                   fhic0=fhic0,
                                   fhic1=fhic1,
                                   min_nbin=min_nbin,f=f)
    result = pd.DataFrame(result)
    return result


def diffdomain_multiple(fhic0,fhic1,fbed,sep='\t',hicnorrm = 'KR',
                       chrn='ALL',reso = 100000,ncore=10,
                       min_nbin = 10, f = 0.5):
    """
    Test if one domain is significantly rewried in another biological condition
    Usage:
    diffdomain_multiple(<fhic0>,<fhic1>,<fbed>,[options])

    --fhic0 the filepath of the first hic file
    --fhic1 the filepath of the second hic file
    --fbed the filepath of TADs' list that you want to test,it usually is the tadlist of hic0(the first hic file)

    Options:
    --sep deliminator for hicfile  [default: '\t']
    --hicnorm hic matrix normalization method  [default: 'KR']
    --chrn chromosome number  [default: 'ALL']
    --reso resolution for hicfile  [default: 100000]
    --ncore number of parallel process  [default: 10]
    --min_nbin effective number of bin  [default: 10]
    --f parameters for filtering the null values of the matrix[0~1)  [default: 0.5]

    """
    # Parallel implementation
    def comp2domins_by_twtest_parallel(i):
        reso = int(reso)
        tadb = loadtads(fbed, sep=sep, chrnum=chrn,
                        min_nbin=int(min_nbin), reso=int(reso))
        tmp_res = comp2domins_by_twtest(chrn=tadb.iloc[i, 0],
                                        start=tadb.iloc[i, 1],
                                        end=tadb.iloc[i, 2],
                                        reso=int(reso),
                                        hicnorm=hicnorm,
                                        fhic0=fhic0,
                                        fhic1=fhic1,
                                        min_nbin=int(min_nbin),
                                        f=f)
        return tmp_res

    P = Pool(int(ncore))
    result = P.map(comp2domins_by_twtest_parallel,
                   range(tadb.shape[0]))
    P.terminate()
    result = pd.DataFrame(result)

    return result


def adjustment(inputdf,alpha=0.05,Filter=False,method='fdr_bh',
               sort=False):
    """
    This func is to adjust the result of diffdomain_multiple
    Usage : 
    adjustment(<inputdf>,<Filter>,[options])
    
    --inputdf the result of diffdomain_multiple (pd.DataFrame)
    
    Options:
    --alpha the threshold of adjusted pvalue [default: 0.05]
    --Filter As long as the pvalue of TADs is less than alpha after adjustment if argument is true [True/False, default: False]
    --method adjustment method you want to use [default: 'fdr_bh']
    --sort wheter to sort the result [default: False]
    """
    
    import statsmodels.stats.multitest as smm
    def read_compared(inputfile, method, alpha):
        df = inputfile.copy()
        df.columns=['chr','start','end','region','stat','pvalue','bins']
        pvalue = df['pvalue']
        pvalue = df['pvalue'].fillna(pvalue.median())
        rej, pval_corr = smm.multipletests(pvalue, is_sorted=sort, alpha=alpha, method=method)[:2]
        df['adj_pvalue'] = pval_corr
        return df

    res = read_compared(inputfile=inputdf, method=method,alpha=alpha)
    if Filter == True:
        res = res[res['adj_pvalue']<=alpha]
        return res

    elif Filter == False:
        return res

    else:
        print("Sorry, this is an invalid parameter ")

        
def visualization(chr,start,end,fhic0,fhic1,ofile,reso=100000,hicnorm='KR',):
    """
    The "visualizing" func is to visualize the two hic files.
    The lower triangle of output is the first input hic0.
    The upper triangle of output is the second input hic1.
    
    Usage:
    visualizing(<chr>,<start>,<end>,<fhic0>,<fhic1>,[options])
    
    --chr chromosome number
    --start the start position of the domain visulized by the func
    --end the end position of the domain visulized by the func
    --fhic0 the filepath of the first hic file
    --fhic1 the filepath of the second hic file
    --ofile filepath for output file
   
    Options:
    --reso resolution for hicfile  [default: 100000]
    --hicnorm hic matrix normalization method [default: 'KR']
    
    """
    import matplotlib
    matplotlib.use('Agg')
    from comparedomains import visualization
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    
    def save_fig(fig_id, tight_layout=True, fig_extension="pdf", resolution=600):
        path = os.path.join(fig_id + "." + fig_extension)
        print("Saving figure", fig_id)
        if tight_layout:
            plt.tight_layout()
            plt.savefig(path, format=fig_extension, dpi=resolution)

    def plot(mat,outputfile):
        fig,ax = plt.subplots(figsize=(6,6))
        sns.heatmap(data=mat,
                vmax=int(round(np.nanmax(mat)))/5,
                cmap=cdict,cbar=False,mask=mat<0.1,
                square=True)
        plt.xticks([])
        plt.yticks([])
        ax.tick_params(bottom=False, top=False,
                       left=False,right=False)
        save_fig(outputfile)

    cdict=LinearSegmentedColormap.from_list('mycmap',['w','#FF0000'])
    result = visualization(chrn=chr,
                           start=int(start),
                           end=int(end),
                           reso=int(reso),
                           hicnorm=hicnorm,
                           fhic0=fhic0,
                           fhic1=fhic1)
    
    
    matrix1 = result[0]
    matrix2 = result[1]
    matrix1[np.isnan(matrix1)] = 1
    matrix2[np.isnan(matrix2)] = 1
    tril_mat = (np.tril(matrix1,-1))
    triu_mat = (np.triu(matrix2,1))
    mat = tril_mat + triu_mat
    plot(mat,ofile)


    
def classification(result_adj_df,tadlist_df,alpha=0.05,limit=40000,kpercent=10,remote=1000000,ofile=None):
    """
    This function is to classify the TADs in adjusted result of diffdomain_multiple into six subtypes to aid bilological analysis and interpretations.
    
    Usage:
    changings(<result_adj_df>,<tadlist_df>,[options])
    --result_adj_df the dateframe of adjusted outcome of diffdomain_multiple
    --tadlist_df the tadlist(dataframe) of the second hic file
    
    Options:
    --alpha the threshold of adjusted p-value
    --limit length of bases, within which the boundaries will be judged as common boundaries [default: 40000]
    --kpercent the common boundareis are within max(l*bin,k% TAD's length) [default: 10]
    --remote the limitation of the biggeset region
    --ofile the filepath of output file
    """

    # load the files
    
    data = result_adj_df.copy()
    tad = tadlist_df.copy()
    #preprocessing
    cols = data.columns
    data.rename(columns={cols['V1']:'chr',cols['V2']:'start',cols['V3']:'end'},inplace=True)
    data_diff = data.loc[data['adj_pvalue']<alpha,['chr','start','end']]
    data_diff['significant'] = 1 
    data_diff.reset_index(inplace=True,drop=True)
    tad = tad.iloc[:,0:3]
    tad.columns = ['chr','start','end']
    tad.sort_values(by=['chr','start','end'],inplace=True)
    tad.reset_index(inplace=True,drop = True)
    tad['range'] = list(map(lambda a,b:(a,b) , tad.start,tad.end))

    # preparation
    chrs = list(map(str,list(range(1,23))))+['X']
    colnames = ['chr','start','end','range','type','origin','subtype','significant']
    tad_ = data_main = loss = single = merge = split = multi = pd.DataFrame(columns=colnames)

    tad_ = pd.concat([tad_,tad],axis=0)
    tad = tad_
    data_main = pd.concat([data_main,data.iloc[:,0:3]],axis=0)
    data_main['significant'] = 0
    data_main = pd.concat([data_main,data_diff],axis=0)
    data_main.drop_duplicates(subset=['chr','start','end'],keep='last',inplace=True)
    data_main['range'] = list(map(lambda a,b:(a,b) , data_main.start,data_main.end))
    data_main['origin'] = 'diffdomain'
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
            if (int(main['end'])-limit > int(vise.loc[i,'start']) and int(main['start'])+limit < int(vise.loc[i,'end']) ):
                note=pd.concat([note,pd.DataFrame(vise.loc[i,:].values.reshape(1,-1),columns=colnames)],axis=0)
        return note

    def n_of_region(outcome):
        # to count the number of regions in the dataframe 
        n_region = 0
        if len(outcome) != 0 :
            n_region = 1
            for i in range(2,len(outcome)):
                if outcome['origin'].values[i]=='diffdomain' and outcome['origin'].values[i-1]=='the other tadlist':
                    n_region = n_region+1
        return n_region

    def n_diffdomain(outcome):
        n_diff = outcome.loc[outcome['origin']=='diffdomain',:].shape[0]
        return n_diff


    # the 4th virsion+ bin
    try:
        for c in chrs:
            temp = data_main.loc[data_main['chr']==c,:].copy()
            tadlist = tad.loc[tad['chr']==c,:].copy()
            tadlist['origin'] = 'the other tadlist'
            temp.reset_index(inplace=True,drop=True)
            tadlist.reset_index(inplace=True,drop=True)
            temp = temp[colnames]
            tadlist = tadlist[colnames]
            temp['start'] = temp['start'].astype(int)
            temp['end'] = temp['end'].astype(int)
            tadlist['start'] = tadlist['start'].astype(int)
            tadlist['end'] = tadlist['end'].astype(int)
            # filter the strength-change diffdomains and other non-significantly differentail tads with common boudaries in vise tadlist
            tad_index = []
            cross_index = []
            for i in range(temp.shape[0]): 
                # the i th TADs in the result of DiffDomain
                # to filter the TADs with common boundaries in different conditions

                # initialize the variables 
                note_tad = note_cross = pd.DataFrame(columns=colnames)
                # set the "limit" for judging the common boundaries
                limit = max(limit,kpercent*0.01*(temp['end'][i]-temp['start'][i]))
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
                note_cross.drop_duplicates(subset=['chr','start','end'],inplace=True)
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

                limit = max(limit,(temp['end'][i]-temp['start'][i])*kpercent*0.01)
                note_cross = pd.DataFrame(columns=colnames)
                note_tad = pd.DataFrame(columns=colnames)

                # to find the tads in tadlist related to the significantly differential tad
                note_cross = pd.concat([note_cross,cross(temp.iloc[i,:],tadlist)],axis=0,
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

                    note_tad['start'] = note_tad['start'].astype(int)
                    note_tad['end'] = note_tad['end'].astype(int)
                    note_cross['start'] = note_cross['start'].astype(int)
                    note_cross['end'] = note_cross['end'].astype(int)
                    while (flag == 1) and (max(note_tad['end'])-min(note_tad['start']) <= int(remote)):
                        for p in range(note_cross.shape[0]):
                            # to find TADs related to the tads found in vise tadlist
                            note_tad = pd.concat([note_tad,cross(note_cross.iloc[p,:],temp)],axis=0,ignore_index = True)
                        for q in range(note_tad.shape[0]):
                            # to find TADs in the tadlist related to the TADs in the result of DiffDomain
                            note_cross = pd.concat([note_cross,cross(note_tad.iloc[q,:],tadlist)],axis=0,ignore_index = True)

                        first_tad =  note_tad.loc[note_tad.start == min(note_tad.start),:]
                        last_tad = note_tad.loc[note_tad.end == max(note_tad.end),:]
                        first_cross = note_cross.loc[note_cross.start == min(note_cross.start),:]
                        last_cross = note_cross.loc[note_cross.end == max(note_cross.end),:] 

                        thres1 = pd.concat([cross(first_tad.iloc[0,:],tadlist),cross(last_tad.iloc[0,:],tadlist)],axis=0)
                        thres2 = pd.concat([cross(first_cross.iloc[0,:],temp),cross(last_cross.iloc[0,:],temp)],axis=0)

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

        if ofile != None:
            result.to_csv('{}_types.txt'.format(filename),sep='\t',index=False)
        print(filename)

        # whether all of the significantly defferential tads have been selected
        len_sig = result.loc[result['significant']==1,:].shape[0]
        if len_sig == data_diff.shape[0]:
            print('Complete.There are %d reorganized TADs.'%data_diff.shape[0])
        else:
            print('Error!')
            print('there are %d reorganized TADs,but the length of result is %d.'%(data_diff.shape[0],len_sig))
