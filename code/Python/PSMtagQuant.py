#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import gridspec
import os
import re
from itertools import combinations,combinations_with_replacement, product
import sys
sys.path.insert(0, "/Users/kevinmcdonnell/Programming/Python/PTI/jmod")
from src.utils.io.read_output import names, dtypes
import src.utils.misc_functions as mf
# from jmod.src.mass_tags import tag_library
# import config
from scipy.stats import gaussian_kde

# from upsetplot import generate_counts
# from upsetplot import plot
# from upsetplot import from_contents
# from upsetplot import UpSet
from matplotlib_venn import venn2, venn3
from sklearn.metrics import auc
from pathlib import Path

from functools import reduce

import seaborn as sb

import tqdm


from src.utils.io.load_fasta import all_fasta_seqs,  all_protein_names


logic_and = np.logical_and.reduce
def intersect_all(*arrs): return reduce(np.intersect1d,arrs)


labels = ["Human",
          "Ecoli",
          "Yeast"]

orgs =np.array(['Human', 'Ecoli', 'Yeast'])
    
org_colors = ["tab:blue","tab:orange","tab:green"]
org_colors = ["tab:blue","tab:orange","tab:red"]

colours = ["tab:blue","tab:orange","tab:green","tab:red",
'tab:purple',
'tab:brown',
'tab:pink',
'tab:gray',
'tab:olive',
'tab:cyan']

plt.rcParams['figure.dpi'] = 500
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['axes.titlesize'] = 15


labels = ["Human",
          "Ecoli",
          "Yeast"]

orgs =np.array(['Human', 'Ecoli', 'Yeast'])
    
org_colors = ["tab:blue","tab:orange","tab:green"]
org_colors = ["tab:blue","tab:orange","tab:red"]

def remove_end(x,end="_y"):
    return x[:-2] if x.endswith(end) else x

def load_jmod_files(fdx_files1,channels,channel_mix,get_dcs=False,file_name =""):
    
    fdxs = []
    dcs = []
    for idx,fdx_file in enumerate(fdx_files1):
        print(f"file {idx+1} of {len(fdx_files1)}")
        fdx_n = pd.read_parquet(fdx_file+file_name)
        if "channel" not in fdx_n.columns:
            fdx_n["channel"] = idx
        fdx_n["file_idx"] = idx
        
        if "MS1_Int_x" in fdx_n.columns:
            col_to_drop = [i for i in list(fdx_n.columns) if i.endswith("_x")]
            fdx_n = fdx_n.drop(col_to_drop,axis=1)
        fdx_n = fdx_n.rename(remove_end,axis="columns")  
        fdxs.append(fdx_n)
        
        if get_dcs:
            dc_path = os.path.dirname(fdx_file)+"/decoylibsearch_coeffs.csv"
            dc_n = pd.read_csv(dc_path,header=None,names=names,dtype=dtypes)
            if "channel" not in dc_n.columns:
                dc_n["channel"] = idx
            dc_n["file_idx"] = idx
            dcs.append(dc_n)
    
    
    # ### Assume they are all of the same type (Tag)
    # params = {}
    # with open(os.path.dirname(fdx_file)+"/params.txt","r") as read_file:
        
    #     line = read_file.readline()
    #     while line:
    #         if ": " in line and not re.match("default_dict",line):
    #             if re.match("tag: TagName",line):
    #                 line = re.sub("tag: ","",line)
    #             try:
    #                 i,j = re.split(": ",line)
    #             except:
    #                 print(line)
    #             params[i]=j.strip()
    #         line = read_file.readline()
    
    # if "TagName" in params:
    #     # if params["TagName"]=="mTRAQ_678":
    #     #     mass_tag1 = mTRAQ_678
    #     #     tag_name1 = "mTRAQ"
            
    #     # elif params["TagName"]=="mTRAQ_02468":
    #     #     mass_tag1 = mTRAQ_02468
    #     #     tag_name1 = "mTRAQ"
            
    #     # elif params["TagName"]=="mTRAQ":
    #     #     mass_tag1 = mTRAQ
    #     #     tag_name1 = "mTRAQ"
        
    #     # elif params["TagName"]=="diethyl_6plex":
    #     #     mass_tag1 = diethyl_6plex
    #     #     tag_name1 = "DE"
    #     mass_tag1 = None
    #     tag_name1 = "plexDIA"
    # else:
    #     # raise ValueError
    #     mass_tag1 = None
    #     tag_name1 = "LF"
    
    fdx1 = pd.concat(fdxs,ignore_index=True)
    dc1 = None
    if get_dcs:
        dc1 = pd.concat(dcs,ignore_index=True)
    # fdx1["tag"] = tag_name1
    
    fdx1["org"] = np.array([";".join(orgs[[i in all_fasta_seqs[j] for j in range(3)]]) for i in fdx1["stripped_seq"]])
    fdx1["Run"] = [Path(i).stem for i in fdx1.file_name]
    fdx1["BestProtein_Qvalue"] = fdx1.groupby(["protein"])["Protein_Qvalue"].transform("min")
    # fdx2["Run"] = [Path(i).stem for i in fdx2.file_name]
    fdx1["mix"] = [channel_mix[channels.index(i)] for i in fdx1.channel]
    # if "num_b" not in fdx1.columns:
    #     fdx1["num_b"] = fdx1.frag_names.str.count("b")
    
    return fdx1,dc1



def merge_and_norm(orig_quant_name,fdx,file_names,file_groups,plex_filter=-np.inf):
    merged_data = []
    for idx,(win,files) in enumerate(zip(file_names,file_groups)):
        # break
        quant_greater_zero = (fdx[orig_quant_name]>0).to_numpy() if type(orig_quant_name)==str else np.logical_and.reduce([(fdx[qn]>0).to_numpy() for qn in orig_quant_name])
        if plex_filter>0:
            subset = fdx[np.logical_and.reduce([fdx['Run'].isin(files).to_numpy(),(fdx["plexfitMS1_p"]>plex_filter).to_numpy(),quant_greater_zero])].copy()  # Filter by file group
        else:
            subset = fdx[np.logical_and.reduce([fdx['Run'].isin(files).to_numpy(),quant_greater_zero])].copy()  # Filter by file group
        subset['WindowNum'] = win  # Add group name
        subset["file_idx"] = idx
        if type(orig_quant_name)==str:
            quant_names = [orig_quant_name]
        else:
            quant_names = orig_quant_name
        for qn in quant_names:
            subset["channel_median_"+qn] = subset.groupby(["file_idx","channel"])[qn].transform(lambda x:  np.power(2,np.nanmedian(np.log2(x))))
            subset["Normalized_"+qn] = subset.groupby(["file_idx","channel"])[qn].transform(lambda x:np.log2(x) - np.nanmedian(np.log2(x)))
            subset["merged_"+qn] = subset.groupby(["untag_prec","channel"])["Normalized_"+qn].transform(lambda x: np.median(x))
        subset = subset.drop_duplicates(["untag_prec","channel"],keep="last").reset_index(drop=True)
        # median_values = subset.groupby(["untag_prec","channel","last_aa","org"])[orig_quant_name].median().reset_index()
        
        merged_data.append(subset)
        
    return pd.concat(merged_data)



def get_intersect(  fdxA,
                    fdxB,
                    quant_name,
                    conditions = [],
                    file_idx1 = None,
                    file_idx2 = None,
                    channels=None,
                    mix=None,
                    at_least_b = 3,
                    Q_value_limit = 0.05,
                    BestChannel_Qvalue_limit = 0.01,
                    last_aa = "R"):
    
    if channels is not None:
        channelA,channelB = channels
    elif mix is not None:
        mixA,mixB = mix
    else:
        raise ValueError("No channel specified")
    
    intersect1 = intersect_all(fdxA["untag_prec"][logic_and([fdxA["channel"]==channelA if channels is not None else np.ones_like(fdxA["channel"],dtype=bool),
                                                             fdxA["mix"]==mixA if mix is not None else np.ones_like(fdxA["channel"],dtype=bool),
                                                               fdxA["file_idx"]==file_idx1 if file_idx1 is not None else np.ones_like(fdxA["channel"],dtype=bool),
                                                                fdxA["plexfitMS1_p"]>.75 if "pearson_fit" in conditions and "plexfitMS1_p" in fdxA.columns else np.ones_like(fdxA["channel"],dtype=bool),
                                                                fdxA["last_aa"]==last_aa if "last_aa" in conditions else np.ones_like(fdxA["channel"],dtype=bool),
                                                                fdxA.Qvalue<Q_value_limit if "Qvalue" in conditions else np.ones_like(fdxA["channel"],dtype=bool),
                                                                fdxA.BestChannel_Qvalue<BestChannel_Qvalue_limit if "BestChannel" in conditions else np.ones_like(fdxA["channel"],dtype=bool),
                                                                fdxA.num_b>=at_least_b if "num_b" in conditions else np.ones_like(fdxA["channel"],dtype=bool),
                                                                # fdxA[quant_name]>0
                                                                
                                                              ])],
                                 fdxB["untag_prec"][logic_and([fdxB["channel"]==channelB if channels is not None else np.ones_like(fdxA["channel"],dtype=bool),
                                                               fdxB["mix"]==mixA if mix is not None else np.ones_like(fdxB["channel"],dtype=bool),
                                                               fdxB["file_idx"]==file_idx2 if file_idx2 is not None else np.ones_like(fdxB["channel"],dtype=bool),
                                                                fdxB["plexfitMS1_p"]>.75 if "pearson_fit" in conditions  and "plexfitMS1_p" in fdxB.columns else np.ones_like(fdxB["channel"],dtype=bool),
                                                                fdxB["last_aa"]==last_aa if "last_aa" in conditions else np.ones_like(fdxB["channel"],dtype=bool),
                                                                fdxB.Qvalue<Q_value_limit if "Qvalue" in conditions else np.ones_like(fdxB["channel"],dtype=bool),
                                                                fdxB.BestChannel_Qvalue<BestChannel_Qvalue_limit if "BestChannel" in conditions else np.ones_like(fdxB["channel"],dtype=bool),
                                                                fdxB.num_b>=at_least_b if "num_b" in conditions else np.ones_like(fdxB["channel"],dtype=bool),
                                                                # fdxB[quant_name]>0
                                                              ])])
    
    return intersect1



def compute_all_comparisons(fdx,
                            quant_name,
                            channels,
                            channel_mix,
                            theoretical_yeast_amounts,
                            plot_type = "all",
                            conditions=['Qvalue', 'BestChannel'],
                            at_least_b = 3,
                            Q_value_limit = 0.05,
                            BestChannel_Qvalue_limit = 0.01,
                            name="",
                            fdx2=None,
                            fdx3=None,
                            bin_num=60,
                            counter_size=1.1,
                            y_lim = 4,
                            scale=1.3,
                            only_A = False,
                            figsize=(7,1),
                            ax=None,
                            ycolor="tab:red",
                            return_vals = False,
                            return_prec = False
                            ):

    
    mix_values={i:[] for i in list(combinations_with_replacement(np.unique(channel_mix), 2))}
    mix_values_h={i:[] for i in list(combinations_with_replacement(np.unique(channel_mix), 2))}
    prec_values={i:[] for i in list(combinations_with_replacement(np.unique(channel_mix), 2))}
    prec_values_h={i:[] for i in list(combinations_with_replacement(np.unique(channel_mix), 2))}
    int_values={i:[] for i in list(combinations_with_replacement(np.unique(channel_mix), 2))}
    int_values_h={i:[] for i in list(combinations_with_replacement(np.unique(channel_mix), 2))}
    # mix_values={i:[] for i in list(product(np.unique(channel_mix), repeat=2))}
    
    # scale=1.3
    if plot_type=="all" and ax is None:
        fig,ax = plt.subplots(figsize=(7*scale,1*scale))
    hist_labels = []
    for idx,chan in enumerate(list(combinations(np.unique(channels), 2))):
        # break
        if fdx2 is None:
            intersect1 = get_intersect(fdx, fdx,
                                           conditions=conditions, 
                                           # file_idx1=file_idx1, file_idx2=file_idx2,
                                            channels=chan,
                                            quant_name=quant_name,
                                           # mix=["A","B"]
                                            at_least_b = at_least_b,
                                            Q_value_limit = Q_value_limit,
                                            BestChannel_Qvalue_limit = BestChannel_Qvalue_limit
                                           )
        elif fdx3 is None:
            intersect1 = intersect_all(get_intersect(fdx, fdx,
                                                       conditions=conditions, 
                                                       # file_idx1=file_idx1, file_idx2=file_idx2,
                                                        channels=chan,
                                                        quant_name=quant_name,
                                                       # mix=["A","B"]
                                                        at_least_b = at_least_b,
                                                        Q_value_limit = Q_value_limit,
                                                        BestChannel_Qvalue_limit = BestChannel_Qvalue_limit
                                                       ),
                                        get_intersect(fdx2, fdx2,
                                                    conditions=conditions, 
                                                    # file_idx1=file_idx1, file_idx2=file_idx2,
                                                     channels=chan,
                                                     quant_name=quant_name,
                                                    # mix=["A","B"]
                                                     at_least_b = at_least_b,
                                                     Q_value_limit = Q_value_limit,
                                                     BestChannel_Qvalue_limit = BestChannel_Qvalue_limit
                                                    )
                                                            )
        else:
            intersect1 = intersect_all(get_intersect(fdx, fdx,
                                                       conditions=conditions, 
                                                       # file_idx1=file_idx1, file_idx2=file_idx2,
                                                        channels=chan,
                                                        quant_name=quant_name,
                                                       # mix=["A","B"]
                                                        at_least_b = at_least_b,
                                                        Q_value_limit = Q_value_limit,
                                                        BestChannel_Qvalue_limit = BestChannel_Qvalue_limit
                                                       ),
                                        get_intersect(fdx2, fdx2,
                                                    conditions=conditions, 
                                                    # file_idx1=file_idx1, file_idx2=file_idx2,
                                                     channels=chan,
                                                     quant_name=quant_name,
                                                    # mix=["A","B"]
                                                     at_least_b = at_least_b,
                                                     Q_value_limit = Q_value_limit,
                                                     BestChannel_Qvalue_limit = BestChannel_Qvalue_limit
                                                    ),
                                        get_intersect(fdx3, fdx3,
                                                    conditions=conditions, 
                                                    # file_idx1=file_idx1, file_idx2=file_idx2,
                                                     channels=chan,
                                                     quant_name=quant_name,
                                                    # mix=["A","B"]
                                                     at_least_b = at_least_b,
                                                     Q_value_limit = Q_value_limit,
                                                     BestChannel_Qvalue_limit = BestChannel_Qvalue_limit
                                                    )
                                                            )
        shared_keys_set = set(intersect1)
        sk_bool = np.array([i in shared_keys_set for i in fdx["untag_prec"]]) ## faster
        
        org_offset = 0
        offset=0
        
        for o_i in [0,2]:
            # break
            
            
            i1 = fdx[logic_and([fdx["channel"]==chan[0],sk_bool,fdx["org"]==orgs[o_i]])]
            i2 = fdx[logic_and([fdx["channel"]==chan[1],sk_bool,fdx["org"]==orgs[o_i]])]
            
            i1 = i1.sort_values(by="untag_prec")
            i2 = i2.sort_values(by="untag_prec")
            # if np.min(i1[quant_name])<0:## log scale
            if "merged" in quant_name:
                # print("log")
                log_ratios = i1[quant_name].to_numpy()-i2[quant_name].to_numpy()
            else:
                # print("no log")
                log_ratios = np.log2(np.divide(i1[quant_name].to_numpy(),i2[quant_name].to_numpy()))
            # log_ratios = i1[quant_name].to_numpy()-i2[quant_name].to_numpy()
            i1["log_ratios"] = log_ratios
            i2["log_ratios"] = log_ratios
            
            if o_i==0:
                offset = np.nanmedian(log_ratios)
                sign = 1 if tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))==tuple(([i1.mix.iloc[0],i2.mix.iloc[0]])) else -1
                mix_values_h[tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))]+=list(sign*(log_ratios-offset))
                prec_values_h[tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))]+=list(i1.prec)
                int_values_h[tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))]+= list((i1[quant_name].to_numpy()+i2[quant_name].to_numpy())/2)
            
            if o_i==2:
                
                sign = 1 if tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))==tuple(([i1.mix.iloc[0],i2.mix.iloc[0]])) else -1
                mix_values[tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))]+=list(sign*(log_ratios-offset))
                # mix_values[tuple([i1.mix.iloc[0],i2.mix.iloc[0]])]+=list(sign*(log_ratios-offset))
                prec_values[tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))]+=list(i1.prec)
                int_values[tuple(sorted([i1.mix.iloc[0],i2.mix.iloc[0]]))]+= list((i1[quant_name].to_numpy()+i2[quant_name].to_numpy())/2)
                
                if plot_type == "all":
                    kde=gaussian_kde(log_ratios-offset)
                    bins = np.linspace(-y_lim,y_lim,bin_num)
                    # ax_hist.plot(kde(bins),bins,color=org_colors[o_i])
                    # plt.fill_betweenx(bins,idx,idx+kde(bins),alpha=.5,color=org_colors[o_i],label=orgs[o_i])
                    _,__,bars =ax.hist(log_ratios-offset,bins,alpha=.5,density=True,color=ycolor,label=orgs[o_i], orientation="horizontal")
                    for bar in bars:
                            bar.set_x(idx)
                # org_offset+=1
                    
                    t_ratio = np.log2(theoretical_yeast_amounts[i1.mix.iloc[0]]/theoretical_yeast_amounts[i2.mix.iloc[0]])
                    ax.hlines(t_ratio,idx,idx+1)
                    # plt.hlines(np.mean(log_ratios)-offset,idx,idx+1,colors="r")
                    hist_labels.append(str(chan[0])+"/"+str(chan[1])+"\n"+i1.mix.iloc[0]+"/"+i2.mix.iloc[0])
                    # hist_labels.append(str(chan[0])+"/"+str(chan[1]))
      
    if plot_type=="all":
        ax.set_xticks(np.arange(len(hist_labels)),hist_labels)
        # ax.set_xlabel("Mixture Comparison",)
        # ax.set_title(f"Jmod Quant {name} (all)")
        # plt.hlines(np.unique(all_ratios),0,len(all_ratios))

    else:
        # counter_size = 1.1#0.9
        
        # y_lim = 4
        tick_offset = 0#.1
        count = 0
        counts = []
        hist_labels = []
        if ax is None:
            fig,ax = plt.subplots(figsize = (figsize[0]*scale,figsize[1]*scale))
        for idx,comb in enumerate(mix_values):
            if mix_values[comb]==[]:
                continue
            if only_A:
                if comb[0]!="A":
                    continue
            kde=gaussian_kde(mix_values[comb])
            bins = np.linspace(-y_lim,y_lim,bin_num)
            # ax_hist.plot(kde(bins),bins,color=org_colors[o_i])
            # plt.fill_betweenx(bins,count,count+ (kde(bins)/max(kde(bins))),alpha=.5,color=org_colors[o_i],label=orgs[o_i]+" ratio" if count==0 else "")
            t_ratio = np.log2(theoretical_yeast_amounts[comb[0]]/theoretical_yeast_amounts[comb[1]])
            vals,__,bars =ax.hist(mix_values[comb],bins,alpha=.5,density=True,color=ycolor,label=orgs[o_i]+" ratio" if count==0 else "", orientation="horizontal")
            for bar in bars:
                    bar.set_x(count)
            # end_of_line = min(counter_size*.95,(max(kde(bins)))*1.2)
            # end_of_line = max(counter_size*.95,(max(vals))*1.2)
            end_of_line = counter_size*.9#
            ax.hlines(t_ratio,count,count+ end_of_line,label="Mixing Ratio" if count==0 else "",linestyles="--",linewidths=2.5,colors="black")
            hist_labels.append(comb[0]+"/"+comb[1])
            counts.append(count)
            # count+=counter_size#
            count+=end_of_line*1.1
        # plt.xticks(range(len(hist_labels)),hist_labels)
        ax.set_xticks(np.array(counts)+tick_offset,hist_labels)
        ax.set_yticks(np.array([-4,-2,0,2]),np.array([-4,-2,0,2]))
        ax.set_xlim(0,count)
        # plt.legend(bbox_to_anchor=(1.1, 1.1))
        # ax.set_ylabel("Log$_2$ Ratios")
        ax.spines[['right', 'top']].set_visible(False)
    if return_vals and return_prec:
        return ax,[mix_values,mix_values_h], [prec_values,prec_values_h],[int_values,int_values_h]
    elif return_vals:
        return ax,[mix_values,mix_values_h]
    else:
        return ax



# fdx_files_5p = ["/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-05-12_5plex_200pg_20win_24nce_1_results/",
# "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-05-12_5plex_200pg_20win_24nce_2_results/",
# "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-05-12_5plex_200pg_20win_24nce_3_results/"]


# fdx_files_9p = [ "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-05-12_9plex_200pg_20win_24nce_1_results/",
# "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-05-12_9plex_200pg_20win_24nce_2_results/",
# "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-05-12_9plex_200pg_20win_24nce_3_results/"]
    
    
# fdx_files_5p = ["/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep11_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep02_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep03_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep04_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep05_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep06_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep07_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep08_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep09_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/5plex_200pg_RA3_rep10_results/"]

# fdx_files_9p = ["/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep01_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep02_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep03_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep04_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep05_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep06_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep07_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep08_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep09_results/",
# "/Volumes/Lab/HS/Revision/kelly_lab/JModResults/9plex_200pg_RA4_rep10_results/"]


# fdx_files_lf = [str(i)+"/" for i in Path("/Volumes/Lab/HS/Revision/kelly_lab/JModResults").glob("LF*_results")]
# fdx_files_lf = ["/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-11-03_LF-A-200pg_40win_1_results/",
#                 "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-11-03_LF-D-200pg_40win_1_results/",
#                 "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-11-03_LF-B-200pg_40win_1_results/",
#                 "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-11-03_LF-C-200pg_40win_1_results/",
#                 "/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/JModResults/2025-11-03_LF-E-200pg_40win_1_results/"]

fdx_files_lf = [str(i)+"/" for i in Path("/Volumes/Lab/HS/Revision/pti_new/Astral/lf_hy_200pg_mixAD/JModResults").glob("*_results")]
fdx_files_5p = [str(i)+"/" for i in Path("/Volumes/Lab/HS/Revision/kelly_lab/JModResults").glob("5plex*_results")]
fdx_files_9p = [str(i)+"/" for i in Path("/Volumes/Lab/HS/Revision/kelly_lab/JModResults").glob("9plex*_results")]


theoretical_ratios = [65/65,20/5,15/30]
theoretical_yeast_amounts = {"A":6,"B":12,"C":18,"D":24,"E":8}
# theoretical_yeast_amounts = {"A":1,"B":4,"C":8,"D":16,"E":32}
# theoretical_human_amounts = {"A":50,"B":50,"C":50,"D":50}
theoretical_human_amounts = {"A":100,"B":100,"C":100,"D":100,"E":100}

channels = [0,2,4,6,8,10,12,14,16]
channel_mix = ["D","A","B","C","D","C","A","B","C"]
tag_name_tag6 = "plexDIA"#"diethyl_6plex"
channel_mix_key = {i:j for i,j in zip(channels,channel_mix)}
lf_mix_channel_key = {i:j for i,j in zip(["A","B","C","D","E"],range(0,17,4))}


fdxs = [load_jmod_files(fdx_files,channels,channel_mix,file_name="filtered_IDs.parquet")[0] for fdx_files in [
                                                                     fdx_files_lf,
                                                                     fdx_files_5p,
                                                                     fdx_files_9p]]

fdxs[0]["prec"] = fdxs[0]["seq"].astype(str) + fdxs[0]["z"].astype(int).astype(str)
fdxs[0].loc[:,"mix"] = [re.findall("LF-([ABCDE])",i)[0] for i in fdxs[0]["Run"]]
fdxs[0].loc[:,"channel"] = [lf_mix_channel_key[i] for i in fdxs[0]["mix"]]

fdxs[1]["prec"] = fdxs[1]["seq"].astype(str) + fdxs[1]["z"].astype(int).astype(str)
fdxs[2]["prec"] = fdxs[2]["seq"].astype(str) + fdxs[2]["z"].astype(int).astype(str)

# diann_5p = pd.read_parquet("//Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/5-plex/report.parquet")
# diann_lf = pd.read_parquet("/Volumes/Lab/HS/Revision/preprint_re-search/DIA/singlecelldata-astral/lf/mixABCDE/report.parquet")
diann_lf = pd.read_parquet("/Volumes/Lab/HS/Revision/pti_new/Astral/lf_hy_200pg_mixAD/report.parquet")
diann_lf = diann_lf.query("~Run.str.contains('blank')").reset_index()

diann_5p = pd.read_parquet("/Volumes/Lab/HS/Revision/kelly_lab/5-plex/report.parquet")
diann_5p = diann_5p.query("~Run.str.contains('blank')").reset_index()


# from FigureFunctions import orgs, all_fasta_seqs
# diann_report=diann_tag6
# for diann_report in [diann_lf,diann_tag6]:
for diann_report in [diann_lf,diann_5p]:
# for diann_report in [diann_lf]:
    # diann_report = diann_report.query("~Run.str.contains('blank')").reset_index()
    d_tag_name = "tag"
    diann_report["untag_prec"] = [re.sub(r"(\d+)",r"_\g<1>",re.sub(f"\({d_tag_name}\)","",i)) for i in diann_report['Precursor.Id']]
    diann_report["last_aa"] = [i[-1] for i in diann_report["Stripped.Sequence"]]
    if "Channel" in diann_report.columns and len(np.unique(diann_report["Channel"]))>1:
        diann_report["channel"] = [int(re.findall("(\d+)",i)[0]) for i in diann_report['Channel']]
    else:
        # diann_report["mix"] = [re.findall("LF-([AD])-200",i)[0] for i in diann_report["Run"]]
        # diann_report["channel"] = [0 if i=="D" else 4 for i in diann_report["mix"]]
        diann_report["mix"] = [re.findall("LF-([ABCDE])_",i)[0] for i in diann_report["Run"]]
        diann_report["channel"] = [lf_mix_channel_key[i] for i in diann_report["mix"]]

        # diann_report["channel"] = 0
    diann_report["org"] = np.array([";".join(orgs[[i in all_fasta_seqs[j] for j in range(3)]]) for i in diann_report["Stripped.Sequence"]])
    # diann_report["protein"] = diann_report['Protein.Names']
    diann_report["protein"] = diann_report['Protein.Group']
    diann_report["Qvalue"]=diann_report["Translated.Q.Value"]
    diann_report["Protein_Qvalue"]=diann_report['PG.Q.Value']
    diann_report["Lib_Protein_Qvalue"]=diann_report['Lib.PG.Q.Value']
    diann_report["BestChannel_Qvalue"]=diann_report.groupby(['Precursor.Id',"Run"])["Q.Value"].transform("min")
    diann_report["BestProtein_Qvalue"] = diann_report.groupby(["protein"])["Protein_Qvalue"].transform("min")

    diann_report["seq"] = diann_report["Modified.Sequence"]
    diann_report["z"] = diann_report["Precursor.Charge"]
    diann_report["prec"] = diann_report["seq"].astype(str) + diann_report["z"].astype(int).astype(str)







################################################################################
################################################################################
#### filter data 
### Keep good data and exclude problem LF data
D_runs = set([i for i in set(fdxs[0].Run) if "-D_" in i and int(i[-1])%2==1])

## keep random subset of A
all_A_runs = sorted([i for i in set(fdxs[0].Run) if "-A_" in i])
np.random.seed(0)
A_runs = set(np.random.choice(all_A_runs,7,replace=False))

### Keep random subset of 5plex and 9plex runs
all_fiveplex_runs = sorted(set(fdxs[1].Run))
np.random.seed(0)
fiveplex_runs = set(np.random.choice(all_fiveplex_runs,7,replace=False))

all_noneplex_runs = sorted(set(fdxs[2].Run))
np.random.seed(0)
nineplex_runs = set(np.random.choice(all_noneplex_runs,7,replace=False))

fdxs[0] = fdxs[0].query("Run in @D_runs or Run in @A_runs")
fdxs[1] = fdxs[1].query("Run in @fiveplex_runs")
fdxs[2] = fdxs[2].query("Run in @nineplex_runs")

diann_lf = diann_lf.query("Run in @D_runs or Run in @A_runs")
diann_5p = diann_5p.query("Run in @fiveplex_runs")

################################################################################
################

################################################################




group_keys = ["LF", "5plex","9plex"]


groups = [[i for i in list(np.unique(fdx.Run))] for fdx in fdxs]

start_idx = 0
fdxs_merged = [merge_and_norm(["plex_Area"], fdx,file_names, [file_groups]) for fdx,file_names,file_groups in zip(fdxs[start_idx:],group_keys[start_idx:],groups[start_idx:])]
fdxs_merged[0]['Precursor.Normalised'] = fdxs_merged[0]["merged_plex_Area"]
fdxs_merged[1]['Precursor.Normalised'] = fdxs_merged[1]["merged_plex_Area"]
fdxs_merged[2]['Precursor.Normalised'] = fdxs_merged[2]["merged_plex_Area"]


diann_merged = [merge_and_norm(['Precursor.Normalised'],diann,file_names, [file_groups]) for diann,file_names,file_groups in zip([diann_lf,diann_5p][start_idx:],
                                                                                                                                 group_keys[start_idx:],
                                                                                                                                 groups[start_idx:])]
diann_merged[1]["mix"] = [channel_mix[channels.index(i)] for i in diann_merged[1].channel]
diann_merged[0]["merged_plex_Area"] = diann_merged[0]['merged_Precursor.Normalised']
diann_merged[1]["merged_plex_Area"] = diann_merged[1]['merged_Precursor.Normalised']



###########################################################################
###########################################################################
###########################################################################
###########################################################################


diann_base = np.array((255,99,71))/255
jmod_base = np.array((100,149,237))/255


plt.rcParams['axes.labelsize'] = 25
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['ytick.labelsize'] = 18##default ="medium"
plt.rcParams['xtick.labelsize'] = 15##default ="medium"


quant_name = 'merged_plex_Area'

max_y = 8
bin_num=25
counter_size = 2
plot_width = 15
scale=2.5
ylims = (-7,2)#(-4.1,1.)#
n=3
fig = plt.figure(figsize=(4*scale,1*scale*n))
gs = gridspec.GridSpec(n,1,wspace=0.05,hspace=.55)
axs = []
for g in gs:
    
    ax = fig.add_subplot(g)
    if len(axs)>0:
        ax.sharex(axs[0])
    rect = plt.Rectangle(
            (0, 1.0),  # (x, y) position relative to the subplot
            1,  # Full width
            0.25,  # Height of the box ## prev=0.3
            edgecolor="black",
            transform=ax.transAxes, facecolor='lightgrey', alpha=0.6, clip_on=False, zorder=1
        )
    ax.add_patch(rect)
    axs.append(ax)



ax,mix_ratios,prec_ratios,int_vals =compute_all_comparisons(
                                            fdx3 =fdxs_merged[1][fdxs_merged[1].Qvalue<.01],
                                            fdx2 = fdxs_merged[2][fdxs_merged[2].Qvalue<.01],
                                            fdx = diann_merged[1][diann_merged[1]["Channel.Q.Value"]<.01],
                                            quant_name=quant_name,
                                                # quant_name='merged_Precursor.Quantity',
                                                    # quant_name='merged_new_plexfitMS1_p',
                                                # quant_name='merged_plex_Area',
                                                channels=channels[::2],
                                                channel_mix=channel_mix,
                                                theoretical_yeast_amounts=theoretical_yeast_amounts,
                                                plot_type = "allm",
                                                Q_value_limit=.01,
                                                # BestChannel_Qvalue_limit=0.005,
                                                at_least_b=0,
                                                bin_num=bin_num,
                                                counter_size=counter_size,
                                                scale=2,
                                                ax=axs[0],
                                                only_A=True,
                                                ycolor=diann_base,
                                                return_vals=True,
                                                return_prec=True,
                                                y_lim=max_y
                                                )

ax.grid()
ax.set_ylim(*ylims)
# plt.title("DIA-NN 5-plex (4 Da)")
ax.set_title("DIA-NN PSMtag (5-plex; 4 Da)")
ax.set_title("DIA-NN 5-plex (4 Da)")

d1 = pd.concat([pd.concat([pd.DataFrame({"Alg":["DIANN"]*len(mix_ratios[0][i]),
                                         "Org":["Yeast"]*len(mix_ratios[0][i]),
                                         "plex":["5plex"]*len(mix_ratios[0][i]),
                                         "prec":prec_ratios[0][i],
                                         "intensity":int_vals[0][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                         "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                pd.concat([pd.DataFrame({"Alg":["DIANN"]*len(mix_ratios[1][i]),
                                         "Org":["Human"]*len(mix_ratios[1][i]),
                                         "plex":["5plex"]*len(mix_ratios[1][i]),
                                         "prec":prec_ratios[1][i],
                                         "intensity":int_vals[1][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                         "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])

ax,mix_ratios,prec_ratios,int_vals =compute_all_comparisons(
                                            fdx =fdxs_merged[1][fdxs_merged[1].Qvalue<.01],
                                            fdx2 = fdxs_merged[2][fdxs_merged[2].Qvalue<.01],
                                            fdx3 = diann_merged[1][diann_merged[1]["Channel.Q.Value"]<.01],
                                               # quant_name='merged_Precursor.Quantity',
                                                    # quant_name='merged_new_plexfitMS1_p',
                                                          # quant_name='merged_plexfitMS1',
                                                quant_name=quant_name,
                                               channels=channels[::2],
                                               channel_mix=channel_mix,
                                               theoretical_yeast_amounts=theoretical_yeast_amounts,
                                               plot_type = "allm",
                                               Q_value_limit=.01,
                                               # BestChannel_Qvalue_limit=0.005,
                                               at_least_b=0,
                                               bin_num=bin_num,
                                                counter_size=counter_size,
                                                scale=2,
                                                ax=axs[1],
                                                only_A=True,
                                                ycolor=jmod_base,
                                                return_vals=True,
                                                return_prec=True,
                                                y_lim=max_y
                                               )
ax.grid()
ax.set_ylim(*ylims)
ax.set_title("JMod PSMtag (5-plex; 4 Da)")
ax.set_title("JMod 5-plex (4 Da)")

# d2 = pd.concat([pd.DataFrame({"Alg":["JMod"]*len(ax_mix_ratios[1][i]),"plex":["5plex"]*len(ax_mix_ratios[1][i]),"Combination":["/".join(i)]*len(ax_mix_ratios[1][i]),"Ratio":ax_mix_ratios[1][i]}) for i in list(ax_mix_ratios[1])])

d2 = pd.concat([pd.concat([pd.DataFrame({"Alg":["JMod"]*len(mix_ratios[0][i]),
                                         "Org":["Yeast"]*len(mix_ratios[0][i]),
                                         "plex":["5plex"]*len(mix_ratios[0][i]),
                                         "prec":prec_ratios[0][i],
                                         "intensity":int_vals[0][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                         "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                pd.concat([pd.DataFrame({"Alg":["JMod"]*len(mix_ratios[1][i]),
                                         "Org":["Human"]*len(mix_ratios[1][i]),
                                         "plex":["5plex"]*len(mix_ratios[1][i]),
                                         "prec":prec_ratios[1][i],
                                         "intensity":int_vals[1][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                         "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])

ax,mix_ratios,prec_ratios,int_vals =compute_all_comparisons(
                                            fdx2 =fdxs_merged[1][fdxs_merged[1].Qvalue<.01],
                                            fdx = fdxs_merged[2][fdxs_merged[2].Qvalue<.01],
                                            fdx3 = diann_merged[1][diann_merged[1]["Channel.Q.Value"]<.01],
                                               # quant_name='merged_Precursor.Quantity',
                                                    # quant_name='merged_new_plexfitMS1_p',
                                                          # quant_name='merged_plexfitMS1',
                                                quant_name=quant_name,
                                               channels=channels[::2],
                                               channel_mix=channel_mix,
                                               theoretical_yeast_amounts=theoretical_yeast_amounts,
                                               plot_type = "allm",
                                               Q_value_limit=.01,
                                               # BestChannel_Qvalue_limit=0.005,
                                               at_least_b=0,
                                               bin_num=bin_num,
                                                counter_size=counter_size,
                                                scale=2,
                                                ax=axs[2],
                                                only_A=True,
                                                ycolor=jmod_base,
                                                return_vals=True,
                                                return_prec=True,
                                                y_lim=max_y
                                               )
ax.grid()
ax.set_ylim(*ylims)
ax.set_title("JMod PSMtag (9-plex; 2 Da)")
ax.set_title("JMod 9-plex (2 Da)")

fig.text(.05, 0.5, "log$_2$ Yeast Precursor Ratios", va='center', rotation='vertical',fontsize=20)

d3 = pd.concat([pd.concat([pd.DataFrame({"Alg":["JMod"]*len(mix_ratios[0][i]),
                                         "Org":["Yeast"]*len(mix_ratios[0][i]),
                                         "plex":["9plex"]*len(mix_ratios[0][i]),
                                         "prec":prec_ratios[0][i],
                                         "intensity":int_vals[0][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                         "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                pd.concat([pd.DataFrame({"Alg":["JMod"]*len(mix_ratios[1][i]),
                                         "Org":["Human"]*len(mix_ratios[1][i]),
                                         "plex":["9plex"]*len(mix_ratios[1][i]),
                                         "prec":prec_ratios[1][i],
                                         "intensity":int_vals[1][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                         "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])



pd.concat([d1,d2,d3],ignore_index=True).to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/Combined_Ratios_PSMtag_RKdata_QuantUMS.csv",index=False)



###### compute all comparisons for all 9plex runs individually

all_runs = []
for a,plex in [["JMod","9plex"],["JMod","5plex"],["DIANN","5plex"]]:
    if a=="JMod" and plex=="9plex":
        fdx_all = fdxs[2].query("Qvalue<.01 and plex_Area>0")
        runs = set(fdx_all.Run)
        fdx_n = fdxs[2].query("Qvalue<.01 and Run==@run and plex_Area>0")
        quant_name = "plex_Area"
        
    elif a=="JMod" and plex=="5plex":
        fdx_all = fdxs[1].query("Qvalue<.01 and plex_Area>0")
        quant_name = "plex_Area"
        
# for run in set(diann_5p.Run):
    elif  a=="DIANN" and plex=="5plex":
        fdx_all = diann_5p.query("`Q.Value`<.01 and `Precursor.Normalised`>0")
        quant_name = "Precursor.Normalised"

        fdx_all["mix"] = [channel_mix[channels.index(i)] for i in fdx_all.channel]

    all_runs = []
    runs = set(fdx_all.Run)
    for run in runs:
        fdx_n = fdx_all.query("Run==@run")
        # fig,ax = plt.subplots()
        ax,mix_ratios, prec_ratios,int_vals = compute_all_comparisons(fdx=fdx_n,
                                                                       # quant_name='merged_Precursor.Quantity',
                                                                            # quant_name='merged_new_plexfitMS1_p',
                                                                                  # quant_name='merged_plexfitMS1',
                                                                        quant_name=quant_name,
                                                                       channels=channels[::2],
                                                                       channel_mix=channel_mix,
                                                                       theoretical_yeast_amounts=theoretical_yeast_amounts,
                                                                       plot_type = "allm",
                                                                       Q_value_limit=.01,
                                                                       # BestChannel_Qvalue_limit=0.005,
                                                                       at_least_b=0,
                                                                       bin_num=bin_num,
                                                                        counter_size=counter_size,
                                                                        scale=2,
                                                                        ax=ax,
                                                                        only_A=True,
                                                                        ycolor=jmod_base,
                                                                        return_vals=True,
                                                                        return_prec=True,
                                                                        y_lim=max_y
                                                                       )
        dn = pd.concat([pd.concat([pd.DataFrame({"Alg":[a]*len(mix_ratios[0][i]),
                                                 "Org":["Yeast"]*len(mix_ratios[0][i]),
                                                 "plex":[plex]*len(mix_ratios[0][i]),
                                                 "prec":prec_ratios[0][i],
                                                 "intensity":int_vals[0][i],
                                                 "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                                 "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                        pd.concat([pd.DataFrame({"Alg":[a]*len(mix_ratios[1][i]),
                                                 "Org":["Human"]*len(mix_ratios[1][i]),
                                                 "plex":[plex]*len(mix_ratios[1][i]),
                                                 "prec":prec_ratios[1][i],
                                                 "intensity":int_vals[1][i],
                                                 "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                                 "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])
        dn["Run"] = run
        
        all_runs.append(dn)
    all_runs_df = pd.concat(all_runs)
    all_runs_df.prec
    
    plt.subplots()
    sb.boxplot(dn,x="Combination",y="Ratio",hue="Org",showfliers=False)
    plt.grid()

    all_runs_df.to_csv(f"/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/{a}{plex}_All_Run_Ratios_RKdata.csv")

sb.boxplot(all_runs_df.query("Combination=='A/D'"),x="Run",y="Ratio",hue="Org",showfliers=False)
plt.grid()


# all_runs_df.to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/JMod5plex_All_Run_Ratios_RKdata.csv")
# all_runs_df.to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/JMod9plex_All_Run_Ratios_RKdata.csv")
# all_runs_df.to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/DIANN5plex_All_Run_Ratios_RKdata.csv")






##########################################################################################

diann_base = np.array((255,99,71))/255
jmod_base = np.array((100,149,237))/255


plt.rcParams['axes.labelsize'] = 25
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['ytick.labelsize'] = 18##default ="medium"
plt.rcParams['xtick.labelsize'] = 15##default ="medium"


quant_name = 'merged_plex_Area'
channels_subset = [0,12]
max_y = 8
bin_num=25
counter_size = 2
plot_width = 15
scale=2.5
ylims = (-7,2)#(-4.1,1.)#
n=3
fig = plt.figure(figsize=(4*scale,1*scale*n))
gs = gridspec.GridSpec(n,1,wspace=0.05,hspace=.55)
axs = []
for g in gs:
    
    ax = fig.add_subplot(g)
    if len(axs)>0:
        ax.sharex(axs[0])
    rect = plt.Rectangle(
            (0, 1.0),  # (x, y) position relative to the subplot
            1,  # Full width
            0.25,  # Height of the box ## prev=0.3
            edgecolor="black",
            transform=ax.transAxes, facecolor='lightgrey', alpha=0.6, clip_on=False, zorder=1
        )
    ax.add_patch(rect)
    axs.append(ax)



ax,mix_ratios,prec_ratios,int_vals =compute_all_comparisons(
                                            fdx2 = fdxs_merged[0][fdxs_merged[0].Qvalue<.01],
                                            fdx = diann_merged[0][diann_merged[0]["Channel.Q.Value"]<.01],
                                            quant_name=quant_name,
                                                # quant_name='merged_Precursor.Quantity',
                                                    # quant_name='merged_new_plexfitMS1_p',
                                                # quant_name='merged_plex_Area',
                                                channels=channels_subset,
                                                channel_mix=channel_mix,
                                                theoretical_yeast_amounts=theoretical_yeast_amounts,
                                                plot_type = "allm",
                                                Q_value_limit=.01,
                                                # BestChannel_Qvalue_limit=0.005,
                                                at_least_b=0,
                                                bin_num=bin_num,
                                                counter_size=counter_size,
                                                scale=2,
                                                ax=axs[0],
                                                only_A=True,
                                                ycolor=diann_base,
                                                return_vals=True,
                                                return_prec=True,
                                                y_lim=max_y
                                                )

ax.grid()
ax.set_ylim(*ylims)
# plt.title("DIA-NN 5-plex (4 Da)")
ax.set_title("DIA-NN PSMtag (5-plex; 4 Da)")
ax.set_title("DIA-NN 5-plex (4 Da)")

d1 = pd.concat([pd.concat([pd.DataFrame({"Alg":["DIANN"]*len(mix_ratios[0][i]),
                                         "Org":["Yeast"]*len(mix_ratios[0][i]),
                                         "plex":["5plex"]*len(mix_ratios[0][i]),
                                         "prec":prec_ratios[0][i],
                                         "intensity":int_vals[0][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                         "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                pd.concat([pd.DataFrame({"Alg":["DIANN"]*len(mix_ratios[1][i]),
                                         "Org":["Human"]*len(mix_ratios[1][i]),
                                         "plex":["5plex"]*len(mix_ratios[1][i]),
                                         "prec":prec_ratios[1][i],
                                         "intensity":int_vals[1][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                         "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])

ax,mix_ratios,prec_ratios,int_vals =compute_all_comparisons(
                                            fdx =fdxs_merged[0][fdxs_merged[0].Qvalue<.01],
                                            fdx2 = diann_merged[0][diann_merged[0]["Channel.Q.Value"]<.01],
                                               # quant_name='merged_Precursor.Quantity',
                                                    # quant_name='merged_new_plexfitMS1_p',
                                                          # quant_name='merged_plexfitMS1',
                                                quant_name=quant_name,
                                               channels=channels_subset,
                                               channel_mix=channel_mix,
                                               theoretical_yeast_amounts=theoretical_yeast_amounts,
                                               plot_type = "allm",
                                               Q_value_limit=.01,
                                               # BestChannel_Qvalue_limit=0.005,
                                               at_least_b=0,
                                               bin_num=bin_num,
                                                counter_size=counter_size,
                                                scale=2,
                                                ax=axs[1],
                                                only_A=True,
                                                ycolor=jmod_base,
                                                return_vals=True,
                                                return_prec=True,
                                                y_lim=max_y
                                               )
ax.grid()
ax.set_ylim(*ylims)
ax.set_title("JMod PSMtag (5-plex; 4 Da)")
ax.set_title("JMod 5-plex (4 Da)")

# d2 = pd.concat([pd.DataFrame({"Alg":["JMod"]*len(ax_mix_ratios[1][i]),"plex":["5plex"]*len(ax_mix_ratios[1][i]),"Combination":["/".join(i)]*len(ax_mix_ratios[1][i]),"Ratio":ax_mix_ratios[1][i]}) for i in list(ax_mix_ratios[1])])

d2 = pd.concat([pd.concat([pd.DataFrame({"Alg":["JMod"]*len(mix_ratios[0][i]),
                                         "Org":["Yeast"]*len(mix_ratios[0][i]),
                                         "plex":["5plex"]*len(mix_ratios[0][i]),
                                         "prec":prec_ratios[0][i],
                                         "intensity":int_vals[0][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                         "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                pd.concat([pd.DataFrame({"Alg":["JMod"]*len(mix_ratios[1][i]),
                                         "Org":["Human"]*len(mix_ratios[1][i]),
                                         "plex":["5plex"]*len(mix_ratios[1][i]),
                                         "prec":prec_ratios[1][i],
                                         "intensity":int_vals[1][i],
                                         "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                         "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])


combined = pd.concat([d1,d2])
combined.prec

sb.boxplot(combined,x="Combination",y="Ratio",hue="Org",showfliers=False)
plt.grid()

sb.relplot(combined.query("Alg=='DIANN'").merge(diann_report.query("Run=='2026-10-19_LF-A_200pg_1'")[["prec","RT"]]),
           x="RT",y="Ratio",s=1)
plt.grid()

sb.boxplot(combined.query("Combination=='A/D'"),x="Alg",y="Ratio",hue="Org",showfliers=False)
plt.grid()

combined.to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/Combined_Ratios_internalLF_AD.csv",index=False)





#### LF All runs

all_lf_runs = set(fdxs[0].Run)
A_runs = set([i for i in all_lf_runs if "-A_" in i])
D_runs = set([i for i in all_lf_runs if "-D_" in i])


import warnings
import pandas as pd



channels_subset = [0,12]

for a in ["JMod","DIANN"]:
    all_runs = []
    for comb in tqdm.tqdm(list(product(A_runs,D_runs))):
        # continue
    
        run = "_".join(["".join(*re.findall("([ABCDE])_200pg_(\d+)",i)) for i in comb])
        
        if a=="JMod":
            fdx_n = fdxs[0].query("Run in @comb").query("plex_Area>0").copy()
            quant_name = "plex_Area"
            plex="LF";
        
        if a=="DIANN":
            fdx_n = diann_lf.query("Run in @comb").query("`Precursor.Normalised`>0")
            quant_name = "Precursor.Normalised"
            plex="5plex";
        
        fdx_n["Run"] = 1
        
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
        
        # fig,ax = plt.subplots()
            ax,mix_ratios, prec_ratios,int_vals = compute_all_comparisons(fdx=fdx_n,
                                                                           # quant_name='merged_Precursor.Quantity',
                                                                                # quant_name='merged_new_plexfitMS1_p',
                                                                                      # quant_name='merged_plexfitMS1',
                                                                            quant_name=quant_name,
                                                                           channels=channels_subset,
                                                                           channel_mix=channel_mix,
                                                                           theoretical_yeast_amounts=theoretical_yeast_amounts,
                                                                           plot_type = "allm",
                                                                           Q_value_limit=.01,
                                                                           # BestChannel_Qvalue_limit=0.005,
                                                                           at_least_b=0,
                                                                           bin_num=bin_num,
                                                                            counter_size=counter_size,
                                                                            scale=2,
                                                                            ax=ax,
                                                                            only_A=True,
                                                                            ycolor=jmod_base,
                                                                            return_vals=True,
                                                                            return_prec=True,
                                                                            y_lim=max_y
                                                                           )
        dn = pd.concat([pd.concat([pd.DataFrame({"Alg":[a]*len(mix_ratios[0][i]),
                                                 "Org":["Yeast"]*len(mix_ratios[0][i]),
                                                 "plex":[plex]*len(mix_ratios[0][i]),
                                                 "prec":prec_ratios[0][i],
                                                 "intensity":int_vals[0][i],
                                                 "Combination":["/".join(i)]*len(mix_ratios[0][i]),
                                                 "Ratio":mix_ratios[0][i]}) for i in list(mix_ratios[0])]),
                        pd.concat([pd.DataFrame({"Alg":[a]*len(mix_ratios[1][i]),
                                                 "Org":["Human"]*len(mix_ratios[1][i]),
                                                 "plex":[plex]*len(mix_ratios[1][i]),
                                                 "prec":prec_ratios[1][i],
                                                 "intensity":int_vals[1][i],
                                                 "Combination":["/".join(i)]*len(mix_ratios[1][i]),
                                                 "Ratio":mix_ratios[1][i]}) for i in list(mix_ratios[1])])])
        dn["Run"] = run
        
        all_runs.append(dn)
        
        
    all_runs_df = pd.concat(all_runs)
    all_runs_df.prec
    
    sb.boxplot(all_runs_df,x="Combination",y="Ratio",hue="Org",showfliers=False)
    plt.grid()
    all_runs_df.to_csv(f"/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/{a}_LF_All_Run_Ratios_internalLF_AD.csv",index=False)


# all_runs_df.to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/JMod_LF_All_Run_Ratios_internalLF_AD.csv",index=False)
# all_runs_df.to_csv("/Volumes/Lab/KMD/For_HS/PSMtag_Revision/SevenReps/DIANN_LF_All_Run_Ratios_internalLF_AD.csv",index=False)



sb.relplot(all_runs_df.query("Alg=='DIANN'").query("Run=='A5_D5'").merge(diann_report.query("Run=='2026-10-19_LF-A_200pg_1'")[["prec","RT","last_aa"]]),
           x="RT",y="Ratio",s=1,hue="Org")
plt.grid()



sb.relplot(combined.query("Alg=='JMod'").merge(diann_report.query("Run=='2026-10-19_LF-A_200pg_1'")[["prec","RT","last_aa"]]),
           x="RT",y="Ratio",s=1,hue="Org")
plt.grid()


