import numpy as np
import pandas as pd



def split_col(a,type=float):
    return([np.array(i.split(';')).astype(float).astype(type) for i in a])


#a = period_dat.donor_copies_date.values
def split_pad_col(a,type=float):
    split = [np.array(i.split(';')).astype(float).astype(type) for i in a]
    q = max([i.shape[0] for i in split])
    padded = [np.pad(i, (0, q-i.shape[0]), mode='constant', constant_values=np.nan) for i in split]
    return(np.vstack(padded))


def get_trans_form(fp):
    with open(fp+'.log', 'r') as f:
        for line in f.readlines():
            if line[:9] == 'Namespace':
                trans_form = [j.split('=')[1].replace("'", "") for j in line.split(', ') if j[:9] == "transForm"][0]
                return(trans_form)
                break

                
def get_ind_n_trans_bd(dr, t=1, v=np.linspace(0,7,1000), form='ve'):
    if form == 've':
        # number of expected transmisisons over a given time period from someone with a given viral load
        # asusmes no "saturation" of transmission e.g. fully susceptible, well-mixed population
        # thus, 1 person over 100 years is equivalent to 100 people over 1 year
        return(t * 108 * dr.use_param1.values * (1 - (1 - dr.use_param2.values)**(dr.use_param3.values * 10**v[:,np.newaxis])))
    elif form == 'sigmoid':
        return(t * dr.use_param1.values / (1 + 10**(-dr.use_param2.values*(v[:,np.newaxis] - dr.use_param3.values))))


def get_couple_prob_trans_bd(dr, t=1, v=np.linspace(0,7,1000), form='ve'):
    if form == 've':
        # probability of intra-couple transmission occurring over time t 
        # account for saturation (i.e. transmission can only occur once)
        return(1 - np.exp(-t * 108 * dr.use_param1.values * (1 - (1 - dr.use_param2.values)**(dr.use_param3.values * 10**v[:,np.newaxis]))))
    elif form == 'sigmoid':
        return(1 - np.exp(-t * dr.use_param1.values / (1 + 10**(-dr.use_param2.values*(v[:,np.newaxis] - dr.use_param3.values)))))



def summarize_draws(dr_to_sum):
    return(pd.DataFrame(np.quantile(dr_to_sum,
        [0.025, 0.25, 0.5, 0.75, 0.975],
        axis=1).T,
        columns=[0.025, 0.25, 0.5, 0.75, 0.975]))


def import_config(config_path):
    import pandas as pd
    return(
        {i[0]:[float(j) for j in i[1].split(';')] if ';' in i[1] else i[1] for idx, i in 
            pd.read_csv(config_path, header=None).iterrows()})


def format_list_col(x):
    fmt_x = [np.array(i.replace('[', '').lstrip().replace(']', '').split()).astype(np.dtypes.StringDType()) for i in x]
    return(np.array([np.pad(i, (0, 2-len(i)), constant_values=np.nan) for i in fmt_x]))



def map_arr(a, d):
    u,inv = np.unique(a,return_inverse = True)
    a_map = np.array([d[x] for x in u])[inv].reshape(a.shape)
    return(a_map)   

    
def format_map_dict(map_files):
    map_dfs = [pd.read_csv(i, sep='\t', index_col=0) for i in map_files]
    map_dict = {}
    for seg_map in map_dfs:
        for name, name_map in seg_map.iteritems():
            map_dict[name] = {int(i):idx for idx, i in name_map.iteritems() if not pd.isna(i)} 
    return(map_dict)


def convert_to_datetime(dt, fmt):
    # pandas to_datetime infers the day and month when it's missing
    # we want to return NAN when that is the case
    from datetime import datetime
    import numpy as np
    try: 
        return(datetime.strptime(dt, fmt))
    except:
        return(np.nan)

        
def datetime_from_numeric(numdate):
    # borrowed from the treetime utilities
    # https://github.com/neherlab/treetime/blob/de6947685fbddc758e36fc4008ddd5f9d696c6d3/treetime/utils.py
    """convert a numeric decimal date to a python datetime object
    Note that this only works for AD dates since the range of datetime objects
    is restricted to year>1.
    Parameters
    ----------
    numdate : float
        numeric date as in 2018.23
    Returns
    -------
    datetime.datetime
        datetime object
    """
    from calendar import isleap
    import datetime
    days_in_year = 366 if isleap(int(numdate)) else 365
    # add a small number of the time elapsed in a year to avoid
    # unexpected behavior for values 1/365, 2/365, etc
    days_elapsed = int(((numdate%1)+1e-10)*days_in_year) -1
    date = datetime.datetime(int(numdate),1,1) + datetime.timedelta(days=days_elapsed)
    return date


# move to utils
# from treetime
# https://github.com/neherlab/treetime/blob/1e378bfbbb98451ce4d6dd309d3699323102ee64/treetime/utils.py
def numeric_from_datetime(dt):
    from calendar import isleap
    import datetime
    if dt is None:
        dt = datetime.datetime.now()
    days = 366 if isleap(dt.year) else 365
    # removing the 0.5 adjustment in the standard treetime code 
    # to align with tajimas D inference method
    #res = dt.year + (dt.timetuple().tm_yday-0.5) / days
    res =  dt.year + (dt.timetuple().tm_yday) / days
    return(res)



def filter_dvg_dat(dat, only_shared=False, min_reads=0, min_del=0, filter_plasmid=False):
    # replace nan support values with 0
    # if we want only shared DVGs
    # note: this *only* applies to samples with technical rpelicates
    # otherwise we just take the identified reads
    replicated_dat = dat[(~dat['Run_y'].isnull()) & 
        (~dat['Run_x'].isnull())].copy()
    unreplicated_dat = dat[(dat['Run_y'].isnull()) | 
        (dat['Run_x'].isnull())].copy()
    if only_shared:
        replicated_dat = replicated_dat[(~replicated_dat['Total_support_x'].isnull())
            & (~replicated_dat['Total_support_y'].isnull())]
    # get average read support for eplicated data
    replicated_dat['Total_support'] = \
        replicated_dat[['Total_support_x', 'Total_support_y']].fillna(0).mean(axis=1)
    replicated_dat['Rel_support'] = \
        replicated_dat[['Rel_support_x', 'Rel_support_y']].fillna(0).mean(axis=1)
    # for unreplicated dat, we just take raw support values
    unreplicated_dat['Total_support'] = \
        unreplicated_dat[['Total_support_x', 'Total_support_y']].fillna(0).max(axis=1)
    unreplicated_dat['Rel_support'] = \
        unreplicated_dat[['Rel_support_x', 'Rel_support_y']].fillna(0).max(axis=1)
    # now combine and filter
    filtered_dat = pd.concat([replicated_dat, unreplicated_dat])
    filtered_dat = \
        filtered_dat[(filtered_dat['Total_support'] >= min_reads)
        & (filtered_dat['Stop'] - filtered_dat['Start'] >= min_del)]
    if filter_plasmid == True: 
        filtered_dat = filtered_dat[(filtered_dat['plasmid_x'] != True) & 
            (filtered_dat['plasmid_y'] != True)]
    return(filtered_dat)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def import_aln(fh):
    s_names = []
    all_s = ''
    #fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            all_s += s
    fh.close()
    n_seqs = len(s_names)
    size = int(len(all_s)/n_seqs)
    s_arr = np.array(list(all_s.lower())).reshape((n_seqs, size))
    return(np.array(s_names), s_arr)


def import_seqs(fh):
    s_names = []
    s_list = []
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            s_list.append(s)
    fh.close()    
    return(np.array(s_names), s_list)


def format_seqs_arr(s, n_seqs):
    n_seqs = int(n_seqs)
    size = int(len(s)/n_seqs)
    seqs_arr = \
        np.frombuffer(s.lower().encode(), dtype=np.int8)
    seqs_arr = np.copy(seqs_arr)
    #[(97, 'A'), (114, 'R'), (119, 'W'), (109, 'M'), (100, 'D'), (104, 'H'), (118, 'V'), 
    #(110, 'N'), (99, 'C'), (121, 'Y'), (115, 'S'), (109, 'M'), (98, 'B'), (104, 'H'), 
    #(118, 'V'), (110, 'N'), (117, 'U'), (121, 'Y'), (119, 'W'), (107, 'K'), (98, 'B'), 
    #(100, 'D'), (104, 'H'), (110, 'N'), (103, 'G'), (114, 'R'), (115, 'S'), (107, 'K'), 
    #(98, 'B'), (100, 'D'), (118, 'V'), (110, 'N'), (116, 'T')]
    seqs_arr = \
        seqs_arr.reshape((n_seqs, int(seqs_arr.shape[0]/n_seqs)))
    return(seqs_arr)



def plot_style(grey='#333333'):
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'arial'
    mpl.rcParams['text.color'] = grey
    mpl.rcParams['axes.labelcolor'] = grey
    mpl.rcParams['xtick.color'] = '#707070'
    mpl.rcParams['ytick.color'] = '#707070'
    # Font sizes
    mpl.rcParams['figure.titlesize'] = 14
    mpl.rcParams['axes.titlesize'] = 14
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['xtick.labelsize'] = 14
    mpl.rcParams['ytick.labelsize'] = 14
    # Border colors
    mpl.rcParams['axes.edgecolor'] = grey
    mpl.rcParams['grid.color'] = '#eaeaea'
    #mpl.rcParams['grid.zorder'] = 0
    # Legend
    mpl.rcParams['legend.fontsize'] = 14
    mpl.rcParams['legend.frameon'] = False