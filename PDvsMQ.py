import numpy as npy

nan = npy.nan
import pandas as pds
import matplotlib.pyplot as plt
%matplotlib qt
import re
from itertools import groupby

#com = pds.read_table(r'D:\aaaaaWork\Python\Compare_exampledata.txt')
com = pds.read_table(r'D:\aaaaaWork\Data\PD_MQ_compare\20161007_Compare_PD_MQ_wt_inh_only.txt')

def mj_givehead( yourfile ):
    head = list()
    for entry in yourfile.iteritems():
        head.append(entry[0])
    return head

head = mj_givehead(com)
head


#C:\Anaconda\Scripts\ipython profile create mj 
#C:\Anaconda\Scripts\ipython --profile=mj

def mj_sgdict( fastafile ):

    '''Reads a Saccharomyces Genome Database (SGD) fasta file and
    enables protein summary information access by trivial or systematic name

    :param fastafile: file path of fasta file. 

    :returns: A dict with a protein's trivial and systematic name (uppercase!)
    as key to access SGD summary information
    (systematic name, trivial name, gene & protein description, sequence)
    '''

    reSys = '(?P<sys>(>\S+\s))'
    reTriv = '(?P<triv>(\S+\s))'
    reDescr = '(?P<descr>SGDID.+\n)'
    reSeq = '(?P<seq>[A-Z\*\n]+)'
    fastapattern = reSys + reTriv + reDescr + reSeq
    _line = lambda s: s.replace('\n', '').replace('\r', '').strip('*').strip()
    with open(fastafile, 'rb') as openfile: 
        #with: closes file at end of 'loop', also in case of error 
        #(when file.close() is not processed)
        readfile = openfile.read() 
        #reads the whole file character by charactrer, incl. \n
        
        sgdict = dict()

        reFastapattern = re.compile(fastapattern, re.VERBOSE) 
        reIters = reFastapattern.finditer(str(readfile)) 

        for iterable in reIters:
            sys = _line(iterable.group('sys').replace('>', ''))
            triv = _line(iterable.group('triv'))
            descr = _line(iterable.group('descr'))
            seq = _line(iterable.group('seq'))
            prot_info = {'sys': sys, 'triv': triv, 'descr': descr, 'seq': seq}
            sgdict[sys] = prot_info
            sgdict[triv] = prot_info
    return sgdict

fasta = r'D:\aaaaaWork\orf_trans_all_R64-1-1_20110203.fasta'
sgd = mj_sgdict(fasta)


#com.prot_site is a list of all found protsites
#next step separates protein and P-site into two lists
prots = list()
sites = list()
sitepattern = '(?P<pr>(\S+\s))(?P<site>\d.+)'
regex_sitepattern = re.compile(sitepattern, re.VERBOSE)
for entry in com.prot_site:
    sitepattern_iters = regex_sitepattern.finditer(entry)
    for iterable in sitepattern_iters:
        pr = iterable.group('pr').replace(',','').strip()
        site = iterable.group('site')
        prots.append(pr)
        sites.append(site)

#separates P-sites into site-position(s) and residue(s)
pos =list()
aa = list()
numbers = '(?P<number>\d{1,4})(?P<aa>\W\D\W)'
regex_numbers = re.compile(numbers, re.VERBOSE)
for entry in sites:
    number_iters = regex_numbers.finditer(entry)
    s = list()
    a = list()
    for iterable in number_iters:
        s.append(iterable.group('number'))
        sfloat = [int(x) for x in s]
        a.append(iterable.group('aa')[1])
    pos.append(sfloat)
    aa.append(a)

#makes a dict out of that info
#accessible by com.prot_site entry
#gives protname, siteposition and residue
sitedict = dict()
for index, entry in enumerate(com.prot_site):
    sitedict[entry] = {'prot': prots[index], 'sites': pos[index], 'aa': aa[index]}

#changes a few entries in the sgd file that have analog names
sgd['GSG1'] = sgd['TRS85']
sgd['KRE11'] = sgd['TRS65']
sgd['KEM1'] = sgd['XRN1']
sgd['DUR12'] = sgd['DUR1,2']
sgd['ARG56'] = sgd['ARG5,6']



def mj_findMotif(sites, sequence, motiflist, limlist, names):
    
    '''Searches an amino acid motif in a translated protein sequence,
    starting from a known amino acid position of a known protein

    :param sites:       iterable ints, list of sites, e.g. [264, 268]
    :param sequence:    str of single-letter-aa-code, e.g. 'MDKLSLAKDPE...'
    :param motiflist:   list of regex-able str describing the motif, e.g. '[ST]P'
    :param limlist:     list of tuples consisting of exactly 2 ints, 
                        motif boundaries relative to site (site pos = -1), 
                        e.g. xSp: (-1, 0), SxxxL: (-1, +3)
    :param names:       list of strings describing the motif (multiple nominations possible)
                        e.g. 'stp', 'stp', 'baso'

    :returns:           list of bools with len(bool-list) = len(site-list)
                        --> can be used with any()
    '''
    
    m = len(motiflist)
    l = len(limlist)
    s = len(sites)
    if m != l:
        print 'Match Error: #motif != #lim'
    else:
        pass
    allMotifs = list()
    rePatterns = [re.compile(motif, re.VERBOSE) for motif in motiflist]
    allSeqs = list()
    for site in sites:
        siteSeqs = [sequence[site + entry[0]:site + entry[1]+1] for entry in limlist]
        allSeqs.append(siteSeqs)
    allMotifs = list()
    for entry in allSeqs:
        whichMotifs = list()
        for number in range(m):
            if rePatterns[number].search(entry[number]) != None:
                mo = names[number]
            else:
                continue
            if mo in whichMotifs:
                continue
            else: 
                whichMotifs.append(mo)
        allMotifs.append(whichMotifs)
    
    return allMotifs


#use mj_findMotif for STP and RXXS
def mj_findSTP(sites, sequence):
    return mj_findMotif(sites, sequence, ['[ST]P'], [(-1, 0)], ['stp'])
#   exampledata for stp:
#   motif = '[ST]P'
#   sequence = sgd['ACA1']['seq']
#   sites = [9]
#   lim = [-1, +1]

def mj_findCAMC(sites, sequence):
    return mj_findMotif(sites, sequence, ['[RK]\D{2}[ST]'], [(-4, -1)], ['camc'])
#   exampledata for camc:
#   motif = 'R\D{2}S'
#   sequence = sgd['ACC1']['seq']
#   sites = [1157]
#   lim = [-3, 0]

#looks for STP or RXXS motifs
#notFound --> in case there are prot missing from the sgd-dict
#appends information about motif to sitedict
notFound = list()

PB = ['[ST]P', '(SS|TT)P', '[RK]\D{2}[ST]', '[ST]\D{3}L']
limPB = [(-1, 0), (-1, 1), (-4, -1), (-1, 3)]
namPB = ['stp', 'stp?','baso', 'baso']
for key, value in sitedict.items():
    try:
        isSTP_isBaso = mj_findMotif(value['sites'], sgd[value['prot']]['seq'], PB, limPB, namPB)
    except KeyError:
        if value['prot'] not in notFound:
            notFound.append(value['prot'])
        else:
            pass
    if any(isSTP_isBaso) == False:
        value['STP_Baso'] = npy.nan
    elif len(isSTP_isBaso) == 1:
        value['STP_Baso'] = isSTP_isBaso[0]
    else:
        value['STP_Baso'] = [entry for entry in isSTP_isBaso]

com['stp_baso'] = [sitedict[x]['STP_Baso'] for x in com.prot_site]

head = mj_givehead(com)
head


#TAKES 3 MINUTES TO EXACUTE!
#allprots = [value['prot'] for key, value in sitedict.items()]
#allprots.sort()
#prots = [x[0] for x in groupby(allprots)]
#
#protdict = dict()
#for entry in prots:
#    l = list()
#    for key, value in sitedict.items():
#        if entry in str(key):
#            for sites in value['sites']:
#                l.append(sites)
#                l.sort()
#    s = [x[0] for x in groupby(l)]
#    protdict[entry] = s


#testdict['FAR8, 132(T)'] = sitedict['FAR8, 132(T)']
#
#testdict = dict()
#
#testdict['FAR8, 132(T)'] = sitedict['FAR8, 132(T)']
#testdict['PMD1, 1641(S)'] = sitedict['PMD1, 1641(S)']
#testdict['PMD1, 1642(S)'] = sitedict['PMD1, 1642(S)']
#testdict['PMD1, 1355(S) / 1356(S)'] = sitedict['PMD1, 1355(S) / 1356(S)']
#testdict['PMD1, 1355(S) / 1358(T)'] = sitedict['PMD1, 1355(S) / 1358(T)']

#samefieldoverlap = 60%
F1 = list()
for x in com.index:
    if 'Field 1' in str(com.FieldMQ[x]):
        F1.append('Field 1')
    else: 
        F1.append(com.FieldMQ[x])

com['F1'] = F1

graph = dict()

cols = ['quantified P-Sites', 
        'quantified Proteins', 
        'P-Sites assigned to Fields',
        'S/TP Sites assigned to Fields', 
        'Basophilic Kinase Sites assigned to Fields']

#quantified sites --> entries only identified in rck2/kss1 shotguns are dismissed from PD
#just for counting
dismissed = com[[com.ProteinPD[x] is nan and com.ProteinMQ[x] is nan for x in com.index]]


#overlaps excluding rck2/kss1 only data
both = com[com.ProteinPD == com.ProteinMQ]
mixed = com[~(com.ProteinMQ == com.ProteinPD)]
pd = mixed[[x is not nan for x in mixed.ProteinPD]]
mq = mixed[[x is not nan for x in mixed.ProteinMQ]]

B = len(both)
PD = len(pd)
MQ = len(mq)
quant = {'group': cols[0], 'both': B, 'PD': PD, 'tot': B+PD+MQ, 'MQ': MQ}
graph[cols[0]] = quant

#proteins
Pr_both = [x[0] for x in groupby(both.ProteinPD)]
Pr_pd = [x[0] for x in groupby(pd.ProteinPD)]
Pr_mq = [x[0] for x in groupby(mq.ProteinMQ)]

Pr_pdonly= [pr for pr in Pr_pd if not any([x == pr for x in Pr_both])]
Pr_mqonly = [pr for pr in Pr_mq if not any([x == pr for x in Pr_both])]

B = len(Pr_both)
PD = len(Pr_pdonly)
MQ = len(Pr_mqonly)
Pr_quant = {'group': cols[1], 'both': B, 'PD': PD, 'tot': B+PD+MQ, 'MQ': MQ}
graph[cols[1]] = Pr_quant

#sites assigned to fields
F = 'Field' 
m = 'middle'
allF_MQ = [any([F in str(x), m in str(x)]) for x in com.FieldMQ]
allF_PD = [any([F in str(x), m in str(x)]) for x in com.FieldPD]
bothass = [all([allF_MQ[x], allF_PD[x]]) for x in com.index]
F_both = com[bothass]
mqass = [all([allF_MQ[x], allF_PD[x] == False]) for x in com.index]
pdass = [all([allF_PD[x], allF_MQ[x] == False]) for x in com.index]
F_pdonly = com[pdass]
F_mqonly = com[mqass]

B = len(F_both)
PD = len(F_pdonly)
MQ = len(F_mqonly)
assigned = {'group': cols[2], 'both': B, 'PD': PD, 'tot': B+PD+MQ, 'MQ': MQ}
graph[cols[2]] = assigned

#stp sites
stp_both = both[['stp' in str(x) for x in both.stp_baso]]
stp_pd = pd[['stp' in str(x) for x in pd.stp_baso]]
stp_mq = mq[['stp' in str(x) for x in mq.stp_baso]]

#baso sites
baso_both = both[['baso' in str(x) for x in both.stp_baso]]
baso_pd = pd[['baso' in str(x) for x in pd.stp_baso]]
baso_mq = mq[['baso' in str(x) for x in mq.stp_baso]]

#stp sites assigned to fields
stp_F_both = F_both[['stp' in str(x) for x in F_both.stp_baso]]
stp_F_pd = F_pdonly[['stp' in str(x) for x in F_pdonly.stp_baso]]
stp_F_mq = F_mqonly[['stp' in str(x) for x in F_mqonly.stp_baso]]

B = len(stp_F_both)
PD = len(stp_F_pd)
MQ = len(stp_F_mq)
assigned_stp = {'group': cols[3], 'both': B, 'PD': PD, 'tot': B+PD+MQ, 'MQ': MQ}
graph[cols[3]] = assigned_stp

#baso sites assigned to fields
baso_F_both = F_both[['baso' in str(x) for x in F_both.stp_baso]]
baso_F_pd = F_pdonly[['baso' in str(x) for x in F_pdonly.stp_baso]]
baso_F_mq = F_mqonly[['baso' in str(x) for x in F_mqonly.stp_baso]]

B = len(baso_F_both)
PD = len(baso_F_pd)
MQ = len(baso_F_mq)
assigned_baso = {'group': cols[4], 'both': B, 'PD': PD, 'tot': B+PD+MQ, 'MQ': MQ}
graph[cols[4]] = assigned_baso

perc = dict()

def MJ_normalize (list_of_fractions, total_number):

    perced = list()
    for entry in list_of_fractions:
        equa = round(100 / (float(total_number) / float(entry)),2)
        perced.append(equa)
    return perced

nums = ['MQ', 'PD', 'both']

perc = dict()
for group in graph.keys():
    fracs = [graph[group][entry] for entry in nums]
    tot = graph[group]['tot']
    perced = MJ_normalize(fracs, tot)
    perc[group] = {nums[0]: perced[0], nums[1]: perced[1],nums[2]: perced[2],
                   'group': group, 'tot': graph[group]['tot']}


stacks = dict()
for colName in ['group', 'both', 'PD', 'MQ', 'tot']:
    stacks[colName] = [graph[name][colName] for name in cols]

stacksperc = dict()
for colName in ['group', 'both', 'PD', 'MQ', 'tot']:
    stacksperc[colName] = [perc[name][colName] for name in cols]

df1 = pds.DataFrame(data = stacks)
df2 = pds.DataFrame(data = stacksperc)

gData = df1[['both', 'PD', 'MQ']]

#gData.plot(kind = 'barh', stacked=True, 
#           x = df1.group, figsize=(10,5), 
#           color = ['#f6ea9c','#e3d1de','#9586b2'])
#plt.title('Comparison of sites/proteins found in MQ and PD')
#plt.xlabel('number of quantified sites/proteins')
#plt.ylabel('')
#plt.gca().invert_yaxis()
#plt.legend(loc='lower right')
#plt.tight_layout()

percData = df2[['both', 'PD', 'MQ']]


from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')

from matplotlib import gridspec as gs

plt.figure(figsize = (10,10))
grid = gs.GridSpec(2, 1, left=0.40, right=0.88, top=0.92, bottom=0.08, hspace=0.25) 
plt.suptitle('Comparison of sites/proteins quantified in MQ and PD', fontsize=16)
plt.subplot(grid[1,0])
percData.plot(kind = 'barh', stacked=True, 
              x = df2.group, ax = plt.gca(), 
              color = ['#f6ea9c','#e3d1de','#9586b2'])
plt.xlabel('normalised number of quantified sites/proteins [%]')
plt.ylabel('')
plt.xlim(-5,105)
plt.gca().invert_yaxis()
plt.box('off')
plt.legend(loc=3, bbox_to_anchor=(1, 0), prop = fontP, frameon = False)

plt.subplot(grid[0,0])
gData.plot(kind = 'barh', stacked=True, 
           x = df1.group, ax = plt.gca(), 
           color = ['#f6ea9c','#e3d1de','#9586b2'])
plt.xlabel('number of quantified sites/proteins')
plt.ylabel('')
plt.xlim(xmax = (max(df1.tot)*1.05))
plt.gca().invert_yaxis()
plt.tick_params(axis= 'x', direction = 'out')
plt.box('off')
plt.legend(loc=3, bbox_to_anchor=(1, 0), prop = fontP, frameon = False)


#field 1 sites
same = both[both.FieldPD == both.F1]
#3569 (includes not assigned)
F_same = F_both[F_both.FieldPD == F_both.F1]
#2587

allF1 = com[[com.FieldPD[x] == 'Field 1' or com.F1[x] == 'Field 1' for x in com.index]]
# 376
allF1_bothquant = both[[both.FieldPD[x] == 'Field 1' or both.F1[x] == 'Field 1' for x in both.index]]
#260
allF1_bothfields = F_both[[F_both.FieldPD[x] == 'Field 1' or F_both.F1[x] == 'Field 1' for x in F_both.index]]
#147

F1_same = F_same[F_same.F1 == 'Field 1']
#101
PrF1_same = [x[0] for x in groupby(F1_same.ProteinPD)]
#82


#pd
F1_pdtotal = com[[com.FieldPD[x] == 'Field 1' for x in com.index]]
len(F1_pdtotal) #204
PrF1_pdtotal = [x[0] for x in groupby(F1_pdtotal.ProteinPD)]
len(PrF1_pdtotal) #151

F1_pdonly = com[[com.FieldPD[x] == 'Field 1' and com.F1[x] != 'Field 1' for x in com.index]]
len(F1_pdonly) #103
PrF1_pdonly = [x[0] for x in groupby(F1_pdonly.ProteinPD)]
len(PrF1_pdonly) #88
PrF1_pd_ONLY = [pr for pr in PrF1_pdonly if not any([x == pr for x in PrF1_same])]
len(PrF1_pd_ONLY) #69


#mq
F1_mqtotal = com[[com.F1[x] == 'Field 1' for x in com.index]]
len(F1_mqtotal) #273
PrF1_mqtotal = [x[0] for x in groupby(F1_mqtotal.ProteinMQ)]
len(PrF1_mqtotal) #194

F1_mqonly = com[[com.F1[x] == 'Field 1' and com.FieldPD[x] != 'Field 1' for x in com.index]]
len(F1_mqonly) #172
PrF1_mqonly = [x[0] for x in groupby(F1_mqonly.ProteinMQ)]
len(PrF1_mqonly) #135
PrF1_mq_ONLY = [pr for pr in PrF1_mqonly if not any([x == pr for x in PrF1_same])]
len(PrF1_mq_ONLY) #112

#pd&mq
PrF1_pd_notmq = [pr for pr in PrF1_pd_ONLY if not any([x == pr for x in PrF1_mq_ONLY])]
len(PrF1_pd_notmq) #56
PrF1_mq_notpd = [pr for pr in PrF1_mq_ONLY if not any([x == pr for x in PrF1_pd_ONLY])]
len(PrF1_mq_notpd) #99

com.to_csv(r'D:\aaaaaWork\Data\PD_MQ_compare\ois.txt', header=True, index=False, sep='\t')

F1stp = com[[com.F1[x] == 'Field 1' and 'stp' in str(com.stp_baso[x]) for x in com.index]]
F1stp_mq = F1stp[[F1stp.ProteinMQ[x] != 'PKH1' and F1stp.ProteinMQ[x] != 'PMD1' for x in F1stp.index]]
F1stp_mqonly = F1stp_mq[F1stp_mq.FieldPD != 'Field 1']

PrF1stp_mq = [x[0] for x in groupby(F1stp_mq.ProteinMQ)]
PrF1stp_mqonly = [x[0] for x in groupby(F1stp_mqonly.ProteinMQ)]

mots = list()
for entry in F1stp_mqonly.prot_site:
    site = sitedict[entry]
    seq = sgd[site['prot']]['seq']
    pos = site['sites'][0]
    motif = seq[pos-5:pos+5]
    mots.append((motif, entry))

oldF1 = com[[com.FieldPD[x] == 'Field 1' and 'stp' in str(com.stp_baso[x]) for x in com.index]]
oldPrF1 = [x[0] for x in groupby(oldF1.ProteinPD)]

PrF1stp_mqnew = [pr for pr in PrF1stp_mqonly if not any([x == pr for x in oldPrF1])]
new = PrF1stp_mqnew

F1stpboth = [pr for pr in PrF1stp_mqonly if not any([x == pr for x in PrF1stp_mqnew])]
