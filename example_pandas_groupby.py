import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Path = r'D:\aaaaaWork\Python\exampledata_MJ.txt'

rc_evi = 'Leading Proteins'


#function for deleting lines containing rev and cont
def delete_revcon ( yourpath, yourcolumn ):
    data = pd.read_table(yourpath)
    rev = data[~data[yourcolumn].str.contains('REV_')]
    revcon = rev[~rev[yourcolumn].str.contains('[cont]')]
    return revcon 


# use function (usually done on various files at once)
evi = delete_revcon(Path, rc_evi)


#setting values for the needed columns
ratio = 'Ratio H/L'
pcount = 'Phospho (STY)'
pfacs = 'Prolinefactors'
logratios = 'Log2_Ratios'
normrats = 'Pcorr_Norm_Ratios'


#append new column containing log2-transformed ratios
evi[logratios] = np.log2(evi[ratio])


#calculation the mean of ratios with pcount = 0
M0 = evi[logratios][evi[pcount] == 0].median()

#opening new df for corrected ratios
lines = len(evi)
evi[normrats] = np.empty(lines)

#grouping all lines by their pcount-value (here: 0, 1, 2)
for c, subset in evi.groupby(pcount): 
    M = subset[logratios].median()	#calculate median of subset
    PF = M0 - M #calaculate p-factor / shift of the bellcurve in reference to p0-subset
    evi[normrats].loc[subset.index] = (evi[logratios].loc[subset.index] + PF) / M0 # works nicely w/o loop

#loop version
#for c, subset in evi.groupby(pcount): 
#    M = subset[logratios].median()	#calculate median of subset
#    PF = M0 - M #calaculate p-factor / shift of the bellcurve in reference to p0-subset
#    for entry in subset.index: #works without loop as well
#        evi[normrats].loc[entry] = (evi[logratios].loc[entry] + PF) / M0


#to check if the groupby-command really does what I want

evi['test1'] = np.empty(lines)
for c, subset in evi.groupby(pcount): 
    for entry in subset.index:
        evi['test1'].loc[entry] = (str(evi[pcount].loc[entry]) + 'P in pcount')

evi['test2'] = np.empty(lines)
for c, subset in evi.groupby(pcount): 
    evi['test2'].loc[subset.index] = str(evi[pcount].loc[subset.index]) + ' P in pcount' 

evi['test3'] = np.empty(lines)
for c, subset in evi.groupby(pcount):
    evi['test3'].loc[subset.index] = evi[pcount].loc[subset.index] 

print evi['test1'] #works fine
print evi['test2'] #does sth weird...
print evi['test3'] #works again. not fine with str?