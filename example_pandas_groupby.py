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
# lines = len(evi)
# evi[normrats] = np.empty(lines)

#grouping all lines by their pcount-value (here: 0, 1, 2)
# for c, subset in evi.groupby(pcount):
#     M = subset[logratios].median()	#calculate median of subset
#     PF = M0 - M #calaculate p-factor / shift of the bellcurve in reference to p0-subset
#     evi[normrats].loc[subset.index] = (evi[logratios].loc[subset.index] + PF) / M0 # works nicely w/o loop
'''
This doesn't work, because the groupby command actually does not generate a separate subset *object*
Manipulating the original dataframe in this conext is not really possible, even with the lines above.
Groupby does allow iteration, but not manipulation of the original data.

I have previously gotten around that in a rather clunky manner by doing something like the following
'''
# fullNormRats = None # Set this to the special python object `None`, so we can check against it later on
# for c, subset in evi.groupby(pcount):
#     M = subset[logratios].median()
#     PF = M0 - M
#     normratSeries = (subset[lograts] + PF)/M0 # Calculate the normRatio series
#     if fullNormRats is None: # Essentially this is only for the first run
#         fullNormRats = normratSeries
#     else: # For every other run, just concatenate the new series with the old
#         fullNormRats = pd.concat([fullNormRats, normratSeries])
'''
This results in a pandas Series, which should have maintained indexes, in reference to the input dataframe,
but I haven't tried it with this code, as there is a much more elegant and 'panda-onic' way.
'''

#loop version
#for c, subset in evi.groupby(pcount):
#    M = subset[logratios].median()	#calculate median of subset
#    PF = M0 - M #calaculate p-factor / shift of the bellcurve in reference to p0-subset
#    for entry in subset.index: #works without loop as well
#        evi[normrats].loc[entry] = (evi[logratios].loc[entry] + PF) / M0

# This function does exactly what the for loop would do on the subset of the data
# The first argument of the function will always be the subset of the data.
# Any other variables have to be passed through apply
def normRatio(df, newcol, normcol, normVal):
    M = df[normcol].median()
    PF = normVal - M
    df[newcol] = (df[normcol] + PF)/M0
    return df

'''
The pandas function apply can be used to send the current dataframe (or subset)
to a function.
One example is to use it to calculate sums or averages and output those, but
another application is to generate a new column.
Generally, you only do `apply(func)`, but you can also pass keyword arguments,
such as here.
See here for some ideas:
http://pandas.pydata.org/pandas-docs/stable/groupby.html
'''
evi = evi.groupby(pcount).apply(normRatio, newcol=normrats, normcol=logratios, normVal=M0)

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
