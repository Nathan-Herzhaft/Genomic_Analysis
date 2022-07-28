#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import Bio




#%%
root = 'data/genes-'
augustus = pd.read_csv(root + 'augustus.csv')
genscan = pd.read_csv(root + 'genscan.csv')
refseq = pd.read_csv(root + 'refseq.csv')
xeno_refseq = pd.read_csv(root + 'xeno-refseq.csv')
ensembl = pd.read_csv(root + 'ensembl.csv')




# %%
augustus.head()




# %%
genscan.head()



# %%

methods = {'Augustus' : augustus, 'Genscan' : genscan, 'RefSeq' : refseq, 'Xeno-RefSeq' : xeno_refseq, 'Ensembl' : ensembl}
for method, dataset in methods.items() :
    print(method + ' number of predicted genes : ' +str(dataset.shape[0]) + '\n')
    
    
    
    
# %%
refseq = pd.read_csv(root + 'refseq.csv')
refseq = refseq.groupby('chrom')['chrom'].count()
refseq= refseq.loc[refseq.values > 10]

ensembl = pd.read_csv(root + 'ensembl.csv')
ensembl = ensembl.groupby('chrom')['chrom'].count()
ensembl = ensembl.loc[ensembl.values > 10]

counts = pd.DataFrame.from_dict(data={'RefSeq' : refseq, 'Ensembl' : ensembl})
counts
# %%
