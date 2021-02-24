#!/usr/bin/env python
# coding: utf-8

# In[1]:


# modules required for handling dataframes
import os
import pandas as pd


# In[2]:


basedir = '/home/yiheng/MinION_data' # the directory where all the documents of each sequencing run are stored.
k2_output_db = pd.read_csv(os.path.join(basedir, 'genome_reference_mock.standardkraken2_output'), sep='\t', header=None)


# In[3]:


# you can name whatever you like, or just use numbers.
headers = ['classification', 'accession', 'taxid', 'length', 'kmer_profile']
k2_output_db.columns = headers


# In[8]:


k2_output_db_classified = k2_output_db[k2_output_db.classification == 'C']


# In[9]:


k2_output_db_classified = k2_output_db_classified.reset_index(drop=True)


# In[10]:


def generate_final_bed(input_df, row_index):
    a = input_df.iloc[row_index,-1]
    b = a.split(' ')
# important filter function
    b = list(filter(None, b))
    end = [int(c.split(':')[-1]) for c in b]
    
    taxid_list = [c.split(':')[0] for c in b]
    
    for x in range(0,len(end)):  
        if x == 0:
            end[x] = end[x]
        elif x == len(end) - 1:
            end[x] = int(input_df.iloc[row_index,3]) - 33

        else:
            end[x] = end[x] + end[x-1]



    start = end.copy()
    for y in range(0,len(start)):
        if y == 0:
            start[y] = 0
        else:
            start[y] = end[y-1]

    bed_df = pd.DataFrame({'taxid':taxid_list, 'start':start, 'end':end})
    bed_df_trim = bed_df[(bed_df.taxid != '0') | ((bed_df.taxid == '0') & (bed_df.end - bed_df.start < 100))]
    bed_df_trim = bed_df_trim.reset_index(drop=True)
    bed_df_drop_taxid = bed_df_trim.drop(columns='taxid')

    bed_df_drop_taxid['accession'] = input_df.iloc[row_index,1]
    return bed_df_drop_taxid
    


# In[11]:


# These following lines need to be sticked together otherwise the bed_df_final will keep accumulating
columns = ['start','end', 'accession']
bed_df_final = pd.DataFrame(columns=columns)
for index in k2_output_db_classified.index:
    
    bed_df_final = bed_df_final.append(generate_final_bed(k2_output_db_classified, index), ignore_index = True)


# In[12]:


bed_df_final = bed_df_final[['accession', 'start', 'end']]


# In[13]:


bed_df_final.to_csv(os.path.join(basedir, 'genome_reference_mock_final.bed'), sep='\t', header=False, index=False)


# In[ ]:




