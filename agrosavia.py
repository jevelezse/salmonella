#!/usr/bin/env python
# coding: utf-8

# 2 Data parsing

# In[1]:


import pandas as pd


# In[2]:


ariba = pd.read_csv(r'~/Documents/agrosavia/coding_assessment_exercise-master/exercise_2/ariba_amr_output.csv')
ncbi = pd.read_csv(r'~/Documents/agrosavia/coding_assessment_exercise-master/exercise_2/ncbi_acquired_genes_metadata.csv')


# In[3]:


ariba.head(3)


# In[4]:


ncbi.head()


# In[5]:


resistance = ['CEPHALOSPORIN','CARBAPENEM']


# In[6]:


ncbi_filtered = ncbi[ncbi['subclass'].str.contains('|'.join(resistance))]


# In[7]:


#gene_family = list(ncbi_filtered['gene_family'].str.replace('-','_'))


# In[8]:


gene_family = ncbi_filtered['gene_family'].tolist()


# In[9]:


refseq = ncbi_filtered['refseq_nucleotide_accession'].tolist()


# In[10]:


gene_family[0:10]


# In[11]:


gene_filtred = list(set(gene_family))
gene_filtred[0:10]


# In[12]:


x = list(ariba)


# In[13]:


list_columns = [i for i in x if any(i for j in gene_filtred if str(j) in i)]


# In[14]:


list_columns.append('name')


# In[15]:


ariba1 = ariba[list_columns]


# In[16]:


ariba2 = ariba1.filter(regex=r'(assembled|ctg_cov|name|ref_seq)')


# In[17]:


ariba.head(2)


# In[18]:


ariba2 = ariba2.replace('yes_nonunique','yes')


# In[19]:


x1 = list(ariba2)


# In[20]:


list(ariba2)


# In[21]:


blaAFM_1 = x1[0:3]
blaAFM_1.append('name')


# In[22]:


blaAFM_1 = ariba2[blaAFM_1]
blaAFM_1.head(2)


# In[23]:


blaAFM_1 = blaAFM_1[(blaAFM_1['blaAFM_1.ctg_cov'] > 10.0) & (blaAFM_1['blaAFM_1.assembled'] == 'yes')]
blaAFM_1


# In[24]:


sample_blaAFM_1 = blaAFM_1['name'].tolist()


# In[25]:


blaEC = x1[3:6]
blaEC.append('name')


# In[26]:


blaEC = ariba2[blaEC]
blaEC.head(2)


# In[27]:


blaEC = blaEC[(blaEC['blaEC-.ctg_cov'] > 10.0) & (blaEC['blaEC-.assembled'] == 'yes')]
blaEC.head(2)


# In[28]:


sample_blaEC = blaEC['name'].tolist()


# In[32]:


blaNDM = x1[6:9]
blaNDM.append('name')


# In[33]:


blaNDM = ariba2[blaNDM]
blaNDM.head(2)


# In[34]:


blaNDM = blaNDM[(blaNDM['blaNDM_-.ctg_cov'] > 10.0) & (blaNDM['blaNDM_-.assembled'] == 'yes')]
blaNDM.head(2)


# In[35]:


sample_blaNDM = blaNDM['name'].tolist()


# In[36]:


blaOXA_10 = x1[9:12]
blaOXA_10.append('name')


# In[37]:


blaOXA_10 = ariba2[blaOXA_10]
blaOXA_10.head(2)


# In[38]:


blaOXA_10  = blaOXA_10[(blaOXA_10 ['blaOXA_-_10.ctg_cov'] > 10.0) & (blaOXA_10['blaOXA_-_10.assembled'] == 'yes')]
blaOXA_10.head(2)


# In[39]:


sample_blaOXA_10 = blaOXA_10['name'].tolist()


# In[40]:


blaOXA_31 = x1[12:15]
blaOXA_31.append('name')


# In[41]:


blaOXA_31 = ariba2[blaOXA_31]
blaOXA_31.head(2)


# In[42]:


blaOXA_31 = blaOXA_31[(blaOXA_31['blaOXA_-_31.ctg_cov'] > 10.0) & (blaOXA_31['blaOXA_-_31.assembled'] == 'yes')]
blaOXA_31.head(2)


# In[43]:


sample_blaOXA_31 = blaOXA_31['name'].tolist()


# In[47]:


blaOXA_35 = x1[15:18]
blaOXA_35.append('name')


# In[48]:


blaOXA_35 = ariba2[blaOXA_35]
blaOXA_35.head(2)


# In[49]:


blaOXA_35 = blaOXA_35[(blaOXA_35['blaOXA_-_35.ctg_cov'] > 10.0) & (blaOXA_35['blaOXA_-_35.assembled'] == 'yes')]
blaOXA_35.head(2)


# In[50]:


sample_blaOXA_35 = blaOXA_35['name'].tolist()


# In[51]:


blaTEM = x1[18:-1]
blaTEM.append('name')


# In[52]:


blaTEM = ariba2[blaTEM]
blaTEM.head(2)


# In[53]:


blaTEM= blaTEM[(blaTEM['blaTEM_-.ctg_cov'] > 10.0) & (blaTEM['blaTEM_-.assembled'] == 'yes')]
blaTEM.head(2)


# In[54]:


sample_blaTEM = blaTEM['name'].tolist()


# In[55]:


sample_blaOXA = sample_blaOXA_10 + sample_blaOXA_31 + sample_blaOXA_35


# In[56]:


total_samples = sample_blaAFM_1 + sample_blaEC + sample_blaNDM + sample_blaOXA + sample_blaTEM


# In[57]:


len(total_samples)


# In[58]:


samples_drop = list(set(total_samples))


# In[59]:


len(samples_drop)


# In[63]:


print('Total samples with cephalosporin or carbapenem resistance:', len(samples_drop))


# In[60]:


samples = pd.DataFrame({'blaAFM_1': pd.Series(sample_blaAFM_1,dtype=str), 'blaEC': pd.Series(sample_blaEC,dtype=str), 'blaNDM': pd.Series(sample_blaNDM,dtype=str), 'blaOXA': pd.Series(sample_blaOXA,dtype=str),
                        'blaTEM': pd.Series(sample_blaTEM,dtype=str)})


# In[61]:


samples.to_csv('samples.csv')

