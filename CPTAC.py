
import cptac
import cptac.utils as ut
import os

os.getcwd()
os.chdir("/Users/naturexv/Datasets/CPTAC/")
os.getcwd()

import fireducks.pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cptac.get_cancer_info()
cptac.get_source_options()

### 读取数据 ###
brca = cptac.Brca()
brca.list_data_sources()
brca_pro = brca.get_dataframe('proteomics', 'umich')
brca_clinical = brca.get_clinical('mssm')   # brca_clinical.type_of_analyzed_samples = NaN

ccrcc = cptac.Ccrcc()
ccrcc.list_data_sources()
ccrcc_pro = ccrcc.get_dataframe('proteomics', 'umich')
ccrcc_clinical = ccrcc.get_clinical('mssm')
ccrcc_clinical.type_of_analyzed_samples

coad = cptac.Coad()
coad.list_data_sources()
coad_pro = coad.get_dataframe('proteomics', 'umich')
coad_clinical = coad.get_clinical('mssm')
coad_clinical.type_of_analyzed_samples   # coad_clinical.type_of_analyzed_samples = NaN


gbm = cptac.Gbm()
gbm.list_data_sources()
gbm_pro = gbm.get_dataframe('proteomics', 'umich')
gbm_clinical = gbm.get_clinical('mssm')
gbm_clinical.type_of_analyzed_samples    # gbm_clinical.type_of_analyzed_samples = tumor

hnscc = cptac.Hnscc()
hnscc.list_data_sources()
hnscc_pro = hnscc.get_dataframe('proteomics', 'umich')
hnscc_clinical = hnscc.get_clinical('mssm')
hnscc_clinical.type_of_analyzed_samples

lscc = cptac.Lscc()
lscc.list_data_sources()
lscc_pro = lscc.get_dataframe('proteomics', 'umich')
lscc_clinical = lscc.get_clinical('mssm')
lscc_clinical.type_of_analyzed_samples

luad = cptac.Luad()
luad.list_data_sources()

ov = cptac.Ov()
ov.list_data_sources()

pdac = cptac.Pdac()
pdac.list_data_sources()

ucec = cptac.Ucec()
ucec.list_data_sources()


