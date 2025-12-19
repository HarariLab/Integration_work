# Jai Babe Di
# Jai Guru Maa Ji

# Jai Babe Di
# Jai Guru Maa Ji


import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import pandas as pd
pd.set_option('display.max_rows', None)

os.chdir('/home/garga7/Anjali_data/int_data/Dolan_data/')

sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(8, 8))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

condition_key = 'batch'
cell_type_key = 'Integrated_cell_type'
target_conditions = ['Dolan']


adata_all  = sc.read('/home/garga7/Anjali_data/int_data/Dolan_data/int_data_both_upgene.h5ad')
adata = adata_all.raw.to_adata()

source_adata = adata[~adata.obs[condition_key].isin(target_conditions)].copy()
target_adata = adata[adata.obs[condition_key].isin(target_conditions)].copy()

sca.models.SCVI.setup_anndata(source_adata, labels_key=cell_type_key)
vae = sca.models.SCVI(source_adata,n_layers=2,encode_covariates=True,deeply_inject_covariates=False,use_layer_norm="both",use_batch_norm="none",)
vae.train()

scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")
print("Labelled Indices: ", len(scanvae._labeled_indices))
print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))
scanvae.train(max_epochs=20)

#Create anndata file of latent representation and compute UMAP

reference_latent = sc.AnnData(scanvae.get_latent_representation())
reference_latent.obs["Integrated_cell_type"] = source_adata.obs[cell_type_key].tolist()
reference_latent.obs["batch"] = source_adata.obs[condition_key].tolist()
sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
sc.pl.umap(reference_latent,color=['Integrated_cell_type', 'batch'],frameon=False,wspace=0.6,)
plt.savefig('umap_reference.png')

reference_latent.obs['predictions'] = scanvae.predict()
print("Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.Integrated_cell_type)))

ref_path = 'ref_model/'
scanvae.save(ref_path, overwrite=True)

#Perform surgery on reference model and train on query dataset without cell type labels

#################
model = sca.models.SCANVI.load_query_data(target_adata,ref_path,freeze_dropout = True,)
model._unlabeled_indices = np.arange(target_adata.n_obs)
model._labeled_indices = []
print("Labelled Indices: ", len(model._labeled_indices))
print("Unlabelled Indices: ", len(model._unlabeled_indices))


model.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['Integrated_cell_type'] = target_adata.obs[cell_type_key].tolist()
query_latent.obs['batch'] = target_adata.obs[condition_key].tolist()
  
sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
plt.figure()
sc.pl.umap(query_latent,color=["batch","Integrated_cell_type"], frameon=False,wspace=0.6,)
plt.savefig('umap_query.png')
surgery_path = 'surgery_model'
model.save(surgery_path, overwrite=True)

#Compute Accuracy of model classifier for query dataset and compare predicted and observed cell types

query_latent.obs['predictions'] = model.predict()
print("Acc: {}".format(np.mean(query_latent.obs.predictions == query_latent.obs.Integrated_cell_type)))

f= open('query_latent.obs-all.txt', 'a') 
print(query_latent.obs, file=f)

df = query_latent.obs.groupby(["Integrated_cell_type", "predictions"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

plt.figure(figsize=(8, 8))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed")
plt.savefig('heatmap_all.png')
