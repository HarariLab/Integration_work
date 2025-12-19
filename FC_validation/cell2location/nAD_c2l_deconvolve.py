import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

import scanpy as sc
import cell2location
from cell2location.models import RegressionModel
import pandas as pd
import numpy as np
import squidpy as sq
import jax
import torch
import pickle
import gc

print(jax.devices(), flush=True)
print(f"PyTorch: {torch.cuda.is_available()}")

ref_adata = sc.read_h5ad("cell2location_ann/sea_10k_ref_adata.h5ad")
ref_model = cell2location.models.RegressionModel.load("cell2location_ann/model_sea_10k", ref_adata)

ref_adata = ref_model.export_posterior(ref_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})

# Export estimated gene expression levels
if 'means_per_cluster_mu_fg' in ref_adata.varm.keys():
    inf_aver = ref_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in ref_adata.uns['mod']['factor_names']]].copy()
else:
    inf_aver = ref_adata.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in ref_adata.uns['mod']['factor_names']]].copy()
inf_aver.columns = ref_adata.uns['mod']['factor_names']

####
base_path = "data/ST_AN1792_GSE263034"
sample_ids = [
    "GSM8184312",
    "GSM8184313",
    "GSM8184314",
    "GSM8184315",
    "GSM8184316",
    "GSM8184317"
]

print(f"Processing {len(sample_ids)} samples\n", flush=True)

for sample_id in sample_ids:

    print(f"Importing data for sample: {sample_id}", flush=True)
    adata_vis = sq.read.visium(f"{base_path}/{sample_id}", library_id=sample_id)
    adata_vis.var_names_make_unique()
    n_initial = len(adata_vis)
    
    # QC with dynamic thresholds
    sc.pp.calculate_qc_metrics(adata_vis, inplace=True)
    qc_counts = adata_vis.obs['total_counts'].quantile(0.05)
    qc_genes = adata_vis.obs['n_genes_by_counts'].quantile(0.05)
    
    sc.pp.filter_cells(adata_vis, min_counts=qc_counts)
    sc.pp.filter_cells(adata_vis, min_genes=qc_genes)
    sc.pp.filter_genes(adata_vis, min_cells=5)
    n_qc = len(adata_vis)
    
    # Normalize
    adata_vis.obs['sample'] = sample_id

    # Remove MT genes from spatial data
    adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]
    
    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis.var_names_make_unique()
    adata_vis = adata_vis[:, list(intersect)].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    
    # prepare anndata for cell2location model
    print("Prepare anndata for cell2location model", flush=True)
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

    torch.cuda.empty_cache()

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=30,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )
    print(f"Train model to estimate cell abundance, sample {sample_id}", flush=True)
    mod.train(max_epochs=10000,
              # train using full data (batch_size=None)
              batch_size=None,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
              train_size=1
    )

    torch.cuda.empty_cache()

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )
    
    # Save model
    print(f"Save models for sample {sample_id}")
    mod.save(f"c2l_res/nAD/{sample_id}", overwrite=True)
    
    # Save anndata object with results
    print(f"Save anndata for sample {sample_id}")
    adata_vis.write(f"c2l_res/nAD/{sample_id}.h5ad")

    del mod
    del adata_vis
    gc.collect()
