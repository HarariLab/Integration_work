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

torch.set_float32_matmul_precision('high')

ref_adata = sc.read_h5ad("cell2location_ann/ref_adata_wo_cellstates.h5ad")

# Create regression model
print(f"Setting up regression model", flush=True)
RegressionModel.setup_anndata(ref_adata, labels_key='celltype', batch_key="external_donor_name_label")

print(f"Creating regression model", flush=True)
ref_model = RegressionModel(ref_adata)

torch.cuda.empty_cache()

print(f"Training regression model", flush=True)
n_epochs = 10000
print(f"Number of epochs: {n_epochs}", flush=True)
ref_model.train(max_epochs=n_epochs,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1)

torch.cuda.empty_cache()

print(f"Exporting posterior", flush=True)
ref_adata = ref_model.export_posterior(ref_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})

# Save model
print(f"Saving models")
model_inf = "cell2location_ann/model_sea_10k"
ref_model.save(model_inf, overwrite=True)
adata_file = "cell2location_ann/sea_10k_ref_adata.h5ad"
ref_adata.write(adata_file)
print(f"Inference model for gene expression saved at: {model_inf}")
print(f"Anndata for inference model for gene expression saved at: {adata_file}")



