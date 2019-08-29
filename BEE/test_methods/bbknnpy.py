#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# @author: Drizzle_Zhang
# @file: bbknnpy.py
# @time: 8/29/19 8:42 AM

import anndata
import scanpy as sc
import bbknn


def bbknn_py(pca_input, batch, pc_num):
    adata = anndata.AnnData(X=pca_input, obs=batch)
    sc.tl.pca(adata, n_comps=int(pc_num))
    adata.obsm['X_pca'] = pca_input
    bbknn.bbknn(adata, batch_key=0)
    sc.tl.umap(adata)

    return adata.obsm['X_umap']
