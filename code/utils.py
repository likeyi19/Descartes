import numpy as np
import scanpy as sc
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import scipy
from sklearn.metrics.cluster import silhouette_score
import anndata
import multiprocessing as mp
import tempfile
import subprocess
from sklearn.metrics import adjusted_rand_score,adjusted_mutual_info_score,homogeneity_score,normalized_mutual_info_score
from sklearn.feature_extraction.text import TfidfTransformer
import itertools
from scipy.io import mmwrite

from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.preprocessing import scale

import sklearn
import os

from scipy.spatial.distance import pdist
from multiprocessing import Pool
import scipy.cluster.hierarchy as sch

def tfidf2(count_mat): 
    tf_mat = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    signac_mat = np.log(1 + np.multiply(1e4*tf_mat,  np.tile((1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1]))))
#     return scipy.sparse.csr_matrix(signac_mat)
    return signac_mat

def tfidf1(count_mat): 
    nfreqs = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    tfidf_mat = np.multiply(nfreqs, np.tile(np.log(1 + 1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1])))
    return tfidf_mat

def tfidf3(count_mat): 
    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(np.transpose(count_mat))
    model.idf_ -= 1
    tf_idf = np.transpose(model.transform(np.transpose(count_mat)))
    return tf_idf.todense()


def silhouette(adata, group_key, embed, metric="euclidean", scale=True):
    """Average silhouette width (ASW)
    Wrapper for sklearn silhouette function values range from [-1, 1] with
        * 1 indicates distinct, compact clusters
        * 0 indicates overlapping clusters
        * -1 indicates core-periphery (non-cluster) structure
    By default, the score is scaled between 0 and 1 (``scale=True``).
    :param group_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    :param scale: default True, scale between 0 (worst) and 1 (best)
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")
    asw = silhouette_score(
        X=adata.obsm[embed], labels=adata.obs[group_key], metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw


def check_adata(adata):
    if type(adata) is not anndata.AnnData:
        raise TypeError("Input is not a valid AnnData object")


def check_batch(batch, obs, verbose=False):
    if batch not in obs:
        raise ValueError(f"column {batch} is not in obs")
    elif verbose:
        print(f"Object contains {obs[batch].nunique()} batches.")


def recompute_knn(adata, type_):
    """Recompute neighbours"""
    if type_ == "embed":
        return sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_emb", copy=True)
    elif type_ == "full":
        if "X_pca" not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver="arpack")
        return sc.pp.neighbors(adata, n_neighbors=15, copy=True)
    else:
        # if knn - do not compute a new neighbourhood graph (it exists already)
        return adata.copy()


def Hbeta(D_row, beta):
    """
    Helper function for simpson index computation
    """
    P = np.exp(-D_row * beta)
    sumP = np.nansum(P)
    if sumP == 0:
        H = 0
        P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.nansum(D_row * P) / sumP
        P /= sumP
    return H, P


def convert_to_one_hot(vector, num_classes=None):
    """
    Converts an input 1-D vector of integers into an output 2-D array of one-hot vectors,
    where an i'th input value of j will set a '1' in the i'th row, j'th column of the
    output array.
    Example:
    .. code-block:: python
        v = np.array((1, 0, 4))
        one_hot_v = convertToOneHot(v)
        print(one_hot_v)
    .. code-block::
        [[0 1 0 0 0]
         [1 0 0 0 0]
         [0 0 0 0 1]]
    """

    # assert isinstance(vector, np.ndarray)
    # assert len(vector) > 0

    if num_classes is None:
        num_classes = np.max(vector) + 1
    # else:
    #    assert num_classes > 0
    #    assert num_classes >= np.max(vector)

    result = np.zeros(shape=(len(vector), num_classes))
    result[np.arange(len(vector)), vector] = 1
    return result.astype(int)


def compute_simpson_index_graph(
    file_prefix=None,
    batch_labels=None,
    n_batches=None,
    n_neighbors=90,
    perplexity=30,
    chunk_no=0,
    tol=1e-5,
):
    """
    Simpson index of batch labels subset by group.
    :param file_prefix: file_path to pre-computed index and distance files
    :param batch_labels: a vector of length n_cells with batch info
    :param n_batches: number of unique batch labels
    :param n_neighbors: number of nearest neighbors
    :param perplexity: effective neighborhood size
    :param chunk_no: for parallelization, chunk id to evaluate
    :param tol: a tolerance for testing effective neighborhood size
    :returns: the simpson index for the neighborhood of each cell
    """
    index_file = file_prefix + "_indices_" + str(chunk_no) + ".txt"
    distance_file = file_prefix + "_distances_" + str(chunk_no) + ".txt"

    # initialize
    P = np.zeros(n_neighbors)
    logU = np.log(perplexity)

    # check if the target file is not empty
    if os.stat(index_file).st_size == 0:
        print("File has no entries. Doing nothing.")
        lists = np.zeros(0)
        return lists

    # read distances and indices with nan value handling
    indices = pd.read_table(index_file, index_col=0, header=None, sep=",")
    indices = indices.T

    distances = pd.read_table(distance_file, index_col=0, header=None, sep=",")
    distances = distances.T

    # get cell ids
    chunk_ids = indices.columns.values.astype("int")

    # define result vector
    simpson = np.zeros(len(chunk_ids))

    # loop over all cells in chunk
    for i, chunk_id in enumerate(chunk_ids):
        # get neighbors and distances
        # read line i from indices matrix
        get_col = indices[chunk_id]

        if get_col.isnull().sum() > 0:
            # not enough neighbors
            print(f"Chunk {chunk_id} does not have enough neighbors. Skipping...")
            simpson[i] = 1  # np.nan #set nan for testing
            continue

        knn_idx = get_col.astype("int") - 1  # get 0-based indexing

        # read line i from distances matrix
        D_act = distances[chunk_id].values.astype("float")

        # start lisi estimation
        beta = 1
        betamin = -np.inf
        betamax = np.inf

        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0

        # first get neighbor probabilities
        while np.logical_and(np.abs(Hdiff) > tol, tries < 50):
            if Hdiff > 0:
                betamin = beta
                if betamax == np.inf:
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if betamin == -np.inf:
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2

            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1

        if H == 0:
            simpson[i] = -1
            continue
            # then compute Simpson's Index
        batch = batch_labels[knn_idx]
        B = convert_to_one_hot(batch, n_batches)
        sumP = np.matmul(P, B)  # sum P per batch
        simpson[i] = np.dot(sumP, sumP)  # sum squares

    return simpson


def lisi_graph_py(
    adata,
    group_key,
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
):
    """
    Function to prepare call of compute_simpson_index
    Compute LISI score on shortes path based on kNN graph provided in the adata object.
    By default, perplexity is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """

    # use no more than the available cores
    n_cores = max(1, min(n_cores, mp.cpu_count()))

    if "neighbors" not in adata.uns:
        raise AttributeError(
            "Key 'neighbors' not found. Please make sure that a kNN graph has been computed"
        )
    elif verbose:
        print("using precomputed kNN graph")

    # get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")

    batch = adata.obs[group_key].cat.codes.values
    n_batches = len(np.unique(adata.obs[group_key]))

    if perplexity is None or perplexity >= n_neighbors:
        # use LISI default
        perplexity = np.floor(n_neighbors / 3)

    # setup subsampling
    subset = 100  # default, no subsampling
    if subsample is not None:
        subset = subsample  # do not use subsampling
        if isinstance(subsample, int) is False:  # need to set as integer
            subset = int(subsample)

    # run LISI in python
    if verbose:
        print("Compute knn on shortest paths")

    # set connectivities to 3e-308 if they are lower than 3e-308 (because cpp can't handle double values smaller than that).
    connectivities = adata.obsp["connectivities"]  # csr matrix format
    large_enough = connectivities.data >= 3e-308
    if verbose:
        n_too_small = np.sum(large_enough is False)
        if n_too_small:
            print(
                f"{n_too_small} connectivities are smaller than 3e-308 and will be set to 3e-308"
            )
            print(connectivities.data[large_enough is False])
    connectivities.data[large_enough is False] = 3e-308

    # temporary file
    tmpdir = tempfile.TemporaryDirectory(prefix="lisi_")
    prefix = tmpdir.name + "/graph_lisi"
    mtx_file_path = prefix + "_input.mtx"

    mmwrite(mtx_file_path, connectivities, symmetry="general")
    # call knn-graph computation in Cpp

    # root = pathlib.Path(scib.__file__).parent  # get current root directory
    
#     cpp_file_path = (
#         root / "knn_graph/knn_graph.o"
#     )  # create POSIX path to file to execute compiled cpp-code
    # comment: POSIX path needs to be converted to string - done below with 'as_posix()'
    # create evenly split chunks if n_obs is divisible by n_chunks (doesn't really make sense on 2nd thought)
    cpp_file_path = './knn_graph.o'
    args_int = [
        cpp_file_path, 
        #cpp_file_path.as_posix(),
        mtx_file_path,
        prefix,
        str(n_neighbors),
        str(n_cores),  # number of splits
        str(subset),
    ]
    if verbose:
        print(f'call {" ".join(args_int)}')
    try:
        subprocess.run(args_int)
    except rinterface.embedded.RRuntimeError as ex:
        print(f"Error computing LISI kNN graph {ex}\nSetting value to np.nan")
        return np.nan

    if verbose:
        print("LISI score estimation")

    if n_cores > 1:

        if verbose:
            print(f"{n_cores} processes started.")
        pool = mp.Pool(processes=n_cores)
        chunk_no = np.arange(0, n_cores)

        # create argument list for each worker
        results = pool.starmap(
            compute_simpson_index_graph,
            zip(
                itertools.repeat(prefix),
                itertools.repeat(batch),
                itertools.repeat(n_batches),
                itertools.repeat(n_neighbors),
                itertools.repeat(perplexity),
                chunk_no,
            ),
        )
        pool.close()
        pool.join()

        simpson_estimate_batch = np.concatenate(results)

    else:
        simpson_estimate_batch = compute_simpson_index_graph(
            file_prefix=prefix,
            batch_labels=batch,
            n_batches=n_batches,
            perplexity=perplexity,
            n_neighbors=n_neighbors,
        )

    tmpdir.cleanup()

    return 1 / simpson_estimate_batch


def clisi_graph(
    adata,
    label_key,
    k0=90,
    type_=None,
    subsample=None,
    scale=True,
    n_cores=1,
    verbose=False,
):
    """Cell-type LISI (cLISI) score
    Local Inverse Simpson’s Index metrics adapted from https://doi.org/10.1038/s41592-019-0619-0 to run on all full
    feature, embedding and kNN integration outputs via shortest path-based distance computation on single-cell kNN
    graphs.
    By default, this function returns a value scaled between 0 and 1 instead of the original LISI range of 0 to the
    number of labels.
    :param adata: adata object to calculate on
    :param group_key: group column name in ``adata.obs``
    :param k0: number of nearest neighbors to compute lisi score
        Please note that the initial neighborhood size that is
        used to compute shortest paths is 15.
    :param `type_`: type of data integration, either knn, full or embed
    :param subsample: Percentage of observations (integer between 0 and 100)
        to which lisi scoring should be subsampled
    :param scale: scale output values between 0 and 1 (True/False)
    :param n_cores: number of cores (i.e. CPUs or CPU cores to use for multiprocessing)
    :return: Median of cLISI scores per cell type labels
    """

    check_adata(adata)
    check_batch(label_key, adata.obs)

    adata_tmp = recompute_knn(adata, type_)

    scores = lisi_graph_py(
        adata=adata_tmp,
        group_key=label_key,
        n_neighbors=k0,
        perplexity=None,
        subsample=subsample,
        n_cores=n_cores,
        verbose=verbose,
    )

    # cLISI: 1 good, nlabs bad
    clisi = np.nanmedian(scores)

    if scale:
        nlabs = adata.obs[label_key].nunique()
        clisi = (nlabs - clisi) / (nlabs - 1)

    return clisi


def run_louvain(
    adata,
    label_key,
    cluster_key,
    range_min=0,
    range_max=3,
    max_steps=30,
    opt_function="NMI",
    resolutions=None,
    use_rep=None,
    inplace=True,
    plot=False,
    force=True,
    verbose=True,
    seed=0,
    **kwargs,
):
    """Optimised Louvain clustering
    Louvain clustering with resolution optimised against a metric
    :param adata: anndata object
    :param label_key: name of column in adata.obs containing biological labels to be
        optimised against
    :param cluster_key: name of column to be added to adata.obs during clustering.
        Will be overwritten if exists and ``force=True``
    :param resolutions: list of resolutions to be optimised over. If ``resolutions=None``,
        default resolutions of 20 values ranging between 0.1 and 2 will be used
    :param use_rep: key of embedding to use only if ``adata.uns['neighbors']`` is not
        defined, otherwise will be ignored
    :returns:
        Tuple of ``(res_max, score_max, score_all)`` or
        ``(res_max, score_max, score_all, clustering)`` if ``inplace=False``.
        ``res_max``: resolution of maximum score;
        ``score_max``: maximum score;
        ``score_all``: ``pd.DataFrame`` containing all scores at resolutions. Can be used to plot the score profile.
        ``clustering``: only if ``inplace=False``, return cluster assignment as ``pd.Series``
    """

    if verbose:
        print("Clustering...")


    if cluster_key in adata.obs.columns:
        if force:
            print(
                f"Warning: cluster key {cluster_key} already exists "
                "in adata.obs and will be overwritten"
            )
        else:
            raise ValueError(
                f"cluster key {cluster_key} already exists in "
                + "adata, please remove the key or choose a different name."
                + "If you want to force overwriting the key, specify `force=True`"
            )

    if resolutions is None:
        n = 20
        resolutions = [2 * x / n for x in range(1, n + 1)]

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    NMI_score_all = []
    ARI_score_all = []
    Homo_score_all = []
    AMI_score_all = []
    
    try:
        adata.uns["neighbors"]
    except KeyError:
        if verbose:
            print("computing neighbours for opt_cluster")
        sc.pp.neighbors(adata, use_rep=use_rep)
    # find optimal cluster parameter and metric
    n_cluster = np.unique(adata.obs[label_key]).shape[0]
    for res in resolutions:
        sc.tl.louvain(adata, resolution=res, key_added=cluster_key, random_state=seed)
        NMI_score = normalized_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
        ARI_score = adjusted_rand_score(adata.obs[label_key], adata.obs[cluster_key])
        Homo_score = homogeneity_score(adata.obs[label_key], adata.obs[cluster_key])
        AMI_score = adjusted_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
        #(adata, label_key, cluster_key, **kwargs)
        if verbose:
            print(f"resolution: {res}, NMI: {NMI_score}, ARI: {ARI_score}, Homo: {Homo_score}, AMI: {AMI_score}")
        if opt_function=='NMI':
            score = NMI_score
        elif opt_function=='ARI':
            score = ARI_score
        elif opt_function=='Homo':
            score = Homo_score
        elif opt_function=='AMI':
            score = AMI_score
#         score_all.append(score)
        NMI_score_all.append(NMI_score)
        ARI_score_all.append(ARI_score)
        Homo_score_all.append(Homo_score)
        AMI_score_all.append(AMI_score)
        if score_max < score:
            score_max = score
            [NMI_score_max, ARI_score_max, Homo_score_max, AMI_score_max] = [NMI_score, ARI_score, Homo_score, AMI_score]
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]

    if verbose:
        print(f"optimised clustering against {label_key}")
        print(f"optimal cluster resolution: {res_max}")
        print(f"selected optimal cluster metrics: {opt_function} {score_max}")
        print(f"NMI: {NMI_score_max}, ARI: {ARI_score_max}, Homo: {Homo_score_max}, AMI: {AMI_score_max}")
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        sc.tl.louvain(adata,resolution=this_resolution,key_added=cluster_key, random_state=seed)
        this_clusters = adata.obs[cluster_key].nunique()
        if this_clusters > n_cluster:
            this_max = this_resolution
        elif this_clusters < n_cluster:
            this_min = this_resolution
        else:
            NMI_score = normalized_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
            ARI_score = adjusted_rand_score(adata.obs[label_key], adata.obs[cluster_key])
            Homo_score = homogeneity_score(adata.obs[label_key], adata.obs[cluster_key])
            AMI_score = adjusted_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
            print("louvain clustering with a binary search")
            print(f"NMI: {NMI_score}, ARI: {ARI_score}, Homo: {Homo_score}, AMI: {AMI_score}")
            resolutions = [this_resolution] + resolutions
            NMI_score_all = [NMI_score] + NMI_score_all
            ARI_score_all = [ARI_score] + ARI_score_all
            Homo_score_all = [Homo_score] + Homo_score_all
            AMI_score_all = [AMI_score] + AMI_score_all
            break
        this_step += 1
    resolutions = [res_max] + resolutions
    NMI_score_all = [NMI_score_max] + NMI_score_all
    ARI_score_all = [ARI_score_max] + ARI_score_all
    Homo_score_all = [Homo_score_max] + Homo_score_all
    AMI_score_all = [AMI_score_max] + AMI_score_all
    score_all = pd.DataFrame(
        zip(resolutions, NMI_score_all, ARI_score_all, Homo_score_all, AMI_score_all), columns=("resolution", "NMI", "ARI", "Homo", "AMI")
    )
    if plot:
        # score vs. resolution profile
        sns.lineplot(data=score_all, x="resolution", y="NMI").set_title(
            "Optimal cluster resolution profile"
        )
        plt.show()

    if inplace:
        adata.obs[cluster_key] = clustering
        return score_all
    else:
        return score_all, clustering

def run_leiden(
    adata,
    label_key,
    cluster_key,
    range_min=0,
    range_max=3,
    max_steps=30,
    opt_function="NMI",
    resolutions=None,
    use_rep=None,
    inplace=True,
    plot=False,
    force=True,
    verbose=True,
    seed=0,
    **kwargs,
):
    """Optimised Louvain clustering
    Louvain clustering with resolution optimised against a metric
    :param adata: anndata object
    :param label_key: name of column in adata.obs containing biological labels to be
        optimised against
    :param cluster_key: name of column to be added to adata.obs during clustering.
        Will be overwritten if exists and ``force=True``
    :param resolutions: list of resolutions to be optimised over. If ``resolutions=None``,
        default resolutions of 20 values ranging between 0.1 and 2 will be used
    :param use_rep: key of embedding to use only if ``adata.uns['neighbors']`` is not
        defined, otherwise will be ignored
    :returns:
        Tuple of ``(res_max, score_max, score_all)`` or
        ``(res_max, score_max, score_all, clustering)`` if ``inplace=False``.
        ``res_max``: resolution of maximum score;
        ``score_max``: maximum score;
        ``score_all``: ``pd.DataFrame`` containing all scores at resolutions. Can be used to plot the score profile.
        ``clustering``: only if ``inplace=False``, return cluster assignment as ``pd.Series``
    """

    if verbose:
        print("Clustering...")


    if cluster_key in adata.obs.columns:
        if force:
            print(
                f"Warning: cluster key {cluster_key} already exists "
                "in adata.obs and will be overwritten"
            )
        else:
            raise ValueError(
                f"cluster key {cluster_key} already exists in "
                + "adata, please remove the key or choose a different name."
                + "If you want to force overwriting the key, specify `force=True`"
            )

    if resolutions is None:
        n = 20
        resolutions = [2 * x / n for x in range(1, n + 1)]

    score_max = 0
    res_max = resolutions[0]
    clustering = None
    NMI_score_all = []
    ARI_score_all = []
    Homo_score_all = []
    AMI_score_all = []
    
    try:
        adata.uns["neighbors"]
    except KeyError:
        if verbose:
            print("computing neighbours for opt_cluster")
        sc.pp.neighbors(adata, use_rep=use_rep)
    # find optimal cluster parameter and metric
    n_cluster = np.unique(adata.obs[label_key]).shape[0]
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=cluster_key, random_state=seed)
        NMI_score = normalized_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
        ARI_score = adjusted_rand_score(adata.obs[label_key], adata.obs[cluster_key])
        Homo_score = homogeneity_score(adata.obs[label_key], adata.obs[cluster_key])
        AMI_score = adjusted_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
        #(adata, label_key, cluster_key, **kwargs)
        if verbose:
            print(f"resolution: {res}, NMI: {NMI_score}, ARI: {ARI_score}, Homo: {Homo_score}, AMI: {AMI_score}")
        if opt_function=='NMI':
            score = NMI_score
        elif opt_function=='ARI':
            score = ARI_score
        elif opt_function=='Homo':
            score = Homo_score
        elif opt_function=='AMI':
            score = AMI_score
#         score_all.append(score)
        NMI_score_all.append(NMI_score)
        ARI_score_all.append(ARI_score)
        Homo_score_all.append(Homo_score)
        AMI_score_all.append(AMI_score)
        if score_max < score:
            score_max = score
            [NMI_score_max, ARI_score_max, Homo_score_max, AMI_score_max] = [NMI_score, ARI_score, Homo_score, AMI_score]
            res_max = res
            clustering = adata.obs[cluster_key]
        del adata.obs[cluster_key]

    if verbose:
        print(f"optimised clustering against {label_key}")
        print(f"optimal cluster resolution: {res_max}")
        print(f"selected optimal cluster metrics: {opt_function} {score_max}")
        print(f"NMI: {NMI_score_max}, ARI: {ARI_score_max}, Homo: {Homo_score_max}, AMI: {AMI_score_max}")
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        sc.tl.leiden(adata,resolution=this_resolution,key_added=cluster_key, random_state=seed)
        this_clusters = adata.obs[cluster_key].nunique()
        if this_clusters > n_cluster:
            this_max = this_resolution
        elif this_clusters < n_cluster:
            this_min = this_resolution
        else:
            NMI_score = normalized_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
            ARI_score = adjusted_rand_score(adata.obs[label_key], adata.obs[cluster_key])
            Homo_score = homogeneity_score(adata.obs[label_key], adata.obs[cluster_key])
            AMI_score = adjusted_mutual_info_score(adata.obs[label_key], adata.obs[cluster_key])
            print("leiden clustering with a binary search")
            print(f"NMI: {NMI_score}, ARI: {ARI_score}, Homo: {Homo_score}, AMI: {AMI_score}")
            resolutions = [this_resolution] + resolutions
            NMI_score_all = [NMI_score] + NMI_score_all
            ARI_score_all = [ARI_score] + ARI_score_all
            Homo_score_all = [Homo_score] + Homo_score_all
            AMI_score_all = [AMI_score] + AMI_score_all
            break
        this_step += 1
    resolutions = [res_max] + resolutions
    NMI_score_all = [NMI_score_max] + NMI_score_all
    ARI_score_all = [ARI_score_max] + ARI_score_all
    Homo_score_all = [Homo_score_max] + Homo_score_all
    AMI_score_all = [AMI_score_max] + AMI_score_all
    score_all = pd.DataFrame(
        zip(resolutions, NMI_score_all, ARI_score_all, Homo_score_all, AMI_score_all), columns=("resolution", "NMI", "ARI", "Homo", "AMI")
    )
    if plot:
        # score vs. resolution profile
        sns.lineplot(data=score_all, x="resolution", y="NMI").set_title(
            "Optimal cluster resolution profile"
        )
        plt.show()

    if inplace:
        adata.obs[cluster_key] = clustering
        return score_all
    else:
        return score_all, clustering

def run_descartes(anndatas,
                select_num,              # 最后要选多少个peak
                seed_base,               # 随机数种子
                tfidf="tfidf2",          # tfidf的方法
                ifPCA=True,              # 是否进行pca
                pc=10,
                k=20,
                similarity="cosine",
                iters=4,
                PCA_select=False,        # 是否去除第一个pc
                spmethod='threshold',          # 空间数据用什么方法构建相似矩阵
                neighbor=5,
                sp_dist = 'recip',    # 距离的度量，const表示常数，recip表示倒数，recip_sq表述平方倒数
                pre_select='highest',         # 预训练选择peak的方法
                peaks_num=None,          # 预训练选择peak的数量
                distance='euclidean',
                cosine_filter=None,
                r=0.4):
    print(anndatas)
    selected_peaks_all = []
    anndatas.X = scipy.sparse.csc_matrix(anndatas.X)
    if tfidf == "tfidf1":
        count_mat_test = anndatas.X.todense().T
        nfreqs = 1.0 * count_mat_test / np.tile(np.sum(count_mat_test,axis=0), (count_mat_test.shape[0],1))
        ATAC_count_test = np.multiply(nfreqs, np.tile(np.log(1 + 1.0 * count_mat_test.shape[1] / np.sum(count_mat_test,axis=1)).reshape(-1,1), (1,count_mat_test.shape[1])))
    elif tfidf == "tfidf2":
        count_mat_test = anndatas.X.todense().T
        tf_mat = 1.0 * count_mat_test / np.tile(np.sum(count_mat_test,axis=0), (count_mat_test.shape[0],1))
        ATAC_count_test = np.log(1 + np.multiply(1e4*tf_mat,  np.tile((1.0 * count_mat_test.shape[1] / np.sum(count_mat_test,axis=1)).reshape(-1,1), (1,count_mat_test.shape[1]))))
    else:
        count_mat_test = anndatas.X.todense().T
        model = TfidfTransformer(smooth_idf=False, norm="l2")
        model = model.fit(np.transpose(count_mat_test))
        model.idf_ -= 1
        tf_idf = np.transpose(model.transform(np.transpose(count_mat_test)))
        ATAC_count_test = tf_idf.todense()
    ATAC_count_test = np.array(ATAC_count_test)
    ATAC_count_test = ATAC_count_test.T
    count = ATAC_count_test.copy()
    print(count.shape)
    similarity_matrix_acb = np.zeros([count.shape[0], count.shape[0]])
    similarity_matrix_spatial = np.zeros([count.shape[0], count.shape[0]])
    for item in range(iters):
        if item != 0:
            if peaks_num != None:
                count = ATAC_count_test[:,sorted_index[-peaks_num:]]
            else:
                count = ATAC_count_test[:,sorted_index[-select_num:]]
        else:
            if pre_select == 'highest':
                peak_sum = np.sum(anndatas.X, axis=0)
                peak_sum = np.array(peak_sum).reshape(-1)
                idx = np.argsort(peak_sum)[::-1]
                if peaks_num != None:
                    idx_highest = idx[:peaks_num]
                else:
                    idx_highest = idx[:peaks_num]
                count = ATAC_count_test[:,idx_highest]
            elif pre_select == 'random':
                if peaks_num != None:
                    idx = np.random.choice(ATAC_count_test.shape[1], size=peaks_num, replace=False)
                else:
                    idx = np.random.choice(ATAC_count_test.shape[1], size=select_num, replace=False)
                count = ATAC_count_test[:,idx]
                
        distance_matrix = np.zeros([count.shape[0], count.shape[0]])
        print(distance_matrix.shape)
    
        if ifPCA:
            count = PCA(n_components=pc,random_state=int(seed_base*1000)).fit_transform(count)
            if PCA_select:
                count = count[:,1:]
            print(count.shape)
            #KNN算法默认使用欧几里得距离
            if similarity == 'Jaccard':
                if distance == 'euclidean':
                    diff = count[:, np.newaxis, :] - count[np.newaxis, :, :]
                    distance_matrix = np.linalg.norm(diff, axis=2)
                    np.fill_diagonal(distance_matrix, np.inf)
                elif distance == 'cosine':
                    count_norm = np.linalg.norm(count, axis=1, keepdims=True)
                    dot_product_matrix = np.dot(count, count.T)
                    distance_matrix = 1 - dot_product_matrix / (count_norm * count_norm.T)
                    np.fill_diagonal(distance_matrix, np.inf)
        else:
            distance_matrix = scipy.spatial.distance.pdist(count, metric=distance)
            distance_matrix = scipy.spatial.distance.squareform(distance_matrix)
            for i in range(count.shape[0]):
                distance_matrix[i][i] = float('inf')

        distance_matrix_acb = np.copy(distance_matrix)
        
        spatial_data = anndatas.obsm['spatial']
        diff = spatial_data[:, np.newaxis, :] - spatial_data[np.newaxis, :, :]
        distance_matrix = np.linalg.norm(diff, axis=2)
        np.fill_diagonal(distance_matrix, np.inf)
        distance_matrix_spatial = np.copy(distance_matrix)

        similarity_matrix_acb = np.zeros([count.shape[0], count.shape[0]])
        similarity_matrix_spatial = np.zeros([count.shape[0], count.shape[0]])
        if spmethod == 'threshold':
            min_dist = np.min(distance_matrix_spatial)
            if sp_dist == 'const':
                similarity_matrix_spatial = np.array(distance_matrix_spatial <= neighbor * min_dist,dtype=int)
            elif sp_dist == 'recip':
                similarity_matrix_spatial = min_dist / distance_matrix_spatial
                similarity_matrix_spatial = similarity_matrix_spatial * (distance_matrix_spatial <= neighbor * min_dist)
            else:
                similarity_matrix_spatial = (min_dist / distance_matrix_spatial)**2
                similarity_matrix_spatial = similarity_matrix_spatial * (distance_matrix_spatial <= neighbor * min_dist)
        elif spmethod == 'SNN':
            neighbor_index = np.argsort(distance_matrix_spatial, axis=1)[:,0:k]
            if similarity == "Jaccard":
                for i in range(count.shape[0]):
                    for j in range(i):
                        intersect_num = len(np.intersect1d(neighbor_index[i,:], neighbor_index[j,:]))
                        similarity_matrix_spatial[i][j] = 1.0*intersect_num/(2*k-intersect_num)
                        similarity_matrix_spatial[j][i] = 1.0*intersect_num/(2*k-intersect_num)
        if similarity == "Jaccard":
            neighbor_index = np.argsort(distance_matrix_acb, axis=1)[:,0:k]
            for i in range(count.shape[0]):
                for j in range(i):
                    intersect_num = len(np.intersect1d(neighbor_index[i,:], neighbor_index[j,:]))
                    similarity_matrix_acb[i][j] = 1.0*intersect_num/(2*k-intersect_num)
                    similarity_matrix_acb[j][i] = 1.0*intersect_num/(2*k-intersect_num)
        elif similarity == "cosine":
            similarity_matrix_acb = sklearn.metrics.pairwise.cosine_similarity(count)
            for i in range(count.shape[0]):
                similarity_matrix_acb[i][i] = -float('inf')
            neighbor_index = np.argsort(similarity_matrix_acb, axis=1)[:,0:(count.shape[0]-k)]
            for i in range(count.shape[0]):
                similarity_matrix_acb[i,neighbor_index[i,:]] = 0
            print(similarity_matrix_acb[0,0])
            if cosine_filter != None:
                similarity_matrix_acb[similarity_matrix_acb < cosine_filter] = 0
        similarity_matrix = (1 - r) * similarity_matrix_acb + r*similarity_matrix_spatial
                    
        #对所有peaks进行计算，这里是否考虑提前筛掉一些peaks
        print('compute scores')
        scores = np.zeros(anndatas.n_vars)
        #z-score
        X_processed = scale(ATAC_count_test, axis=0)
        temp_matrix = np.matmul(similarity_matrix, X_processed)
        scores = np.sum(X_processed * temp_matrix, axis=0)
        sorted_index = np.argsort(scores)
        selected_peaks = sorted_index[-select_num:]
        selected_peaks_all.append(selected_peaks)
        print('scores:')
        print(scores[selected_peaks])

    return selected_peaks, sorted_index, similarity_matrix, selected_peaks_all, scores, X_processed[:,selected_peaks], similarity_matrix_acb, similarity_matrix_spatial

# 获取当前进程的 PID
def get_peak_memory_usage(pid):
    """
    读取 /proc/[PID]/status 文件来获取进程的峰值内存占用（VmHWM）
    """
    with open(f"/proc/{pid}/status") as f:
        for line in f:
            if "VmHWM" in line:
                # 提取内存使用量并转换为 MB
                return int(line.split()[1]) / 1024
            

def fx_PAS(clusterlabel, location):
    matched_location = location.copy()
    NAs = np.where(np.isnan(clusterlabel))[0]
    
    if len(NAs) > 0:
        clusterlabel = np.delete(clusterlabel, NAs)
        matched_location = np.delete(matched_location, NAs, axis=0)

    with Pool() as pool:
        results = pool.starmap(fx_kNN, zip(range(matched_location.shape[0]), [matched_location]*matched_location.shape[0], [10]*matched_location.shape[0], [clusterlabel]*matched_location.shape[0]))

    return np.sum(results) / len(clusterlabel)

def fx_kNN(i, location_in, k, cluster_in):
    line_i = np.zeros(location_in.shape[0])
    
    for j in range(location_in.shape[0]):
        line_i[i] = float('inf')
        if i != j:
            pairwise_distances = pdist(np.vstack([location_in[i, :], location_in[j, :]]))
            line_i[j] = pairwise_distances[0]

    ind = np.argsort(line_i)[:k]
    cluster_use = np.delete(cluster_in, i)
    ind = ind - (ind > i)

    if np.sum(cluster_use[ind] != cluster_in[i]) > (k / 2):
        return 1
    else:
        return 0
    

def fx_1NN(args):
    """
    Helper function to calculate 1NN distances for a specific spot.
    """
    i, location_in = args
    line_i = np.zeros(location_in.shape[0])

    for j in range(location_in.shape[0]):
        line_i[i] = float('inf')
        if i != j:
            pairwise_distances = pdist(np.vstack([location_in[i, :], location_in[j, :]]))
            line_i[j] = pairwise_distances[0]
    return np.min(line_i)

def fx_CHAOS(clusterlabel, location):
    """
    Calculate CHAOS metric based on cluster labels and spatial locations.

    Parameters:
    - clusterlabel: 1D array, cluster labels for each spot.
    - location: 2D array, spatial locations of spots.

    Returns:
    - chaos_value: float, calculated CHAOS metric.
    """
    matched_location = location.copy()
    NAs = np.where(np.isnan(clusterlabel))[0]

    if len(NAs) > 0:
        clusterlabel = np.delete(clusterlabel, NAs)
        matched_location = np.delete(matched_location, NAs, axis=0)

    matched_location = scale(matched_location)
    dist_val = np.zeros(len(np.unique(clusterlabel)))
    count = 0

    for k in np.unique(clusterlabel):
        count += 1
        location_cluster = matched_location[clusterlabel == k, :]

        if len(location_cluster) == 2:
            continue

        with Pool() as pool:
            results = pool.map(fx_1NN, [(i, location_cluster) for i in range(len(location_cluster))])

        dist_val[count - 1] = np.sum(results)

    dist_val = dist_val[~np.isnan(dist_val)]
    return np.sum(dist_val) / len(clusterlabel)


def run_leiden_(
    adata,
    label_key,
    cluster_key,
    range_min=0,
    range_max=3,
    max_steps=30,
    opt_function="NMI",
    resolutions=None,
    use_rep=None,
    inplace=False,
    plot=False,
    force=True,
    verbose=True,
    seed=0,
    n_cluster=None,
    **kwargs,
):
    """Optimised Louvain clustering
    Louvain clustering with resolution optimised against a metric
    :param adata: anndata object
    :param label_key: name of column in adata.obs containing biological labels to be
        optimised against
    :param cluster_key: name of column to be added to adata.obs during clustering.
        Will be overwritten if exists and ``force=True``
    :param resolutions: list of resolutions to be optimised over. If ``resolutions=None``,
        default resolutions of 20 values ranging between 0.1 and 2 will be used
    :param use_rep: key of embedding to use only if ``adata.uns['neighbors']`` is not
        defined, otherwise will be ignored
    :returns:
        Tuple of ``(res_max, score_max, score_all)`` or
        ``(res_max, score_max, score_all, clustering)`` if ``inplace=False``.
        ``res_max``: resolution of maximum score;
        ``score_max``: maximum score;
        ``score_all``: ``pd.DataFrame`` containing all scores at resolutions. Can be used to plot the score profile.
        ``clustering``: only if ``inplace=False``, return cluster assignment as ``pd.Series``
    """

    if verbose:
        print("Clustering...")


    if cluster_key in adata.obs.columns:
        if force:
            print(
                f"Warning: cluster key {cluster_key} already exists "
                "in adata.obs and will be overwritten"
            )
        else:
            raise ValueError(
                f"cluster key {cluster_key} already exists in "
                + "adata, please remove the key or choose a different name."
                + "If you want to force overwriting the key, specify `force=True`"
            )

    
    try:
        adata.uns["neighbors"]
    except KeyError:
        if verbose:
            print("computing neighbours for opt_cluster")
        sc.pp.neighbors(adata, use_rep=use_rep)
    # find optimal cluster parameter and metric
    # n_cluster = np.unique(adata.obs[label_key]).shape[0]
    # n_cluster = 10
    
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        sc.tl.leiden(adata,resolution=this_resolution,key_added=cluster_key, random_state=seed)
        this_clusters = adata.obs[cluster_key].nunique()
        if this_clusters > n_cluster:
            this_max = this_resolution
        elif this_clusters < n_cluster:
            this_min = this_resolution
        else:
            break
        this_step += 1
        
    return


def peak_modules_(selected_peaks_data, similarity_matrix, method):
    pdist = np.matmul(np.matmul(selected_peaks_data.T, similarity_matrix), selected_peaks_data) \
    +np.matmul(np.matmul(selected_peaks_data.T, similarity_matrix), selected_peaks_data)
    return pdist


def peak_module_cluster(adata,idx,peak_distance,method,num_clusters = 10):

    Z = sch.linkage(peak_distance, method=method) 
    dendrogram = sch.dendrogram(Z)
    
    labels = sch.fcluster(Z, t=num_clusters, criterion='maxclust')
    print("Cluster labels:", labels)

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Peaks')
    plt.ylabel('Distance')
    plt.show()

    unique, counts = np.unique(labels, return_counts=True)

    total_samples = len(labels)
    proportions = counts / total_samples

    for num, count, prop in zip(unique, counts, proportions):
        print(f"Number {num}: Count = {count}, Proportion = {prop:.2%}")

    plt.bar(unique, proportions)
    plt.xlabel('Number')
    plt.ylabel('Proportion')
    plt.title('Proportion of Each Number')
    plt.xticks(unique)
    plt.show()
    sc.pl.spatial(adata, color=['label'], spot_size=0.5)
    ATAC_matrix = scale(np.array(adata.X.todense()))
    coord_x = np.array(adata.obsm['spatial'][:,0])
    coord_y = -np.array(adata.obsm['spatial'][:,1])
    rand_num = 9999
    for i in np.unique(labels):
        print("Number of peaks in module %d: %d"%(i,counts[i-1]))
        print("Visualize peak module of label %d:"%i)
        index = (np.array(labels)==i)
        color = np.mean(ATAC_matrix[:,idx[index]],axis = 1).squeeze()
        print(rand_num)
        plt.scatter(coord_x, coord_y, c=color, s=3, cmap='Blues')
        plt.colorbar()
        plt.show()
    return Z, labels
