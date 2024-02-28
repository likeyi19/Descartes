from utils import *
import argparse
import episcanpy.api as epi
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fp', '--file_path', type=str, default=None)
    parser.add_argument('-sp', '--save_path', type=str, default=None)

    parser.add_argument('-n', '--num_select_peak', type=int, default=10000)
    parser.add_argument('-s', '--seed_base', type=int, default=1)
    parser.add_argument('-tf', '--TF_IDF', type=str, default='tfidf2')
    parser.add_argument('-pc', '--pc_number', type=int, default=10)
    parser.add_argument('-k', '--k_number', type=int, default=20)
    parser.add_argument('-s', '--similarity', type=str, default='cosine')
    parser.add_argument('-iter', '--iter_time', type=int, default=4)
    parser.add_argument('-spm', '--sp_method', type=str, default='threshold')
    parser.add_argument('-n', '--neighbor', type=int, default=5)
    parser.add_argument('-spd', '--sp_dist', type=str, default='recip')
    parser.add_argument('-ps', '--pre_select', type=str, default='highest')
    parser.add_argument('-pn', '--peaks_num', type=int, default=50000)
    parser.add_argument('-d', '--distance', type=str, default='euclidean')
    parser.add_argument('-r', '--ratio', type=float, default=0.4)

    opt = parser.parse_args()
    file_path = opt.file_path
    save_path = opt.save_path
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    num_select_peak = opt.num_select_peak
    seed_base = opt.seed_base
    tf = opt.TF_IDF
    pc = opt.pc_number
    k = opt.k_number
    similarity = opt.similarity
    iter_time = opt.iter_time
    spmethod = opt.sp_method
    neighbor = opt.neighbor
    sp_dist = opt.sp_dist
    pre_select = opt.pre_select
    peaks_num = opt.peaks_num
    distance = opt.distance
    r = opt.ratio

    adata_raw = sc.read_h5ad(file_path)
    print('load data: ', adata_raw)

    epi.pp.filter_features(adata_raw, min_cells=1)
    epi.pp.filter_cells(adata_raw, min_features=1)
    num_all_peak = adata_raw.n_vars
    print('data after pre-filtering: ', adata_raw)

    start_time = time.time()
    
    adata = sc.AnnData(adata_raw.X,dtype = 'float32')
    idx, sorted_index, simi_matrix, idx_all, scores, selected_peaks_data,similarity_matrix_acb, similarity_matrix_spatial = run_descartes(adata_raw, num_select_peak, seed_base=seed_base, tfidf=tf, ifPCA=True, pc=pc, k=k, similarity=similarity, iters=iter_time, spmethod=spmethod,neighbor=neighbor,sp_dist=sp_dist, pre_select=pre_select, peaks_num=peaks_num, distance=distance,r=r)
    
    end_time = time.time()

    result = pd.DataFrame({'idx':idx})
    filename = save_path + 'result_idx.csv'
    result.to_csv(filename,header=True)

    current_pid = os.getpid()
    peak_memory = get_peak_memory_usage(current_pid)
    run_time = end_time - start_time
    meta = pd.DataFrame({'Peak memory':[peak_memory], 'run time':[run_time]})
    meta.to_csv(save_path + 'result_meta.csv',index=None)
