#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
优化版 HSP 并行计算脚本：
- Numba JIT 加速延伸核心
- 流式任务生成，减少峰值内存
- 半核进程池，动态 chunksize
- orjson 序列化
"""

import os
import sys
import time
import glob
import math
import shutil
from multiprocessing import Pool, cpu_count, set_start_method

import pandas as pd
from Bio import SeqIO, Align
import numpy as np
try:
    import orjson
    _dumps = orjson.dumps
except ImportError:
    import json
    _dumps = lambda obj: json.dumps(obj, separators=(",", ":"), ensure_ascii=False).encode()

from numba import njit

# ======== 用户配置区，请根据实际路径修改 ========
QUERY_FASTA      = "/Users/liweihang/BT1051/BLAST_Algorithm/TestFile/collagen.fasta"
DB_FASTA         = "/Users/liweihang/BT1051/BLAST_Algorithm/TestFile/uniprot_sprot.fasta"
VARIANT_CSV      = "/Users/liweihang/BT1051/BLAST_Algorithm/Database/all_3mer_variants.csv"
HASH_FOLDER      = "/Users/liweihang/BT1051/BLAST_Algorithm/Database/Final_hash_chunks_Table"
OUTPUT_DIR       = "/Users/liweihang/BT1051/BLAST_Algorithm/Database/HSP"

K                = 3
DROPOFF_RATIO    = 0.30
MATRIX_NAME      = "BLOSUM62"
LAMBDA_VAL       = 0.323
K_COEFF          = 0.138
E_THRESHOLD      = 1e-5
# ================================================

# 强制 fork，以继承全局数据
try:
    set_start_method('fork')
except RuntimeError:
    pass

# 构建氨基酸索引
AA_LIST = list("ARNDCQEGHILKMFPSTWYV")
AA2IDX = {aa: i for i, aa in enumerate(AA_LIST)}

# 加载替代矩阵并转换为 numpy 数组
_blosum = Align.substitution_matrices.load(MATRIX_NAME)
blosum_arr = np.full((20,20), -4, dtype=np.int16)
for (a,b), v in _blosum.items():
    if a in AA2IDX and b in AA2IDX:
        blosum_arr[AA2IDX[a], AA2IDX[b]] = v

# Numba 加速核心
@njit
def extend_hit_jit(q_idx, d_idx, qs, ds, k, init_s, dropoff):
    max_s = init_s
    cur_s = init_s
    ql, qr = qs, qs + k - 1
    dl, dr = ds, ds + k - 1
    best = (ql, qr, dl, dr)
    # 右延伸
    while qr + 1 < q_idx.shape[0] and dr + 1 < d_idx.shape[0]:
        qr += 1; dr += 1
        cur_s += blosum_arr[q_idx[qr], d_idx[dr]]
        if cur_s > max_s:
            max_s = cur_s
            best = (ql, qr, dl, dr)
        elif cur_s < max_s * dropoff:
            break
    # 左延伸
    cur_s = init_s
    ql, qr = qs, qs + k - 1
    dl, dr = ds, ds + k - 1
    while ql > 0 and dl > 0:
        ql -= 1; dl -= 1
        cur_s += blosum_arr[q_idx[ql], d_idx[dl]]
        if cur_s > max_s:
            max_s = cur_s
            best = (ql, qr, dl, dr)
        elif cur_s < max_s * dropoff:
            break
    return max_s, best[0], best[1], best[2], best[3]

# 1) 预加载 variant-scores
def load_variant_scores(path):
    df = pd.read_csv(path, usecols=['Original_3mer','Variant','BLOSUM62_Score'])
    valid = set(AA2IDX.keys())
    df = df[df['Original_3mer'].str.len().eq(3) &
            df['Variant'].str.len().eq(3) &
            df['Original_3mer'].apply(lambda s: all(c in valid for c in s)) &
            df['Variant'].apply(lambda s: all(c in valid for c in s))]
    return {(r.Original_3mer, r.Variant): r.BLOSUM62_Score
            for r in df.itertuples()}

VARIANT_SCORES = load_variant_scores(VARIANT_CSV)

# 2) 预加载数据库序列，并转为索引数组
DB_SEQS = []
for rec in SeqIO.parse(DB_FASTA, 'fasta'):
    seq = str(rec.seq)
    idx = np.array([AA2IDX.get(c, -1) for c in seq], dtype=np.int8)
    DB_SEQS.append((rec.id, seq, idx))

# 3) 预加载 query 序列及 kmer 位置，并转为索引
_qrec = next(SeqIO.parse(QUERY_FASTA, 'fasta'))
QUERY_ID = _qrec.id
QUERY_SEQ = str(_qrec.seq)
QUERY_IDX = np.array([AA2IDX.get(c, -1) for c in QUERY_SEQ], dtype=np.int8)
QUERY_KMER_POS = {}
for i in range(len(QUERY_SEQ) - K + 1):
    kmer = QUERY_SEQ[i:i+K]
    QUERY_KMER_POS.setdefault(kmer, []).append(i)

# 4) 子进程初始化（无需额外对象）
def init_worker():
    pass

# 5) 单个任务：JIT 延伸 + E-value 过滤
def process_task(task):
    kmer, qpos, dbidx, dpos = task
    if dbidx >= len(DB_SEQS):
        return None
    db_id, db_seq, db_idx = DB_SEQS[dbidx]
    if dpos + K > db_idx.shape[0]:
        return None
    frag_idx = db_idx[dpos:dpos+K]
    if np.any(frag_idx < 0):
        return None
    sc0 = VARIANT_SCORES.get((kmer, db_seq[dpos:dpos+K]))
    if sc0 is None or sc0 < 0:
        return None
    # 调用 JIT 核心
    max_s, ql, qr, dl, dr = extend_hit_jit(QUERY_IDX, db_idx, qpos, dpos, K, sc0, DROPOFF_RATIO)
    E = K_COEFF * QUERY_IDX.shape[0] * db_idx.shape[0] * math.exp(-LAMBDA_VAL * max_s)
    if E > E_THRESHOLD:
        return None
    return {
        'query_id':    QUERY_ID,
        'query_range': (int(ql), int(qr)),
        'query_frag':  QUERY_SEQ[ql:qr+1],
        'db_id':       db_id,
        'db_index':    dbidx,
        'db_range':    (int(dl), int(dr)),
        'db_frag':     db_seq[dl:dr+1],
        'score':       float(max_s),
        'e_value':     float(E),
        'seed_kmer':   kmer
    }

# 流式生成任务
def iter_tasks(hash_path):
    with open(hash_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            kmer, locs = line.split("\t")
            if kmer not in QUERY_KMER_POS:
                continue
            hits = [tuple(map(int, x.split(","))) for x in locs.split(";") if x]
            for qpos in QUERY_KMER_POS[kmer]:
                for dbidx, dpos in hits:
                    yield (kmer, qpos, dbidx, dpos)

# 6) 处理单个 hash chunk
def process_hash_chunk(hash_path):
    chunk = os.path.splitext(os.path.basename(hash_path))[0]
    out_file = os.path.join(OUTPUT_DIR, f"{chunk}_hsp.json")
    results = []
    n_proc = max(1, cpu_count() // 2)
    chunksize = 1024
    with Pool(n_proc, initializer=init_worker) as pool:
        for r in pool.imap_unordered(process_task, iter_tasks(hash_path), chunksize):
            if r:
                results.append(r)
    # 序列化并写文件
    with open(out_file, "wb") as fw:
        fw.write(_dumps(results))
    return len(results)

def main():
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    hash_files = sorted(glob.glob(os.path.join(HASH_FOLDER, "hash_chunk_*_withVariants.txt")))
    print(f"🔍 检测到 {len(hash_files)} 个哈希分块")
    start = time.time()
    total = 0
    for hp in hash_files:
        cnt = process_hash_chunk(hp)
        total += cnt
        print(f"  ▸ {os.path.basename(hp)} → {cnt} 条 HSP")
    elapsed = time.time() - start
    print(f"\n🎉 完成，总 HSP 数: {total}")
    print(f"⏱️ 用时: {elapsed:.2f} 秒")
    print(f"📁 输出: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
