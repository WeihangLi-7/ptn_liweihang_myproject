#!/usr/bin/env python3
"""
Standalone optimized pipeline:
- Optionally skip query-based variant filtering to build a full hash DB.
- Parallel hash construction, in-memory index, and final deduplication.
Steps:
1. (可选) Filter & simplify variants  
2. Build hash chunks  
3. Load index & test queries  
4. Deduplicate "withVariants" chunk files
"""
import os
import shutil
import csv
from glob import glob
from collections import defaultdict
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters
from tqdm import tqdm

# ---------------- Configuration ----------------
k = 3
T_THRESHOLD = 13
FULL_DB = True   # True: 忽略 query, 构建完整哈希数据库

# 文件路径配置
QUERY_FASTA = '/Users/liweihang/BT1051/BLAST_Algorithm/TestFile/8HHU.fa'
TARGET_FASTA = '/Users/liweihang/BT1051/BLAST_Algorithm/TestFile/uniprot_sprot.fasta'
VARIANT_CSV = '/Users/liweihang/BT1051/BLAST_Algorithm/Database/all_3mer_variants.csv'
FILTERED_VARIANTS = '/Users/liweihang/BT1051/BLAST_Algorithm/Database/filtered_variants.csv'
SIMPLIFIED_VARIANTS = '/Users/liweihang/BT1051/BLAST_Algorithm/Database/simplified_variants.csv'
HASH_DIR = '/Users/liweihang/BT1051/BLAST_Algorithm/Database/hash_chunks'
DEDUP_DIR = '/Users/liweihang/BT1051/BLAST_Algorithm/Database/Final_hash_chunks_Table'
# ------------------------------------------------

# ---------------- Part1_Hash & Variant Filter ----------------
STANDARD_AA = set(protein_letters)
EXTENDED_AA_MAP = {'B':'DN','Z':'EQ','J':'IL','X':'','U':'C','O':'K','*':''}
VALID_AA = STANDARD_AA.union(set(EXTENDED_AA_MAP.keys())) - {'X','*'}

class KmerVariantFilter:
    def __init__(self, fasta_path, variant_csv_path, k=3, T_threshold=10):
        self.fasta_path = fasta_path
        self.variant_csv_path = variant_csv_path
        self.k = k
        self.T_threshold = T_threshold
        self.query_sequence = ''
        self.query_kmers = []

    def read_protein_fasta(self):
        rec = next(SeqIO.parse(self.fasta_path, 'fasta'))
        seq = str(rec.seq).upper()
        invalid = set(seq) - VALID_AA
        if invalid:
            print(f"⚠️ 非法字符 {invalid} 被跳过")
        self.query_sequence = seq

    def generate_kmers(self):
        seq = self.query_sequence
        self.query_kmers = [seq[i:i+self.k] for i in range(len(seq)-self.k+1)
                             if all(aa in VALID_AA for aa in seq[i:i+self.k])]

    def filter_variants(self, out_csv):
        qset = set(self.query_kmers)
        count = 0
        with open(self.variant_csv_path) as inf, open(out_csv, 'w', newline='') as outf:
            reader = csv.DictReader(inf)
            writer = csv.writer(outf)
            writer.writerow(['Original_3mer','Variant','BLOSUM62_Score'])
            for row in reader:
                if row['Original_3mer'] in qset and float(row['BLOSUM62_Score'])>=self.T_threshold:
                    writer.writerow([row['Original_3mer'],row['Variant'],row['BLOSUM62_Score']])
                    count+=1
        print(f"✅ 筛选保留 {count} 条变体")

def simplify_variant_csv(inf, outf):
    merged = defaultdict(set)
    with open(inf) as f:
        for row in csv.DictReader(f): merged[row['Original_3mer']].add(row['Variant'])
    with open(outf,'w',newline='') as f:
        writer = csv.writer(f); writer.writerow(['Original_3mer','Variants'])
        for o,vs in merged.items(): writer.writerow([o,';'.join(sorted(vs))])
    print('✅ 变体简化完成')

# ---------------- Hash Database ----------------
class ProteinHashDatabase:
    def __init__(self, k=3, variants_file=None, output_dir='hash_chunks'):
        self.k = k
        self.variant_map = {}
        if variants_file:
            for row in csv.DictReader(open(variants_file)):
                for v in row['Variants'].split(';'):
                    self.variant_map[v] = row['Original_3mer']
        self.out_dir = output_dir
        self.chunk_idx = 0
        self.chunk_size = 30000
        os.makedirs(self.out_dir, exist_ok=True)

    def _save(self, table):
        f1 = os.path.join(self.out_dir, f'hash_chunk_{self.chunk_idx}.txt')
        f2 = os.path.join(self.out_dir, f'hash_chunk_{self.chunk_idx}_withVariants.txt')
        with open(f1,'w') as o1, open(f2,'w') as o2:
            for kmer, coords in sorted(table.items()):
                pos=';'.join(f"{i},{p}" for i,p in coords)
                o1.write(f"{kmer}\t{pos}\n")
                o2.write(f"{kmer}\t{pos}\n")
                if kmer in self.variant_map:
                    o2.write(f"{self.variant_map[kmer]}\t{pos}\n")
        print(f"💾 保存分块 {self.chunk_idx}")
        self.chunk_idx+=1

    def build(self, fasta):
        recs=list(SeqIO.parse(fasta,'fasta'))
        for start in range(0,len(recs),self.chunk_size):
            table=defaultdict(list)
            for idx,rec in enumerate(recs[start:start+self.chunk_size], start):
                seq=str(rec.seq).upper()
                for p in range(len(seq)-self.k+1): table[seq[p:p+self.k]].append((idx,p))
            self._save(table)
        print('✅ 哈希数据库构建完成')

# ---------------- Utility ----------------
def load_index(input_dir):
    idx={}
    for fn in os.listdir(input_dir):
        if fn.endswith('_withVariants.txt'):
            for line in open(os.path.join(input_dir,fn)):
                k,p=line.strip().split('\t')
                idx.setdefault(k,[]).extend(p.split(';'))
    return idx

def dedup(input_file,out_dir):
    merged=defaultdict(set)
    for ln in open(input_file):
        k,ps=ln.strip().split('\t')
        merged[k].update(ps.split(';'))
    with open(os.path.join(out_dir,os.path.basename(input_file)),'w') as o:
        for k in sorted(merged): o.write(f"{k}\t{';'.join(sorted(merged[k]))}\n")

def dedup_all(in_dir,out_dir):
    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    for fn in glob(os.path.join(in_dir,'*_withVariants.txt')): dedup(fn,out_dir)
    print(f"✅ 去重完成, 输出在 {out_dir}")

# ---------------- Main ----------------
def main():
    if os.path.exists(HASH_DIR): shutil.rmtree(HASH_DIR)
    # 1. 可选变体过滤
    if not FULL_DB:
        print('=== STEP1: 过滤 & 简化变体 ===')
        kvf=KmerVariantFilter(QUERY_FASTA, VARIANT_CSV, k, T_THRESHOLD)
        kvf.read_protein_fasta(); kvf.generate_kmers()
        kvf.filter_variants(FILTERED_VARIANTS)
        simplify_variant_csv(FILTERED_VARIANTS, SIMPLIFIED_VARIANTS)
        var_file = SIMPLIFIED_VARIANTS
    else:
        var_file = None
    # 2. 构建哈希
    print('=== STEP2: 构建哈希分块 ===')
    ph=ProteinHashDatabase(k, variants_file=var_file, output_dir=HASH_DIR)
    ph.build(TARGET_FASTA)
    # 3. 测试查询
    print('=== STEP3: 测试查询 ===')
    idx=load_index(HASH_DIR)
    for km in ['AGH','XBB']: print(f"{km}: {len(idx.get(km,[]))} 个位置")
    # 4. 去重
    print('=== STEP4: 去重 withVariants ===')
    dedup_all(HASH_DIR, DEDUP_DIR)

if __name__=='__main__': main()
