# BLAST_Alogorithm

## Part 1 Hash Database Build

首先先通过generate_kmer_variants.py，快速生成自定义的k中所有得分大于T的kmer，将文件（CSV format）保存于本地。

Note: k的提高会使得数据库的大小和计算时间呈指数级上增。选择K=3, T=10(T可以设置的稍微小点，以备不时之需)

这个部分通过两个py文件实现

1. Part1_Hash.py（存放了两个class和一个函数）

   - class **KmerVariantFilter**: 读取query序列，打成k-mer,去根据query的k-mer过滤generate_kmer_variants.py生成的csv文件,只保留那些存在query中的完全匹配的kmer以及其变体，减少后续的工作量，避免不必要的哈希库添加。

   - **simplify_variant_csv**：将上一步生成的csv文件，进一步整理，简化变体CSV文件，合并相同原始k-mer的变体：

     输出如下：

     ```
     Original_3mer,Variants
     AAM,AAM
     AAQ,AAQ
     ADM,ADI;ADL;ADM;ADV;AEM;CDM;GDM;SDM;TDM;VDM
     AEG,ADG;AEG;AKG;AQG;CEG;GEG;SEG;TEG;VEG
     ...
     ```

   - **ProteinHashDatabase**: Hash化数据库，输入是一个存放多条序列的fasta文件。针对每一条序列，都生成可能的匹配，最终是输出哈希化的数据库，会以TXT的形式保存。

     ```
     MKK	0,0;1,278;2,278
     KKA	0,1;0,1;1,114;1,114;2,114;2,114;3,4;3,4
     KKS	0,1;1,114;2,114;3,4
     KRA	0,1;1,114;2,114;3,4
     RKA	0,1;1,114;2,114;3,4
     KAT	0,2;3,5
     ATV	0,3;0,3;1,105;1,105;2,105;2,105
     ...
     ```

2. runProteinHashDatabase.py

运行上述程序，需要输入query，database，设置（filtered_variants.csv，simplified_variants.csv，complete_hash_table.txt）的输出的路径。

大概结果如下：

```bash
(base) liweihang@Mac-3 BLAST_Algorithm % /Users/liweihang/miniconda3/bin/python /Users/liweihang/BT1051/BLAST_Algorithm/Part1_HASH/run
ProteinHashDatabase.py
=== 步骤1: 基于查询序列过滤变体 ===
✅ 筛选完成，保留了query中 6258 个K-mer（含有Variant）。输出文件：/Users/liweihang/BT1051/BLAST_Algorithm/Database/filtered_variants.csv

=== 步骤2: 简化变体文件 ===
✅ filtered_variants.csv简化完成！输出路径：/Users/liweihang/BT1051/BLAST_Algorithm/Database/simplified_variants.csv

=== 步骤3: 构建完整哈希数据库 ===
已加载 320 个k-mer的变体信息
哈希数据库构建完成，包含 134 条序列和 5295 个唯一k-mer
哈希表已保存到 /Users/liweihang/BT1051/BLAST_Algorithm/Database/complete_hash_table.txt

=== 总运行时间: 0.75秒 ===

=== 数据库统计信息 ===
ProteinHashDatabase(k=3) 包含 134 条序列, 5295 个唯一k-mer

示例查询结果:
查询 'AAM': [(0, 131), (1, 133), (2, 133), (3, 133), (4, 133), (5, 133), (6, 133), (7, 133), (8, 133), (34, 133), (35, 133), (36, 133), (41, 134), (42, 133), (96, 136), (102, 343), (105, 131), (109, 134), (111, 343), (121, 131), (122, 131), (123, 131), (127, 133), (128, 383), (129, 143)]
查询 'VHQ': [(1, 307), (2, 307), (3, 307), (4, 307), (5, 307), (6, 307), (7, 307), (8, 307), (34, 307), (35, 307), (36, 307), (42, 307), (85, 189), (108, 106), (127, 307)]
查询 'AHQ': [(1, 112), (1, 307), (2, 112), (2, 307), (3, 112), (3, 307), (4, 112), (4, 307), (5, 112), (5, 307), (6, 112), (6, 307), (7, 112), (7, 307), (8, 112), (8, 307), (34, 112), (34, 307), (35, 112), (35, 307), (36, 112), (36, 307), (42, 112), (42, 307), (85, 189), (108, 106), (124, 117), (125, 117), (127, 112), (127, 307), (128, 126), (131, 111)]
```

​	

## Part 2 

blastp -query 8HHU.fa  -db uniprot_sprot  -matrix BLOSUM62 -gapopen 11 -gapextend 1 -word_size 4 -outfmt 6 -threshold 11        

blastp -query GZMK.fasta  -db uniprot_sprot  -matrix BLOSUM62 -gapopen 11 -gapextend 1 -word_size 4 -outfmt 6 -threshold 11 -out gzmk.txt -max_target_seqs 2000 -evalue 5

blastp -query GZMK.fasta  -db uniprot_sprot  -matrix BLOSUM62 -gapopen 11 -gapextend 1 -word_size 4 -outfmt 6 -threshold 11 -out gzmk.txt -max_target_seqs 2000 -evalue 5

(BT1051) liweihang@Mac-3 TestFile % blastp -task blastp-short \                                                                                
       -query 8HHU.fa \
       -db uniprot_sprot \
       -matrix BLOSUM62 \
       -gapopen 11 \
       -gapextend 1 \
       2>&1 | grep -A 3 "Lambda"
Lambda      K        H        a         alpha
   0.323    0.138    0.433    0.792     4.96 

Gapped
Lambda      K        H        a         alpha    sigma
   0.267   0.0410    0.140     1.90     42.6     43.6 

Effective search space used: 26336991252
(BT1051) liweihang@Mac-3 TestFile % 
