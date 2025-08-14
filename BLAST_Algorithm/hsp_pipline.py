#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hsp_pipeline.py

How to use:
  ./hsp_pipeline.py -i /path/to/HSP -o /path/to/Results

/Users/liweihang/miniconda3/bin/python /Users/liweihang/BT1051/BLAST_Algorithm/Part3_Output/hsp_pipline.py -i HSP -o Results
"""

import os
import sys
import json
import argparse
import shutil
import glob
from collections import defaultdict

import pandas as pd

def is_overlap(r1, r2):
    return not (r1[1] < r2[0] or r2[1] < r1[0])

def process_one(input_json, output_dir):
    base = os.path.basename(input_json).rsplit('.', 1)[0]
    # 1) 读取 JSON
    try:
        with open(input_json) as f:
            data = json.load(f)
    except Exception as e:
        print(f"❌ {base}: JSON 解析失败：{e}", file=sys.stderr)
        return 0

    if not data:
        print(f"⚠️ {base}: 无 HSP 数据，跳过", file=sys.stderr)
        return 0

    # 2) 规范化 e_value 字段
    sample_keys = set(data[0].keys())
    if 'e_value' not in sample_keys:
        if 'evalue' in sample_keys:
            for h in data:
                h['e_value'] = h.pop('evalue')
        elif 'ev' in sample_keys:
            for h in data:
                h['e_value'] = h.pop('ev')
        else:
            print(f"❌ {base}: 找不到 e_value 字段，现有字段：{list(sample_keys)}", file=sys.stderr)
            return 0

    # 3) 排序
    sorted_data = sorted(data, key=lambda x: x['e_value'])
    out_sorted = os.path.join(output_dir, f"{base}_sorted.json")
    with open(out_sorted, 'w') as fo:
        json.dump(sorted_data, fo, indent=2)

    # 4) 按 db_id 分组 & 去除 query_range 重叠
    by_db = defaultdict(list)
    for h in sorted_data:
        by_db[h['db_id']].append(h)

    full_results = []
    for db_id, hits in by_db.items():
        sel = []
        for h in hits:
            if any(is_overlap(h['query_range'], x['query_range']) for x in sel):
                continue
            sel.append(h)
        full_results.extend(sel)

    # 5) 写出完整 CSV
    df_full = pd.DataFrame(full_results)
    out_full = os.path.join(output_dir, f"{base}_full.csv")
    df_full.to_csv(out_full, index=False)

    # 6) 写出每个 db_id 的最佳一条
    df_top1 = (
        df_full
        .sort_values('e_value')
        .groupby('db_id', as_index=False)
        .first()
        .sort_values('e_value')
    )
    out_top1 = os.path.join(output_dir, f"{base}_top1.csv")
    df_top1.to_csv(out_top1, index=False)

    print(f"📄 {base}: 原始条数={len(data)}, 去重后={len(full_results)}, top1={len(df_top1)}")
    return len(full_results)

def main():
    p = argparse.ArgumentParser(description="一键式 HSP 处理管道")
    p.add_argument('-i','--input-dir', required=True,
                   help="HSP JSON 文件夹路径")
    p.add_argument('-o','--output-dir', required=True,
                   help="结果输出目录")
    args = p.parse_args()

    # 准备输出目录
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    # 找 JSON
    json_files = sorted(
        os.path.join(args.input_dir, f)
        for f in os.listdir(args.input_dir)
        if f.endswith('.json')
    )
    if not json_files:
        print("⚠️ 未找到任何 .json 文件，请检查输入目录。", file=sys.stderr)
        sys.exit(1)

    total = 0
    for fn in json_files:
        total += process_one(fn, args.output_dir)

    print(f"\n✅ 全部完成，总保留 HSP 条数: {total}")

    # 7) 合并所有 top1 并按 e_value 排序
    top1_pattern = os.path.join(args.output_dir, '*_top1.csv')
    top1_files = sorted(glob.glob(top1_pattern))
    if top1_files:
        df_list = [pd.read_csv(f) for f in top1_files]
        df_all_top1 = pd.concat(df_list, ignore_index=True)
        df_all_top1 = df_all_top1.sort_values('e_value')
        out_all = os.path.join(args.output_dir, 'all_top1.csv')
        df_all_top1.to_csv(out_all, index=False)
        print(f"📑 合并所有 top1，输出: {out_all}")
    else:
        print("⚠️ 未找到任何 top1 文件以合并。")

if __name__ == '__main__':
    main()
