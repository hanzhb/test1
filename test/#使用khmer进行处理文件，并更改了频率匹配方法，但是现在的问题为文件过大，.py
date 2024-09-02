#使用khmer进行处理文件，并更改了频率匹配方法，但是现在的问题为文件过大，
# 使用khmer.ReadParser(input_fasta)函数无法处理，会报错，保存作为中间版本，
import csv
import hashlib
import math
import random

import khmer
import numpy as np
from datasketch import MinHash
import pandas as pd
from scipy.stats import pearsonr


# 第一步：使用khmer计算k-mer频率
def count_kmers_with_khmer(input_fasta, k, mincount):

    ht = khmer.Countgraph(k, 1e7, 4)
    def is_valid_kmer(kmer):
        # 只允许标准碱基（A、T、C、G）组成的k-mer
        return all(base in 'ATCG' for base in kmer)

    # 过滤掉出现次数少于两次的k-mers
    def filter_kmers(kmer_counts, min_count):
        filtered_kmers = {kmer: count for kmer, count in kmer_counts.items() if
                          count > min_count }#and is_valid_kmer(kmer)}
        return filtered_kmers

    # 使用ReadParser读取FASTA文件
    parser = khmer.ReadParser(input_fasta)
    a = parser.num_reads
    print(parser,a)
    kmer_counts = {}
    sequences = ""
    # 对每条序列进行处理
    for record in parser:
        sequence = record.sequence.upper()
        sequences += sequence
        #print(sequences,type(sequences))
    ht.consume(sequences)#consume方法将整个序列添加到Countgraph中，并对序列中的每个k-mer进行计数
        #print(sequence)
    c =ht.get_kmers(sequences)
    #print(c)
    for kmer in ht.get_kmers(sequences):
            count = ht.get(kmer)
            if is_valid_kmer(kmer):  # 检查kmer是否有效
                kmer_counts[kmer] = count
    kmer_counts = filter_kmers(kmer_counts,mincount)
    #print(type(kmer_counts))
    return kmer_counts


# 第二步：生成MinHash签名
def compute_minhash(kmer_counts, num_perm=1000):
    m = MinHash(num_perm=num_perm)
    for kmer, count in kmer_counts.items():
        #for _ in range(count):
            m.update(kmer.encode('utf-8'))
            #print("c",m.digest())
    return m

def jaccard_similarity(minhash1, minhash2,kmer_size):
    intersection_count = np.sum(minhash1.hashvalues == minhash2.hashvalues)
    total_count = len(minhash1.hashvalues)
    similarity1 = intersection_count / total_count
    #similarity2 = minhash1.jaccard(minhash2)
    print("pipei:",intersection_count,"/",total_count)
    print("jaccard:",similarity1)
    if intersection_count == 0:
        jaccard = 0
    else:
        jaccard = intersection_count / total_count
    mash_distance = -1 / kmer_size * math.log(2 * jaccard / (1 + jaccard))
    print("mashdistance",mash_distance)
    return similarity1


# 第三步：找到MinHash签名对应的k-mer频率---pass
def find_kmer_frequencies(minhash, kmer_counts):
    kmer_frequencies = {}
    for kmer in kmer_counts:
        if minhash.is_member(kmer.encode('utf-8')):  # 判断kmer是否在MinHash中
            kmer_counts.append(kmer)
            kmer_frequencies.append(kmer_counts[kmer])
    print("kmer_frequence\n",kmer_frequencies)
    return kmer_frequencies

#得到频率列表
def count_list3(list1,list2):
    list1 = list(list1.items())
    list2 = list(list2.items())
    #print(type(list1),list1)
    filename = 'data.csv'
    #df1 = pd.DataFrame(list1,list2)
    #df1.to_csv(filename, index=False, mode='a', header='list1')

    # 使用这种方法集合 保证每次抽取的值都是相同的
    def stable_hash(item):
        # 将项目转换为字符串，并计算其稳定的哈希值
        item_str = str(item).encode('utf-8')
        return int(hashlib.md5(item_str).hexdigest(), 16)

    def hash_sample(lst, n):

        # 对列表中的每个元素计算稳定的哈希值，并将哈希值与元素对进行排序
        hashed_items = sorted((stable_hash(item), item) for item in lst)
        # 选择排序后的前 n 个元素
        sampled_items = [item for _, item in hashed_items[:n]]
        return sampled_items
    list1_hash = hash_sample(list1,min(500, len(list1)))
    #list2_hash = hash_sample(list2,500)
    random.seed(42)
    list1 = random.sample(list1, min(500, len(list1)))
    #list2 = random.sample(list2, 1000)
    #print(list1,"\nlist2",list2)

    # 将两个列表中的 k-mer 提取并合并成一个集合，去重.  使用字典集,保证排序一致性
    #all_kmers = sorted(set(kmer for kmer, _ in list1).union(kmer for kmer, _ in list2))
    all_kmers = sorted(set(kmer for kmer, _ in list1_hash))
    #print("\nall",all_kmers)
    hebing = pd.DataFrame(all_kmers)
    hebing.to_csv(filename, index=False, mode='a', header="hebing")
   # 创建两个空的频率列表
    freq_list1 = []
    freq_list2 = []

    # 构建文档一的频率列表
    kmer_dict1 = dict(list1_hash)
    print("k1\n",kmer_dict1)

    for kmer in all_kmers:
        freq_list1.append(kmer_dict1.get(kmer, 0))

    # 构建文档二的频率列表
    kmer_dict2 = dict(list2)
    for kmer in all_kmers:
        freq_list2.append(kmer_dict2.get(kmer, 0))
    print("k2\n",kmer_dict2,"list1\n",freq_list1,"list2\n",freq_list2)
    # randoml1 = pd.DataFrame(list1)
    # randoml2 = pd.DataFrame(list2)
    # randoml1.to_csv(filename,index=False,mode='a',header="randoml1")
    # randoml2.to_csv(filename, index=False, mode='a', header="randoml2")
    return freq_list1, freq_list2

def cosine_similarity(vector1, vector2):
    # 计算向量的点积
    dot_product = np.dot(vector1, vector2)
    # 计算向量的范数（即向量的长度）
    norm_vector1 = np.linalg.norm(vector1)
    norm_vector2 = np.linalg.norm(vector2)
    # 计算余弦相似性
    cosine_sim = dot_product / (norm_vector1 * norm_vector2)
    print('余弦相似度:',cosine_sim)
    return cosine_sim

#欧几里得距离
def euclidean_distance(vec1, vec2):
    a = np.array(vec1)
    b = np.array(vec2)
    c = np.linalg.norm(a-b)
    print("欧几里得距离:",c)
    return c

#pearson corrlation
def pearson(freq_list1,freq_list2):
    correlation, p_value = pearsonr(freq_list1, freq_list2)
    #corr,p = stats.spearmanr(freq_list1,freq_list2)
    # 输出结果
    print("Pearson 相关系数:", correlation)
    print("p-value:", p_value)
    #print("spearman 相关系数:", correlation)
    #print("p-value:", p_value)
    return correlation

def save_to_file(kmer_counts, output_file):
    # 将k-mer频率保存到文件
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        # 写入数据行
        for row in kmer_counts:
            writer.writerow(row)


# 示例代码
input_fasta1 = '/home/hzb/dna/YC_con_3_BDDP210000934-1A_1.clean.fq'
input_fasta2 = '/home/hzb/dna/YC_con_1_BDDP210000932-1A_1.clean.fq'
output_file = '/home/hzb/dna/test1.csv'
kmer_size = 16
min_count = 5
kmer_counts1 = count_kmers_with_khmer(input_fasta1,kmer_size,min_count)
kmer_counts2 = count_kmers_with_khmer(input_fasta2,kmer_size,min_count)
print("kmercount2",kmer_counts2)
minhash1 = compute_minhash(kmer_counts1)
minhash2 = compute_minhash(kmer_counts2)
jaccard = jaccard_similarity(minhash1,minhash2,kmer_size)
list1,list2 = count_list3(kmer_counts1,kmer_counts2)
print("list1",list1,'list2',list2)
cos = cosine_similarity(list1,list2)
pearson = pearson(list1,list2)
ou = euclidean_distance(list1,list2)
#kmer_frequencies = find_kmer_frequencies(minhash1, kmer_counts1)
#print("MinHash对应的k-mer频率:")
#for kmer, freq in kmer_frequencies.items():
    #print(f"{kmer}: {freq}")

