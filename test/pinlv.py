import csv
import hashlib
import re
import subprocess
import numpy as np
from datasketch import MinHash
import os
import scipy.stats as stats
import random
#测试各种计算方法所得到的相似性结果
#该代码通过jellyfish包计算每一个kmer的出现次数，然后生成每个kmer与频率的list列表，根据出现频率添加到minhash草图中，计算jaccard距离
# 定义一个函数运行 Jellyfish count
def run_jellyfish_count(input_fastq, output_file, kmer_size):
    run_jellyfish_count_cmd = [
        'jellyfish', 'count',
        '-m', str(kmer_size),
        '-s', '1G',
        '-C',  # count canonical k-mers (to consider reverse complements)
        '-o', output_file,
        input_fastq
    ]
    subprocess.run(run_jellyfish_count_cmd, check=True)
    print(f"Jellyfish count completed: {run_jellyfish_count_cmd}")

# 定义一个函数运行 Jellyfish query
def run_jellyfish_query(output_file):
    jellyfish_query_cmd = [
        'jellyfish', 'dump',
        output_file
    ]
    result = subprocess.run(jellyfish_query_cmd, stdout=subprocess.PIPE, check=True)
    return result.stdout.decode('utf-8')

# 从 Jellyfish 输出中加载 k-mer 计数
def load_kmer_counts_from_jf(output_file):
    kmer_counts = []
    output = run_jellyfish_query(output_file)
    lines = output.splitlines()
    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):  # Ensure there are two lines to process
            kmer_line = lines[i].strip().split()
            count_line = lines[i + 1].strip().split()
            if kmer_line and count_line:
                count_str = ' '.join(map(str, kmer_line))
                # 使用正则表达式匹配数字
                kmer_line= re.search(r'\d+', count_str)
                #print(kmer_line)
                kmer = int(kmer_line.group())
                count = count_line[0] # Assuming count is in the first part

                kmer_counts.append(( count,kmer))
                #print(kmer,"--@@",count)
    # for line in output.splitlines():
    #     parts = line.strip().split()
    #     print(parts,type(parts),len(parts))
    #     #if len(parts) == 1:
    #     kmer = parts[0]
    #     count = parts[0]
    #     print(kmer,"~@@@",count)
    #     kmer_counts.append((kmer, count))
    return kmer_counts


#传入两个kmer计数，进行合并输出kmer频率列表
def count_list1(list1,list2):
    #list1 = [('GCATGCAGAAGGGAAA', 1), ('GCAGTGCCAGGTAGCC', 2), ('ATCAGCGGACAAAGTG', 3), ('GGGCGCTTAGCATGAC', 4), ('ACTTAGGTCGGGGTCA', 5), ('GGGCCTGCTCGAATGA', 6), ('ATAGTTTAAAAGCACC', 7)]
    #list2 = [('GCATGCAGAAGGGAAAT', 1), ('GCAGTGCCAGGTAGCCC', 2), ('ATCAGCGGACAAAGTGA', 3), ('GGGCGCTTAGCATGAC', 4), ('ACTTAGGTCGGGGTCA', 5), ('GGGCCTGCTCGAATGA', 6), ('ATAGTTTAAAAGCACC', 7)]

    # 将两个列表中的 k-mer 提取并合并成一个集合，去重
    all_kmers = set(kmer for kmer, _ in list1).union(kmer for kmer, _ in list2)

    # 创建两个空的频率列表
    freq_list1 = []
    freq_list2 = []

    # 构建文档一的频率列表
    kmer_dict1 = dict(list1)
    for kmer in all_kmers:
        freq_list1.append(kmer_dict1.get(kmer, 0))

    # 构建文档二的频率列表
    kmer_dict2 = dict(list2)
    for kmer in all_kmers:
        freq_list2.append(kmer_dict2.get(kmer, 0))

    return freq_list1, freq_list2

def count_list2(list1,list2):
    # 将两个列表中的 k-mer 提取并合并成一个集合，去重.  使用字典集,保证排序一致性
    all_kmers = sorted(set(kmer for kmer, _ in list1).union(kmer for kmer, _ in list2))
    #all_kmers = random.sample(all_kmers,1000)

    #使用这种方法集合 保证每次抽取的值都是相同的
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
    all_kmers = hash_sample(all_kmers,1000)

    # 创建两个空的频率列表
    freq_list1 = []
    freq_list2 = []

    # 构建文档一的频率列表
    kmer_dict1 = dict(list1)
    for kmer in all_kmers:
        freq_list1.append(kmer_dict1.get(kmer, 0))

    # 构建文档二的频率列表
    kmer_dict2 = dict(list2)
    for kmer in all_kmers:
        freq_list2.append(kmer_dict2.get(kmer, 0))

    return freq_list1, freq_list2


def count_list3(list1,list2):
    #print(list2)
    # 使用这种方法集合 保证每次抽取的值都是相同的
    def stable_hash(item):
        # 将项目转换为字符串，并计算其稳定的哈希值
        item_str = str(item).encode('utf-8')
        return int(hashlib.md5(item_str).hexdigest(), 16)

    def hash_sample(lst, n):
        # 过滤出频率值大于 2 的元素
        filtered_items = [item for item in lst if item[1] > 10]
        # 对列表中的每个元素计算稳定的哈希值，并将哈希值与元素对进行排序
        hashed_items = sorted((stable_hash(item), item) for item in filtered_items)
        # 选择排序后的前 n 个元素
        sampled_items = [item for _, item in hashed_items[:n]]
        return sampled_items

    list1 = hash_sample(list1,500)
    list2 = hash_sample(list2,500)
    #print(list1,"\nlist2",list2)
    # 将两个列表中的 k-mer 提取并合并成一个集合，去重.  使用字典集,保证排序一致性
    all_kmers = sorted(set(kmer for kmer, _ in list1).union(kmer for kmer, _ in list2))
    #print("\nall",all_kmers)
   # 创建两个空的频率列表
    freq_list1 = []
    freq_list2 = []

    # 构建文档一的频率列表
    kmer_dict1 = dict(list1)
    print("k1\n",kmer_dict1)
    for kmer in all_kmers:
        freq_list1.append(kmer_dict1.get(kmer, 0))

    # 构建文档二的频率列表
    kmer_dict2 = dict(list2)
    for kmer in all_kmers:
        freq_list2.append(kmer_dict2.get(kmer, 0))
    print("k2\n",kmer_dict2,"list1\n",freq_list1,"list2\n",freq_list2)
    return freq_list1, freq_list2

def count_list4(list1,list2):
    def hash_sample(lst, n):
        # 过滤出频率值大于 2 的元素
        filtered_items = [item for item in lst if item[1] > 2]
        all_kmers = random.sample(filtered_items, n)
        return all_kmers
    list1 = hash_sample(list1,500)
    list2 = hash_sample(list2,500)
    #print(list1,"\nlist2",list2)
    # 将两个列表中的 k-mer 提取并合并成一个集合，去重.  使用字典集,保证排序一致性
    all_kmers = sorted(set(kmer for kmer, _ in list1).union(kmer for kmer, _ in list2))
    #print("\nall",all_kmers)
   # 创建两个空的频率列表
    freq_list1 = []
    freq_list2 = []

    # 构建文档一的频率列表
    kmer_dict1 = dict(list1)
    print("k1\n",kmer_dict1)
    for kmer in all_kmers:
        freq_list1.append(kmer_dict1.get(kmer, 0))

    # 构建文档二的频率列表
    kmer_dict2 = dict(list2)
    for kmer in all_kmers:
        freq_list2.append(kmer_dict2.get(kmer, 0))
    print("k2\n",kmer_dict2,"list1\n",freq_list1,"list2\n",freq_list2)
    return freq_list1, freq_list2


# 计算 MinHash
#将根据count计数将每一个kmer上传到草图中
def compute_minhash(kmer_counts, num_perm):
    #print(kmer_counts)
    m = MinHash(num_perm=num_perm)
    for kmer, count in kmer_counts:
        for _ in range(count):
            #print(kmer)
            m.update(kmer.encode('utf-8'))
    return m

# 计算两个 MinHash 之间的 Jaccard 距离
def compute_jaccard_distance(minhash1, minhash2):
    return minhash1.jaccard(minhash2)

#pearson corrlation
def pearson(freq_list1,freq_list2):
    correlation, p_value = stats.pearsonr(freq_list1, freq_list2)
    corr,p = stats.spearmanr(freq_list1,freq_list2)
    # 输出结果
    print("Pearson 相关系数:", correlation)
    print("p-value:", p_value)
    #print("spearman 相关系数:", correlation)
    #print("p-value:", p_value)
    return correlation
#余弦相似度
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


# 主函数
def main(input_fastq1,input_fastq2,kmer_size):#
    # 文件路径和参数
    input_fastq1 = '/home/hzb/dna/test1.fq'
    input_fastq2 = '/home/hzb/dna/test5.fq'
    kmer_size = 12
    num_perm = 12
    output_file1 = 'output1.jf'
    output_file2 = 'output2.jf'

    # 运行 Jellyfish count,对DNA序列进行kmer分割
    run_jellyfish_count(input_fastq1, output_file1, kmer_size)
    run_jellyfish_count(input_fastq2, output_file2, kmer_size)
    #print("--------")
    # 从 Jellyfish 文件中加载 k-mer 计数
    kmer_counts1 = load_kmer_counts_from_jf(output_file1)
    kmer_counts2 = load_kmer_counts_from_jf(output_file2)
    #print("~~~~~~\n",kmer_counts1[:100],"~~~~~~\n",kmer_counts2[:100])

    #计算得到两个头部序列的频率list
    freq_list1,freq_list2 = count_list4(kmer_counts1,kmer_counts2)
    freq_list1 = freq_list1
    freq_list2 = freq_list2
    #print(freq_list1)
    #print(freq_list2)
    #pearson
    corr = pearson(freq_list1,freq_list2)

    # 计算余弦相似度
    cos = cosine_similarity(freq_list1,freq_list2)
    # 计算欧几里得距离
    ed = euclidean_distance(freq_list1,freq_list2)
    return corr,cos,ed

def main1():
    result = []
    kmersize = 21
    base_string = '/mnt/00.zhibo/6MnPCR/chuan/YC_con_2_BDDP210000933-1A_1.clean.fq'#dizhi27
    for j in range(1):
        new_string1 = base_string[:34] + str(j + 1) + base_string[35:47] + str(j + 32) + base_string[49:]

        for i in range(j+1,2):
            a = i + 48
            new_string2 = base_string[:34] + str(i + 1) + base_string[35:47] + str(i + 32) + base_string[49:]
            #print(new_string1,new_string2)
            #output = [new_string1,new_string2]
            output = main(new_string1,new_string2,kmersize)
            result.append(output)
    output_file = "output.csv"
    # 将 result 列表中的数据写入文件
    with open(output_file, "w",newline="") as file:
        writer = csv.writer(file)
        writer.writerows(result)
    print(result)
    # 计算 MinHash
    #minhash1 = compute_minhash(kmer_counts1, num_perm)
    #minhash2 = compute_minhash(kmer_counts2, num_perm)
    #print(minhash1,"~~~~~~~~~~~~\n",minhash2)
    # 计算 Jaccard 距离
    #jaccard_distance = compute_jaccard_distance(minhash1, minhash2)
    #print(f"Jaccard Distance: {jaccard_distance}")

def main2():
    result = []
    kmersize = 21
    base_string = '/mnt/00.zhibo/6MnPCR/chuan/YC_con_2_BDDP210000933-1A_1.clean.fq'#dizhi27
    for j in range(2):
        new_string1 = base_string[:34] + str(j + 1) + base_string[35:47] + str(j + 32) + base_string[49:]


        #a = i + 48
        new_string2 = base_string[:34] + str(2) + base_string[35:47] + str(33) + base_string[49:]
        new_string3 = base_string[:34] + str(3) + base_string[35:47] + str(34) + base_string[49:]
        new_string4 = base_string[:34] + str(4) + base_string[35:47] + str(35) + base_string[49:]
        new_string5 = base_string[:34] + str(5) + base_string[35:47] + str(36) + base_string[49:]
        new_string6 = base_string[:34] + str(6) + base_string[35:47] + str(37) + base_string[49:]
        #print(new_string1,new_string2,new_string3,new_string4,new_string5,new_string6)
            #output = [new_string1,new_string2]
        output = main(new_string1,new_string2,kmersize)
        result.append(output)
        output = main(new_string1, new_string3, kmersize)
        result.append(output)
        output = main(new_string1, new_string4, kmersize)
        result.append(output)
        output = main(new_string1, new_string5, kmersize)
        result.append(output)
        output = main(new_string1, new_string6, kmersize)
        result.append(output)
    output_file = "output.csv"
    # 将 result 列表中的数据写入文件
    with open(output_file, "w",newline="") as file:
        writer = csv.writer(file)
        writer.writerows(result)
    print(result)
    # 计算 MinHash
    #minhash1 = compute_minhash(kmer_counts1, num_perm)
    #minhash2 = compute_minhash(kmer_counts2, num_perm)
    #print(minhash1,"~~~~~~~~~~~~\n",minhash2)
    # 计算 Jaccard 距离
    #jaccard_distance = compute_jaccard_distance(minhash1, minhash2)
    #print(f"Jaccard Distance: {jaccard_distance}")
if __name__ == "__main__":
    main1()
