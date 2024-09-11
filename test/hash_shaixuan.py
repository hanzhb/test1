#oringin test 使用mur32 得到最小hash值和频率。问题在于双层for循环
import khmer

import time

def murmur3_32(key: bytes, seed: int = 0) -> int:
    length = len(key)
    nblocks = length // 4
    h1 = seed
    c1 = 0xcc9e2d51
    c2 = 0x1b873593
    r1 = 15
    r2 = 13
    m = 5
    n = 0xe6546b64

    # Process each 4-byte block of the key
    for i in range(nblocks):
        k1 = int.from_bytes(key[i * 4:(i + 1) * 4], byteorder='little')
        k1 = (k1 * c1) & 0xffffffff
        k1 = (k1 << r1) | (k1 >> (32 - r1))
        k1 = (k1 * c2) & 0xffffffff
        h1 ^= k1
        h1 = (h1 << r2) | (h1 >> (32 - r2))
        h1 = (h1 * m + n) & 0xffffffff

    # Handle remaining bytes
    tail = key[nblocks * 4:]
    k1 = 0
    if len(tail) >= 3:
        k1 ^= tail[2] << 16
    if len(tail) >= 2:
        k1 ^= tail[1] << 8
    if len(tail) >= 1:
        k1 ^= tail[0]
    if len(tail) > 0:
        k1 = (k1 * c1) & 0xffffffff
        k1 = (k1 << r1) | (k1 >> (32 - r1))
        k1 = (k1 * c2) & 0xffffffff
        h1 ^= k1

    # Finalization
    h1 ^= length
    h1 = (h1 ^ (h1 >> 16)) & 0xffffffff
    h1 = (h1 * 0x85ebca6b) & 0xffffffff
    h1 = (h1 ^ (h1 >> 13)) & 0xffffffff
    h1 = (h1 * 0xc2b2ae35) & 0xffffffff
    h1 = (h1 ^ (h1 >> 16)) & 0xffffffff

    return h1
class ListDict:
    def __init__(self, max_size=10):
        self.max_size = max_size
        self.data = []

    def add(self, value):
        if len(self.data) < self.max_size:
            self.data.append(value)
        else:
            if value < max(self.data):
                self.data[self.data.index(max(self.data))] = value

    def get_data(self):
        return self.data


def main(file,k):
    counts = khmer.Countgraph(k, 1e7, 4)
    # 处理文件中的序列，计算 k-mer 频率
    counts.consume_seqfile(file)
    list_dict = ListDict()
    # 循环获取 k-mer 并进行 hash 计算和筛选
    filtered_kmer_freqs = {}  # 存储筛选后的 k-mer 和其频率,正在想办法优化掉这个列表
    # 使用 khmer 提供的 iterator 一次性遍历 k-mer 和频率
    for record in khmer.ReadParser(file):
        sequence = record.sequence
        # 遍历并输出每个 k-mer 和频率
        for kmer in counts.get_kmers(sequence):#对于这两层循环也在想办法去改
            frequency = counts.get(kmer)
            #print(f"K-mer: {kmer}, Frequency: {frequency}")
            hashvalue1 = murmur3_32(kmer.encode("utf-8"),42)
            filtered_kmer_freqs[hashvalue1] = frequency
            list_dict.add(hashvalue1)
    #print(list_dict.get_data())  # 打印最终的数据列表
    c = list_dict.get_data()
    print(c)
    # 遍历 c 中的每个 最小hash集合,查找它在 filtered_kmer_freqs 中的频率
    for hash in c:
        if hash in filtered_kmer_freqs:
            freq = filtered_kmer_freqs[hash]  # 获取该 k-mer 的频率
            print(f"K-mer: {hash}, Frequency: {freq}")
        else:
            print(f"K-mer: {hash} not found in filtered_kmer_freqs")



# 创建计数表对象，参数根据你的数据调整
k = 21  # k-mer 的长度

# 输入的 FASTA 文件
input_fasta1 = '/home/hzb/dna/test5.fq'
input_fasta1 = '/mnt/00.zhibo/6M100Retrs/bing/YC_con_6_BDDP210000937-1A_1.clean.fq'
main(input_fasta1,k)