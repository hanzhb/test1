#添加了mur32和thomas——wang两个hash函数，其中mur32更加复杂一些，thomas可以进行反向函数，获得kmer
#main方法中的双层for循环还没有解决，由于其中拆分包含其他字符，导致thomas方法无法正常运行。
#任务，找到双层for循环替代方法。
import khmer

import time

def murmur3_32(kmer) -> int:
    key = kmer.encode("utf-8")
    h1 = seed = 42
    length = len(key)
    nblocks = length // 4

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

class thomas_wang1:
    mask = 0xFFFFFFFFFFFFFFFF

    def __init__(self, kmer):
        self.kmer = kmer
        self.k = len(self.kmer)

    def hash64(self):
        def dna_to_binary(dna_sequence):
            # 映射碱基到二进制
            base_map = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
            # 将碱基转换为二进制
            binary_sequence = ''.join(base_map[base] for base in dna_sequence)
            return int(binary_sequence, 2)

        key = dna_to_binary(self.kmer)
        key = (~key + (key << 21)) & self.mask
        key = key ^ (key >> 24)
        key = ((key + (key << 3)) + (key << 8)) & self.mask
        key = key ^ (key >> 14)
        key = ((key + (key << 2)) + (key << 4)) & self.mask
        key = key ^ (key >> 28)
        key = (key + (key << 31)) & self.mask

        return key

    def hash64i(self, key):
        def binary_to_dna(binary_sequence, k):
            # 映射二进制到碱基
            reverse_base_map = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
            # 将二进制按每两位分割并转换为碱基
            binary_str = format(key, '0{}b'.format(2 * k))
            sequence = ''.join(reverse_base_map[binary_str[i:i + 2]] for i in range(0, len(binary_str), 2))
            return sequence

        tmp = (key - (key << 31)) & self.mask
        key = (key - (tmp << 31)) & self.mask
        tmp = key ^ (key >> 28)
        key = key ^ (tmp >> 28)
        key = (key * 14933078535860113213) & self.mask
        tmp = key ^ (key >> 14)
        tmp = key ^ (tmp >> 14)
        tmp = key ^ (tmp >> 14)
        key = key ^ (tmp >> 14)
        key = (key * 15244667743933553977) & self.mask
        tmp = key ^ (key >> 24)
        key = key ^ (tmp >> 24)
        tmp = ~key
        tmp = ~(key - (tmp << 21))
        tmp = ~(key - (tmp << 21))
        key = ~(key - (tmp << 21)) & self.mask

        kmer1 = binary_to_dna(key, self.k)
        return kmer1, key


class ListDict:
    def __init__(self, max_size):
        self.max_size = max_size
        self.data = {}  # 用于存储值和频率的字典

    def add(self, value, frequency):
        # 如果字典的大小还没达到最大值，直接添加
        if len(self.data) < self.max_size:
            self.data[value] = frequency
        else:
            # 如果新传入的值比当前字典中最大的值还要小
            max_value = max(self.data.keys())  # 获取字典中最大的值
            if value < max_value:
                # 替换掉最大值，并保留频率
                del self.data[max_value]
                self.data[value] = frequency

    def get_data(self):
        return self.data  # 返回包含值和频率的字典


def main(file,k):
    counts = khmer.Countgraph(k, 1e7, 4)
    # 处理文件中的序列，计算 k-mer 频率
    print("time",time.time())
    counts.consume_seqfile(file)
    list_dict = ListDict(20)
    # 循环获取 k-mer 并进行 hash 计算和筛选
    # 使用 khmer 提供的 iterator 一次性遍历 k-mer 和频率
    print("time2", time.time())
    for record in khmer.ReadParser(file):
        sequence = record.sequence
        # 遍历并输出每个 k-mer 和频率
        for kmer in counts.get_kmers(sequence):#对于这两层循环也在想办法去改
            frequency = counts.get(kmer)
            #print(kmer)
            #print(f"K-mer: {kmer}, Frequency: {frequency}")
            hashvalue1 = murmur3_32(kmer)
            #hashvalue2 = thomas_wang1(kmer)
            #hashvalue2 = hashvalue2.hash64()

            list_dict.add(hashvalue1,frequency)
    #print(list_dict.get_data())  # 打印最终的数据列表
    print("time", time.time())
    c = list_dict.get_data()
    print(c)




# 创建计数表对象，参数根据你的数据调整
k = 16  # k-mer 的长度

# 输入的 FASTA 文件
input_fasta1 = '/home/hzb/dna/test5.fq'
# = '/mnt/00.zhibo/6M100Retrs/bing/YC_con_6_BDDP210000937-1A_1.clean.fq'
main(input_fasta1,k)