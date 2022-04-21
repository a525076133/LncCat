import sys
from FrameKmer import kmer_freq_file

def coding_nocoding_potential(input_file):
    coding = {}
    noncoding = {}
    for line in open(input_file).readlines():
        fields = line.split()
        if fields[0] == 'hexamer':
            continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])
    return coding, noncoding

def maker_hexmer(nc_path,pc_path,hexmer_path,ff="w"):
    f = open(hexmer_path,ff)
    cod = kmer_freq_file(fastafile=pc_path, word_size=6, step_size=3, frame=0)
    noncod = kmer_freq_file(fastafile=nc_path, word_size=6, step_size=1, frame=0)
    cod_sum = 0.0
    cod_sum += sum(cod.values())
    noncod_sum = 0.0
    noncod_sum += sum(noncod.values())
    print('hexamer' + '\t' + 'coding' + '\t' + 'noncoding',file=f)
    for kmer in cod:
        if 'N' in kmer:
            continue
        print(kmer + '\t' + str(float(cod[kmer] / cod_sum)) + '\t' + str(float(noncod[kmer] / noncod_sum)),file=f)
    f.close()

# path_p_train = "/Users/bbbojack/Downloads/文献/datasets/Human/lnc.train.human.B.fa"
# path_n_train = "/Users/bbbojack/Downloads/文献/datasets/Human/pct.train.human.B.fa"
# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Human/lnc.test.human.B.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Human/pct.test.human.B.fa"

# path_p_train = "/Users/bbbojack/Downloads/文献/datasets/Mouse/lnc.train.Mouse.fa"
# path_n_train = "/Users/bbbojack/Downloads/文献/datasets/Mouse/pct.train.Mouse.fa"
# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Mouse/lnc.test.Mouse.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Mouse/pct.test.Mouse.fa"

path_p_train = "/Users/bbbojack/Downloads/文献/datasets/Wheat/lnc.train.Wheat.fa"
path_n_train = "/Users/bbbojack/Downloads/文献/datasets/Wheat/pct.train.Wheat.fa"
path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Wheat/lnc.test.Wheat.fa"
path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Wheat/pct.test.Wheat.fa"

# path_p_train = "/Users/bbbojack/Downloads/文献/datasets/Chicken/lnc.train.chicken.fa"
# path_n_train = "/Users/bbbojack/Downloads/文献/datasets/Chicken/pct.train.chicken.fa"
# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Chicken/lnc.test.chicken.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Chicken/pct.test.chicken.fa"

# path_p_train = "/Users/bbbojack/Downloads/文献/datasets/Zebrafish/lnc.train.zebrafish.fa"
# path_n_train = "/Users/bbbojack/Downloads/文献/datasets/Zebrafish/pct.train.zebrafish.fa"
# path_p_test = "/Users/bbbojack/Downloads/文献/datasets/Zebrafish/lnc.test.zebrafish.fa"
# path_n_test = "/Users/bbbojack/Downloads/文献/datasets/Zebrafish/pct.test.zebrafish.fa"



nc_path_train = path_p_train
pc_path_train = path_n_train

nc_path_test = path_p_test
pc_path_test = path_n_test

hexmer_path ="hexmer.txt"


import os

os.system("cat {path1} {path2} > nc.fa".format(path1=nc_path_test,path2=nc_path_train))
os.system("cat {path1} {path2} > pc.fa".format(path1=pc_path_test,path2=pc_path_train))

maker_hexmer("nc.fa","pc.fa",hexmer_path,ff="w")

os.system("rm -rf ./nc.fa")
os.system("rm -rf ./pc.fa")
