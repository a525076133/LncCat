import warnings
warnings.filterwarnings("ignore")
import sys
sys.path.append("/")
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import FrameKmer
from Bio.SeqUtils import ProtParam
import re
import orfipy_core
import numpy as np
import math
import random
top = 3

class Fickett:
    def __init__(self):
        self.content_parameter = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]
        self.position_parameter = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
        '''
        newly calculated lookup table for RNA full length
        '''
        self.position_probability = {
            "A": [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
            "C": [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
            "G": [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
            "T": [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24],
        }
        self.position_weight = {"A": 0.062, "C": 0.093, "G": 0.205, "T": 0.154}
        self.content_probability = {
            "A": [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
            "C": [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
            "G": [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
            "T": [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51]
        }
        self.content_weight = {"A": 0.084, "C": 0.076, "G": 0.081, "T": 0.055}

    def look_up_position_probability(self, value, base):
        '''
        look up positional probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx, val in enumerate(self.position_parameter):
            if (float(value) >= val):
                return float(self.position_probability[base][idx]) * float(self.position_weight[base])

    def look_up_content_probability(self, value, base):
        '''
        look up content probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx, val in enumerate(self.content_parameter):
            if (float(value) >= val):
                return float(self.content_probability[base][idx]) * float(self.content_weight[base])

    def fickett_value(self, dna):
        '''
        calculate Fickett value from full RNA transcript sequence
        '''
        if len(dna) < 2:
            return 0
        fickett_score = 0
        dna = dna
        total_base = len(dna)
        A_content = float(dna.count("A")) / total_base
        C_content = float(dna.count("C")) / total_base
        G_content = float(dna.count("G")) / total_base
        T_content = float(dna.count("T")) / total_base

        phase_0 = dna[::3]
        phase_1 = dna[1::3]
        phase_2 = dna[2::3]

        phase_0_A = phase_0.count("A")
        phase_1_A = phase_1.count("A")
        phase_2_A = phase_2.count("A")
        phase_0_C = phase_0.count("C")
        phase_1_C = phase_1.count("C")
        phase_2_C = phase_2.count("C")
        phase_0_G = phase_0.count("G")
        phase_1_G = phase_1.count("G")
        phase_2_G = phase_2.count("G")
        phase_0_T = phase_0.count("T")
        phase_1_T = phase_1.count("T")
        phase_2_T = phase_2.count("T")

        A_content = float(phase_0_A + phase_1_A + phase_2_A) / total_base
        C_content = float(phase_0_C + phase_1_C + phase_2_C) / total_base
        G_content = float(phase_0_G + phase_1_G + phase_2_G) / total_base
        T_content = float(phase_0_T + phase_1_T + phase_2_T) / total_base
        A_position = np.max([phase_0_A, phase_1_A, phase_2_A]) / (np.min([phase_0_A, phase_1_A, phase_2_A]) + 1.0)
        C_position = np.max([phase_0_C, phase_1_C, phase_2_C]) / (np.min([phase_0_C, phase_1_C, phase_2_C]) + 1.0)
        G_position = np.max([phase_0_G, phase_1_G, phase_2_G]) / (np.min([phase_0_G, phase_1_G, phase_2_G]) + 1.0)
        T_position = np.max([phase_0_T, phase_1_T, phase_2_T]) / (np.min([phase_0_T, phase_1_T, phase_2_T]) + 1.0)

        fickett_score += self.look_up_content_probability(A_content, "A")
        fickett_score += self.look_up_content_probability(C_content, "C")
        fickett_score += self.look_up_content_probability(G_content, "G")
        fickett_score += self.look_up_content_probability(T_content, "T")

        fickett_score += self.look_up_position_probability(A_position, "A")
        fickett_score += self.look_up_position_probability(C_position, "C")
        fickett_score += self.look_up_position_probability(G_position, "G")
        fickett_score += self.look_up_position_probability(T_position, "T")

        return fickett_score

class FindCDS:
    '''
    Find the most like CDS in a given sequence
    The most like CDS is the longest ORF found in the sequence
    When having same length, the upstream ORF is printed
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
    '''

    def __init__(self, seq):
        self.seq = seq
        self.result = (0, 0, 0, 0, 0)
        self.longest = 0
        self.basepair = {"A": "T", "T": "A", "U": "A", "C": "G", "G": "C", "N": "N", "X": "X"}

    def _reversecompliment(self):
        return "".join(self.basepair[base] for base in self.seq)[::-1]

    def get_codons(self, frame_number):
        coordinate = frame_number
        while coordinate + 3 <= len(self.seq):
            yield (self.seq[coordinate:coordinate + 3], coordinate)
            coordinate += 3

    def find_longest_in_one(self, myframe, direction, start_codon, stop_codon):
        triplet_got = self.get_codons(myframe)
        starts = start_codon
        stops = stop_codon
        while True:
            try:
                codon, index = triplet_got.__next__()
            except StopIteration:
                break
            if codon in starts and codon not in stops:
                orf_start = index
                end_extension = False
                while True:
                    try:
                        codon, index = triplet_got.__next__()
                    except StopIteration:
                        end_extension = True
                        integrity = -1
                    if codon in stops:
                        integrity = 1
                        end_extension = True
                    if end_extension:
                        orf_end = index + 3
                        Length = (orf_end - orf_start)
                        if Length > self.longest:
                            self.longest = Length
                            self.result = [direction, orf_start, orf_end, Length, integrity]
                        if Length == self.longest and orf_start < self.result[1]:
                            self.result = [direction, orf_start, orf_end, Length, integrity]
                        break

    def longest_orf(self, direction, start_codon={"ATG": None}, stop_codon={"TAG": None, "TAA": None, "TGA": None}):
        return_orf = ""
        for frame in range(3):
            self.find_longest_in_one(frame, "+", start_codon, stop_codon)
        return_orf = self.seq[self.result[1]:self.result[2]][:]
        start_coordinate = self.result[1]
        strand_direction = "+"
        orf_integrity = self.result[4]
        if direction == "-":
            self.seq = self._reversecompliment()
            for frame in range(3):
                self.find_longest_in_one(frame, "-", start_codon, stop_codon)
            if self.result[0] == "-":
                return_orf = self.seq[self.result[1]:self.result[2]][:]
                start_coordinate = self.result[1]
                strand_direction = "-"
                orf_integrity = self.result[4]
        return return_orf, orf_integrity

# len = 30
def CTD_atcg(seq):
	n=len(seq)-1
	n=float(n)
	num_A,num_T,num_G,num_C=0,0,0,0
	AT_trans,AG_trans,AC_trans,TG_trans,TC_trans,GC_trans=0,0,0,0,0,0
	for i in range(len(seq)-1):
		if seq[i]=="A":
			num_A=num_A+1
		if seq[i]=="T":
			num_T=num_T+1
		if seq[i]=="G":
			num_G=num_G+1
		if seq[i]=="C":
			num_C=num_C+1
		if (seq[i]=="A" and seq[i+1]=="T") or (seq[i]=="T" and seq[i+1]=="A"):
			AT_trans=AT_trans+1
		if (seq[i]=="A" and seq[i+1]=="G") or (seq[i]=="G" and seq[i+1]=="A"):
			AG_trans=AG_trans+1
		if (seq[i]=="A" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="A"):
			AC_trans=AC_trans+1
		if (seq[i]=="T" and seq[i+1]=="G") or (seq[i]=="G" and seq[i+1]=="T"):
			TG_trans=TG_trans+1
		if (seq[i]=="T" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="T"):
			TC_trans=TC_trans+1
		if (seq[i]=="G" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="G"):
			GC_trans=GC_trans+1

	a,t,g,c=0,0,0,0
	A0_dis,A1_dis,A2_dis,A3_dis,A4_dis=0.0,0.0,0.0,0.0,0.0
	T0_dis,T1_dis,T2_dis,T3_dis,T4_dis=0.0,0.0,0.0,0.0,0.0
	G0_dis,G1_dis,G2_dis,G3_dis,G4_dis=0.0,0.0,0.0,0.0,0.0
	C0_dis,C1_dis,C2_dis,C3_dis,C4_dis=0.0,0.0,0.0,0.0,0.0
	for i in range(len(seq)-1):
		if seq[i]=="A":
			a=a+1
			if a == 1:
				A0_dis=((i*1.0)+1)/n
			if a == int(round(num_A/4.0)):
				A1_dis=((i*1.0)+1)/n
			if a == int(round(num_A/2.0)):
				A2_dis=((i*1.0)+1)/n
			if a == int(round((num_A*3/4.0))):
				A3_dis=((i*1.0)+1)/n
			if a == num_A:
				A4_dis=((i*1.0)+1)/n
		if seq[i]=="T":
			t=t+1
			if t == 1:
				T0_dis=((i*1.0)+1)/n
			if t == int(round(num_T/4.0)):
				T1_dis=((i*1.0)+1)/n
			if t == int(round((num_T/2.0))):
				T2_dis=((i*1.0)+1)/n
			if t == int(round((num_T*3/4.0))):
				T3_dis=((i*1.0)+1)/n
			if t == num_T:
				T4_dis=((i*1.0)+1)/n
		if seq[i]=="G":
			g=g+1
			if g == 1:
				G0_dis=((i*1.0)+1)/n
			if g == int(round(num_G/4.0)):
				G1_dis=((i*1.0)+1)/n
			if g == int(round(num_G/2.0)):
				G2_dis=((i*1.0)+1)/n
			if g == int(round(num_G*3/4.0)):
				G3_dis=((i*1.0)+1)/n
			if g == num_G:
				G4_dis=((i*1.0)+1)/n
		if seq[i]=="C":
			c=c+1
			if c == 1:
				C0_dis=((i*1.0)+1)/n
			if c == int(round(num_C/4.0)):
				C1_dis=((i*1.0)+1)/n
			if c == int(round(num_C/2.0)):
				C2_dis=((i*1.0)+1)/n
			if c == int(round(num_C*3/4.0)):
				C3_dis=((i*1.0)+1)/n
			if c == num_C:
				C4_dis=((i*1.0)+1)/n
	return [num_A/n,num_T/n,
             num_G/n,num_C/n,
             AT_trans/(n-1),AG_trans/(n-1),
             AC_trans/(n-1),TG_trans/(n-1),
             TC_trans/(n-1),GC_trans/(n-1),
             A0_dis,A1_dis,A2_dis,
             A3_dis,A4_dis,
             T0_dis,T1_dis,
             T2_dis,T3_dis,
             T4_dis,G0_dis,
             G1_dis,G2_dis,
             G3_dis,G4_dis,
             C0_dis,C1_dis,
             C2_dis,C3_dis,C4_dis]


# 1 float
def get_length(seq):
    return np.log(len(seq)+1)

# 1 float
def get_fickett_value(seq):
    fickett = Fickett()
    return fickett.fickett_value(seq)

# 1 int
def get_stop_codon_num(seq):
    translate_prot = Seq(seq).translate()
    stop_num = translate_prot.count("*")
    return stop_num

# 1 float
def get_stop_codon_frequency(seq):
    stop_num = get_stop_codon_num(seq)
    transript_length = get_length(seq)
    stop_freq = float(stop_num) / transript_length
    return stop_freq

def get_orf(seq):
    findCDS = FindCDS(seq)
    return_orf, orf_integrity = findCDS.longest_orf(seq)
    return return_orf, orf_integrity

def get_orf_coverge(seq):
    transript_length = get_length(seq)
    orf, _ = get_orf(seq)
    orf_length = len(orf)
    ORF_coverage = float(orf_length) / transript_length
    return ORF_coverage

def get_orf_frame_score(seq):
    ORF_length_in_frame1, _ = get_orf(seq)
    ORF_length_in_frame2, _ = get_orf(seq[1:])
    ORF_length_in_frame3, _ = get_orf(seq[2:])

    ORF_length_in_frame1 = len(ORF_length_in_frame1)
    ORF_length_in_frame2 = len(ORF_length_in_frame2)
    ORF_length_in_frame3 = len(ORF_length_in_frame3)

    ORF_len = [ORF_length_in_frame1, ORF_length_in_frame2, ORF_length_in_frame3]
    ORF_frame = ((ORF_len[0] - ORF_len[1]) ** 2 + (ORF_len[0] - ORF_len[2]) ** 2 + (ORF_len[1] - ORF_len[2]) ** 2) / 2
    return ORF_frame

def get_GC1(mRNA):
    if len(mRNA) < 3:
        numGC = 0
        mRNA = 'ATG'
    else:
        numGC = mRNA[0::3].count("C") + mRNA[0::3].count("G")
    return numGC * 1.0 / len(mRNA) * 3

# 0
def get_GC2(mRNA):
    if len(mRNA) < 3:
        numGC = 0
        mRNA = 'ATG'
    else:
        numGC = mRNA[1::3].count("C") + mRNA[1::3].count("G")
    return numGC * 1.0 / len(mRNA) * 3

def get_GC3(mRNA):
    if len(mRNA) < 3:
        numGC = 0
        mRNA = 'ATG'
    else:
        numGC = mRNA[2::3].count("C") + mRNA[2::3].count("G")
    return numGC * 1.0 / len(mRNA) * 3

def get_gc1_frame_score(seq):
    GC1_in_frame1 = get_GC1(seq)
    GC1_in_frame2 = get_GC1(seq[1:])
    GC1_in_frame3 = get_GC1(seq[2:])
    GC1_all = [GC1_in_frame1, GC1_in_frame2, GC1_in_frame3]
    GC1_frame = ((GC1_all[0] - GC1_all[1]) ** 2 + (GC1_all[0] - GC1_all[2]) ** 2 + (GC1_all[1] - GC1_all[2]) ** 2) / 2
    return GC1_frame

def get_gc2_frame_score(seq):
    GC2_in_frame1 = get_GC2(seq)
    GC2_in_frame2 = get_GC2(seq[1:])
    GC2_in_frame3 = get_GC2(seq[2:])
    GC2_all = [GC2_in_frame1, GC2_in_frame2, GC2_in_frame3]
    GC2_frame = ((GC2_all[0] - GC2_all[1]) ** 2 + (GC2_all[0] - GC2_all[2]) ** 2 + (GC2_all[1] - GC2_all[2]) ** 2) / 2
    return GC2_frame

def get_gc3_frame_score(seq):
    GC3_in_frame1 = get_GC3(seq)
    GC3_in_frame2 = get_GC3(seq[1:])
    GC3_in_frame3 = get_GC3(seq[2:])
    GC3_all = [GC3_in_frame1, GC3_in_frame2, GC3_in_frame3]
    GC3_frame = ((GC3_all[0] - GC3_all[1]) ** 2 + (GC3_all[0] - GC3_all[2]) ** 2 + (GC3_all[1] - GC3_all[2]) ** 2) / 2
    return GC3_frame

def get_stop_frame_score(seq):
    stop_num_in_frame1 = get_stop_codon_num(seq)
    stop_num_in_frame2 = get_stop_codon_num(seq[1:])
    stop_num_in_frame3 = get_stop_codon_num(seq[2:])
    stop_num_all = [stop_num_in_frame1, stop_num_in_frame2, stop_num_in_frame3]
    stop_num_frame = ((stop_num_all[0] - stop_num_all[1]) ** 2 + (stop_num_all[0] - stop_num_all[2]) ** 2 + (
            stop_num_all[1] - stop_num_all[2]) ** 2) / 2
    return stop_num_frame

# a string e.g. "IASGSFG*DQM*CLMFLSVYLIRLR*KGFILKGL*NGALWVAALWPMGRWEDCGAVSCRGAVPHRGA"
def mRNA_translate(mRNA):
    return Seq(mRNA).translate()


def get_Mw(seq):
    seqprot = mRNA_translate(seq)
    strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
    newseqprot = strinfoAmbiguous.sub("", str(seqprot))
    protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot).replace("*",""))
    mw = protparam_obj.molecular_weight()
    return mw

def get_pI(seq):
    seqprot = mRNA_translate(seq)
    strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
    newseqprot = strinfoAmbiguous.sub("", str(seqprot))
    protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot).replace("*", ""))
    pI = protparam_obj.isoelectric_point()
    return pI

def get_gravy(seq):
    seqprot = mRNA_translate(seq)
    strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
    newseqprot = strinfoAmbiguous.sub("", str(seqprot))
    protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot).replace("*", ""))
    Gravy = protparam_obj.gravy()
    return Gravy

def get_instablility_index(seq):
    seqprot = mRNA_translate(seq)
    strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
    newseqprot = strinfoAmbiguous.sub("", str(seqprot))
    protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot).replace("*", ""))
    instablility_index = protparam_obj.instability_index()
    return instablility_index

def coding_nocoding_potential():
    coding = {}
    noncoding = {}
    for line in open("hexmer.txt").readlines():
        fields = line.split()
        if fields[0] == 'hexamer': continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])
    return coding, noncoding

start_codons = 'ATG'
stop_codons = 'TAG,TAA,TGA'
Coverage = 0

maxlen = 100000

def find_max_orf(seq):
    seq = seq.upper()
    res = []
    for start,stop,strand,description in orfipy_core.orfs(seq,minlen=10,maxlen=maxlen):
        tmp = [int(start),int(stop),int(description.split(";")[2].split("=")[-1])]
        res.append(tmp)
    if res:
        max_index = max(res,key=lambda x:x[2])
    else:
        max_index = [0,0,0]
    return max_index

def find_top_orf(seq,top=2):
    seq = seq.upper()
    res = []
    for start,stop,strand,description in orfipy_core.orfs(seq,minlen=10,maxlen=maxlen):
        tmp = [int(start),int(stop),int(description.split(";")[2].split("=")[-1])]
        res.append(tmp)
    if res:
        max_index = list(sorted(res,key=lambda x:x[2],reverse=True))[0:top]
        if len(max_index)<top:
            max_index = max_index+[[0,0,0]]*(top-len(max_index))
    else:
        max_index = [[0,0,0]]*top
    return max_index

def kmer_count(seq):
    # std = [x+y+z+p+q for x in "ATCG" for y in "ATCG" for z in "ATCG" for p in "ATCG" for q in "ATCG"]
    std = [x+y+z for x in "ATCG" for y in "ATCG" for z in "ATCG"]
    cc = [seq.count(mer)/len(seq) for mer in std]
    return cc

def kmer_count2(seq):
    std = [x+y for x in "ATCG" for y in "ATCG"]
    cc = [seq.count(mer)/len(seq) for mer in std]
    return cc

def gc(seq):
    seq0 = seq[0::3]
    seq1 = seq[1::3]
    seq2 = seq[2::3]
    res = []
    for seqx in [seq0,seq1,seq2]:
        gcx = (seqx.count("G")+seqx.count("C"))/len(seq)
        res.append(gcx)
    return res

def orf_num(seq):
    num = len(orfipy_core.orfs(seq, minlen=10, maxlen=len(seq)))/len(seq)
    return num

def eiip(seq):
    eiip_dict = {
        "A": 0.12601,
        "T": 0.13400,
        "C": 0.08060,
        "G": 0.13350
    }
    res = []
    for x in seq:
        res.append(eiip_dict[x])
    res = np.abs(np.fft.fft(np.array(res)))
    res = np.percentile(res,q=[0,25,50,75,100]).tolist()
    return res

def get_avgList(path_p,path_n):
    kmers_list = [x+y+z for x in "ATCG" for y in "ATCG" for z in "ATCG"]

    p_seqs = []
    tmp_list = list(read_fa(path_p).values())
    p_seqs.extend(tmp_list)

    n_seqs = []
    tmp_list = list(read_fa(path_n).values())
    n_seqs.extend(tmp_list)

    lnc_avg = np.zeros(64)
    for seq in p_seqs:
        for i in range(len(kmers_list)):
            lnc_avg[i] += seq.count(kmers_list[i])

    pct_avg = np.zeros(64)
    for seq in n_seqs:
        for i in range(len(kmers_list)):
            pct_avg[i] += seq.count(kmers_list[i])

    return lnc_avg/len(p_seqs), pct_avg/len(n_seqs)

def DistFeature(seq):
    kmers_list = [x+y+z for x in "ATCG" for y in "ATCG" for z in "ATCG"]
    fre = np.zeros(64)
    for i in range(len(kmers_list)):
        fre[i] += seq.count(kmers_list[i])

    EucDist_LNC = 0
    LogDist_LNC = 0
    EucDist_PCT = 0
    LogDist_PCT = 0

    for i in range(64):
        EucDist_LNC += pow(fre[i]-freq_lnc[i],2)
        EucDist_PCT += pow(fre[i]-freq_pct[i],2)
        if fre[i] == 0:
            continue
        else:
            LogDist_LNC += math.log(fre[i]/freq_lnc[i])
            LogDist_PCT += math.log(fre[i]/freq_pct[i])

    EucDist_LNC = math.sqrt(EucDist_LNC)
    LogDist_LNC = LogDist_LNC/64
    EucDist_PCT = math.sqrt(EucDist_PCT )
    LogDist_PCT = LogDist_PCT/64

    EucDist_Ratio = EucDist_LNC/EucDist_PCT
    LogDist_Ratio = LogDist_LNC/LogDist_PCT

    return [EucDist_LNC, LogDist_LNC, EucDist_PCT, LogDist_PCT, EucDist_Ratio, LogDist_Ratio]

def calc_fs(seq):
    seq = seq.upper()
    codon_num = str(Seq(seq).translate()).count("*")
    fickett_val = Fickett().fickett_value(seq)
    gcc = gc(seq)
    conda_ratio = codon_num/len(seq)

    coding, noncoding = coding_nocoding_potential()
    hexamer = FrameKmer.kmer_ratio(seq, 6, 3, coding, noncoding)

    gc1f = get_gc1_frame_score(seq)
    gc2f = get_gc2_frame_score(seq)
    gc3f = get_gc3_frame_score(seq)

    mw = get_Mw(seq)
    pi = get_pI(seq)
    gravy = get_gravy(seq)
    instability_index = get_instablility_index(seq)
    gc_list = [gc1f,gc2f,gc3f,mw,pi,gravy,instability_index,conda_ratio]

    start,end,length = find_max_orf(seq)
    cover = length/len(seq)
    ctd_list = CTD_atcg(seq)
    df_list = DistFeature(seq)
    k_mer = kmer_count(seq)
    k_mer2 = kmer_count2(seq)
    res = [codon_num,fickett_val,cover,length,hexamer] + gcc  + ctd_list + df_list + gc_list + k_mer + k_mer2
    return res

def orf_fs(seq):
    orfs = []
    seq = seq.upper()
    s_e_l_s = find_top_orf(seq,top=3)

    orf_seq = ""
    for start,end,length in s_e_l_s:
        if length:
            orf_seq = seq[start:end]
        else:
            orf_seq = seq
        orfs.append(orf_seq)

    ress = []
    for orf_seq in orfs:
        seq = orf_seq

        seq_pro = str(Seq(seq).translate()).replace("*","")
        mw = Bio.SeqUtils.ProtParam.ProteinAnalysis(seq_pro).molecular_weight()
        gravy = Bio.SeqUtils.ProtParam.ProteinAnalysis(seq_pro).gravy()

        len_p = Bio.SeqUtils.ProtParam.ProteinAnalysis(seq_pro).length
        ss = list(Bio.SeqUtils.ProtParam.ProteinAnalysis(seq_pro).secondary_structure_fraction())
        gc1f = get_gc1_frame_score(seq)
        gc2f = get_gc2_frame_score(seq)
        gc3f = get_gc3_frame_score(seq)
        gc_list = [gc1f,gc2f,gc3f]
        fickett_val = Fickett().fickett_value(seq)
        coding, noncoding = coding_nocoding_potential()
        hexamer = FrameKmer.kmer_ratio(seq, 6, 3, coding, noncoding)
        ctd_list = CTD_atcg(seq)
        k_mer = kmer_count(seq)
        gcc = gc(seq)
        eiip_list = eiip(seq)
        tmp = [mw,gravy,fickett_val,hexamer,len_p] + gcc + ss + ctd_list + gc_list + k_mer + eiip_list
        ress.extend(tmp)
    return ress

def pip_fs(seq):
    res = calc_fs(seq) + orf_fs(seq)
    return res

def read_fa(path):
    rx = list(SeqIO.parse(path, format="fasta"))
    res = {}
    for x in rx:
        id = str(x.id)
        seq = str(x.seq).upper()
        if set(list(seq)) == set(list("ATCG")):
            res[id] = seq
    return res

def upseq(seq,k=1):
    new_seq = list(seq)
    for i in range(k):
        index = random.choice(list(range(len(seq))))
        new_seq[index] = random.choice(["A", "T", "C", "G"])
    return "".join(new_seq)

def aug(seqs:list):
    aseqs =  seqs + [upseq(seq, k=int(len(seq)*0.01)) for seq in seqs] + [upseq(seq, k=int(len(seq)*0.03)) for seq in seqs] + [upseq(seq, k=int(len(seq)*0.05)) for seq in seqs]+[upseq(seq, k=int(len(seq)*0.07)) for seq in seqs]
    return aseqs


# path_p_train = "./Mouse/lnc.train.mouse.fa"
# path_n_train = "./Mouse/pct.train.mouse.fa"
# path_p_test = "./Mouse/lnc.test.mouse.fa"
# path_n_test = "./Mouse/pct.test.mouse.fa"

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

freq_lnc, freq_pct = get_avgList(path_p=path_p_train,path_n=path_n_train)
