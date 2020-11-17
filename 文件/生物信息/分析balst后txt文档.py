from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import re, argparse


# 酸性磷酸酶蛋白序列.fas
class Select():
    def __init__(self):
        self.hmmer_txt = '/users/yanrui/desktop/nature/result'
        par = argparse.ArgumentParser()
        par.add_argument('-db', '--hmmerdb')
        par.add_argument('-f', '--file')
        arg = par.parse_args()
        self.hmmdb = '/users/yanrui/desktop/nature/target_db/' + str(arg.hmmerdb) + '.hmm'
        self.file = arg.file
        txt = pd.read_table(self.file, names=['条带编号', '起始位点', '终止位点', '酶编号', '得分'])
        txt.drop_duplicates(['条带编号'])
        self.seqname = [str(i) for i in txt['条带编号']]
        if '担子' in str(self.file):
            self.query = '/users/yanrui/desktop/nature/担子菌门/担子菌prot.fasta'
        else:
            self.query = '/users/yanrui/desktop/nature/子囊菌门/子囊菌门prot.fasta'

    def hmmer(self):
        txt_path = '/users/yanrui/desktop/nature/result/hmm.txt'
        os.system(f'hmmsearch {self.hmmdb} {self.query} > {txt_path}')
        with open(txt_path) as file:
            wd = file.read()
            a = re.findall(r'jgi.*? ', wd)
            data = list(set([i.strip() for i in a]))
            self.sec_selcted_seq = [i for i in data if i in self.seqname]
            print('blast，hmmer筛选后有%s条目标序列序列' % len(self.sec_selcted_seq))

    def get_target_gene(self):
        target_prot = []
        for i in SeqIO.parse(self.query, 'fasta'):
            if str(i.id) in self.sec_selcted_seq:
                record = SeqRecord(i.seq[0:-1], id=str(i.id).replace('%', '').replace('#', ''), description='')
                target_prot.append(record)
        outname = str(self.file).replace('txt', 'fasta')
        SeqIO.write(target_prot, outname, 'fasta')
        print('目标序列已保存至', outname)


yr = Select()
yr.hmmer()
yr.get_target_gene()
