import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
par = argparse.ArgumentParser()
par.add_argument('-f', '--filename')
par.add_argument('-r', '--resultname')
arg = par.parse_args()
gs=['Amamu1','Hebcy2','Lacam2','Lacbi2','Paxin1','Paxru2',
    'Rhice_1','Rhidi1','RhiirA_1','Rhives1','Rhivi1']
newfile=[]
for seq in SeqIO.parse(arg.filename,'fasta'):
    if str(seq.id).split('|')[1] in gs :
        record=SeqRecord(seq.seq[0:-1], id=seq.id, description='')
        newfile.append(record)
for seq in SeqIO.parse(arg.filename,'fasta'):
    if str(seq.id).split('|')[1]  not in gs :
        record=SeqRecord(seq.seq[0:-1], id=seq.id, description='')
        newfile.append(record)
SeqIO.write(newfile,arg.resultname,'fasta')