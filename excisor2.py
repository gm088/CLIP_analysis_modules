import multiprocessing
import re
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fname', type=str)
parser.add_argument('--direction', type=str)
parser.add_argument('--outname', type=str)
parser.add_argument('--outdir', type=str)
parser.add_argument('--minlen', type=int)
args = parser.parse_args()

def process_read(nom,seeq,somm,quali, minlen,direc):
    sequence = Seq(seeq)
    if(direc=="read1"):
        umi = str(Seq(nom.split("_")[1].split()[0][-10:]).reverse_complement())
    if(direc=="read2"):
        umi = str(Seq(nom.split("_")[1].split()[0][:12]).reverse_complement())
    z = re.search("(.*)({0})".format(umi[:5]),str(sequence))
    try:
        if len(z.group(1))>minlen:
            newseeq = z.group(1)
            newqual = quali[:len(newseeq)]
            print("read trimmed and kept")
        else:
            print("read discarded")
            return
    except AttributeError:
        if(len(seeq)>minlen):
            newseeq = seeq
            newqual = quali
            print("read not trimmed; kept because long enough")
        else:
            print("read discarded")
            return
    newread = nom+'\n'+newseeq+'\n'+somm+'\n'+newqual+'\n'
    with open("{0}/{1}".format(args.outdir, args.outname), "a") as f:
        f.write(newread)

mini_len = args.minlen
dirn = args.direction

with gzip.open(args.fname,'rt') as file:
    while True:
        try:
            name = next(file).rstrip('\n')
            seq = next(file).rstrip('\n')
            som = next(file).rstrip('\n')
            qual = next(file).rstrip('\n')
            process_read(name,seq,som,qual,mini_len,dirn)

        except StopIteration:
            break

