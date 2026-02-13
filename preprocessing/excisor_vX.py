
import threading
import multiprocessing
import re
from Bio.Seq import Seq
import gzip
import time
import argparse
from math import ceil

class Read():

    def __init__(self, name, seq, som, qual):
        self.name = name
        self.seq = seq
        self.som = som
        self.qual = qual

class Statistics(object):

	def __init__(self) -> None:
		self.discarded_reads = 0
		self.trimmed = 0
		self.not_trimmed = 0

def write_to_file(lock, reads_obj, fname):

	#acquire lock
	lock.acquire(blocking=True)

	## write to file
	with open(fname, 'a') as f:
		for i in reads_obj:
			try:
				f.write(i)
			except TypeError:
				pass

	##release lock
	lock.release()

def process_read(read, minlen, direc):
    
    sequence = Seq(read.seq)
    if(direc=="read1"):
        umi = str(Seq(read.name.split("_")[1].split()[0][-10:]).reverse_complement())
    if(direc=="read2"):
        umi = str(Seq(read.name.split("_")[1].split()[0][:12]).reverse_complement())
    z = re.search("(.*)({0})".format(umi[:5]),str(sequence))
    try:
        if len(z.group(1))>minlen:
            newseeq = z.group(1)
            newqual = read.qual[:len(newseeq)]
            #print("read trimmed and kept")
        else:
            #print("read discarded")
            return
    except AttributeError:
        if(len(read.seq)>minlen):
            newseeq = read.seq
            newqual = read.qual
            #print("read not trimmed; kept because long enough")
        else:
            #print("read discarded")
            return
    newread = read.name+'\n'+newseeq+'\n'+read.som+'\n'+newqual+'\n'
    return(newread)

def chunked(iterable, n):
	# n is the number of batches
	num_elements_per_chunk = ceil(len(iterable)/n)
	for i in range(0, len(iterable), num_elements_per_chunk):
		yield iterable[i:i+num_elements_per_chunk]

def batch_process(reads, n_procs, n_threads, outfname):
	
	##multiprocess for processing
	with multiprocessing.Pool( n_procs ) as pool:
		newreads = pool.starmap(process_read, reads, len(reads)//4)
			
		#### split reads into chunks, one chunk for each thread
		read_chunks = list(chunked(newreads, n_threads))
		
		threads = []
		for chunk in read_chunks:
			threads.append(
				threading.Thread(target=write_to_file, 
				kwargs={"lock": thread_lock, "reads_obj": chunk, "fname": outfname})
			)
		# Start each thread
		for t in threads:
			t.start()

		# Wait for all threads to finish
		for t in threads:
			t.join()

####
# TO DO - stats
####


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('--fname', type=str)
	parser.add_argument('--direction', type=str)
	parser.add_argument('--outname', type=str)
	parser.add_argument('--minlen', type=int)
	parser.add_argument('--cores', type=int)
	args = parser.parse_args()

	start_time = time.time()
	thread_lock = threading.Lock()
	reporting_stats = Statistics()

	### dirn must be either 'read1' or 'read2'
	if not(args.direction == "read1" or args.direction == "read2"):
		raise AssertionError("direction must be read1 or read2")

	fname = args.fname
	#fname = "/Users/IEO5559/enhancedClip/test_fastq/test_R1.fastq.gz"
	mini_len = args.minlen
	dirn = args.direction
	outfname = args.outname

	batch_size = 50000
	num_threads = 10
	num_processes = args.cores

	### read 50k reads at a time
	reads = []

	with gzip.open(fname,'rt') as file:
		while True:
			try:
				name = next(file).rstrip('\n')
				seq = next(file).rstrip('\n')
				som = next(file).rstrip('\n')
				qual = next(file).rstrip('\n')

				# each element is a tuple, with the read obj and the function args
				reads.append( (Read(name, seq, som, qual), mini_len, dirn ) )
				if(len(reads) == batch_size):
					# hand off to a batch processing function
					batch_process(reads, num_processes, num_threads, outfname)
					print("wrote 50k reads: ", " --- %s seconds ---" % (time.time() - start_time))
					reads = []
		
			except StopIteration:
				if(len(reads) > 0):
					## process the remaining reads, then exit the while loop
					print("processing last batch of reads")
					batch_process(reads, num_processes, num_threads, outfname)
					break
				else:
					#exit the while loop
					break



