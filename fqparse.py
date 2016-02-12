import gzip
import argparse
from subprocess import Popen,PIPE

from itertools import izip,izip_longest
from difflib import SequenceMatcher as SM

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def mismatch_eval(bcodes,bcode_seq,mismatches=1):
    candidates=[]
    tolerance=(len(bcode_seq)-mismatches)/float(len(bcode_seq))
    for bcode in bcodes:
        ratio=SM(None,bcode,bcode_seq).ratio()
        if ratio>=tolerance:
            candidates.append(barcode)
    if len(candidates)==1:
        return candidates[0]

def match_to_sample(bcodelist,seq,idx,stats_dict,idmatch_dict,notfound_dict,barcodelength=5):
    barcode=seq1[:barcodelength]
    countflag=0;retainflag=0
    if (barcode,idx) in stats_dict and (barcode,idx) in idmatch_dict:
        stats_dict[(barcode,idx)]+=1
        countflag=1
        retainflag=1
    elif (barcode,idx) in idmatch_dict:
        stats_dict[(barcode,idx)]=1
        countflag=1
        retainflag=1
    else: 
        fuzzy_barcode=mismatch_eval(bcodelist,barcode)
        if fuzzy_barcode!=None:
            if (fuzzy_barcode,idx) in stats_dict and (fuzzy_barcode,idx) in idmatch_dict:
                stats_dict[(fuzzy_barcode,idx)]+=1
                countflag=1
                retainflag=1
            elif (fuzzy_barcode,idx) in idmatch_dict:
                stats_dict[(fuzzy_barcode,idx)]=1
                countflag=1
                retainflag=1 
            barcode=fuzzy_barcode              
        else:
            if (barcode,idx) in notfound_dict:
                notfound_dict[(barcode,idx)]+=1
            else:
                notfound_dict[(barcode,idx)]=1    
            countflag=1

    if retainflag>0:
        test_pass=True
    else:
        test_pass=False
            
    return test_pass,barcode,countflag,retainflag

def barcode_check(read1,read2,overhang1,overhang2,barcodelength=5):
    if read1[barcodelength:barcodelength+len(overhang1)]==overhang1 and read2[:len(overhang2)]==overhang2:
        return True
    
def get_illumina_run_id(read1file):
    cmdline="zcat %s |head -1 |awk -F: '{print $1\":\"$2\":\"$3}'" % read1file
    stdout,stderr=Popen(cmdline,shell=True,stderr=PIPE,stdout=PIPE).communicate()
    return stdout.split()[0]
        
if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description='Arguments for ddRAD demultiplexing of inline barcode - Illumina indexed samples')
    parser.add_argument('-R1','--read1_file',dest='r1',type=str,help='read 1 fastq.gz file')
    parser.add_argument('-R2','--read2_file',dest='r2',type=str,help='read 2 fastq.gz file')
    parser.add_argument('-c','--indices_barcode_id_file',dest='sampleinfo',type=str,help='space or tab separated file without header')
    #parser.add_argument('-z','--gz',action='store_true',destination='gzipped',help='infiles are gzipped fastq')
    parser.add_argument('-o1','--R1_overhang_sequence',type=str,dest='overseq1',help='overhang sequence for rad tag on R1')
    parser.add_argument('-o2','--R2_overhang_sequence',type=str,dest='overseq2',help='overhang sequence for rad tag on R2')
    parser.add_argument('-l','--logfile',dest='logout',type=str,help='log of binning for rad data')
    #parser.add_argument('-p','--outfile_prefix',dest='outprefix',type=str,help='prefix to demultiplexed fastq so multiple processescontaining same index-barcode-sampleID combination do not overwrite each other')
    
    opts = parser.parse_args()
    
    barcodes=open(opts.sampleinfo,'r')
    barcode_index_dict={}
    counter_dict={}
    unfound_dict={}
    index_to_barcodes={}
    R1_handles={}
    R2_handles={}
    barcode_drops={}
    
    runid=get_illumina_run_id(opts.r1)
    #print 'runid is', runid

    

    for line in barcodes:
        illumina_idx,barcode,sample_id=line.strip().split()
        barcode_index_dict[(barcode,illumina_idx)]=sample_id

    for bc,idx in barcode_index_dict.keys():
        if idx in index_to_barcodes:
            if bc not in index_to_barcodes[idx]:
                index_to_barcodes[idx].append(bc)        
        else:
            index_to_barcodes[idx]=[bc]
            
    #print 'test of barcode_index_dict',barcode_index_dict[('GTCGA', 'CTTGTA')]  
    #print 'test of index to barcodes','GTCGA' in index_to_barcodes['CTTGTA']     

#print 'index to barcodes', index_to_barcodes

    total_reads=0;reads_retained=0    
 
    with gzip.open(opts.r1,'rb') as f1, gzip.open(opts.r2,'rb') as f2:
        R1=grouper(f1,4)
        R2=grouper(f2,4)
        counter=0
        for entry1 in R1:
            counter+=1
            if counter%100000==0:
                print "%s reads processed" % counter
            head1,seq1,placeholder1,qual1=[i.strip() for i in entry1]
            head2,seq2,placeholder2,qual2=[j.strip() for j in R2.next()]
           

            index_from_read1=head1.split()[-1].split(':')[-1]
            index_from_read2=head2.split()[-1].split(':')[-1]
    
            if index_from_read1!=index_from_read2:
                print "Illumina indices don't match in R1 and R2"
 
            else:
                found,barcode_found,countadd,retainadd=match_to_sample(index_to_barcodes[index_from_read1],seq1,index_from_read1,counter_dict,barcode_index_dict,unfound_dict)
                #print 'stats ',found,barcode_found,countadd,retainadd
                #print 'id',barcode_index_dict[(barcode_found,index_from_read1)]
                #print 'index from R1',index_from_read1
                                              
                total_reads+=countadd
                #reads_retained+=retainadd
                
                if found==True:
                    read_handle="_".join([barcode_index_dict[(barcode_found,index_from_read1)],barcode_found,index_from_read1,runid])
                    
                    if barcode_check(seq1,seq2,opts.overseq1,opts.overseq2)==True:
                    
                        if read_handle in R1_handles:
                            R1_handles[read_handle].write("%s\n" % "\n".join([head1,seq1,placeholder1,qual1]))
                            R2_handles[read_handle].write("%s\n" % "\n".join([head2,seq2,placeholder2,qual2]))
                        else:
                            R1_handles[read_handle]=open('%s_R1.fastq' % read_handle,'w')     
                            R2_handles[read_handle]=open('%s_R2.fastq' % read_handle,'w')
                            R1_handles[read_handle].write("%s\n" % "\n".join([head1,seq1,placeholder1,qual1]))
                            R2_handles[read_handle].write("%s\n" % "\n".join([head2,seq2,placeholder2,qual2]))
                        
                            reads_retained+=retainadd     
                    else:
                        if (barcode_index_dict[(barcode_found,index_from_read1)],barcode_found,index_from_read1) in barcode_drops:
                            barcode_drops[(barcode_index_dict[(barcode_found,index_from_read1)],barcode_found,index_from_read1)]+=1
                        else:
                            barcode_drops[(barcode_index_dict[(barcode_found,index_from_read1)],barcode_found,index_from_read1)]=1              


logfile=open(opts.logout,'w')
logsummary=open('recovered_summary_'+opts.logout,'w')
logsummary.write('sampleID,IlluminaIdx,Barcode,Count\n')

logfile.write("total read pairs: %s\n" %total_reads)
logfile.write("total reads retained (barcoded and good RAD-tags): %s\n " % reads_retained)

logfile.write("\n\n\n### SAMPLES RECOVERED GOOD+BAD RADS ###\n")
counter_keys=counter_dict.keys()
counter_keys.sort()
for key in counter_keys:
    logfile.write("%s\t%s\t%s\n" % (str(barcode_index_dict[key]),key,counter_dict[key]))

unfound_keys=unfound_dict.keys()
unfound_keys.sort()
logfile.write("\n\n\n### UNMATCHING BARCODES ###\n")
for key in unfound_keys:
    logfile.write("%s\t%s\n" % (key,unfound_dict[key]))
    
logfile.write("\n\n### SAMPLES WITH BAD RADTAGS###\n")
drop_keys=barcode_drops.keys()
drop_keys.sort()

for key in drop_keys:
    logfile.write("%s,%s\n" % (str(key),barcode_drops[key])) 

logfile.write("total bad barcode read pairs: %s\n" % sum(barcode_drops.values()))   

logfile.write("\n\n###RECOVERED READS WITH GOOD RAD-TAGS###\n")

for dropkey in drop_keys:
    counter_dict[tuple(dropkey[1:3])]=counter_dict[tuple(dropkey[1:3])]-barcode_drops[dropkey]    

for retainkey in counter_keys:
    logfile.write('%s,%s,%s\n' % (str(barcode_index_dict[retainkey]),retainkey,counter_dict[retainkey]))
    logsummary.write('%s,%s,%s,%s\n' % (str(barcode_index_dict[retainkey]),retainkey[1],retainkey[0],counter_dict[retainkey]))           

percent_passing=100*sum(counter_dict.values())/float(total_reads)


logfile.write("%s of %s = %s \n" % (sum(counter_dict.values()),total_reads,percent_passing))

for r1_file in R1_handles.keys():
    R1_handles[r1_file].close()
    
for r2_file in R2_handles.keys():
    R2_handles[r2_file].close()
    
    
    
    
    
    
