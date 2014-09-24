# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

from shared_imports import *
from shared_functions import *

# <codecell>

experiments = OrderedDict([
    ('expt01_workflow_1_8_4', OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt01_workflow_1_8_4'), ('ylim_hist', 1000)])),
    ('expt01_workflow_1_9',   OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt01_workflow_1_9'),   ('ylim_hist', 1000)])),
    ('expt01_workflow_1_9_1', OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt01_workflow_1_9_1'), ('ylim_hist', 1000)])),
    ('expt02_workflow_1_9_1', OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt02_workflow_1_9_1'), ('ylim_hist', 100)])),
    ('expt03_workflow_1_9_1', OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt03_workflow_1_9_1'), ('ylim_hist', 100)]))
])
processed_reads_dir_format = '/data/minion/PROCESSED_DATA/{machine}/{directory}'
output_dir_format = '/data/minion/work/{directory}'
fast5stats_output_format = output_dir_format + '/fast5stats/{sample_name}.fast5stats'

# <codecell>

reference_genomes = OrderedDict([
    ('ont_lambda', '/data/minion/reference_genomes/lambda/lambda_ref.fasta'),
    ('cs', '/data/minion/reference_genomes/cs/dna_cs.fasta'),
    ('3D7_V3', '/data/minion/reference_genomes/3D7_V3/3D7_V3.fasta'),
    ('AgamP3', '/data/minion/reference_genomes/AgamP3/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fasta')
])
    

# <codecell>

read_types = ['template', 'complement', '2D']

# <codecell>

import brewer2mpl
set1 = brewer2mpl.get_map('Set1', 'qualitative', 8).mpl_colors
set3 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors

# <headingcell level=1>

# Plan

# <markdowncell>

# - For each read get
#     - 
# - Summary of logs from https://nanoporetech.atlassian.net/wiki/download/attachments/2983444/summarise_ont_logs.py.txt?version=1&modificationDate=1402616251617&api=v2
# - Extract template, complement and 2D reads to fastq
# - Ensure can get last running on fastq and outputting base qualities. Also, any way of getting mapping qualities? And does last give multiple mappings for any reads?
# - Create report of number of reads mapped for template, complement and 2D
# - Figure out how to get a sequence identity metric for each read
# - Run last with multiple different parameters
# - For each set of parameters, calculate number of fast5 files, number base called, %2D, %template only, %reads mapping, mean identity of mapped reads, total read length, number of variants (SNPs and indels) after GATK calling (UG and HC), best perfect kmer, median best perfect kmer
# - Best perfect kmer per read distribution
# - Also try using BLAST via Biopython
# - 
# - Is there anything in big file that is not in individual reads files?
# - Are there time stamps in fast5 reads files? If so, change fast5stats to write out times. Also, maybe include num of events and possibly other data
# 
# 

# <headingcell level=1>

# Functions

# <codecell>

ids = OrderedDict([
        ('2D basecalling successful'                  , 'Basecalling 2D hairpin data'),
        ('Missing Complement data'                    , 'Splitter indicates no complement data. Workflow cannot continue'),
        ('Missing Template data'                      , 'Splitter indicates no template data. Workflow cannot continue'),
        ('Too few events'                             , 'Read has too few events. Workflow cannot continue'),
        ('Poor ratio of template to complement data'  , 'Ratio of template to complement data is no good. Workflow cannot continue'),
        ('Inconsistent template vs complement events' , 'Number of template events out of range aligning template'),
        ('Inconsistent complement vs template events' , 'Number of complement events out of range aligning'),
        ('Failed template to complement alignment'    , 'Template to complement alignment has failed'),
        ('Base calling script score scaling error'    , 'Basecall score increased during scaling'),
        ('Base calling script IndexError'             , 'IndexError: index out of bounds'),
        ('Base calling script ValueError'             , 'ValueError')
])

def classify_read_from_log(log_filename='/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.log'):
    if log_filename.endswith('.log') and os.path.exists(log_filename):
        f = open(log_filename,'r')
        filecontents = f.read()
        flag=0
        for i,v in ids.iteritems():
            if v in filecontents:
                f.close()
                return i
        f.close()
        return None
    else:
        return None
    
def classify_reads_from_logs(processed_reads_dir):
    freqs=OrderedDict()
    for i,v in ids.iteritems():
            freqs[i]=0
    
    number_of_files=0
    for filename in os.listdir(processed_reads_dir):
        filename = os.path.join(processed_reads_dir, filename)
        if filename.endswith('.log'):
            number_of_files = number_of_files + 1
            f = open(filename,'r')
            filecontents = f.read()
            flag=0
            for i,v in ids.iteritems():
                if v in filecontents:
                            freqs[i]=freqs[i]+1
                            flag=1
            if flag==0:
                print filename
    total=0
    for i,v in freqs.iteritems():
            total = total + v
    for i,v in freqs.iteritems():
            percent = round(100*float(v)/float(total),2)
            print "%50s %4i (%2.0f percent)" % (i,v,percent)
    
    print "%50s %i" % ('Total', total)
    print "%50s %i" % ('Number of files', number_of_files)
    
    return freqs

classify_read_from_log()

# <codecell>

def blast_summary(fastq_record, hit_def='Enterobacteria phage lambda, complete genome'):
    result_handle = NCBIWWW.qblast("blastn", "nt", fastq_record.seq)
    blast_record = NCBIXML.read(result_handle)
    
    query_cover = 0
    hsp_identities = 0
    hsp_query_lengths = 0
    
    for alignment in blast_record.alignments:
        if alignment.hit_def == hit_def:
            for hsp in alignment.hsps:
                query_cover = query_cover + (hsp.query_end - hsp.query_start)
                hsp_identities = hsp_identities + hsp.identities
                hsp_query_lengths = hsp_query_lengths + len(hsp.query)
    if hsp_query_lengths == 0:
        return (0.0, 0.0)
    else:
        return (query_cover*1.0/len(fastq_record.seq), hsp_identities*1.0/hsp_query_lengths)

# hdf = h5py.File('/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5', 'r')
# fastq_record = SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'][()]), 'fastq')
# blast_summary(fastq_record)

# <codecell>

def read_summary_from_fast5(fast5_filename='/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5'):
#     print '.',
    read_summary_dtype = [
        ('fast5_filename', 'a200'),
        ('read_type_from_log', 'a100'),
        ('channel_number', np.int),
        ('read_number', np.int),
        ('device_id', 'S100'),
        ('asic_id', 'S100'),
        ('flow_cell_id', 'S100'),
        ('asic_temp', np.float),
        ('heatsink_temp', np.float),
        ('exp_start_time', np.int),
        ('version_name', 'S100'),
        ('run_id', 'S100'),
        ('is_analysed', np.bool),
        ('is_base_called_2d', np.bool),
        ('has_template_read', np.bool),
        ('has_complement_read', np.bool),
        ('has_2D_read', np.bool),
        ('read_length_template', np.int),
        ('read_length_complement', np.int),
        ('read_length_2D', np.int),
        ('start_time_read', np.float),
        ('duration_read', np.float),
        ('start_time_template', np.float),
        ('duration_template', np.float),
        ('start_time_complement', np.float),
        ('duration_complement', np.float),
        ('basecall_called_events_template', np.int),
        ('basecall_num_events_template', np.int),
        ('basecall_num_skips_template', np.int),
        ('basecall_num_stays_template', np.int),
        ('basecall_mean_qscore_template', np.float),
        ('basecall_strand_score_template', np.float),
        ('basecall_called_events_complement', np.int),
        ('basecall_num_events_complement', np.int),
        ('basecall_num_skips_complement', np.int),
        ('basecall_num_stays_complement', np.int),
        ('basecall_mean_qscore_complement', np.float),
        ('basecall_strand_score_complement', np.float),
        ('basecall_mean_qscore_2D', np.float),
        ('basecall_sequence_length_2D', np.int),
        ('basecall_num_raw_events_template', np.int),
        ('basecall_num_merged_events_template', np.int),
        ('basecall_num_raw_events_complement', np.int),
        ('basecall_num_merged_events_complement', np.int),
    ]
    hdf = h5py.File(fast5_filename, 'r')
    
    channel_number =  hdf['UniqueGlobalKey/read_id'].attrs['channel_number']
    read_number =  hdf['UniqueGlobalKey/read_id'].attrs['read_number']
    device_id =  hdf['UniqueGlobalKey/tracking_id'].attrs['device_id']
    asic_id =  hdf['UniqueGlobalKey/tracking_id'].attrs['asic_id']
    flow_cell_id =  hdf['UniqueGlobalKey/tracking_id'].attrs['flow_cell_id']
    asic_temp =  hdf['UniqueGlobalKey/tracking_id'].attrs['asic_temp']
    heatsink_temp =  hdf['UniqueGlobalKey/tracking_id'].attrs['heatsink_temp']
    exp_start_time =  hdf['UniqueGlobalKey/tracking_id'].attrs['exp_start_time']
    version_name =  hdf['UniqueGlobalKey/tracking_id'].attrs['version_name']
    run_id =  hdf['UniqueGlobalKey/tracking_id'].attrs['run_id']
    log_filename = fast5_filename.replace('.fast5', '.log')
    read_type_from_log = classify_read_from_log(log_filename)
    is_analysed = 'Analyses' in hdf.keys()
    is_base_called_2d = False if is_analysed==False else 'Basecall_2D_000' in hdf['Analyses'].keys()
    has_template_read = False if is_base_called_2d==False else 'BaseCalled_template' in hdf['Analyses/Basecall_2D_000'].keys()
    has_complement_read = False if is_base_called_2d==False else 'BaseCalled_complement' in hdf['Analyses/Basecall_2D_000'].keys()
    has_2D_read = False if is_base_called_2d==False else 'BaseCalled_2D' in hdf['Analyses/Basecall_2D_000'].keys()
    read_length_template = 0 if has_template_read==False else len(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_template/Fastq'][()]), 'fastq').seq)
    read_length_complement = 0 if has_complement_read==False else len(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Fastq'][()]), 'fastq').seq)
    read_length_2D = 0 if has_2D_read==False else len(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'][()]), 'fastq').seq)
    start_time_read = hdf['Analyses/EventDetection_000/Reads/Read_%s/Events' % read_number].attrs['start_time']
    duration_read = hdf['Analyses/EventDetection_000/Reads/Read_%s/Events' % read_number].attrs['duration']
    start_time_template = 0.0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/BaseCalled_template/Events'].attrs['start_time']
    duration_template = 0.0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/BaseCalled_template/Events'].attrs['duration']
    start_time_complement = 0.0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Events'].attrs['start_time']
    duration_complement = 0.0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Events'].attrs['duration']

    basecall_called_events_template = 0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['called_events']
    basecall_num_events_template = 0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['num_events']
    basecall_num_skips_template = 0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['num_skips']
    basecall_num_stays_template = 0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['num_stays']
    basecall_mean_qscore_template = 0.0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['mean_qscore']
    basecall_strand_score_template = 0.0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_template'].attrs['strand_score']
    basecall_called_events_complement = 0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['called_events']
    basecall_num_events_complement = 0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['num_events']
    basecall_num_skips_complement = 0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['num_skips']
    basecall_num_stays_complement = 0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['num_stays']
    basecall_mean_qscore_complement = 0.0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['mean_qscore']
    basecall_strand_score_complement = 0.0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_1d_complement'].attrs['strand_score']
    basecall_mean_qscore_2D = 0.0 if has_2D_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_2d'].attrs['mean_qscore']
    basecall_sequence_length_2D = 0 if has_2D_read==False else hdf['Analyses/Basecall_2D_000/Summary/basecall_2d'].attrs['sequence_length']
    basecall_num_raw_events_template = 0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/post_processing_template'].attrs['num_raw_events']
    basecall_num_merged_events_template = 0 if has_template_read==False else hdf['Analyses/Basecall_2D_000/Summary/post_processing_template'].attrs['num_merged_events']
    basecall_num_raw_events_complement = 0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/post_processing_complement'].attrs['num_raw_events']
    basecall_num_merged_events_complement = 0 if has_complement_read==False else hdf['Analyses/Basecall_2D_000/Summary/post_processing_complement'].attrs['num_merged_events']
    
    read_summary = np.array(
        (
            fast5_filename,
            read_type_from_log,
            channel_number,
            read_number,
            device_id,
            asic_id,
            flow_cell_id,
            asic_temp,
            heatsink_temp,
            exp_start_time,
            version_name,
            run_id,
            is_analysed,
            is_base_called_2d,
            has_template_read,
            has_complement_read,
            has_2D_read,
            read_length_template,
            read_length_complement,
            read_length_2D,
            start_time_read,
            duration_read,
            start_time_template,
            duration_template,
            start_time_complement,
            duration_complement,
            basecall_called_events_template,
            basecall_num_events_template,
            basecall_num_skips_template,
            basecall_num_stays_template,
            basecall_mean_qscore_template,
            basecall_strand_score_template,
            basecall_called_events_complement,
            basecall_num_events_complement,
            basecall_num_skips_complement,
            basecall_num_stays_complement,
            basecall_mean_qscore_complement,
            basecall_strand_score_complement,
            basecall_mean_qscore_2D,
            basecall_sequence_length_2D,
            basecall_num_raw_events_template,
            basecall_num_merged_events_template,
            basecall_num_raw_events_complement,
            basecall_num_merged_events_complement,
        ),
        dtype=read_summary_dtype
    )

    hdf.close()
    
    return read_summary

print read_summary_from_fast5()
# No 2D data
print read_summary_from_fast5('/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch101_file13_strand.fast5')
# RAW fast5 file
# print read_summary_from_fast5('/data/minion/RAW_DATA/KWIAT50/E01/reads_expt01/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5')
    

# <codecell>

def blast_summary_from_fast5(fast5_filename='/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5'):
    print '.',
    blast_summary_dtype = [
        ('blast_query_cover_template', np.float),
        ('blast_identity_template', np.float),
        ('blast_query_cover_complement', np.float),
        ('blast_identity_complement', np.float),
        ('blast_query_cover_2D', np.float),
        ('blast_identity_2D', np.float),
    ]
    hdf = h5py.File(fast5_filename, 'r')
    
    (blast_query_cover_template, blast_identity_template) = (0.0, 0.0) if has_template_read==False else blast_summary(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_template/Fastq'][()]), 'fastq'))
    (blast_query_cover_complement, blast_identity_complement) = (0.0, 0.0) if has_complement_read==False else blast_summary(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Fastq'][()]), 'fastq'))
    (blast_query_cover_2D, blast_identity_2D) = (0.0, 0.0) if has_2D_read==False else blast_summary(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'][()]), 'fastq'))
    
    blast_summary = np.array(
        (
            blast_query_cover_template,
            blast_identity_template,
            blast_query_cover_complement,
            blast_identity_complement,
            blast_query_cover_2D,
            blast_identity_2D,
        ),
        dtype=blast_summary_dtype
    )

    hdf.close()
    
    return blast_summary

# print read_summary_from_fast5()
# # No 2D data
# print read_summary_from_fast5('/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch101_file13_strand.fast5')
# # RAW fast5 file
# # print read_summary_from_fast5('/data/minion/RAW_DATA/KWIAT50/E01/reads_expt01/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5')
    

# <codecell>

def determine_read_summaries(processed_reads_dir='/data/minion/PROCESSED_DATA/KWIAT50/expt03_workflow_1_9_1'):
    filenames = list()
    for filename in os.listdir(processed_reads_dir):
        filename = os.path.join(processed_reads_dir, filename)
        if filename.endswith('.fast5'):
            filenames.append(filename)
    print len(filenames)
    fast5_summaries_list = [read_summary_from_fast5(filename) for filename in filenames]
    fast5_summaries_array = np.hstack(fast5_summaries_list).view(recarray)
    return fast5_summaries_array

# temp = determine_read_summaries()

# <codecell>

temp['read_length_template'][20]

# <codecell>

hdf = h5py.File('/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5', 'r')

# <codecell>

hdf['UniqueGlobalKey/read_id'].attrs['read_number']

# <headingcell level=1>

# QC fast5 files

# <codecell>

processed_reads_dirs = OrderedDict()
read_summaries = OrderedDict()

for experiment in experiments:
    print experiment
    processed_reads_dirs[experiment] = processed_reads_dir_format.format(
        machine=experiments[experiment]['machine'],
        directory=experiments[experiment]['directory']
    )
    read_summaries[experiment] = determine_read_summaries(processed_reads_dirs[experiment])

# <codecell>


# <headingcell level=1>

# Create fastq files

# <codecell>

for experiment in experiments:
    if not os.path.exists(os.path.join(output_dirs[experiment], 'fastq')):
        !mkdir {output_dirs[experiment]}/fastq
    if not os.path.exists(os.path.join(output_dirs[experiment], 'fasta')):
        !mkdir {output_dirs[experiment]}/fasta
    for read_type in read_types:
        print experiment, read_type
        if not os.path.exists(os.path.join(output_dirs[experiment], 'fastq', 'reads_%s.fastq' % read_type)):
            !../alignment/fast5extractbytype.py {processed_reads_dirs[experiment]} {read_type} fastq {output_dirs[experiment]}/fastq
        if not os.path.exists(os.path.join(output_dirs[experiment], 'fasta', 'reads_%s.fasta' % read_type)):
            !../alignment/fast5extractbytype.py {processed_reads_dirs[experiment]} {read_type} fasta {output_dirs[experiment]}/fasta

# <headingcell level=1>

# Map reads with LAST

# <codecell>

for experiment in experiments:
    for reference_genome in reference_genomes:
        if not os.path.exists(os.path.join(output_dirs[experiment], 'last', reference_genome)):
            !mkdir {output_dirs[experiment]}/last/{reference_genome}
        for read_type in read_types:
            print experiment, read_type, reference_genome
            !../alignment/ont_mapwithlast.sh {output_dirs[experiment]}/fasta/reads_{read_type}.fasta {reference_genomes[reference_genome]} {output_dirs[experiment]}/last/{reference_genome}

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>

read_summaries

# <codecell>

processed_reads_dirs = OrderedDict()
output_dirs = OrderedDict()
sample_names = OrderedDict()
fast5stats_output_filenames = OrderedDict()
fast5stats_tables = OrderedDict()
fast5stats_arrays = OrderedDict()

for experiment in experiments:
    print experiment
    processed_reads_dirs[experiment] = processed_reads_dir_format.format(
        machine=experiments[experiment]['machine'],
        directory=experiments[experiment]['directory']
    )
    output_dirs[experiment] = output_dir_format.format(
        directory=experiments[experiment]['directory']
    )
    sample_names[experiment] = experiments[experiment]['directory']
    fast5stats_output_filenames[experiment] = fast5stats_output_format.format(
        directory=experiments[experiment]['directory'],
        sample_name = sample_names[experiment]
    )
    !mkdir -p {processed_reads_dirs[experiment]}
    !mkdir -p {output_dirs[experiment]}
    !mkdir -p {os.path.dirname(fast5stats_output_filenames[experiment])}
    
    if not os.path.exists(fast5stats_output_filenames[experiment]):
        !../QC/fast5stats.py {processed_reads_dirs[experiment]} {sample_names[experiment]} {fast5stats_output_filenames[experiment]}

    fast5stats_tables[experiment] = fromtsv(fast5stats_output_filenames[experiment]).convertnumbers()
    fast5stats_arrays[experiment] = toarray(fast5stats_tables[experiment])
    print shape(fast5stats_arrays[experiment])

# sample_name = os.path.basename(processed_reads_dir)
# fast5stats_output_fn = '%s/fast5stats/%s.fast5stats' % (outputs_dir, sample_name)

# <headingcell level=1>

# BLAST exploration

# <codecell>

from Bio.Blast import NCBIWWW
from Bio import SeqIO

hdf = h5py.File('/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5', 'r')
fastq_record = SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'][()]), 'fastq')
print fastq_record

# record = SeqIO.read("m_cold.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", fastq_record.seq)

# <codecell>

result_handle = NCBIWWW.qblast("blastn", "nt", fastq_record.seq)

# <codecell>

str(fastq_record.seq)

# <codecell>

result_handle

# <codecell>

from Bio.Blast import NCBIXML
blast_record = NCBIXML.read(result_handle)

# <codecell>

blast_record

# <codecell>

print blast_record.alignments[0].title
print blast_record.alignments[0].hit_id
print blast_record.alignments[0].hit_def
print len(blast_record.alignments[0].hsps)
print blast_record.alignments[0].length

# <codecell>

help(NCBIWWW.qblast)

# <codecell>

def blast

# <codecell>

len(blast_record.alignments)

# <codecell>

blast_record.alignments[0].length

# <codecell>

blast_summary(fastq_record)

# <codecell>

>>> for alignment in blast_record.alignments:
...     print len(alignment.hsps)
...     for hsp in alignment.hsps:
...         if hsp.expect < 0.01:
...             print('****Alignment****')
...             print('sequence:', alignment.title)
# ...             print('length:', alignment.length)
# ...             print('e value:', hsp.expect)
# ...             print(hsp.query[0:75] + '...')
# ...             print(hsp.match[0:75] + '...')
# ...             print(hsp.sbjct[0:75] + '...')
...             print(hsp.identities)
# ...             print(hsp.score)
# ...             print(hsp.bits)
# ...             print(hsp.expect)
# ...             print(hsp.num_alignments)
# ...             print(hsp.positives)
# ...             print(hsp.gaps)
# ...             print(hsp.strand)
# ...             print(hsp.frame)
...             print(len(hsp.query))
# ...             print(hsp.query)
...             print(hsp.query_start)
...             print(hsp.query_end)
# ...             print(hsp.match)
...             print
...             print(hsp)
...             print



# <codecell>


# <codecell>

fast5stats_tables['expt01_workflow_1_9_1']

# <codecell>

reads_from_logs_breakdowns = OrderedDict()

for experiment in experiments:
    print experiment
    reads_from_logs_breakdowns[experiment] = classify_reads_from_logs(processed_reads_dirs[experiment])
    print

# <codecell>

reads_from_logs_breakdowns['expt01_workflow_1_8_4']

# <codecell>

for experiment in reads_from_logs_breakdowns:
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(1, 1, 1)
    for i, basecall_outcome in enumerate(reads_from_logs_breakdowns[experiment]):
        ax.bar(i, reads_from_logs_breakdowns[experiment][basecall_outcome], color=set3[i], label=basecall_outcome, align='center')
    ax.set_xticks(range(len(reads_from_logs_breakdowns[experiment])))
    ax.set_xticklabels(reads_from_logs_breakdowns[experiment].keys(), rotation='vertical', horizontalalignment='center')
    ax.set_xlim(-1, 11)
    ax.set_title(experiment)
    fig.savefig(os.path.join(
        output_dir_format.format(
            directory=experiments[experiment]['directory']
        ),
        "reads_from_logs_breakdown_%s.pdf" % experiment
    ))

# <codecell>

for experiment in experiments:
    fig = plt.figure(figsize=(16, 8))
    for i, read_type in enumerate(read_types):
        print experiment, read_type, np.count_nonzero(fast5stats_arrays[experiment]['readtype'] == read_type), "reads, ", np.sum(fast5stats_arrays[experiment]['seqlen'][(fast5stats_arrays[experiment]['readtype'] == read_type)]), " bases"
        ax = fig.add_subplot(len(read_types), 1, i + 1)
        ax.hist(
                fast5stats_arrays[experiment]['seqlen'][(fast5stats_arrays[experiment]['readtype'] == read_type)],
                color=set1[i],
                bins=arange(0, 18000, 500),
                label=read_type)
        ax.set_xlim(0, 18000)
        ax.set_ylim(0, experiments[experiment]['ylim_hist'])
        if i < (len(read_types) - 1):
            ax.set_xticks([])
        ax.set_title(experiment + ' - ' + read_type)
    ax.set_xlabel('Sequence length')
# fig.savefig(PLOT_DIR + '/fig_s8_human_contamination_by_wbc_depletion.pdf')


# <codecell>

seqlen_per_channel = OrderedDict()

for experiment in experiments:
    seqlen_per_channel[experiment] = np.zeros(512)
    for i in arange(512):
        seqlen_per_channel[experiment][i] = np.sum(
            fast5stats_arrays[experiment]['seqlen'][(fast5stats_arrays[experiment]['channel'] == i+1)]
        )

# <codecell>

for experiment in experiments:
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(
            seqlen_per_channel[experiment],
            bins=arange(0, 400000, 10000),
            label=read_type)
    ax.set_xlim(0, 400000)
    ax.set_title(experiment)
    ax.set_xlabel('Total sequence length from pore')
    ax.set_ylabel('Frequency (number of pores)')

# <codecell>

hist(seqlen_per_channel, bins=range(0, 400000, 10000))

# <headingcell level=1>

# Create fastq files

# <codecell>

for experiment in experiments:
    if not os.path.exists(os.path.join(output_dirs[experiment], 'fastq')):
        !mkdir {output_dirs[experiment]}/fastq
    if not os.path.exists(os.path.join(output_dirs[experiment], 'fasta')):
        !mkdir {output_dirs[experiment]}/fasta
    for read_type in read_types:
        print experiment, read_type
        if not os.path.exists(os.path.join(output_dirs[experiment], 'fastq', 'reads_%s.fastq' % read_type)):
            !../alignment/fast5extractbytype.py {processed_reads_dirs[experiment]} {read_type} fastq {output_dirs[experiment]}/fastq
        if not os.path.exists(os.path.join(output_dirs[experiment], 'fasta', 'reads_%s.fasta' % read_type)):
            !../alignment/fast5extractbytype.py {processed_reads_dirs[experiment]} {read_type} fasta {output_dirs[experiment]}/fasta

# <headingcell level=1>

# Map reads with LAST

# <codecell>

for experiment in experiments:
    for reference_genome in reference_genomes:
        if not os.path.exists(os.path.join(output_dirs[experiment], 'last', reference_genome)):
            !mkdir {output_dirs[experiment]}/last/{reference_genome}
        for read_type in read_types:
            print experiment, read_type, reference_genome
            !../alignment/ont_mapwithlast.sh {output_dirs[experiment]}/fasta/reads_{read_type}.fasta {reference_genomes[reference_genome]} {output_dirs[experiment]}/last/{reference_genome}

# <headingcell level=1>

# Exploration of LAST bam files

# <codecell>

import pysam
samfile = pysam.Samfile( "/data/minion/work/expt01_workflow_1_8_4/last/reads_2D.last.sorted.bam", "rb" )

# <codecell>

temp = samfile.next()

# <codecell>

all_reads = list(samfile.fetch())

# <codecell>

mapped_read_summaries = OrderedDict()

for experiment in experiments:
    mapped_read_summaries[experiment] = OrderedDict()
    for reference_genome in reference_genomes:
        mapped_read_summaries[experiment][reference_genome] = OrderedDict()
        for read_type in read_types:
            print experiment, reference_genome, read_type
            samfile = pysam.Samfile("/data/minion/work/%s/last/%s/reads_%s.last.sorted.bam" %(experiment, reference_genome, read_type), "rb" )
            all_reads = list(samfile.fetch())
            samfile.close()
            mapped_read_summaries[experiment][reference_genome][read_type] = summarise_bam_reads(all_reads)

# <codecell>

full_read_summaries = OrderedDict()
read_summary_tables = OrderedDict()
mapped_read_intial_tables = OrderedDict()
read_summary_initial_tables = OrderedDict()
read_summary_intermediate_tables = OrderedDict()
read_summary_intermediate_tables_2 = OrderedDict()

for experiment in experiments:
    print experiment
    read_summary_table = (etl
        .fromarray(read_summaries[experiment])
    )
    read_summary_initial_tables[experiment] = read_summary_table
    mapped_read_intial_tables[experiment] = OrderedDict()
    read_summary_intermediate_tables[experiment] = OrderedDict()
    read_summary_intermediate_tables_2[experiment] = OrderedDict()
    for reference_genome in reference_genomes:
        mapped_read_intial_tables[experiment][reference_genome] = OrderedDict()
        read_summary_intermediate_tables[experiment][reference_genome] = OrderedDict()
        read_summary_intermediate_tables_2[experiment][reference_genome] = OrderedDict()
        for read_type in read_types:
            print experiment, reference_genome, read_type
            mapped_read_table = (etl
                .fromarray(mapped_read_summaries[experiment][reference_genome][read_type])
                .rename(
                    dict(
                        zip(
                            etl.fromarray(mapped_read_summaries[experiment][reference_genome][read_type]).header(),
                            map(
                                lambda x: x + '_%s_%s' % (reference_genome, read_type),
                                etl.fromarray(mapped_read_summaries[experiment][reference_genome][read_type]).header()
                            )
                        )
                    )
                )
            )
            mapped_read_intial_tables[experiment][reference_genome][read_type] = mapped_read_table
#             read_summary_intermediate_tables[experiment][reference_genome][read_type] = mapped_read_table
            read_name_field_name = 'read_names_%s_%s' % (reference_genome, read_type)
            read_name_field_format = 'channel_%%d_read_%%d_%s' % (read_type)
            print read_name_field_name, read_name_field_format
            read_summary_intermediate_tables[experiment][reference_genome][read_type] = (read_summary_table
                .addfield(read_name_field_name, lambda(rec): read_name_field_format % (rec['channel_number'], rec['read_number']))
            )
            read_summary_table = (read_summary_table
                .addfield(read_name_field_name, lambda(rec): read_name_field_format % (rec['channel_number'], rec['read_number']))  
                .leftjoin(mapped_read_table)
            )
            read_summary_intermediate_tables_2[experiment][reference_genome][read_type] = read_summary_table


            print experiment, reference_genome, read_type
            samfile = pysam.Samfile("/data/minion/work/%s/last/%s/reads_%s.last.sorted.bam" %(experiment, reference_genome, read_type), "rb" )
            all_reads = list(samfile.fetch())
            samfile.close()
            mapped_read_summaries[experiment][reference_genome][read_type] = summarise_bam_reads(all_reads)
    read_summary_tables[experiment] = read_summary_table

# <codecell>

read_summary_tables.keys()

# <codecell>

read_summary_initial_tables['expt01_workflow_1_9_1'].head(3)

# <codecell>

mapped_read_intial_tables['expt01_workflow_1_9_1']['ont_lambda']['2D'].head(3)

# <codecell>

read_summary_intermediate_tables['expt01_workflow_1_8_4']['ont_lambda']['template'].header()

# <codecell>

mapped_read_intial_tables['expt01_workflow_1_8_4']['ont_lambda']['template'].header()

# <codecell>

(read_summary_intermediate_tables['expt01_workflow_1_8_4']['ont_lambda']['template']
    .leftjoin(mapped_read_intial_tables['expt01_workflow_1_8_4']['ont_lambda']['template'])
).head(3)

# <codecell>

read_summary_intermediate_tables['expt01_workflow_1_8_4']['ont_lambda']['template'].head(3)

# <codecell>

read_summary_initial_tables['expt01_workflow_1_9_1'].head(3)

# <codecell>

read_summary_tables['expt01_workflow_1_9_1']

# <codecell>

(etl
    .fromarray(mapped_read_summaries['expt01_workflow_1_9_1']['ont_lambda']['2D'])
).head(3)

# <codecell>

read_2D_table = (etl
    .fromarray(mapped_read_summaries['expt01_workflow_1_9_1']['ont_lambda']['2D'])
    .rename(dict(zip(etl.fromarray(mapped_read_summaries['expt01_workflow_1_9_1']['ont_lambda']['2D']).header(), map(lambda x: x + '_ont_lambda_2D', etl.fromarray(mapped_read_summaries['expt01_workflow_1_9_1']['ont_lambda']['2D']).header()))))
#     .rename('read_names', 'read_name_2D')
)
read_2D_table.head(3)

# <codecell>

dict(zip(read_2D_table.header(), map(lambda x: x + '_ont_lambda_2D', read_2D_table.header())))

# <codecell>

read_2D_table.header()

# <codecell>

(etl
    .fromarray(read_summaries['expt01_workflow_1_9_1'])
    .addfield('read_name_ont_lambda_template', lambda(rec): 'channel_%d_read_%d_template' % (rec['channel_number'], rec['read_number']))  
    .addfield('read_name_ont_lambda_complement', lambda(rec): 'channel_%d_read_%d_complements' % (rec['channel_number'], rec['read_number']))  
    .addfield('read_names_ont_lambda_2D', lambda(rec): 'channel_%d_read_%d_2D' % (rec['channel_number'], rec['read_number']))
    .leftjoin(read_2D_table)
    .select(lambda rec: rec['fragments_ont_lambda_2D'] > 0)
).head(3)

# <codecell>

mapped_read_summaries = OrderedDict()

for experiment in experiments:
    print experiment
    processed_reads_dirs[experiment] = processed_reads_dir_format.format(
        machine=experiments[experiment]['machine'],
        directory=experiments[experiment]['directory']
    )
    read_summaries[experiment] = determine_read_summaries(processed_reads_dirs[experiment])

# <codecell>

def summarise_bam_reads(all_reads):
    qnames = np.array(map(lambda x: x.qname, all_reads))
    inferred_lengths = np.array(map(lambda x: x.inferred_length, all_reads))
    mismatches = np.array(map(lambda x: dict(x.tags)['NM'], all_reads))
    Ns = np.array(map(lambda x: np.sum([y == 'N' for y in x.seq]), all_reads))
    deletions = np.array(map(lambda x: np.sum(map(lambda y: y == None, np.array(x.aligned_pairs)[:, 0])), all_reads))
    insertions = np.array(map(lambda x: np.sum(map(lambda y: y == None, np.array(x.aligned_pairs)[:, 1])), all_reads))

    read_names = np.unique( qnames )
    fragments_per_read = np.array([np.sum(qnames == x) for x in read_names])
    total_fragment_length = np.array([np.sum(inferred_lengths[qnames == x]) for x in read_names])
    mean_fragment_length = np.array([np.mean(inferred_lengths[qnames == x]) for x in read_names])
    min_fragment_length = np.array([np.min(inferred_lengths[qnames == x]) for x in read_names])
    max_fragment_length = np.array([np.max(inferred_lengths[qnames == x]) for x in read_names])
    mismatches_per_read = np.array([np.sum(mismatches[qnames == x]) for x in read_names])
    Ns_per_read = np.array([np.sum(Ns[qnames == x]) for x in read_names])
    deletions_per_read = np.array([np.sum(deletions[qnames == x]) for x in read_names])
    insertions_per_read = np.array([np.sum(insertions[qnames == x]) for x in read_names])
    
    read_summaries = np.core.records.fromarrays(
        [
            read_names,
            fragments_per_read,
            total_fragment_length,
            mean_fragment_length,
            min_fragment_length,
            max_fragment_length,
            mismatches_per_read,
            Ns_per_read,
            deletions_per_read,
            insertions_per_read
        ],
        names='read_names,fragments,total_fragment_length,mean_fragment_length,min_fragment_length,max_fragment_length,mismatches,Ns,deletions,insertions'
    )
    
    return read_summaries

temp_read_summaries = summarise_bam_reads(all_reads)
    

# <codecell>

shape(temp_read_summaries)

# <codecell>


# <codecell>

temp_read_summaries['read_names']

# <codecell>

temp_read_table = etl.fromarray(temp_read_summaries)

# <codecell>

temp_read_table.head()

# <codecell>

z = [np.mean(x[exps == i]) for i in np.unique( exps )]

# <codecell>

qnames = np.array(map(lambda x: x.qname, all_reads))
inferred_lengths = np.array(map(lambda x: x.inferred_length, all_reads))
mismatches = np.array(map(lambda x: dict(x.tags)['NM'], all_reads))
Ns = np.array(map(lambda x: np.sum([y == 'N' for y in x.seq]), all_reads))
deletions = np.array(map(lambda x: np.sum(map(lambda y: y == None, np.array(x.aligned_pairs)[:, 0])), all_reads))
insertions = np.array(map(lambda x: np.sum(map(lambda y: y == None, np.array(x.aligned_pairs)[:, 1])), all_reads))

# seq_lengths = np.array(map(lambda x: len(x.seq), all_reads))
print shape(qnames)
print shape(inferred_lengths)
print shape(mismatches)
print shape(seq_lengths)
print shape(Ns)
# print seq_lengths[0:30]
print inferred_lengths[0:30]
print mismatches[0:30]
print Ns[0:30]
print deletions[0:30]
print insertions[0:30]
print (Ns+deletions+insertions)[0:30]
print (mismatches-(Ns+deletions+insertions))[qnames == 'channel_281_read_62_2D']

# <codecell>

read_names = np.unique( qnames )
mean_fragment_length = np.array([np.mean(inferred_lengths[qnames == x]) for x in read_names])
max_fragment_length = np.array([np.max(inferred_lengths[qnames == x]) for x in read_names])
min_fragment_length = np.array([np.min(inferred_lengths[qnames == x]) for x in read_names])
total_fragment_length = np.array([np.sum(inferred_lengths[qnames == x]) for x in read_names])
fragments_per_read = np.array([np.sum(qnames == x) for x in read_names])
mismatches_per_read = np.array([np.sum(mismatches[qnames == x]) for x in read_names])
Ns_per_read = np.array([np.sum(Ns[qnames == x]) for x in read_names])
deletions_per_read = np.array([np.sum(deletions[qnames == x]) for x in read_names])
insertions_per_read = np.array([np.sum(insertions[qnames == x]) for x in read_names])

print shape(inferred_length_per_read)
print shape(fragments_per_read)
print shape(read_names)
print shape(mismatches_per_read)
print total_fragment_length[0:10]
print fragments_per_read[0:10]
print mismatches_per_read[0:10]
print Ns_per_read[0:10]
print deletions_per_read[0:10]
print insertions_per_read[0:10]
print read_names[0:10]
print max_read_length[0:10]
print min_read_length[0:10]
print mean_fragment_length[0:10]

# <codecell>


# <codecell>

inferred_length_per_read

# <codecell>

read_names = sort(np.unique(map(lambda x: x.qname, all_reads)))

# <codecell>

print len(read_names)
print read_names[0:10]

# <codecell>

from pprint import pprint
pprint(temp)
dir(temp)

# <codecell>

# for x in dir(temp):
for x in ['qname', 'aend', 'alen', 'bin', 'inferred_length', 'pos', 'qend', 'qlen', 'rlen', 'tags', 'tlen', 'aligned_pairs']:
    exec 'print \''+x+':\', temp.'+x

# <codecell>

dict(temp.tags)['NM']

# <codecell>

np.sum([x == 'N' for x in temp.seq])

# <codecell>

print np.apply_along_axis(lambda x: np.sum(map(lambda y: y == None, x)), 0, np.array(temp.aligned_pairs))

# <codecell>

print np.sum(map(lambda y: y == None, np.array(temp.aligned_pairs)[:, 0]))
print np.sum(map(lambda y: y == None, np.array(temp.aligned_pairs)[:, 1]))

# <codecell>

temp.qname

# <codecell>

np.array(temp.aligned_pairs)[:, 0]

# <codecell>

def read_summary

# <codecell>

temp = np.array(temp.aligned_pairs)[:, 0]

# <codecell>

np.sum(map(lambda x: x == None, temp))

# <codecell>

print temp.qname
print temp.qname
print temp.qname
print temp.qname
print temp.qname
print temp.qname
print temp.qname
print temp.qname
print temp.qname
print temp.qname

# <codecell>


# <codecell>

print temp

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <headingcell level=1>

# Sandpit

# <codecell>

!../QC/fast5stats.py {processed_reads_dir} {sample_name} {fast5stats_output_fn}

# <codecell>

fast5stats_table = fromtsv(fast5stats_output_fn).convertnumbers()
fast5stats_array = toarray(fast5stats_table)
print shape(fast5stats_array)
fast5stats_table.head()

# <codecell>

tabulate(fast5stats_array['channel'])

# <codecell>

read

# <codecell>


