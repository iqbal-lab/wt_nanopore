# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Setup

# <codecell>

from shared_imports import *
from shared_functions import *

# <codecell>

experiments = OrderedDict([
    ('expt01_workflow_1_8_4',           OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt01_workflow_1_8_4'),           ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('expt01_workflow_1_9',             OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt01_workflow_1_9'),             ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('expt01_workflow_1_9_1',           OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt01_workflow_1_9_1'),           ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('expt02_workflow_1_9_1',           OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt02_workflow_1_9_1'),           ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 100) ])),
    ('expt03_workflow_1_9_1',           OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt03_workflow_1_9_1'),           ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 100) ])),
    ('expt04_workflow_1_9_1',           OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt04_workflow_1_9_1'),           ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 100) ])),
    ('expt05_aliquot_1_workflow_1_9_1', OrderedDict([('machine', 'KWIAT50'), ('directory', 'expt05_aliquot_1_workflow_1_9_1'), ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_3_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT50'), ('directory', 'E05_aliquot_3_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_5_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT50'), ('directory', 'E05_aliquot_5_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_7_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT50'), ('directory', 'E05_aliquot_7_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_2_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT51'), ('directory', 'E05_aliquot_2_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_4_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT51'), ('directory', 'E05_aliquot_4_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_6_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT51'), ('directory', 'E05_aliquot_6_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E05_aliquot_8_workflow_1_9_1',    OrderedDict([('machine', 'KWIAT51'), ('directory', 'E05_aliquot_8_workflow_1_9_1'),    ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E06_ANOPH_workflow_1_9_1',        OrderedDict([('machine', 'KWIAT54'), ('directory', 'E06_ANOPH_workflow_1_9_1'),        ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
#     ('E06_ANOPH_workflow_1_10',         OrderedDict([('machine', 'KWIAT54'), ('directory', 'E06_ANOPH_workflow_1_10'),         ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
    ('E07_ANOPH_workflow_1_9_1.2.1',    OrderedDict([('machine', 'KWIAT53'), ('directory', 'E07_ANOPH_workflow_1_9_1'),        ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
#     ('E07_ANOPH_workflow_1_9_1.2.10',   OrderedDict([('machine', 'KWIAT53'), ('directory', 'E07_ANOPH_workflow_1_9_1'),        ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 10), ('ylim_hist', 1000)])),
#     ('E07_ANOPH_workflow_1_9_1.1.1',    OrderedDict([('machine', 'KWIAT53'), ('directory', 'E07_ANOPH_workflow_1_9_1'),        ('readsfmt', 'fastq'), ('lastal_a', 1), ('lastal_m', 1), ('ylim_hist', 1000)])),
#     ('E07_ANOPH_workflow_1_9_1.3.1',    OrderedDict([('machine', 'KWIAT53'), ('directory', 'E07_ANOPH_workflow_1_9_1'),        ('readsfmt', 'fastq'), ('lastal_a', 3), ('lastal_m', 1), ('ylim_hist', 1000)])),
#     ('E07_ANOPH_workflow_1_10',         OrderedDict([('machine', 'KWIAT53'), ('directory', 'E07_ANOPH_workflow_1_10'),         ('readsfmt', 'fastq'), ('lastal_a', 2), ('lastal_m', 1), ('ylim_hist', 1000)])),
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

# Note I haven't used this as far too slow

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

def gc(seq):
    g_and_c = seq.count('G') + seq.count('C')
    a_and_t = seq.count('A') + seq.count('T')
    return g_and_c * 1.0 / (g_and_c + a_and_t)

# <codecell>

def read_summary_from_fast5(fast5_filename='/data/minion/PROCESSED_DATA/KWIAT50/expt01_workflow_1_9_1/KWIAT50_K50_expt01_std_lambda_DNA_140619_A_4429_1_ch102_file11_strand.fast5'):
#     print fast5_filename
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
        ('gc_template', np.float),
        ('gc_complement', np.float),
        ('gc_2D', np.float),
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
    gc_template = None if has_template_read==False else gc(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_template/Fastq'][()]), 'fastq').seq)
    gc_complement = None if has_complement_read==False else gc(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_complement/Fastq'][()]), 'fastq').seq)
    gc_2D = None if has_2D_read==False else gc(SeqIO.read(StringIO(hdf['Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'][()]), 'fastq').seq)
    start_time_read = 0.0 if is_analysed==False else hdf['Analyses/EventDetection_000/Reads/Read_%s/Events' % read_number].attrs['start_time']
    duration_read = 0.0 if is_analysed==False else hdf['Analyses/EventDetection_000/Reads/Read_%s/Events' % read_number].attrs['duration']
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
            gc_template,
            gc_complement,
            gc_2D,
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

# Note I haven't used this as far too slow

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

def summarise_bam_reads(all_reads):
    qnames = np.array(map(lambda x: x.qname, all_reads))
    inferred_lengths = np.array(map(lambda x: x.inferred_length, all_reads))
    mismatches = np.array(map(lambda x: dict(x.tags)['NM'], all_reads))
    Ns = np.array(map(lambda x: np.sum([y == 'N' for y in x.seq]), all_reads))
    deletions = np.array(map(lambda x: np.sum(map(lambda y: y == None, np.array(x.aligned_pairs)[:, 0])), all_reads))
    insertions = np.array(map(lambda x: np.sum(map(lambda y: y == None, np.array(x.aligned_pairs)[:, 1])), all_reads))
    snps = mismatches - Ns - deletions - insertions
    tids = np.array(map(lambda x: x.tid, all_reads))

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
    snps_per_read = np.array([np.sum(snps[qnames == x]) for x in read_names])
    mapped_chromosomes_per_read = np.array([len(np.unique(tids[qnames == x])) for x in read_names])
    
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
            insertions_per_read,
            snps_per_read,
            mapped_chromosomes_per_read,
        ],
        names='read_names,fragments,total_fragment_length,mean_fragment_length,min_fragment_length,max_fragment_length,mismatches,Ns,deletions,insertions,snps,chromosomes'
    )
    
    return read_summaries

all_reads = list(pysam.Samfile("/data/minion/work/expt01_workflow_1_9_1/last/ont_lambda/reads_2D.2.1.last.sorted.bam", "rb" ).fetch())
temp_read_summaries = summarise_bam_reads(all_reads)

# <codecell>

def convert_time(YOUR_EPOCH_TIME):
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(YOUR_EPOCH_TIME))
convert_time(1403178270)

# <headingcell level=1>

# QC fast5 files

# <codecell>

processed_reads_dirs = OrderedDict()
read_summaries = OrderedDict()
output_dirs = OrderedDict()

for experiment in experiments:
    print experiment
    processed_reads_dirs[experiment] = processed_reads_dir_format.format(
        machine=experiments[experiment]['machine'],
        directory=experiments[experiment]['directory']
    )
    output_dirs[experiment] = output_dir_format.format(
        directory=experiments[experiment]['directory']
    )
    if not os.path.exists(os.path.join(output_dirs[experiment])):
        !mkdir -p {output_dirs[experiment]}
    read_summaries_fn = os.path.join(output_dirs[experiment], "read_summaries.npy")
#     if experiment in read_summaries.keys(): # previously created summary and want to save
#         np.save(read_summaries_fn, read_summaries[experiment])
#     else:
    print read_summaries_fn
    if os.path.exists(read_summaries_fn):
        read_summaries[experiment] = np.load(read_summaries_fn)
    else:
        read_summaries[experiment] = determine_read_summaries(processed_reads_dirs[experiment])
        np.save(read_summaries_fn, read_summaries[experiment])

# <headingcell level=1>

# Create fastq files

# <codecell>

processed_reads_dirs = OrderedDict()
output_dirs = OrderedDict()

for experiment in experiments:
    processed_reads_dirs[experiment] = processed_reads_dir_format.format(
        machine=experiments[experiment]['machine'],
        directory=experiments[experiment]['directory']
    )
    output_dirs[experiment] = output_dir_format.format(
        directory=experiments[experiment]['directory']
    )

    if not os.path.exists(os.path.join(output_dirs[experiment], 'fastq')):
        !mkdir -p {output_dirs[experiment]}/fastq
    if not os.path.exists(os.path.join(output_dirs[experiment], 'fasta')):
        !mkdir -p {output_dirs[experiment]}/fasta
    for read_type in read_types:
        print experiment, read_type
        if not os.path.exists(os.path.join(output_dirs[experiment], 'fastq', 'reads_%s.fastq' % read_type)):
            !../alignment/fast5extractbytype_rdp.py {processed_reads_dirs[experiment]} {read_type} fastq {output_dirs[experiment]}/fastq
        if not os.path.exists(os.path.join(output_dirs[experiment], 'fasta', 'reads_%s.fasta' % read_type)):
            !../alignment/fast5extractbytype_rdp.py {processed_reads_dirs[experiment]} {read_type} fasta {output_dirs[experiment]}/fasta

# <headingcell level=1>

# Map reads with LAST

# <codecell>

for experiment in experiments:
    for reference_genome in reference_genomes:
        if not os.path.exists(os.path.join(output_dirs[experiment], 'last', reference_genome)):
            !mkdir -p {output_dirs[experiment]}/last/{reference_genome}
        for read_type in read_types:
            print experiment, read_type, reference_genome
            if experiments[experiment]['readsfmt'] == 'fasta':
                !../alignment/ont_mapwithlast_rdp.sh {output_dirs[experiment]}/fasta/reads_{read_type}.fasta fasta experiment {reference_genomes[reference_genome]} {output_dirs[experiment]}/last/{reference_genome} {experiments[experiment]['lastal_a']} {experiments[experiment]['lastal_m']}
            elif experiments[experiment]['readsfmt'] == 'fastq':
                !../alignment/ont_mapwithlast_rdp.sh {output_dirs[experiment]}/fastq/reads_{read_type}.fastq fastq experiment {reference_genomes[reference_genome]} {output_dirs[experiment]}/last/{reference_genome} {experiments[experiment]['lastal_a']} {experiments[experiment]['lastal_m']}
            else:
                print "Unknown readsfmt"

# <codecell>

mapped_read_summaries = OrderedDict()

for experiment in experiments:
    mapped_read_summaries[experiment] = OrderedDict()
    for reference_genome in reference_genomes:
        mapped_read_summaries[experiment][reference_genome] = OrderedDict()
        for read_type in read_types:
            lastal_a = experiments[experiment]['lastal_a']
            lastal_m = experiments[experiment]['lastal_m']
            print experiment, reference_genome, read_type, lastal_a, lastal_m
            mapped_read_summaries_fn = os.path.join(output_dirs[experiment], "mapped_read_summaries_%s_%s.%d.%d.npy" % (reference_genome, read_type, lastal_a, lastal_m))
            if os.path.exists(mapped_read_summaries_fn):
                mapped_read_summaries[experiment][reference_genome][read_type] = np.load(mapped_read_summaries_fn)
            else:
                samfile = pysam.Samfile("/data/minion/work/%s/last/%s/reads_%s.%d.%d.last.sorted.bam" %(experiments[experiment]['directory'], reference_genome, read_type, lastal_a, lastal_m), "rb" )
                all_reads = list(samfile.fetch())
                samfile.close()
                mapped_read_summaries[experiment][reference_genome][read_type] = summarise_bam_reads(all_reads)
                np.save(mapped_read_summaries_fn, mapped_read_summaries[experiment][reference_genome][read_type])

# <codecell>

etl.fromarray(mapped_read_summaries['expt01_workflow_1_9_1']['cs']['2D']).head()

# <codecell>

read_summary_tables = OrderedDict()

for experiment in experiments:
    print experiment
    read_summary_table = (etl
        .fromarray(read_summaries[experiment])
    )
    for reference_genome in reference_genomes:
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
            read_name_field_name = 'read_names_%s_%s' % (reference_genome, read_type)
            temp_read_summary_table = (read_summary_table
                .addfield('read_type', read_type)
                .addfield('read_name_temp', lambda(rec): 'run_id_%s_channel_%d_read_%d' % (rec['run_id'], rec['channel_number'], rec['read_number']))
                .addfield('read_name', lambda(rec): rec['read_name_temp'] + '_' + rec['read_type'])
            )
            read_summary_table = (temp_read_summary_table
                .rename('read_name', read_name_field_name)
                .cutout('read_type')
                .cutout('read_name_temp')
                .leftjoin(mapped_read_table, missing=0)
            )
    read_summary_tables[experiment] = read_summary_table

# <codecell>

read_summary_tables['E07_ANOPH_workflow_1_9_1.2.1'].select(lambda rec: rec['fragments_AgamP3_2D'] > 0).head(3)

# <codecell>

def summarise_experiment(experiment = 'expt01_workflow_1_9_1', read_summary_tables=read_summary_tables):
    read_summary_table = read_summary_tables[experiment]
    read_summary_array = toarray(read_summary_table)
    experiment_summary_dict = OrderedDict()
    experiment_summary_dict['non_mapping'] = np.array(
        (
             experiment,
             '_'.join(np.unique(read_summary_array['device_id'])),
             '_'.join(np.unique(read_summary_array['asic_id'])),
             '_'.join(np.unique(read_summary_array['flow_cell_id'])),
             '___'.join(map(convert_time, np.unique(read_summary_array['exp_start_time']))),
             len(read_summary_array),
             np.sum(read_summary_array['is_analysed']),
             np.sum(read_summary_array['is_base_called_2d']),
             np.sum(read_summary_array['has_template_read']),
             np.sum(read_summary_array['has_complement_read']),
             np.sum(read_summary_array['has_2D_read']),
             np.sum(read_summary_array['has_template_read']) * 1.0 / len(read_summary_array),
             np.sum(read_summary_array['has_complement_read']) * 1.0 / np.sum(read_summary_array['has_template_read']),
             np.sum(read_summary_array['has_2D_read']) * 1.0 / np.sum(read_summary_array['has_complement_read']),
             len(np.unique(read_summary_array['channel_number'])),
             len(np.unique(read_summary_array['channel_number'][read_summary_array['has_2D_read']])),
             np.sum(read_summary_array['read_length_template']),
             np.sum(read_summary_array['read_length_complement']),
             np.sum(read_summary_array['read_length_2D']),
             np.mean(read_summary_array['read_length_template'][read_summary_array['read_length_template'] > 0]),
             np.mean(read_summary_array['read_length_complement'][read_summary_array['read_length_complement'] > 0]),
             np.mean(read_summary_array['read_length_2D'][read_summary_array['read_length_2D'] > 0]),
             np.mean(read_summary_array['gc_template'][logical_not(isnan(read_summary_array['gc_template']))]),
             np.mean(read_summary_array['gc_complement'][logical_not(isnan(read_summary_array['gc_complement']))]),
             np.mean(read_summary_array['gc_2D'][logical_not(isnan(read_summary_array['gc_2D']))]),
             np.mean(read_summary_array['duration_read'][(read_summary_array['duration_read'] > 0) & (read_summary_array['duration_template'] > 0)]),
             np.mean(read_summary_array['duration_read'][read_summary_array['duration_read'] > 0]),
             np.mean(read_summary_array['duration_template'][read_summary_array['duration_template'] > 0]),
             np.mean(read_summary_array['duration_complement'][read_summary_array['duration_complement'] > 0]),
             np.mean(read_summary_array['basecall_mean_qscore_template'][read_summary_array['basecall_mean_qscore_template'] > 0]),
             np.mean(read_summary_array['basecall_mean_qscore_complement'][read_summary_array['basecall_mean_qscore_complement'] > 0]),
             np.mean(read_summary_array['basecall_mean_qscore_2D'][read_summary_array['basecall_mean_qscore_2D'] > 0]),
             np.mean(read_summary_array['basecall_strand_score_template'][read_summary_array['basecall_strand_score_template'] > 0]),
             np.mean(read_summary_array['basecall_strand_score_complement'][read_summary_array['basecall_strand_score_complement'] > 0]),
        ),
        dtype  = [
            ('experiment', 'a100'),
            ('device_ids', 'a100'),
            ('asic_ids', 'a100'),
            ('flowcell_ids', 'a100'),
            ('exp_start_times', 'a1000'),
            ('total_number_of_reads', np.int),
            ('number_of_reads_analysed', np.int),
            ('number_of_reads_base_called', np.int),
            ('number_of_reads_with_template_calls', np.int),
            ('number_of_reads_with_complement_calls', np.int),
            ('number_of_reads_with_2D_calls', np.int),
            ('proportion_of_reads_with_template_calls', np.float),
            ('proportion_of_template_called_reads_with_complement_calls', np.float),
            ('proportion_of_complement_called_reads_with_2D_calls', np.float),
            ('number_of_channels_used', np.int),
            ('number_of_channels_with_2D_calls', np.int),
            ('total_read_length_template', np.int),
            ('total_read_length_complement', np.int),
            ('total_read_length_2D', np.int),
            ('mean_read_length_template', np.int),
            ('mean_read_length_complement', np.int),
            ('mean_read_length_2D', np.int),
            ('mean_gc_template', np.float),
            ('mean_gc_complement', np.float),
            ('mean_gc_2D', np.float),
            ('mean_duration_read_where_template_called', np.int),
            ('mean_duration_read', np.int),
            ('mean_duration_template', np.int),
            ('mean_duration_complement', np.int),
            ('mean_qscore_template', np.float),
            ('mean_qscore_complement', np.float),
            ('mean_qscore_2D', np.float),
            ('mean_strand_score_template', np.float),
            ('mean_strand_score_complement', np.float),
        ]
    )
    for reference_genome in reference_genomes:
        for read_type in read_types:
#             print experiment, reference_genome, read_type
            experiment_summary_dict['%s_%s' % (reference_genome, read_type)] = np.array(
                (
                     np.sum(read_summary_array['fragments_%s_%s' % (reference_genome, read_type)] > 0),
                     np.sum(read_summary_array['fragments_%s_%s' % (reference_genome, read_type)] > 0) * 1.0 / experiment_summary_dict['non_mapping']['number_of_reads_with_%s_calls' % (read_type)],
                     0.0 if np.sum(read_summary_array['fragments_%s_%s' % (reference_genome, read_type)] > 0) == 0 else np.mean(read_summary_array['fragments_%s_%s' % (reference_genome, read_type)].astype(np.float)[read_summary_array['fragments_%s_%s' % (reference_genome, read_type)] > 0]),
                     np.sum(read_summary_array['read_length_template'][read_summary_array['fragments_ont_lambda_template'] > 0]),
                     np.sum(read_summary_array['total_fragment_length_%s_%s' % (reference_genome, read_type)]),
                     np.sum(read_summary_array['mismatches_%s_%s' % (reference_genome, read_type)]),
                     np.sum(read_summary_array['mismatches_%s_%s' % (reference_genome, read_type)])*1.0 / np.sum(read_summary_array['total_fragment_length_%s_%s' % (reference_genome, read_type)]),
                     np.sum(read_summary_array['snps_%s_%s' % (reference_genome, read_type)])*1.0 / np.sum(read_summary_array['total_fragment_length_%s_%s' % (reference_genome, read_type)]),
                     (np.sum(read_summary_array['snps_%s_%s' % (reference_genome, read_type)])+np.sum(read_summary_array['Ns_%s_%s' % (reference_genome, read_type)]))*1.0 / np.sum(read_summary_array['total_fragment_length_%s_%s' % (reference_genome, read_type)]),
                     np.max(read_summary_array['total_fragment_length_%s_%s' % (reference_genome, read_type)]),
                ),
                dtype  = [
                    ('number_mapped_to_%s_%s' % (reference_genome, read_type), np.int),
                    ('proportion_mapped_to_%s_%s' % (reference_genome, read_type), np.float),
                    ('mean_number_fragments_%s_%s' % (reference_genome, read_type), np.float),
                    ('total_length_of_reads_mapped_to_%s_%s' % (reference_genome, read_type), np.int),
                    ('mapped_length_of_reads_mapped_to_%s_%s' % (reference_genome, read_type), np.int),
                    ('number_of_mismatches_%s_%s' % (reference_genome, read_type), np.int),
                    ('error_rate_%s_%s' % (reference_genome, read_type), np.float),
                    ('snp_rate_discarding_Ns_%s_%s' % (reference_genome, read_type), np.float),
                    ('snp_rate_including_Ns_%s_%s' % (reference_genome, read_type), np.float),
                    ('longest_match_to_%s_%s' % (reference_genome, read_type), np.int),
                ]
            )
    experiment_summary = numpy.lib.recfunctions.merge_arrays(experiment_summary_dict.values(), flatten=True)
    return experiment_summary

summarise_experiment()['mean_gc_template']

# <codecell>

experiment_summaries_dict = OrderedDict()

for experiment in experiments:
    print experiment
    experiment_summaries_dict[experiment] = summarise_experiment(experiment, read_summary_tables)

experiment_summaries = np.hstack(experiment_summaries_dict.values())

# <codecell>

toxlsx(fromarray(experiment_summaries), "/Users/rpearson/Desktop/experiment_summaries_20140711.xlsx")

# <codecell>

etl.fromarray(experiment_summaries)

# <codecell>


# <codecell>


# <codecell>


