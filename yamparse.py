import csv
import yaml
import os
import time
import glob
import re
import argparse
from decimal import Decimal
parser = argparse.ArgumentParser()

# parser.add_argument("file", type=str )
# parser.add_argument("outdir", type=str)
# args = parser.parse_args()

infile = 'dir.txt'

yaml_out = {}

qc_files= (
'alignment_summary.txt',
'verify_bam_id.selfSM',
'wgs_metric_summary.txt',
'bamutil_stats.txt',
'flagstat.out',
'insert_size_summary.txt',
'mark_dups_metrics.txt',
'GC_bias_summary.txt',
'all_chrom.variant_calling_detail_metrics',
)

wd = os.getcwd()
# path = '/Users/ltrani/Desktop/yamlparser/yamlparser/gscmnt/gc13035/production/compute_158076365/'

def process_alignment_summary(qc_file):
    if not os.path.exists(qc_file):
        print('no file found' + qc_file)
        return None

    with open(qc_file) as readfile:
        for line in readfile:
            if line.startswith('#') or not line.strip():
                continue
            break
        reader = csv.DictReader(readfile, fieldnames=line.split("\t"), delimiter="\t")
        for line in reader:
            if line['CATEGORY'] == 'PAIR':
                yaml_out['ALIGNMENT_RATE'] = str((int(line['PF_READS_ALIGNED']))/int(line['TOTAL_READS']))
                yaml_out['PCT_ADAPTER'] = line['PCT_ADAPTER']
                yaml_out['PF_READS'] = line['PF_READS']
                yaml_out['PF_ALIGNED_BASES'] = line['PF_ALIGNED_BASES']
                yaml_out['TOTAL_READS'] = line['TOTAL_READS']
                yaml_out['PF_HQ_ALIGNED_Q20_BASES'] = line['PF_HQ_ALIGNED_Q20_BASES']

            # ALIGNMENT_RATE	alignment_summary.txt	(PF_READS_ALIGNED,TOTAL_READS) FROM PAIR done
# PCT_ADAPTER	alignment_summary.txt	(PCT_ADAPTER) from PAIR done
# PF_READS	alignment_summary.txt	(PF_READS) from (PAIR) done
# PF_ALIGNED_BASES	alignment_summary.txt	(PF_ALIGNED_BASES) from PAIR  done
# TOTAL_READS	alignment_summary.txt	(TOTAL_READS) from PAIR
# PF_HQ_ALIGNED_Q20_BASES	alignment_summary.txt	(PF_HQ_ALIGNED_Q20_BASES) from PAIR

            if line['CATEGORY'] == 'FIRST_OF_PAIR':
                yaml_out['FIRST_OF_PAIR_MISMATCH_RATE'] = line['PF_MISMATCH_RATE']
# FIRST_OF_PAIR_MISMATCH_RATE	alignment_summary.txt	(PF_MISMATCH_RATE) from FIRST_OF_PAIR

            if line['CATEGORY'] == 'SECOND_OF_PAIR':
                yaml_out['SECOND_OF_PAIR_MISMATCH_RATE'] = line['PF_MISMATCH_RATE']
# SECOND_OF_PAIR_MISMATCH_RATE	alignment_summary.txt	(PF_MISMATCH_RATE) from SECOND_OF_PAIR


def process_verify_bam_id_selfSM(qc_file):
    if not os.path.exists(qc_file):
        print('file not found' + qc_file)
        return None

    with open(qc_file) as readfile:
        for line in readfile:
            reader = csv.DictReader(readfile, fieldnames=line.split("\t"), delimiter="\t")
            for line in reader:
                yaml_out['FREEMIX'] = line['FREEMIX']
        # FREEMIX	verify_bam_id.selfSM	FREEMIX value under this column


def process_wgs_metric_summary(qc_file):
    if not os.path.exists(qc_file):
        print('file not found' + qc_file)
        return None

    with open(qc_file) as readfile:
        for line in readfile:
            if line.startswith('#') or not line.strip():
                continue
            break
        reader = csv.DictReader(readfile, fieldnames=line.strip().split("\t"), delimiter="\t")
        line = next(reader)
        haploid_coverage = float(line['MEAN_COVERAGE']) * ((1 - float(line['PCT_EXC_DUPE'])) / (1 - float(line['PCT_EXC_TOTAL'])))
        yaml_out['PCT_10X'] = line['PCT_10X']
        yaml_out['PCT_20X'] = line['PCT_20X']
        yaml_out['HET_SNP_Q'] = line['HET_SNP_Q']
        yaml_out['HET_SNP_SENSITIVITY'] = line['HET_SNP_SENSITIVITY']
        yaml_out['MEAN_COVERAGE'] = line['MEAN_COVERAGE']
        yaml_out['HAPLOID_COVERAGE'] = str(haploid_coverage)
# HAPLOID_COVERAGE	wgs_metric_summary.txt	MEAN_COVERAGE,PCT_EXC_DUPE,PCT_EXC_TOTAL
# PCT_10X	wgs_metric_summary.txt	PCT_10X
# PCT_20X	wgs_metric_summary.txt	PCT_20X
# HET_SNP_Q	wgs_metric_summary.txt	HET_SNP_Q
# HET_SNP_SENSITIVITY	wgs_metric_summary.txt	HET_SNP_SENSITIVITY
# MEAN_COVERAGE	wgs_metric_summary.txt	MEAN_COVERAGE


def process_bamutil_stats(qc_file):
    if not os.path.exists((qc_file)):
        print('file not found ' + qc_file)
        return None
    with open(qc_file) as readfile:
        phred20 = ''
        phred30 = ''
        for line in readfile:
            if line.startswith('N') or not line.strip():
                continue
            break

        reader = csv.DictReader(readfile, fieldnames=line.strip().split("\t"), delimiter="\t")
        for line in reader:
            if line['Phred'] == '20':
                phred20 = line['Count']
            if line['Phred'] == '30':
                phred30 = line['Count']

        TOTAL_BASES_Q20_OR_MORE = int(phred20) + int(phred30)
        yaml_out['TOTAL_BASES_Q20_OR_MORE'] = str(TOTAL_BASES_Q20_OR_MORE)


# TOTAL_BASES_Q20_OR_MORE	bamutil_stats.txt	column 2 values corresponding to 20,30


def process_flagstat_out(qc_file):
    if not os.path.exists(qc_file):
        print('files not found' + qc_file)
        return None

    with open(qc_file) as readfile:
        reader = csv.reader(readfile)
        mapped = ''
        properly_paired = ''
        with_itself_and_mate_mapped = ''
        for line in reader:
            for match in line:
                if re.search('%', match) and re.search('mapped', match):
                    mapped = re.search("\d+.\d+%", match).group()
                    mapped = re.match("\d+.\d+", mapped).group()
                    mapped = Decimal(mapped)

                if re.search('%', match) and re.search('properly', match):
                    properly_paired = re.search("\d+.\d+%", match).group()
                    yaml_out['reads_mapped_in_proper_pairs_percentage'] = str(properly_paired)
                    properly_paired = re.match("\d+.\d+", properly_paired).group()
                    properly_paired = Decimal(properly_paired)


                if re.search('itself', match):
                    itself = match.split()
                    with_itself_and_mate_mapped = itself[0]

                if re.search('different', match) and not re.search('mapQ', match):
                    mate = match.split()
                    mate_mapped_different_chr = mate[0]

                if re.search('%', match) and re.search('singleton', match):
                    reads_mapped_singleton = re.search("\d+.\d+%", match).group()
                    yaml_out['reads_mapped_as_singleton_percentage'] = str(reads_mapped_singleton)

        discordant_rate = mapped - properly_paired
        yaml_out['discordant_rate'] = str(discordant_rate)

        interchromosomal_rate = int(mate_mapped_different_chr) / int(with_itself_and_mate_mapped)
        yaml_out['interchromosomal_rate'] = str(interchromosomal_rate)

# discordant_rate	flagstat.out	mapped (%), properly paired (%)
# interchromosomal_rate	flagstat.out	with itself and mate mapped,with mate mapped to a different chr
# reads_mapped_as_singleton_percentage	flagstat.out	singletons (%)
# reads_mapped_in_proper_pairs_percentage	flagstat.out	properly paired (%)


def process_insert_size_summary(qc_file):
    if not os.path.exists((qc_file)):
        print('file not found ' + qc_file)
        return None

    with open(qc_file) as readfile:
        for line in readfile:
            if line.startswith('#')or not line.strip():
                continue
            break
        reader = csv.DictReader(readfile, fieldnames=line.split("\t"), delimiter="\t")
        line = next(reader)
        yaml_out['MEAN_INSERT_SIZE'] = line['MEAN_INSERT_SIZE']
        yaml_out['STANDARD_DEVIATION'] = line['STANDARD_DEVIATION']

# MEAN_INSERT_SIZE	insert_size_summary.txt	MEAN_INSERT_SIZE
# STANDARD_DEVIATION	insert_size_summary.txt	STANDARD_DEVIATION

def process_mark_dups_metrics(qc_file):
    if not os.path.exists(qc_file):
        print('file not found ' + qc_file)
        return  None

    with open(qc_file) as readfile:
        for line in readfile:
            if line.startswith('#') or not line.strip():
                continue
            break
        reader = csv.DictReader(readfile, fieldnames=line.strip().split("\t"), delimiter="\t")
        line = next(reader)
        PERCENT_DUPLICATION = ((int(line['UNPAIRED_READ_DUPLICATES']) + int(line['READ_PAIR_DUPLICATES'])) ) / \
                              ((int(line['UNPAIRED_READS_EXAMINED']) + int(line['READ_PAIRS_EXAMINED'])) )
        yaml_out['PERCENT_DUPLICATION'] = str(PERCENT_DUPLICATION)
# PERCENT_DUPLICATION = (UNPAIRED_READ_DUPLICATES + READ_PAIR_DUPLICATES *2) /(double) (UNPAIRED_READS_EXAMINED + READ_PAIRS_EXAMINED *2);
# TOTAL_PERCENT_DUPLICATION	mark_dups_metrics.txt 	UNPAIRED_READS_EXAMINED,READ_PAIR_DUPLICATES,UNPAIRED_READS_EXAMINED,READ_PAIRS_EXAMINED


def process_GC_bias_summary(qc_file):
    if not os.path.exists(qc_file):
        print('file not found ' + qc_file)
        return None

    with open(qc_file) as readfile:
        for line in readfile:
            if line.startswith('#') or not line.strip():
                continue
            break
        reader = csv.DictReader(readfile, fieldnames=line.split("\t"), delimiter="\t")
        for line in reader:
            yaml_out['ALIGNED_READS'] = line['ALIGNED_READS']
    # ALIGNED_READS


def process_all_chrom_variant_calling_detail_metrics(qc_file):
    if not os.path.exists(qc_file):
        print('no file found')
        return None
    # Open the file
    with open( qc_file ) as readfile:
        # Pass file to csv dictreader and skip the header
            # (see http://stackoverflow.com/questions/14158868/python-skip-comment-lines-marked-with-in-csv-dictreader)
        for line in readfile:
            if line.startswith('#') or not line.strip():
                continue
            break

        reader = csv.DictReader(readfile, fieldnames=line.split("\t"), delimiter="\t")
        for line in reader:
            yaml_out['SAMPLE_ALIAS'] = line['SAMPLE_ALIAS']
        # Read the file
        # Return the first element of the reader (a dictionary)

with open(infile) as csvfile:

    reader = csv.reader(csvfile, delimiter="\t")
    for line in reader:
        for path in line:
            path = '/Users/ltrani/Desktop/yamlparser/yamlparser/' + path
            outfile = os.path.basename(os.path.normpath(path))
            outfile = outfile + '.qc.tsv'
            print(outfile)
            out = open(outfile, 'w')

        if not os.path.exists(path):
            print('dir not found')
            continue
            exit()

        for file in qc_files:
            if file == 'all_chrom.variant_calling_detail_metrics':
                qc_file = path + file
                process_all_chrom_variant_calling_detail_metrics(qc_file)

            if file == 'alignment_summary.txt':
                qc_file = path + file
                process_alignment_summary(qc_file)

            if file == 'verify_bam_id.selfSM':
                qc_file = path + file
                process_verify_bam_id_selfSM(qc_file)

            if file == 'wgs_metric_summary.txt':
                qc_file = path + file
                process_wgs_metric_summary(qc_file)

            if file == 'flagstat.out':
                qc_file = path + file
                process_flagstat_out(qc_file)

            if file == 'GC_bias_summary.txt':
                qc_file = path + file
                process_GC_bias_summary(qc_file)

            if file == 'insert_size_summary.txt':
                qc_file = path + file
                process_insert_size_summary(qc_file)

            if file == 'mark_dups_metrics.txt':
                qc_file =  path + file
                process_mark_dups_metrics(qc_file)

            if file == 'bamutil_stats.txt':
                qc_file = path + file
                process_bamutil_stats(qc_file)

        yaml.dump(yaml_out, out, default_flow_style=False)


exit()

