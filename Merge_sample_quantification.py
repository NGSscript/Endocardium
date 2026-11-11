#!/usr/bin/env python
#coding=utf-8

"""
Versions:001
Desc:
"""
import argparse
import os, re

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--samples_info', dest='samples_info', metavar='FILE', type=str, required=True,
                        help="fast pjson")
    parser.add_argument('-w', '--work_dir', dest='work_dir', metavar='DIR', type=str, default="./",
                        help="Working/output directory path")
    parser.add_argument('-a', '--gtf', dest='gtf', metavar='FILE', type=str, required=True,
                        help="GTF annotation file for target species")
    parser.add_argument('-c', '--featureCounts', dest='feature_counts', metavar='STR', type=str, default="gene",
                        help="Feature type for counting: gene or transcript")
    parser.add_argument('-s', '--alignment_soft', dest='alignment_soft', metavar='STR', type=str, default="hisat2",
                        help="Name of read alignment software")
    args = parser.parse_args()
    return args

def get_samples_list(samples_info_file, work_dir):
    samples_list = list()
    with open(samples_info_file) as inf:
        for lines in inf:
            if lines.startswith('sampleid'):
                continue
            tabs = lines.strip().split(",")
            if os.path.exists(f"{work_dir}/{tabs[0]}"):
                samples_list.append(tabs[0])
    return samples_list

def get_read_gtf(gtf_file, feature_counts):
    id_info_dict = dict()
    with open(gtf_file) as inf:
        for lines in inf:
            if lines.startswith('#'):
                continue
            tabs = lines.strip().split("\t")
            if tabs[2] != 'transcript':
                continue
            gene_biotype = re.findall("gene_biotype \"(.+?)\";", tabs[8])[0]
            gene_id = re.findall("gene_id \"(.+?)\";", tabs[8])[0]
            if "gene_name" in tabs[8]:
                gene_name = re.findall("gene_name \"(.+?)\";", tabs[8])[0]
            else:
                gene_name = ""
            trans_id = re.findall("transcript_id \"(.+?)\";", tabs[8])[0] 
            line = f"{tabs[0]}\t{tabs[3]}\t{tabs[4]}\t{tabs[6]}\t{gene_id}\t{gene_name}\t{gene_biotype}"
            if feature_counts == 'transcript':
                if trans_id not in id_info_dict:
                    id_info_dict[trans_id] = line
            else:
                if gene_id not in id_info_dict:
                    id_info_dict[gene_id] = line
    return id_info_dict

def get_file_dict(type_samples_dict, gene_id_list, file_name, sample):
    with open(file_name) as inf:
        for lines in inf:
            if lines.startswith('GeneID'):
                continue
            tabs = lines.strip().split("\t")
            if tabs[0] not in gene_id_list:
                gene_id_list.append(tabs[0])
            if tabs[0] not in type_samples_dict:
                type_samples_dict[tabs[0]] = dict()
            type_samples_dict[tabs[0]][sample] = tabs[1]

def main():
    args = parse_args()
    samples_list = get_samples_list(args.samples_info, args.work_dir)
    gene_id_name_dict = get_read_gtf(args.gtf, args.feature_counts)
    type_list = ["count", "rpkm", "fpkm", "tpm","cpm", "log2cpm"]
    for type_ in type_list:
        gene_id_list = list()
        type_samples_dict = dict()
        flog = 0
        for sample in samples_list:
            file_name = f"{args.work_dir}/{sample}/03_matrix_processing/{sample}.{type_}"
            if args.feature_counts == 'transcript':
                file_name = f"{args.work_dir}/{sample}/03_matrix_processing/{sample}_trans.{type_}"
            if os.path.exists(file_name):
                 flog += 1
                 get_file_dict(type_samples_dict, gene_id_list, file_name, sample)
        if flog != len(samples_list):
            continue
        
        if args.feature_counts == 'transcript':
            write_file_ = f"{args.work_dir}/Quantification/{args.alignment_soft}_{args.feature_counts}.{type_}.tsv"
        else:
            write_file_ = f"{args.work_dir}/Quantification/{args.alignment_soft}.{type_}.tsv"
            
        
        with open(write_file_, 'w') as outf:
            sample_str = '\t'.join(samples_list)
            outf.write("ID\tCHR\tstart\tend\tstrand\tgene_id\tgene_name\tgene_biotype\t%s\n" % sample_str)
            for gene_id in gene_id_list:
                if gene_id in gene_id_name_dict:
                    line_ = f"{gene_id}\t{gene_id_name_dict[gene_id]}"
                else:
                    line_ = f"{gene_id}\t\t\t\t\t\t\t"
                for sample in samples_list:
                    if sample in type_samples_dict[gene_id]:
                        line_ += "\t" + type_samples_dict[gene_id][sample]
                    else:
                        line_ += "\t"
                outf.write("%s\n" % line_)

if __name__ == "__main__":
    main()
