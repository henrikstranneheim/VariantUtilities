#!/usr/bin/env python
# coding: utf-8


#Import standard Library
import os
import sys
import argparse
import csv
from collections import defaultdict

#Import third party library
import vcf
from bx.intervals.intersection import IntervalTree
from pyfasta import Fasta

####User options
parser = argparse.ArgumentParser(description='Read a vcf.')
parser.add_argument('infile',
                    type=str, nargs=1, help='(VCF)'
                    )
parser.add_argument('annotation_file',
                    type=str, nargs=1,
                    help='A annotations file. Default is ref_gene format.'
                    )
parser.add_argument('-at', '--annotation_type',
                    type=str, nargs=1, choices=['ref_gene'],
                    default='ref_gene', help='The format of the annotation file. Default: "ref_gene"'
                    ) 
parser.add_argument('-g', '--genome_reference',
                    type=str, nargs=1,
                    help='The genome reference fasta file.'
                    )
parser.add_argument('-v', '--version',
                    action="store_true", help='Display version'
                    )
args = parser.parse_args()

program_version="v.1.0"

if args.version:
    sys.stdout.write("Splicer ", program_version)
    exit()

#Print version
sys.stderr.write("Splicer " + program_version + '\n')


vcf_file_path = args.infile[0]
vcf_file_name, vcf_file_extension = os.path.splitext(vcf_file_path)
annotation_file_path = args.annotation_file[0]
annotation_file_name, annotation_file_extension = os.path.splitext(annotation_file_path)
annotation_type = args.annotation_type
if args.genome_reference:

    fasta_path = args.genome_reference[0]
    fasta_file_name, fasta_file_extension = os.path.splitext(fasta_path)

def codons():
    """Generate codon table"""
    
    bases = ['t', 'c', 'a', 'g']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, amino_acids))

def translate(seq):
    """From input sequence generate codon unitl you hit stop or unknown codon and return it"""
    
    seq = seq.lower().replace('\n', '').replace(' ', '')
    peptide = ''
    
    for i in xrange(0, len(seq), 3):
        codon = seq[i: i+3]
        amino_acid = codon_table.get(codon, '*')
        if amino_acid != '*':
            peptide += amino_acid
        else:
            break
        
    return peptide
    
def list_to_key(dictionnary, key, array):
    """Inputs a dictionnary, the key to assign to the list and a list"""
    
    dictionnary[key] = dict()
    
    for i, feature in enumerate(array):
        
        dictionnary[key][i] = feature

    
def ref_gene_parser(line):
    """Parse a file in the refGene format. Return gene object with collect info"""
    #dictionary for collecting gene info
    gene = dict()

    #remove the chromosome prefix
    if 'hr' in line[2]:
        gene['chrom'] = line[2][3:]
    else:
        gene['chrom'] = line[2]
    
    #If start and stop is found
    if line[4].isdigit() and line[5].isdigit():
      
        #Gene information
        gene['gene_id'] = line[12]
        gene['feature_id'] = gene['gene_id']
        gene['strand'] = line[3]
        gene['start'] = int(line[4])
        gene['stop'] = int(line[5])
        gene['cds_start'] = int(line[6])
        gene['cds_stop'] = int(line[7])

        #Transcript information
        gene['transcript_id'] = line[1]

        #Exon information
        gene['exon_total_number'] = int(line [8])

        exon_starts = line[9].split(',')
        #remove empty entries from exon_starts from the split
        exon_starts.pop()
        list_to_key(gene, "exon_start", exon_starts)

        exon_stops = line[10].split(',')
        #remove empty entries from exon_stops from the split
        exon_stops.pop()
        list_to_key(gene, "exon_stop", exon_stops)

        return gene

def create_fasta_flat_file(file):
    """Reads a fasta file for fast sequence retrival"""

    fasta_file = Fasta(file, key_fn=lambda key: key.split()[0])

    fasta_headers = set(fasta_file.keys());

    return fasta_file, fasta_headers

def read_fasta_file(file, chromosome, start_position, stop_position):
    """Reads a fasta file for fast sequence retrival"""

    if chromosome in fasta_headers:
        
        return fasta_file[chromosome][start_position:stop_position]

def map_genomic_position_to_mrna_position(coding_region_info, strand, transcript, position_mrna, position_exon, exon_start):
    """Maps the genomic position to the mRNA position. Returns the last exon position to be mapped to the mRNA"""                               

    if "+" in strand:

        for nt in range(len(position_exon)):
        
            coding_region_info[transcript][exon_start + nt] = position_mrna
            position_mrna += 1
            
    return position_mrna

def index_annotation_file(annotation_file_path, annotation_type):
    """"Parses a annotation file and builds an interval tree"""

    #dictionary mapping chromosome names to interval trees, collecting geneID info
    genome = dict()
    #dictionary mapping chromosmoes names to interval trees, collecting transcriptID info
    transcriptome = dict()
    #dictionnary mapping transcript info
    transcripts_info = dict()
    #dictionary mapping chromosmoes names to interval trees, collecting exon number info
    exome = dict()
    #dictionnary mapping coding region info info
    coding_region_info = defaultdict(dict)

    with open(annotation_file_path,  'r') as annotation_file:

        reader = csv.reader(annotation_file, delimiter='\t')
        
        for line in reader:

            #Start with blank tree for each line
            tree_gene = None
            tree_transcript = None
            tree_exon = None

            if annotation_type == 'ref_gene':

                gene = ref_gene_parser(line)

                #one interval tree per chromosome
                if gene['chrom'] in genome:

                    tree_gene = genome[gene['chrom']]
                    tree_transcript = transcriptome[gene['chrom']]
                    tree_exon = exome[gene['chrom']]

                else:
                
                    #Chromosome not seen previously, create interval tree key
                    tree_gene = IntervalTree()
                    tree_transcript = IntervalTree()
                    tree_exon = IntervalTree()
                    genome[gene['chrom']] = tree_gene
                    transcriptome[gene['chrom']] = tree_transcript
                    exome[gene['chrom']] = tree_exon
                
                #index the feature
                tree_gene.add(gene['start'], gene['stop'], gene['gene_id'])
                tree_transcript.add(gene['start'], gene['stop'], gene['transcript_id'])
                
                #Fasta file exists
                if args.genome_reference: 
                    
                    transcripts_info[ gene['transcript_id'] ] = gene['transcript_id']

                    #Collect fasta sequence and coding region 
                    coding_region_info[ gene['transcript_id'] ]['fasta'] = read_fasta_file(fasta_path, gene['chrom'], int(gene['cds_start']), int(gene['cds_stop']))
                    coding_region_info[ gene['transcript_id'] ]['cds_start'] = gene['cds_start']
                    coding_region_info[ gene['transcript_id'] ]['cds_stop'] = gene['cds_stop']
                    coding_region_info[ gene['transcript_id'] ]['strand'] = gene['strand']
                
                mrna_fasta = []
                position_mrna = 0
                for exon in gene['exon_start']:

                    tree_exon.add(int(gene['exon_start'][exon]), int(gene['exon_stop'][exon]), exon) 
                    
                    #print(gene['transcript_id'], exon)
                    if coding_region_info[ gene['transcript_id'] ]['fasta']:

                        start_fasta = 0
                        stop_fasta = 0

                        if "+" in gene['strand']:

                            #Within coding region
                            if (int(gene['exon_start'][exon]) > gene['cds_start']) and (int(gene['exon_stop'][exon]) < gene['cds_stop']):
                                
                                start_fasta = int(gene['exon_start'][exon]) - gene['cds_start']
                                stop_fasta = int(gene['exon_stop'][exon]) - gene['cds_start']
                                position_exon = range(int(gene['exon_start'][exon]), int(gene['exon_stop'][exon]))
                                
                                position_mrna = map_genomic_position_to_mrna_position(coding_region_info, gene['strand'], gene['transcript_id'], position_mrna, position_exon, int(gene['exon_start'][exon]))
                                
                                #Upstream of coding region
                            elif (int(gene['exon_stop'][exon]) < gene['cds_start']):
                                
                                start_fasta = 0
                                stop_fasta = 0
                                
                            #Downstream of coding region
                            elif (int(gene['exon_start'][exon]) > gene['cds_stop']):
                                
                                start_fasta = 0
                                stop_fasta = 0
                            
                            #Start downstream of cds
                            elif int(gene['exon_start'][exon]) < gene['cds_start']:
                            
                                start_fasta = 0
                            
                                #Exon encompasses whole cds
                                if ( (int(gene['exon_stop'][exon]) > gene['cds_start']) and (int(gene['exon_stop'][exon]) > gene['cds_stop']) ):

                                    stop_fasta = gene['cds_stop'] - gene['cds_start']
                                    position_exon = range(gene['cds_start'], gene['cds_stop'])

                                #Finish upstream of cds start, but less than cds stop (handled above)
                                elif int(gene['exon_stop'][exon]) > gene['cds_start']:
                                
                                    stop_fasta = int(gene['exon_stop'][exon]) - gene['cds_start']
                                    position_exon = range(gene['cds_start'], int(gene['exon_stop'][exon]))                                
                                
                                    position_mrna = map_genomic_position_to_mrna_position(coding_region_info, gene['strand'], gene['transcript_id'], position_mrna, position_exon, int(gene['exon_start'][exon]))
                            
                                elif int(gene['exon_stop'][exon]) > gene['cds_stop']:

                                    stop_fasta = gene['cds_stop'] - gene['cds_start']

                                    if int(gene['exon_start'][exon]) < gene['cds_stop']:
                                
                                        start_fasta = int(gene['exon_start'][exon]) - gene['cds_start']
                                        position_exon = range(int(gene['exon_start'][exon]), gene['cds_stop'])
                                        
                                        position_mrna = map_genomic_position_to_mrna_position(coding_region_info, gene['strand'], gene['transcript_id'], position_mrna, position_exon, int(gene['exon_start'][exon]))

                            mrna_fasta.append(coding_region_info[ gene['transcript_id'] ]['fasta'][start_fasta:stop_fasta])

                    coding_region_info[gene['transcript_id'] ]['mRNA'] = ''.join(mrna_fasta)
                
    return genome, transcriptome, exome, transcripts_info, coding_region_info

def collect_annotation(chromosome, position, gene):
    """Collects the gene, transcript and exon information. Lacking a gene feature returns 0"""

    if gene:
        
        previous_gene = None
        
        for annotation in gene:
            
            #Only go through transcripts if previous gene does not equal annotation
            if previous_gene != annotation:
                
                transcripts = transcriptome[chromosome].find(position, position)
                previous_transcript = None

                for transcript_id in transcripts:

                    exons = exome[chromosome].find(position, position)
                    
                    if "+" in coding_region_info[transcript_id]['strand']:
                        
                        if exons:
                            
                            for exon_id in exons:
                                
                                if transcript_id in transcripts_info:
                                    
                                #print(transcript_id)
                                #Within mRNA
                                    if position in coding_region_info[transcript_id]:
                                        
                                        position_mrna = coding_region_info[transcript_id][position]
                                      #exon_id+1 for output reasons only exon numbering usually starts at 1 
                                        sys.stdout.write("\t" + annotation + ":" +  transcript_id + ":exon" + (str(exon_id+1)))
                                        sys.stdout.write(":c." + str(position_mrna))
                                        
                                    else:
                                        
                                        position_utr = (position - coding_region_info[transcript_id]['cds_start']) - len(coding_region_info[transcript_id]['mRNA'])
                                        
                                        if position_utr < 0: 
                                            
                                            sys.stdout.write("\t" + annotation + "\tUTR5")
                                            
                                        else:
                                            
                                            sys.stdout.write("\t" + annotation + "\tUTR3")
                                            
                                else:
                                    
                                    sys.stdout.write("\t" + annotation + "\tintronic")
                                    
                                  #Save the current gene to avoid printing multiple entries of the same gene
                                    previous_gene = annotation                                   

        #Found gene for coordinate
        return 1
    
    else:
        #No gene was found for coordinate
        return 0
        

def parse_vcf_file(vcf_file_path):

    vcf_reader = vcf.Reader(open(vcf_file_path, 'r'))

    for record in vcf_reader:
        
        if record.CHROM in genome:
            
            if genome[record.CHROM]:
                
                sys.stdout.write(str(record.CHROM) + "\t" + str(record.POS))

                found_gene = collect_annotation(record.CHROM, record.POS, genome[record.CHROM].find(record.POS, record.POS)) 

            if found_gene != 1:

                sys.stdout.write("\tIntergenic")
            
            sys.stdout.write("\n")

####MAIN

if args.genome_reference:

    fasta_file, fasta_headers = create_fasta_flat_file(fasta_path)

(genome, transcriptome, exome, transcripts_info, coding_region_info) = index_annotation_file(annotation_file_path, annotation_type)

parse_vcf_file(vcf_file_path)

###Test
#annotation = genome['10'].find(75910939, 77000000)
#print annotation
