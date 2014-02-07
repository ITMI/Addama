import argparse
import csv
import json
import StringIO
from subprocess import check_output

class MultilineTabixResult():
    def __init__(self, chromosome, start, end, values, info=None, snpid=None):
        self.set_chromosome(chromosome)
        self.set_start(start)
        self.set_end(end)
        self.set_values(values)
        self.set_snpid(snpid)
        self.set_info(info)
    
    def get_chromosome(self):
        return self.chromosome

    def set_chromosome(self, chromosome):
        self.chromosome = chromosome

    def get_start(self):
        return self.start
    
    def set_start(self, coordinate):
        self.start = int(coordinate)
    
    def get_end(self):
        return self.end
    
    def set_end(self, coordinate):
        self.end = int(coordinate)
    
    def get_values(self):
        return self.values
    
    def set_values(self, values):
        self.values = values

    def get_snpid(self):
        return self.snpid
    
    def set_snpid(self, snpid):
        if snpid == None:
            self.snpid = ""
        else:
            self.snpid = snpid

    def get_info(self):
        return self.info
    
    def set_info(self, info):
        if info == None:
            self.info = {}
        else:
            self.info = info

    chromosome = property(get_chromosome, set_chromosome)
    start = property(get_start, set_start)
    end = property(get_end, set_end)
    values = property(get_values, set_values)
    snpid = property(get_snpid, set_snpid)
    info = property(get_info, set_info)

def create_tabix_command(tabix_path, vcf_path, chromosome, start, end):
    # Example command for running tabix to fetch one coordinate:
    # tabix file.vcf.gz chr1:12345-12345
    return tabix_path + " " + vcf_path + " chr" + str(chromosome) + ":" + str(start) + "-" + str(end)

def parse_vcf_line(row):
    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT    <values ..... >
    VCF_SNP_ID_COLUMN_INDEX = 2
    VCF_INFO_COLUMN_INDEX = 7
    VCF_VALUE_START_INDEX = 9

    split_row = row.split('\t')
    chromosome, coordinate = split_row[:2]
    snpid = split_row[VCF_SNP_ID_COLUMN_INDEX]
    info_field = split_row[VCF_INFO_COLUMN_INDEX]
    values = map(lambda x: x.strip(), split_row[VCF_VALUE_START_INDEX:])
    return MultilineTabixResult(chromosome, coordinate, coordinate, values, info=info_field, snpid=snpid)

def parse_region_lookup_result(data):
    reader = csv.DictReader(StringIO.StringIO(data), delimiter='\t')
    values = []
    for row in reader:
        result = {}
        # Remove leading '#' characters from field keys.
        for key, v in row.iteritems():
            if key.startswith('#'):
                result[key.lstrip('#')] = v
            else:
                result[key] = v

        values.append(result)

    return values

def vcf_singleline_lookup(tabix_path, vcf_path, chromosome, start_coordinate, end_coordinate):
    coordinate = start_coordinate
    command = create_tabix_command(tabix_path, vcf_path, chromosome, coordinate, coordinate)
    
    tabix_output = check_output(command.split())
    result = parse_vcf_line(tabix_output)
    if result.chromosome[3:] != chromosome or result.start != int(coordinate):
        errmsg = "Asked for " + str(chromosome) + ":" + str(coordinate) + ", got " + str(result.chromosome) + ":" + str(result.coordinate)
        raise Exception(errmsg)

    return result

def parse_triotype_line(row):
    # CHR POS <values ...>
    TRIOTYPE_VALUE_START_INDEX = 2

    split_row = row.split('\t')
    chromosome, coordinate = split_row[:2]
    values = map(lambda x: x.strip(), split_row[TRIOTYPE_VALUE_START_INDEX:])
    return MultilineTabixResult(chromosome, coordinate, coordinate, values)

def triotype_singleline_lookup(tabix_path, tsv_path, chromosome, start_coordinate, end_coordinate):
    coordinate = start_coordinate
    command = create_tabix_command(tabix_path, tsv_path, chromosome, coordinate, coordinate)
    
    tabix_output = check_output(command.split())
    result = parse_triotype_line(tabix_output)
    if result.chromosome[3:] != chromosome or result.start != int(coordinate):
        errmsg = "Asked for " + str(chromosome) + ":" + str(coordinate) + ", got " + str(result.chromosome) + ":" + str(result.coordinate)
        raise Exception(errmsg)

    return result

def tsv_region_lookup(tabix_path, tsv_path, chromosome, start, end):
    command = create_tabix_command(tabix_path, tsv_path, chromosome, start, end)
    
    tabix_output = check_output(command.split())
    values = parse_region_lookup_result(tabix_output)
    
    return MultilineTabixResult(chromosome, start, end, values)
