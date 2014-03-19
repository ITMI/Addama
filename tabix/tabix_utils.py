import csv
import StringIO
import subprocess

# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT    <values ..... >
VCF_SNP_ID_COLUMN_INDEX = 2
VCF_REF_COLUMN_INDEX = 3
VCF_ALT_COLUMN_INDEX = 4
VCF_INFO_COLUMN_INDEX = 7
VCF_VALUE_START_INDEX = 9

def query_description(chromosome, start, end):
    return "chr" + str(chromosome) + ":" + str(start) + "-" + str(end)

class CoordinateRangeEmptyError(Exception):
    def __init__(self, chromosome, start, end, msg):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.msg = msg
    def __str__(self):
        query = query_description(self.chromosome, self.start, self.end)
        return "Tabix - empty result. Query " + query + ". Info: " + repr(self.msg)

class WrongLineFoundError(Exception):
    def __init__(self, chromosome, start, end, msg):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.msg = msg
    def __str__(self):
        query = query_description(self.chromosome, self.start, self.end)
        return "Tabix - wrong line found. Query " + query + ". Info: " + repr(self.msg)

class UnexpectedTabixOutputError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return "Tabix - unexpected output: " + repr(self.msg)

class TabixExecutionError(Exception):
    def __init__(self, error, output=""):
        self.error = error
        self.output = output
    def __str__(self):
        return "Tabix - execution failed: " + repr(self.output)

class MultilineTabixResult():
    def __init__(self, chromosome, start, end, values, info=None, snpid=None, ref=None, alt=None):
        self.set_chromosome(chromosome)
        self.set_start(start)
        self.set_end(end)
        self.set_values(values)
        self.set_snpid(snpid)
        self.set_ref(ref)
        self.set_alt(alt)
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

    def get_ref(self):
        return self.ref
    
    def set_ref(self, ref):
        if ref == None:
            self.ref = ""
        else:
            self.ref = ref

    def get_alt(self):
        return self.alt
    
    def set_alt(self, alt):
        if alt == None:
            self.alt = ""
        else:
            self.alt = alt
            
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
    ref = property(get_ref, set_ref)
    alt = property(get_alt, set_alt)
    info = property(get_info, set_info)

def create_tabix_command(tabix_path, vcf_path, chromosome, start, end):
    # Example command for running tabix to fetch one coordinate:
    # tabix file.vcf.gz chr1:12345-12345
    return tabix_path + " " + vcf_path + " chr" + str(chromosome) + ":" + str(start) + "-" + str(end)

def parse_vcf_line(row, identifiers=None):

    split_row = row.split('\t')
    chromosome, coordinate = split_row[:2]
    snpid = split_row[VCF_SNP_ID_COLUMN_INDEX]
    ref = split_row[VCF_REF_COLUMN_INDEX]
    alt = split_row[VCF_ALT_COLUMN_INDEX]
    info_field = split_row[VCF_INFO_COLUMN_INDEX]
    values = map(lambda x: x.strip(), split_row[VCF_VALUE_START_INDEX:])
    result_values = None

    if identifiers is None:
        result_values = values
    else:
        result_values = {key: value for (key, value) in zip(identifiers, values)}

    return MultilineTabixResult(chromosome,
                                coordinate,
                                coordinate,
                                result_values,
                                info=info_field,
                                snpid=snpid,
                                ref=ref,
                                alt=alt)

def parse_region_lookup_result(data, header_line_identifier='#', line_filter="##"):
    lines = [line for line in data if not line.startswith(line_filter)]

    reader = csv.DictReader(StringIO.StringIO(data), delimiter='\t')
    values = []
    for row in reader:
        result = {}
        # Remove leading header line identifier characters from field keys.
        for key, v in row.iteritems():
            if key.startswith(header_line_identifier):
                result[key.lstrip(header_line_identifier)] = v
            else:
                result[key] = v

        values.append(result)

    return values

def split_and_remove_empty_lines(data):
    lines = data.split('\n')
    return [line for line in lines if line.strip() != '']

def vcf_singleline_lookup(tabix_path, vcf_path, chromosome, start_coordinate, end_coordinate):
    coordinate = start_coordinate
    command = create_tabix_command(tabix_path, vcf_path, chromosome, coordinate, coordinate)
    
    tabix_output = None

    try:
        tabix_output = subprocess.check_output(command.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        raise TabixExecutionError(str(cpe), cpe.output)

    output = split_and_remove_empty_lines(tabix_output)
    if (len(output) == 0):
        raise CoordinateRangeEmptyError(chromosome, start_coordinate, end_coordinate, "VCF lookup")

    if (len(output) > 1):
        raise UnexpectedTabixOutputError("expected 1 line, found " + str(len(output)))

    result = parse_vcf_line(tabix_output)
    if result.chromosome[3:] != chromosome or result.start != int(coordinate):
        errmsg = "got " + str(result.chromosome) + ":" + str(result.start) + ", full line \'" + repr(output) + "\'"
        raise WrongLineFoundError(chromosome, start_coordinate, end_coordinate, errmsg)

    return result

def get_vcf_lookup_identifiers_from_header(header_line):
    split_header_line = [value.strip() for value in header_line.split('\t')]

    identifiers = split_header_line[VCF_VALUE_START_INDEX:]

    return identifiers

def vcf_singleline_lookup_with_header(tabix_path, vcf_path, chromosome, start_coordinate, end_coordinate):
    # Use '-h' flag to include headers in the tabix output.
    tabix_path = tabix_path + " -h"

    coordinate = start_coordinate
    command = create_tabix_command(tabix_path, vcf_path, chromosome, coordinate, coordinate)

    tabix_output = None

    try:
        tabix_output = subprocess.check_output(command.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        raise TabixExecutionError(str(cpe), cpe.output)

    line_filter = "##"
    output = split_and_remove_empty_lines(tabix_output)
    output = [line for line in output  if not line.startswith(line_filter)]

    if (len(output) == 1):
        raise CoordinateRangeEmptyError(chromosome, start_coordinate, end_coordinate, "VCF lookup")

    if (len(output) > 2):
        raise UnexpectedTabixOutputError("expected 2 lines, found " + str(len(output)))

    identifiers = get_vcf_lookup_identifiers_from_header(output[0])
    result = parse_vcf_line(output[1], identifiers)
    if result.chromosome[3:] != chromosome or result.start != int(coordinate):
        errmsg = "got " + str(result.chromosome) + ":" + str(result.start) + ", full line \'" + repr(output) + "\'"
        raise WrongLineFoundError(chromosome, start_coordinate, end_coordinate, errmsg)

    return result

def parse_triotype_line(row, result_dict=False):
    # CHR POS <values ...>
    TRIOTYPE_VALUE_START_INDEX = 2

    split_row = row.split('\t')
    chromosome, coordinate = split_row[:2]
    values = map(lambda x: x.strip(), split_row[TRIOTYPE_VALUE_START_INDEX:])
    return MultilineTabixResult(chromosome, coordinate, coordinate, values)

def triotype_singleline_lookup(tabix_path, tsv_path, chromosome, start_coordinate, end_coordinate):
    coordinate = start_coordinate
    command = create_tabix_command(tabix_path, tsv_path, chromosome, coordinate, coordinate)

    tabix_output = None

    try:
        tabix_output = subprocess.check_output(command.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        raise TabixExecutionError(str(cpe), cpe.output)

    output = split_and_remove_empty_lines(tabix_output)
    if (len(output) == 0):
        raise CoordinateRangeEmptyError(chromosome, start_coordinate, end_coordinate, "triotype lookup")

    if (len(output) > 1):
        raise UnexpectedTabixOutputError("expected 1 line, found " + str(len(output)))

    result = parse_triotype_line(tabix_output)

    if result.chromosome[3:] != chromosome or result.start != int(coordinate):
        errmsg = "got " + str(result.chromosome) + ":" + str(result.start) + ", full line \'" + repr(output) + "\'"
        raise WrongLineFoundError(chromosome, start_coordinate, end_coordinate, errmsg)

    return result

def parse_triotype_lookup_with_headers(data, header_line_identifier='#', line_filter="##"):
    lines = [line for line in data if not line.startswith(line_filter)]

    header_line = lines[0]
    split_header_line = [value.strip() for value in header_line.split('\t')]

    value_line = lines[1]
    split_value_line = [value.strip() for value in value_line.split('\t')]

    chromosome, coordinate = split_value_line[:2]

    identifiers = split_header_line[2:]
    triotype_values = split_value_line[2:]

    values = {key: value for (key, value) in zip(identifiers, triotype_values)}

    return MultilineTabixResult(chromosome, coordinate, coordinate, values)

def triotype_singleline_lookup_with_header(tabix_path, tsv_path, chromosome, start_coordinate, end_coordinate):
    # Use '-h' flag to include headers in the tabix output.
    tabix_path = tabix_path + " -h"

    coordinate = start_coordinate
    command = create_tabix_command(tabix_path, tsv_path, chromosome, coordinate, coordinate)

    tabix_output = None

    try:
        tabix_output = subprocess.check_output(command.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as cpe:
        raise TabixExecutionError(str(cpe), cpe.output)

    output = split_and_remove_empty_lines(tabix_output)
    if (len(output) == 1):
        raise CoordinateRangeEmptyError(chromosome, start_coordinate, end_coordinate, "triotype lookup")

    if (len(output) > 2):
        raise UnexpectedTabixOutputError("expected 2 lines, found " + str(len(output)))

    result = parse_triotype_lookup_with_headers(output, header_line_identifier='#', line_filter='##')

    if result.chromosome[3:] != chromosome or result.start != int(coordinate):
        errmsg = "got " + str(result.chromosome) + ":" + str(result.start) + ", full line \'" + repr(output) + "\'"
        raise WrongLineFoundError(chromosome, start_coordinate, end_coordinate, errmsg)

    return result

def tsv_region_lookup(tabix_path, tsv_path, chromosome, start, end):
    command = create_tabix_command(tabix_path, tsv_path, chromosome, start, end)
    tabix_output = subprocess.check_output(command.split())
    values = parse_region_lookup_result(tabix_output)
    
    return MultilineTabixResult(chromosome, start, end, values)
