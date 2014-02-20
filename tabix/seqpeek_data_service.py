import sys
import argparse
import json
from pprint import pprint
from pymongo import MongoClient
from pymongo import ASCENDING

import tabix_utils as TU

class SeqPeekResult():
    def __init__(self, config):
        self.set_config(config)

    def __repr__(self):
        result = "<SeqPeekResult "
        for (k,v) in self.__dict__.items():
            result+=k+":"+str(v)+" "
        return result+">"

    def __str__(self):
        result = "<SeqPeekResult "
        for (k,v) in self.__dict__.items():
            result+=k+":"+str(v)+" "
        return result+">"

    def get_chromosome(self):
        return self.chromosome

    def set_chromosome(self, chromosome):
        self.chromosome = chromosome

    def get_start(self):
        return self.start
    
    def set_start(self, start):
        self.start = int(start)

    def get_end(self):
        return self.end

    def set_end(self, end):
        self.end = int(end)

    def get_config(self):
        return self.config

    def set_config(self, config):
        self.config = config

    def get_features(self):
        return self.features

    def set_features(self,features):
        self.features = features

    def get_gene_info(self, gene_info):
        self.gene_info = gene_info

    def set_gene_info(self, gene_info):
        return self.gene_info

    def get_full_gene(self, full_gene):
        self.full_gene = full_gene

    def set_full_gene(self, full_gene):
        return self.full_gene

    chromosome = property(get_chromosome, set_chromosome)
    start = property(get_start, set_start)
    end = property(get_end, set_end)
    config = property(get_config, set_config)
    features = property(get_features, set_features)
    gene_info = property(get_gene_info, set_gene_info)
    full_gene = property(get_full_gene, set_full_gene)

mDEBUG = 0

def set_parameters(seqObj, args):
    if 'gene_name' in args and args.gene_name is not None:
        seqObj.gene_name = str(args.gene_name[0])

    if 'chr' in args and args.chr is not None:
        seqObj.chromosome = str(args.chr[0])
        seqObj.start = str(args.start[0])
        seqObj.end = str(args.end[0])

    if 'full_gene' in args and args.full_gene is not None:
        seqObj.full_gene = args.full_gene # true/false
    
    if 'debug' in args and args.debug is not False:
        global mDEBUG
        mDEBUG = 1

# if we want to have this print to log file, can change this function
def control_print(message):
    global mDEBUG
    if mDEBUG:
        print("mdebug is "+str(mDEBUG))
        print(message)

def control_pprint(message):
    global mDEBUG
    if mDEBUG:
        pprint(message)

def default_config():
    control_print("default config called")
    config = {}
    config['tabix_executable'] = 'tabix'
    config['variant_file'] = "/u5/workspace/dmauldin/SeqPeek/SeqPeekData/variants/out3.bgz"
    config['region_data_host'] = "mongodb://bobama.systemsbiology.net"
    config['region_data_database'] = "seqpeek_data_service"
    config['region_data_collection'] = 'region_data'
    config['feature_matrix_host'] = "mongodb://bobama.systemsbiology.net"
    config['feature_matrix_database'] = "seqpeek_data_service"
    config['feature_matrix_collection'] = 'feature_matrix'
    return SeqPeekResult(config)

def parse_config(seqObj,config_file):
    control_print("loading config file "+str(config_file))
    json_data = open(config_file)
    config = json.load(json_data)
    combined_config = seqObj.config
    for key in config:
        combined_config[key] = config[key]
    seqObj.config = combined_config
    return seqObj

def create_mongo_connection(uri):
    c = MongoClient(uri)
    return c


def get_gene_info(seqObj):
    client = create_mongo_connection(seqObj.config['region_data_host'])
    db = client[seqObj.config['region_data_database']]
    region_data_collection = db[seqObj.config['region_data_collection']]
    gene_data = {}

    if seqObj.gene_name and not hasattr(seqObj,'chromosome'):
        control_print("getting information based on gene name")
        gene_results = region_data_collection.find({"gene":str(seqObj.gene_name)}).sort("exonstart",ASCENDING)

        for gene in gene_results:
            if not hasattr(seqObj,'chromosome'):
                seqObj.chromosome = gene['chr']
            if not hasattr(seqObj,'start'):
                seqObj.start = gene['txstart']
            if not hasattr(seqObj,'end'):
                seqObj.end = gene['txend']
            if gene['txstart'] < seqObj.start:
                seqObj.start = gene['txstart']
            if gene['txend'] > seqObj.end:
                seqObj.end = gene['txend']


# needs to create a gene hash:
# hash{gene name}{transcript id}
# keys chromosome, start, stop, strand, exon_starts(array), exon_ends(array)
def get_region_data(seqObj):
    client = create_mongo_connection(seqObj.config['region_data_host'])
    db = client[seqObj.config['region_data_database']]
    region_data_collection = db[seqObj.config['region_data_collection']]
    gene_data = {}

    if hasattr(seqObj,'gene_name') and not hasattr(seqObj,'chromosome'):
        get_gene_info(seqObj)

    search_dict_entire_gene = {
                   "chr":str(seqObj.chromosome),
                   "$or":[
                        {"txstart":{"$gte": int(seqObj.start), "$lt": int(seqObj.end)}},
                        {"txend":  {"$gte": int(seqObj.start),   "$lt": int(seqObj.end)}}
                        ]
                  }
    search_dict_exons = {
                   "chr":str(seqObj.chromosome),
                   "$or":[
                        {"exonstart":{"$gt": int(seqObj.start), "$lt": int(seqObj.end)}},
                        {"exonend":  {"$gt": int(seqObj.start),   "$lt": int(seqObj.end)}}
                        ]
                  }

    control_print("calling region_data for chromosome %s start %s end %s" % (seqObj.chromosome, seqObj.start, seqObj.end))
    if (seqObj.full_gene):
        control_print(search_dict_entire_gene)
        query_results = region_data_collection.find(search_dict_entire_gene).sort("exonstart", ASCENDING)
    else:
        control_print(search_dict_exons)
        query_results = region_data_collection.find(search_dict_exons).sort("exonstart", ASCENDING)
    control_pprint(query_results)
    for region_data in query_results:
        control_print(region_data)
        gene = region_data['gene']
        mrna = region_data['mrna']
        start = region_data['exonstart']
        stop = region_data['exonend']
        chromosome = region_data['chr']
        strand = region_data['strand']
        if gene not in gene_data:
            gene_data[gene] = {}
        if mrna not in gene_data[gene]:
            gene_data[gene][mrna] = {}
        if 'start' not in gene_data[gene][mrna] or gene_data[gene][mrna]['start'] > start:
            gene_data[gene][mrna]['start'] = start
            control_print("set gene_start to region_data "+str(start))
        if 'stop' not in gene_data[gene][mrna] or gene_data[gene][mrna]['stop'] < stop:
            gene_data[gene][mrna]['stop'] = stop
        if 'exon_starts' not in gene_data[gene][mrna]:
            gene_data[gene][mrna]['exon_starts'] = []
        if 'exon_stops' not in gene_data[gene][mrna]:
            gene_data[gene][mrna]['exon_stops'] = []
        if 'strand' not in gene_data[gene][mrna]:
            gene_data[gene][mrna]['strand'] = strand
        if 'chromosome' not in gene_data[gene][mrna]:
            gene_data[gene][mrna]['chromosome'] = chromosome

        gene_data[gene][mrna]['exon_starts'].append(start)
        gene_data[gene][mrna]['exon_stops'].append(stop)
        control_print(gene_data)
    control_print("afterwards gene_data is "+str(gene_data))
    if not gene_data:
        control_print("Information not found for this region: %s:%s-%s" % (seqObj.chromosome, seqObj.start, seqObj.end))
        exit()
    seqObj.gene_info = gene_data

def get_features_from_json(seqObj):
    matrix_datasource_config = seqObj.config['feature_matrix']
    with open(matrix_datasource_config['path']) as json_file:
        json_data = json.load(json_file)
    feature_map = {}
    sample_ids = json_data['sample_id_array']
    feature_value_arrays = json_data['feature_value_array']
    feature_list = json_data['ordered_list']
    for index, (feature_id, text_name) in zip(xrange(len(feature_list)), feature_list):
        feature_map[feature_id] = dict(zip(sample_ids, feature_value_arrays[index]))

    seqObj.features = feature_map

# create a mongo connection to the feature_matrix_uri
# retrieve and store feature results
def get_features_from_mongodb(seqObj):
    matrix_datasource_config = seqObj.config['feature_matrix']

    client = create_mongo_connection(matrix_datasource_config['host'])
    db = client[matrix_datasource_config['database']]
    feature_matrix_collection = db[matrix_datasource_config['collection']]
    all_features = {}
    for feature_matrix in feature_matrix_collection.find():
        all_features[feature_matrix["id"]] = feature_matrix["v"]
    seqObj.features = all_features

def get_features(seqObj):
    matrix_datasource_type = seqObj.config['feature_matrix']['type']
    if matrix_datasource_type  == 'mongodb':
        get_features_from_mongodb(seqObj)
    elif matrix_datasource_type == 'json':
        get_features_from_json(seqObj)

# create the variants dict
# query the variants from the variant file
# calculate the statistics for each gene_name, transcript_id, variant
# results{gene_name}{transcript_id}{coordinate}{variant}
def tabix_query_variant(seqObj):
    tabix_variant_file = seqObj.config['variant_file']
    tabix_exe = seqObj.config['tabix_executable'] + " -h" 
    tabix_output = TU.tsv_region_lookup(tabix_exe, tabix_variant_file, seqObj.chromosome, seqObj.start, seqObj.end)

    features = seqObj.features
    # chromosome position gene_name transcript_id variant family_id variant_type uniprot_id zygosity amino_acid_variant
    #chr1    908929  PLEKHN1 NM_032129       C->T    101-641-NB      SUBSTITUTION    Q494U1-2  ''  heterozygous '' S448L
    counter = {}
    for row in tabix_output.values:
        if not row:
            break
        composite_key_fields = ['chr',
                                'coordinate',
                                'gene',
                                'transcript',
                                'variant',
                                'type',
                                'uniprot_id',
                                'protein_change']
        
        composite_key = "\t".join([(lambda key: row[key])(key) for key in composite_key_fields])
        id = row['sample_id']
        zygosity = row['genotype']
        control_print("processing composite_key %s id %s zygosity %s" % (composite_key, id, zygosity))
        if composite_key not in counter:
            counter[composite_key] = {}
        for feature in features:
            if feature not in counter[composite_key]:
                counter[composite_key][feature] = {}

            # take the member off the end of the id
            family_info = id.split('-')
            member_id = family_info.pop()
            family_id = ('-').join(family_info)
            # get the feature matrix value of this 
            #control_print("looking for feature %s id %s in features[feature] %s" % (feature, family_id, features[feature]))
            feature_value = features[feature][family_id] # should be true, false, or na
            #control_print("found feature_value %s for id %s for feature %s" % (feature_value, family_id, feature))
            #feature_value = choice(possible_values) # mock code for random true/false/na
            if feature_value not in counter[composite_key][feature]:
                counter[composite_key][feature][feature_value] = []
            counter[composite_key][feature][feature_value].append(id)
        
    control_pprint(counter)
    results = {}
    transcripts = {}

    for composite_key in counter:
        split_key = composite_key.split('\t')
        chromosome = split_key[0]
        coordinate = split_key[1]
        gene_name = split_key[2]
        transcript = split_key[3]
        variant = split_key[4]
        variant_type = split_key[5]
        uniprot = split_key[6]
        protein_change = split_key[7]

        
        if gene_name not in results:
            results[gene_name] = {}

        if transcript not in transcripts:
            transcripts[transcript] = {}

        if 'variants' not in transcripts[transcript]:
            transcripts[transcript]['variants'] = []

        # add the global results if they're not already there
        if 'start' not in results[gene_name]:
            gene_info = seqObj.gene_info
            if gene_name in gene_info:
                control_print("gene info")
                control_pprint(gene_info)

                results[gene_name]['chromosome'] = gene_info[gene_name][transcript]['chromosome']
                results[gene_name]['start'] = gene_info[gene_name][transcript]['start']
                results[gene_name]['stop'] = gene_info[gene_name][transcript]['stop']
                results[gene_name]['strand'] = gene_info[gene_name][transcript]['strand']
                results[gene_name]['exon_starts'] = gene_info[gene_name][transcript]['exon_starts']
                results[gene_name]['exon_stops'] = gene_info[gene_name][transcript]['exon_stops']

                control_print("results")
                control_pprint(results)

        local_result = {}
        local_result['base_change'] = variant
        local_result['coordinate'] = int(coordinate)
        local_result['uniprot_id'] = uniprot
        local_result['protein_change'] = protein_change
        local_result['variant_type'] = variant_type
        local_result['statistics'] = {}
        
        for feature in features:
            control_print("feature is "+feature)
            local_statistic = {}
            
            for count_type in counter[composite_key][feature]:
                local_statistic[count_type] = len(counter[composite_key][feature][count_type])
            local_result['statistics'][feature] = local_statistic
        control_pprint(local_result)
        transcripts[transcript]['id'] = transcript
        transcripts[transcript]['variants'].append(local_result)

    results[gene_name]['transcripts'] = transcripts

    control_print("results")
    control_pprint(results)
#    if result.chromosome[3:] != chromosome or result.coordinate != int(coordinate):
#        errmsg = "Asked for " + str(chromosome) + ":" + str(coordinate) + ", got " + str(result.chromosome) + ":" + str(result.coordinate)
#        raise Exception(errmsg)

    return results


def error(self,message):
    sys.stderr.write('Error: %s\n' % message)
    self.print_help()
    sys.exit(2)

def example_object():
    example = {}
    example['OR4F5'] = {}
    example['OR4F5']['NM_001005484'] = {}
    example['OR4F5']['NM_001005484']['chromosome'] = '1'
    example['OR4F5']['NM_001005484']['start'] = '69090'
    example['OR4F5']['NM_001005484']['end'] = '70008'
    example['OR4F5']['NM_001005484']['strand'] = '+'
    example['OR4F5']['NM_001005484']['exon_starts'] = ['69090']
    example['OR4F5']['NM_001005484']['exon_stops'] = ['79501']
    example['OR4F5']['NM_001005484']['variants'] = []
    variants = {}
    variants['base_change'] = 'A->G'
    variants['coordinate'] = '69511'
    variants['type'] = 'SUBSTITUTION'
    variants['uniprot'] = 'Q8NH21'
    variants['protein_change'] = 'T141A'
    variants['statistics'] = {}
    variants['statistics']['feature_id'] = "B:MRGE:Uterine_Related:NB::::"
    variants['statistics']['true'] = 63
    variants['statistics']['false'] = 30
    variants['statistics']['na'] = 90
    example['OR4F5']['NM_001005484']['variants'].append(variants)
    return example


def do_gene_query(seq_obj):
    get_region_data(seq_obj)
    get_features(seq_obj)
    return tabix_query_variant(seq_obj)
    return {}


def main():
    mainparser = argparse.ArgumentParser(description="SeqPeekDataService")

    mainparser.add_argument('--config',nargs=1,help='Config File')
    mainparser.add_argument('--chr', nargs=1, help='Chromosome')
    mainparser.add_argument('--start', nargs=1, type=int, help='Coordinate Start')
    mainparser.add_argument('--end', nargs=1, type=int, help='Coordinate End')
    mainparser.add_argument('--gene_name', nargs=1, help='Gene Name')
    mainparser.add_argument('--full_gene', help='Return entire gene information (not variants) for partial overlap of coordinates',action="store_true")
    mainparser.add_argument('--debug', help='Print debug information',action="store_true")
    mainparser.add_argument('--example',help="Produce an Example object",action="store_true")

    try:
        args = mainparser.parse_args()
    except IOError, msg:
        mainparser.error(str(msg))
        exit()

    if len(sys.argv) < 2:
        error(mainparser, 'Must provide either a gene_name or a chr, start, end')
    if args.gene_name is not None and args.chr is not None:
        error(mainparser, 'Must provide either a gene_name or a chr,start,end but not both')
    if (args.chr is not None and args.start is None) or (args.chr is not None and args.end is None): 
        error(mainparser, 'Must provide chr, start, and end')

    seqObj = default_config()
    set_parameters(seqObj,args) 
    control_pprint("seqObj default is "+str(seqObj))
    control_print(args)
    if 'config' in args and args.config is not None:
        seqObj = parse_config(seqObj,args.config[0])
        control_print("seqObj with configs is"+str(seqObj))

    if 'example' in args and args.example is True:
        control_print("running example")
        tabix_result = example_object()
        control_pprint(tabix_result)
        json.dumps(tabix_result)
    else:
        get_region_data(seqObj)
        get_features(seqObj)
        results = tabix_query_variant(seqObj)
        control_print(seqObj)
        control_print("\n\n")
        control_pprint(results)
        results_json_string = json.dumps(results)
        print(results_json_string)

if __name__ == "__main__":
    main()










