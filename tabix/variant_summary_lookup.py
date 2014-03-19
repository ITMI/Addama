from collections import Counter

import argparse

from tabix_utils import vcf_singleline_lookup_with_header, triotype_singleline_lookup_with_header
import feature_data_source
from tabix_utils import CoordinateRangeEmptyError, WrongLineFoundError, TabixExecutionError, UnexpectedTabixOutputError

CHROMOSOME_SET = frozenset([str(x) for x in xrange(1, 23)] + list(['M', 'X', 'Y']))
TRIO_TYPES = ['12','21','22','31','40','41','42','50','51','60','NA','MIE']

FEATURE_CATEGORIES = ["false", "true"]

def process_triotypes(data, feature):
    plot_data = []
    category_trio_types = {}
    values = [{'id': x[0], 'value': x[1]} for x in data.values.iteritems()]
    for category in FEATURE_CATEGORIES:
        category_values = [x['value'] for x in values if feature[x['id']] == category]
        category_trio_types[category] = Counter(category_values)

    for mtype in TRIO_TYPES:
        cat_data = []

        for category in FEATURE_CATEGORIES:
            counts_by_trio_type = category_trio_types[category]
            num = 0

            if mtype in counts_by_trio_type:
                num = counts_by_trio_type[mtype]

            cat_data.append({
                'name': category,
                'count': num
            })

        plot_data.append({
            'trio_type': mtype,
            'categories': cat_data
        })

    return {
        'plot_data': plot_data
    }

########################
# Data loading methods #
########################

def get_triotype_data(tabix_exe, data_file_path, chromosome, coordinate):
    result = triotype_singleline_lookup_with_header(tabix_exe, data_file_path, chromosome, coordinate, coordinate)
    return result

def get_vcf_data(tabix_exe, data_file_path, chromosome, coordinate):
    result = vcf_singleline_lookup_with_header(tabix_exe, data_file_path, chromosome, coordinate, coordinate)
    return result

def do_query(configuration, chromosome_digit, coordinate, feature_id):
    chromosome = str(chromosome_digit)
    tabix_exe = "tabix"
    triotypes = get_triotype_data(tabix_exe, configuration['triotype_file'], chromosome, coordinate)
    vcf_data = get_vcf_data(tabix_exe, configuration['vcf_file'], chromosome, coordinate)

    feature = feature_data_source.get_feature_by_id(configuration['feature_matrix'], feature_id)
    triotype_response = process_triotypes(triotypes, feature)

def main():
    mainparser = argparse.ArgumentParser(description="Variant summary service")

    subparsers = mainparser.add_subparsers()
    cmd_line_parser = subparsers.add_parser('query', help="Read all parameters from command line")

    cmd_line_parser.add_argument('--fmx', nargs=1, help='Feature matrix JSON path')
    cmd_line_parser.add_argument('--chr', nargs=1, help='Chromosome')
    cmd_line_parser.add_argument('--coord', nargs=1, type=int, help='Coordinate')
