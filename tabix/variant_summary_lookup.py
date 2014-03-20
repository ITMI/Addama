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
        category_values = [x['value'] for x in values if feature.values[x['id']] == category]
        category_trio_types[category] = Counter(category_values)

    for mtype in TRIO_TYPES:
        cat_data = []

        for category in FEATURE_CATEGORIES:
            counts_by_trio_type = category_trio_types[category]
            category_size = feature.aggregate[category]
            num = 0
            percentage = 0

            if mtype in counts_by_trio_type:
                num = counts_by_trio_type[mtype]

            if category_size > 0:
                percentage = float(num) / float(category_size)

            cat_data.append({
                'name': category,
                'value': percentage,
                'count': num
            })

        plot_data.append({
            'trio_type': mtype,
            'categories': cat_data
        })

    return {
        'category_sizes': feature.aggregate,
        'plot_data': plot_data
    }

def process_vcf_sample_data(values, feature, count_map, bar_groups):
    plot_data = []
    category_family_members = {}

    for category in FEATURE_CATEGORIES:
        category_values = []

        for v in values:
            family_id = v['id'].rsplit('-', 1)[0]
            if feature.values[family_id] == category:
                val = v['value']
                if val in count_map:
                    category_values.append(count_map[val])
                else:
                    category_values.append(val)

        category_family_members[category] = Counter(category_values)

    for groupinfo in bar_groups:
        cat_data = []
        grouptype = groupinfo['type']

        for category in FEATURE_CATEGORIES:
            counts_by_trio_type = category_family_members[category]
            num = 0
            category_size = feature.aggregate[category]
            percentage = 0

            if grouptype in counts_by_trio_type:
                num = counts_by_trio_type[grouptype]

            if category_size > 0:
                percentage = float(num) / float(category_size)

            cat_data.append({
                'name': category,
                'value': percentage,
                'count': num
            })

        plot_data.append({
            'trio_type': groupinfo['label'],
            'categories': cat_data
        })

    return {
        'category_sizes': feature.aggregate,
        'plot_data': plot_data
    }

def process_vcf_NB_or_M(filtered_values, feature):
    count_map = {
        '0/1': '0/1 1/0',
        '1/0': '0/1 1/0'
    }

    bar_groups = [
        {
            'type': '0/0',
            'label': '0/0'
        },
        {
            'type': '0/1 1/0',
            'label': '0/1'
        },
        {
            'type': '1/1',
            'label': '1/1'
        },
        {
            'type': './.',
            'label': './.'
        }
    ]

    return process_vcf_sample_data(filtered_values, feature, count_map, bar_groups)

def process_vcf(data, feature):
    values = [{'id': x[0], 'value': x[1]} for x in data.values.iteritems()]

    m_result = process_vcf_NB_or_M([v for v in values if v['id'].split('-')[2] == 'M'], feature)
    nb_result = process_vcf_NB_or_M([v for v in values if v['id'].split('-')[2].startswith('NB')], feature)

    return {
        'm': m_result,
        'nb': nb_result
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

#############
# Query API #
#############

def do_query(configuration, chromosome_digit, coordinate, feature_id):
    chromosome = str(chromosome_digit)
    tabix_exe = "tabix"

    feature = feature_data_source.get_feature_by_id(configuration['feature_matrix'], feature_id)

    triotypes = get_triotype_data(tabix_exe, configuration['triotype_file'], chromosome, coordinate)
    triotype_response = process_triotypes(triotypes, feature)

    vcf_data = get_vcf_data(tabix_exe, configuration['vcf_file'], chromosome, coordinate)
    vcf_response = process_vcf(vcf_data, feature)

def main():
    mainparser = argparse.ArgumentParser(description="Variant summary service")

    subparsers = mainparser.add_subparsers()
    cmd_line_parser = subparsers.add_parser('query', help="Read all parameters from command line")

    cmd_line_parser.add_argument('--fmx', nargs=1, help='Feature matrix JSON path')
    cmd_line_parser.add_argument('--chr', nargs=1, help='Chromosome')
    cmd_line_parser.add_argument('--coord', nargs=1, type=int, help='Coordinate')
