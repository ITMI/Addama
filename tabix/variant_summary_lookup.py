from collections import Counter

import argparse

from tabix_utils import vcf_singleline_lookup_with_header, triotype_singleline_lookup_with_header
import feature_data_source

CHROMOSOME_SET = frozenset([str(x) for x in xrange(1, 23)] + list(['M', 'X', 'Y']))
TRIO_TYPES = ['12','21','22','31','40','41','42','50','51','60','NA','MIE']

FEATURE_CATEGORIES = ["false", "true"]

class BarGroup:
    def __init__(self, type, label):
        self.type = type
        self.label = label

class VCFValue:
    def __init__(self, id, value):
        self.id = id
        self.value = value

DEFAULT_COUNT_MAP = {
    '0/1': '0/1 1/0',
    '1/0': '0/1 1/0'
}

DEFAULT_BAR_GROUPS = [
    BarGroup('0/0', '0/0'),
    BarGroup('0/1 1/0', '0/1'),
    BarGroup('1/1', '1/1'),
    BarGroup('./.', './.')
]

def process_triotypes(data, feature):
    plot_data = []
    category_trio_types = {}
    values = [VCFValue(x[0], x[1]) for x in data.values.iteritems()]
    for category in FEATURE_CATEGORIES:
        category_values = [x.value for x in values if feature.values[x.id] == category]
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
            family_id = v.id.rsplit('-', 1)[0]
            if feature.values[family_id] == category:
                if v.value in count_map:
                    category_values.append(count_map[v.value])
                else:
                    category_values.append(v.value)

        category_family_members[category] = Counter(category_values)

    for groupinfo in bar_groups:
        cat_data = []
        grouptype = groupinfo.type

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
            'trio_type': groupinfo.label,
            'categories': cat_data
        })

    return {
        'category_sizes': feature.aggregate,
        'plot_data': plot_data
    }

def process_vcf_F(filtered_values, feature, chromosome_digit):
    if chromosome_digit in ['X', 'Y', 'M']:
        count_map = {}

        bar_groups = [
            BarGroup('0', 'Reference'),
            BarGroup('1', 'Non-Reference'),
            BarGroup('.', 'Missing data')
        ]

        return process_vcf_sample_data(filtered_values, feature, count_map, bar_groups)
    else:
        return process_vcf_sample_data(filtered_values, feature, DEFAULT_COUNT_MAP, DEFAULT_BAR_GROUPS)


def process_vcf_NB_or_M(filtered_values, feature):
    return process_vcf_sample_data(filtered_values, feature, DEFAULT_COUNT_MAP, DEFAULT_BAR_GROUPS)

def process_vcf(data, feature, chromosome_digit):
    values = [VCFValue(x[0], x[1]) for x in data.values.iteritems()]

    f_result = process_vcf_F([v for v in values if v.id.split('-')[2] == 'F'], feature, chromosome_digit)
    m_result = process_vcf_NB_or_M([v for v in values if v.id.split('-')[2] == 'M'], feature)
    nb_result = process_vcf_NB_or_M([v for v in values if v.id.split('-')[2].startswith('NB')], feature)

    return {
        'f': f_result,
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
    tabix_exe = configuration['tabix_executable']

    feature = feature_data_source.get_feature_by_id(configuration['feature_matrix'], feature_id)

    triotypes = get_triotype_data(tabix_exe, configuration['triotype_file'], chromosome, coordinate)
    triotype_response = process_triotypes(triotypes, feature)

    vcf_data = get_vcf_data(tabix_exe, configuration['vcf_file'], chromosome, coordinate)
    vcf_response = process_vcf(vcf_data, feature, chromosome_digit)

    return {
        'triotypes': triotype_response,
        'vcf': vcf_response
    }

def main():
    mainparser = argparse.ArgumentParser(description="Variant summary service")

    subparsers = mainparser.add_subparsers()
    cmd_line_parser = subparsers.add_parser('query', help="Read all parameters from command line")

    cmd_line_parser.add_argument('--fmx', nargs=1, help='Feature matrix JSON path')
    cmd_line_parser.add_argument('--chr', nargs=1, help='Chromosome')
    cmd_line_parser.add_argument('--coord', nargs=1, type=int, help='Coordinate')
