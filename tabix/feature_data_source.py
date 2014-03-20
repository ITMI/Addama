from collections import Counter

import json

from pymongo import MongoClient

class FeatureData:
    def __init__(self, feature_id, value_dict):
        self.id = feature_id
        self.values = value_dict
        self.value_array = [{'id': x[0], 'value': x[1]} for x in value_dict.iteritems()]
        self.aggregate = Counter(value_dict.values())

def create_mongo_connection(uri):
    c = MongoClient(uri)
    return c

def get_feature_by_id_from_json(config, param_feature_id):
    with open(config['path']) as json_file:
        json_data = json.load(json_file)
    feature_map = {}
    sample_ids = json_data['sample_id_array']
    feature_value_arrays = json_data['feature_value_array']
    feature_list = json_data['ordered_list']
    for index, (feature_id, text_name) in zip(xrange(len(feature_list)), feature_list):
        feature_map[feature_id] = dict(zip(sample_ids, feature_value_arrays[index]))

    return feature_map[param_feature_id]

def get_feature_by_id_from_mongodb(config, feature_id):

    client = create_mongo_connection(config['host'])
    db = client[config['database']]
    feature_matrix_collection = db[config['collection']]
    all_features = {}
    query = {
        "id": feature_id
    }
    for feature in feature_matrix_collection.find(query):
        all_features[feature["id"]] = feature["v"]

    return all_features[feature_id]

def get_feature_by_id(config, feature_id):
    matrix_datasource_type = config['type']
    value_dict = None
    if matrix_datasource_type  == 'mongodb':
        value_dict = get_feature_by_id_from_mongodb(config, feature_id)
    elif matrix_datasource_type == 'json':
        value_dict = get_feature_by_id_from_json(config, feature_id)

    return FeatureData(feature_id, value_dict)
