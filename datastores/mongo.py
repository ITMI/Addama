import logging

from tornado.options import options

from pymongo import MongoClient
from pymongo.errors import ConnectionFailure

import tornado.web
import tornado.escape
import csv
import re

RESERVED_KEYS = ["output", "output_filename", "sort_by", "sort_direction"]

class MongoDbQueryHandler(tornado.web.RequestHandler):
    datastores_config = {}

    def initialize(self):
        self._datastore_map = self.datastores

    def get(self, *uri_path):
        try:
            if options.verbose: logging.info("GET [uri=%s] [arguments=%s]" % (self.request.uri, self.request.arguments))

            sub_path = self.request.path.replace("/datastores", "")
            uri_parts = sub_path.split("/")
            if options.verbose: logging.info("GET [sub_path=%s] [len=%d]" % (sub_path, len(uri_parts)))

            if len(uri_parts) == 1:
                self.list_datastores()
                self.set_status(200)
                return

            datastore_id = uri_parts[1]
            if not datastore_id in self._datastore_map.keys():
                if options.verbose: logging.info("unknown datastore [%s]" % datastore_id)
                raise tornado.web.HTTPError(404)

            if len(uri_parts) == 2:
                self.list_databases(datastore_id)
                self.set_status(200)
                return

            db_name = uri_parts[2]
            if len(uri_parts) == 3:
                self.list_collections(datastore_id, db_name)
                self.set_status(200)
                return

            collection_id = uri_parts[3]
            datastore = self._datastore_map[datastore_id]
            collection = self.open_collection(datastore_id, db_name, collection_id)
            if len(uri_parts) == 4:
                datatypes = self.get_datatypes(datastore_id, db_name, collection_id)
                query = self.transpose_query_arguments(db_name, datastore, datatypes)
                sort_fld = self.get_argument("sort_by", None)
                sort_dir = int(self.get_argument("sort_direction", "1"))
                json_items = self.query_collection(collection, query, sort_fld, sort_dir)

                if self.get_argument("output", "json") == "tsv":
                    self.write_tsv(json_items)
                    self.set_status(200)
                    return

                self.write({"items": json_items})
                self.set_status(200)
                return

            last_part = uri_parts[4]
            if last_part == "fields":
                self.list_fields(collection)
                self.set_status(200)
                return

            raise tornado.web.HTTPError(404)
        except ConnectionFailure as cfe:
            raise tornado.web.HTTPError(500, str(cfe))

    def post(self, *uri_path):
        try:
            if options.verbose: logging.info("POST [uri=%s] [body=%s]" % (self.request.uri, self.request.body))

            sub_path = self.request.path.replace("/datastores", "")
            uri_parts = sub_path.split("/")
            if options.verbose: logging.info("POST [sub_path=%s] [len=%d]" % (sub_path, len(uri_parts)))

            if len(uri_parts) == 1:
                self.list_datastores()
                self.set_status(200)
                return

            datastore_id = uri_parts[1]
            if not datastore_id in self._datastore_map.keys():
                if options.verbose: logging.info("unknown datastore [%s]" % datastore_id)
                raise tornado.web.HTTPError(404)

            if len(uri_parts) == 2:
                self.list_databases(datastore_id)
                self.set_status(200)
                return

            db_name = uri_parts[2]
            if len(uri_parts) == 3:
                self.list_collections(datastore_id, db_name)
                self.set_status(200)
                return
                
            collection_id = uri_parts[3]
            # datastore = self._datastore_map[datastore_id]
            collection = self.open_collection(datastore_id, db_name, collection_id)

            if len(uri_parts) == 4:
                datatypes = self.get_datatypes(datastore_id, db_name, collection_id)
                query = tornado.escape.json_decode(self.request.body) 
                sort_fld = query.pop("sort_by", None)
                sort_dir = query.pop("sort_direction", 1)
                json_items = self.query_collection(collection, query, sort_fld, sort_dir)

                if query.pop("output", "json") == "tsv":
                    self.write_tsv(json_items)
                    self.set_status(200)
                    return

                self.write({"items": json_items})
                self.set_status(200)
                return

                raise tornado.web.HTTPError(404)
        except ConnectionFailure as cfe:
            raise tornado.web.HTTPError(500, str(cfe))

    def list_datastores(self):
        if options.verbose: logging.info("list_datastores [%s]" % self.request.uri)

        items = []
        for datastore_id in self._datastore_map.keys():
            items.append({ "id": datastore_id, "uri": self.request.path + "/" + datastore_id })
        self.write({"items": items, "data_type": "datastores" })

    def list_databases(self, datastore_id):
        if options.verbose: logging.info("list_databases [%s] [%s]" % (self.request.uri, datastore_id))

        mongo_uri = self._datastore_map[datastore_id].uri
        if options.verbose: logging.info("list_databases [%s] [%s] [%s]" % (self.request.uri, datastore_id, mongo_uri))

        mongoClient = MongoClient(mongo_uri)
        items = []
        for database_name in mongoClient.database_names():
            items.append({ "id": database_name, "uri": self.request.path + "/" + database_name })
        self.write({"items": items, "data_type": "databases" })

    def list_collections(self, datastore_id, database_id):
        if options.verbose: logging.info("list_collections [%s] [%s] [%s]" % (self.request.uri, datastore_id, database_id))

        mongo_uri = self._datastore_map[datastore_id].uri
        if options.verbose: logging.info("list_collections [%s] [%s] [%s] [%s]" % (self.request.uri, datastore_id, database_id, mongo_uri))

        mongoClient = MongoClient(mongo_uri)
        database = mongoClient[database_id]

        items = []
        for collection_name in database.collection_names(False):
            items.append({ "id": collection_name, "uri": self.request.path + "/" + collection_name })
        self.write({"items": items, "data_type": "collections" })

    def open_collection(self, datastore_id, db_name, collection_id):
        if options.verbose: logging.info("open_collection [%s] [%s] [%s]" % (datastore_id, db_name, collection_id))

        mongo_uri = self._datastore_map[datastore_id].uri
        mongoClient = MongoClient(mongo_uri)
        database = mongoClient[db_name]
        return database[collection_id]

    def list_fields(self, collection):
        if options.verbose: logging.info("list_fields [%s]" % (collection.name))
        self.write({"items": collection.find_one().keys(), "data_type": "fields" })

    def query_collection(self, collection, query, sort_fld=None, sort_dir=1):
        if options.verbose: logging.info("query_collection [%s] [%s] [%s] [%s]" % (collection.name, query, sort_fld, sort_dir))
        json_items = []
        query_limit = options.mongo_rows_limit

        if sort_fld:
            for idx, item in enumerate(collection.find(query).sort(sort_fld, sort_dir)):
                if idx > query_limit: break
                json_items.append(self.jsonable_item(item))

        else:
            for idx, item in enumerate(collection.find(query)):
                if idx > query_limit: break
                json_items.append(self.jsonable_item(item))

        return json_items

    def get_datatypes(self, datasource_id, db_name, collection_id):
        c_dtypes = {}
        if options.verbose: logging.info("get_datatypes(%s, %s, %s)" % (datasource_id, db_name, collection_id))

        if not self.datastores_config is None and "datastores" in self.datastores_config:
            c_datastores = self.datastores_config["datastores"]
            if not c_datastores is None and datasource_id in c_datastores:
                if db_name in c_datastores[datasource_id]:
                    c_db = c_datastores[datasource_id][db_name]
                    if not c_db is None and "datatypes" in c_db and collection_id in c_db["datatypes"]:
                        c_dtypes = c_db["datatypes"][collection_id]

        if options.verbose: logging.info("get_datatypes(%s, %s, %s): %s" % (datasource_id, db_name, collection_id, str(c_dtypes)))
        return c_dtypes


    def transpose_query_arguments(self, db_name, datasource, datatypes={}):
        # by default, queries are case-insensitive
        normalize_fn = lambda x: re.compile("^" + x + "$", re.IGNORECASE)

        if datasource.is_case_sensitive_database(db_name):
            normalize_fn = lambda x: x

        query = {}
        args = self.request.arguments
        for key in args.keys():
            if not key in RESERVED_KEYS:
                if len(args[key]) == 1:
                    if key in datatypes:
                        if datatypes[key] == "int":
                            query[key] = int(args[key][0])
                        elif datatypes[key] == "float":
                            query[key] = float(args[key][0])
                        else:
                            query[key] = normalize_fn(args[key][0])
                    else:
                        query[key] = normalize_fn(args[key][0])
                else:
                    if key in datatypes:
                        if datatypes[key] == "int":
                            query[key] = {"$in": map(int, args[key])}
                        elif datatypes[key] == "float":
                            query[key] = {"$in": map(float, args[key])}
                        else:
                            query[key] = {"$in": map(normalize_fn, args[key])}
                    else:
                        query[key] = {"$in": map(normalize_fn, args[key])}
        return query

    def jsonable_item(self, item):
        json_item = {}
        for k in item.iterkeys():
            if k == "_id":
                json_item["id"] = str(item["_id"])
            elif "[]" in k:
                json_item[k.replace("[]", "")] = item[k]
            else:
                json_item[str(k)] = item[k]
        return json_item

    def write_tsv(self, items):
        filename = self.get_argument("output_filename", "data_export.tsv")
        attachment = "attachment; filename=\"%s\"" % filename

        if options.verbose: logging.info("write_tsv [%s]" % attachment)

        self.set_header("Content-Type", "text/tab-separated-values")
        self.set_header("Content-Disposition", attachment)

        tsvwriter = csv.writer(self, delimiter="\t")
        excludedheaders = ["id", "values"]

        if len(items) > 0:
            colheaders = [a for a in items[0].keys() if a not in excludedheaders]
            pivotkeys = []
            if "values" in items[0]:
                values_keys = items[0]["values"].keys()
                for value_key in values_keys: pivotkeys.append(str(value_key))

            combinedkeys = []
            combinedkeys.extend(colheaders)
            if pivotkeys: combinedkeys.extend(pivotkeys)
            tsvwriter.writerow(combinedkeys)

            for item in items:
                vals = []
                for colheader in colheaders: vals.append(item[colheader])

                if "values" in item:
                    item_values = item["values"]
                    for pivotkey in pivotkeys:
                        vals.append(item_values[pivotkey])

                tsvwriter.writerow(vals)
