from tornado.options import logging
import tornado.web

import json

import variant_summary_lookup as vsl

REQUIRED_ARGUMENTS = frozenset(['chromosome', 'coordinate', 'feature_id'])

class VariantSummaryHandler(tornado.web.RequestHandler):
    def initialize(self):
        self._config_map = self.data_map

    def _handle_request_exception(self, error):
        logging.error("Uncaught exception: " + str(error))

    def write_error(self, status_code, **kwargs):
        if status_code == 400:
            self.write("Bad request")
        else:
            self.write("Server error occurred")

    def get(self, *uri_path):
        sub_path = self.request.path.replace("/variant_summary", "")
        uri_parts = sub_path.split("/")
        data_id = uri_parts[1]

        if data_id not in self._config_map.keys():
            logging.error("Unknown Variant Summary data lookup ID [%s]" % data_id)
            self.send_error(400)
            return

        config = self._config_map[data_id]

        if not set(self.request.arguments).issubset(REQUIRED_ARGUMENTS):
            logging.error("Variant Summary - invalid arguments: [%s]" % str(self.request.arguments))
            self.send_error(400)
            return

        chromosome = int(self.get_argument("chromosome"))
        coordinate = int(self.get_argument("coordinate"))
        feature_id = self.get_argument("feature_id")

        logging.debug("Querying Variant Summary for \'" + str((chromosome, coordinate, feature_id)) + "\'")

        try:
            result = vsl.do_query(config, chromosome, coordinate, feature_id)
            self.write(json.dumps(result, sort_keys=True))
            self.set_status(200)

        except Exception as e:
            logging.error("Running Variant Summary service failed: " + str(e))
            self.send_error(500)
