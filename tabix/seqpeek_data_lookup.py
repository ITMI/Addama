from tornado.options import options, logging
import tornado.web

import json

import seqpeek_data_service as sds

sds.local_pprint = logging.debug

class SeqPeekDataHandler(tornado.web.RequestHandler):
    def initialize(self):
        sds.mDEBUG = options.verbose

        self._config_map = self.seqpeek_data_map

    def _handle_request_exception(self, error):
        logging.error("Uncaught exception: " + str(error))

    def write_error(self, status_code, **kwargs):
        if status_code == 400:
            self.write("Bad request")
        else:
            self.write("Server error occurred")

    def get(self, *uri_path):
        sub_path = self.request.path.replace("/seqpeek_data", "")
        uri_parts = sub_path.split("/")
        data_id = uri_parts[1]

        if data_id not in self._config_map.keys():
            logging.error("Unknown SeqPeek data loookup ID [%s]" % data_id)
            self.send_error(400)
            return

        config = self._config_map[data_id]

        if "gene" not in self.request.arguments:
            logging.error("Gene label missing in request arguments: [%s]" % str(self.request.arguments))
            self.send_error(400)
            return

        gene_label = self.get_argument("gene")
        logging.debug("Querying SeqPeek data for gene \'" + gene_label + "\'")

        seqObj = sds.SeqPeekResult(config)
        
        seqObj.full_gene = False
        seqObj.gene_name = gene_label

        try:
            result = sds.do_gene_query(seqObj)
            self.write(json.dumps(result, sort_keys=True))
            self.set_status(200)

        except Exception as e:
            logging.error("Running SeqPeek data service failed: " + str(e))
            self.send_error(500)