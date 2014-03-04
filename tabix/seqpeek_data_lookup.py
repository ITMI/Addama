from tornado.options import options, logging
import tornado.web

import json

import seqpeek_data_service as sds

sds.local_pprint = logging.debug

class SeqPeekDataHandler(tornado.web.RequestHandler):
    def initialize(self):
        sds.mDEBUG = options.verbose

        self._config_map = self.seqpeek_data_map

    def get(self, *uri_path):
        sub_path = self.request.path.replace("/seqpeek_data", "")
        uri_parts = sub_path.split("/")
        data_id = uri_parts[1]

        if data_id not in self._config_map.keys():
            if options.verbose:
                logging.info("Unknown SeqPeek data loookup ID [%s]" % data_id)
            raise tornado.web.HTTPError(404)

        config = self._config_map[data_id]

        if "gene" not in self.request.arguments:
            raise tornado.web.HTTPError(404)

        gene_label = self.get_argument("gene")
        logging.debug("Querying SeqPeek data for gene \'" + gene_label + "\'")

        seqObj = sds.SeqPeekResult(config)
        
        seqObj.full_gene = False
        seqObj.gene_name = gene_label

        result = sds.do_gene_query(seqObj)
        self.write(json.dumps(result, sort_keys=True))
        
        self.set_status(200)
