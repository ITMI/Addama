from tornado.options import options, logging
import tornado.web

import json

from tabix_utils import tsv_region_lookup, vcf_singleline_lookup, triotype_singleline_lookup

class TabixLookupHandler(tornado.web.RequestHandler):
    def initialize(self):
        self._config_map = self.tabix_file_map
    
    def get(self, tabix_id, chromosome, start_coordinate, end_coordinate=None):
        if end_coordinate is None:
            end_coordinate = start_coordinate
        
        if tabix_id not in self._config_map.keys():
            if options.verbose:
                logging.info("Unknown tabix loookup ID [%s]" % tabix_id)
            raise tornado.web.HTTPError(404)
        
        file_info = self._config_map[tabix_id]
        file_path = file_info['path']
        lookup_fn = None
        tabix_exe = None
        
        if file_info['type'] == 'vcf':
            lookup_fn = vcf_singleline_lookup
            tabix_exe = options.tabix_executable
        elif file_info['type'] == 'trio':
            lookup_fn = triotype_singleline_lookup
            tabix_exe = options.tabix_executable
        elif file_info['type'] == 'tsv':
            lookup_fn = tsv_region_lookup
            # For parsing the tabix output for TSV files, the header has to be included. Therefore, the "-h"
            # flag has to be included in the command line.
            tabix_exe = options.tabix_executable + " -h"
        else:
            logging.error("Unknown type for file " + file_path)
            raise tornado.web.HTTPError(404)
        
        try:
            result = lookup_fn(tabix_exe, file_path, chromosome, start_coordinate, end_coordinate)
            response = {
                "chr": result.chromosome,
                "start": result.start,
                "end": result.end,
                "values": result.values,
                "snpid": result.snpid,
                "info": result.info
            }
            
            self.write(json.dumps(response, sort_keys=True))
            self.set_status(200)
        except Exception as e:
            logging.error("Running tabix failed:")
            logging.error(e)
            raise tornado.web.HTTPError(404)






