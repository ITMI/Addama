from tornado.options import options, logging
import tornado.web

import json

from tabix_utils import tsv_region_lookup, vcf_singleline_lookup, triotype_singleline_lookup

from tabix_utils import CoordinateRangeEmptyError, WrongLineFoundError, UnexpectedTabixOutputError

class TabixLookupHandler(tornado.web.RequestHandler):
    def initialize(self):
        self._config_map = self.tabix_file_map

    def _handle_request_exception(self, error):
        logging.error("Uncaught exception: " + str(error))

    def write_error(self, status_code, **kwargs):
        if status_code == 400:
            self.write("Bad request")
        else:
            self.write("Server error occurred")

    def build_response_object(self, chromosome, start, end, values=None, snpid="", ref="", alt="", info=""):
        if values == None:
            values = []

        return {
            "chr": chromosome,
            "start": start,
            "end": end,
            "values": values,
            "snpid": snpid,
            "ref": ref,
            "alt": alt,
            "info": info
        }

    def get(self, tabix_id, chromosome, start_coordinate, end_coordinate=None):
        if end_coordinate is None:
            end_coordinate = start_coordinate
        
        if tabix_id not in self._config_map.keys():
            logging.error("Unknown tabix lookup ID [%s]" % tabix_id)
            self.send_error(400)
            return
        
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
            logging.error("Unknown type for file in configuration: " + file_path)
            self.send_error(500)
            return
        
        try:
            result = lookup_fn(tabix_exe, file_path, chromosome, start_coordinate, end_coordinate)
            response = self.build_response_object(result.chromosome,
                                                  result.start,
                                                  result.end,
                                                  values=result.values,
                                                  snpid=result.snpid,
                                                  ref=result.ref,
                                                  alt=result.alt,
                                                  info=result.info)
            
            self.write(json.dumps(response, sort_keys=True))
            self.set_status(200)

        except CoordinateRangeEmptyError as cnf:
            logging.info(cnf)
            response = self.build_response_object(chromosome, start_coordinate, end_coordinate)
            self.write(json.dumps(response, sort_keys=True))
            self.set_status(200)

        except WrongLineFoundError as wlf:
            logging.error(wlf)
            self.send_error(500)

        except UnexpectedTabixOutputError as eto:
            logging.error(eto)
            self.send_error(500)

        except Exception as e:
            logging.error("Running tabix failed: " + str(e))
            self.send_error(500)
