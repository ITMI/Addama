from tornado.options import options, logging, define
import tornado.web

import json
import sys

from tabix_lookup import vcf_lookup, triotype_lookup

server_settings = {
    "xheaders" : True,
    "address" : "0.0.0.0"
}

settings = {
    "debug": True
}

DATA_MAP = None

class TabixLookup(tornado.web.RequestHandler):
    def get(self, vcf_id, chromosome, coordinate):
        global VCF_MAP
        
        if vcf_id not in DATA_MAP:
            logging.info("Unknown VCF id \'" + vcf_id + "\'")
            raise tornado.web.HTTPError(404)
        
        file_info = DATA_MAP[vcf_id]
        file_path = file_info['path']
        lookup_fn = None
        if file_info['type'] == 'vcf':
            lookup_fn = vcf_lookup
        elif file_info['type'] == 'trio':
            lookup_fn = triotype_lookup
        else:
            logging.error("Unknown type for file " + file_path)
            raise tornado.web.HTTPError(404)
        
        try:
            result = lookup_fn(options.tabix_executable, file_path, chromosome, coordinate, coordinate)
            response = {
                "chr": result.chromosome,
                "coordinate": result.coordinate,
                "values": result.values
            }
            
            self.write(json.dumps(response, sort_keys=True))
            self.set_status(200)
        except Exception as e:
            logging.error("Running tabix failed:")
            logging.error(e)
            raise tornado.web.HTTPError(404)


def parse_vcf_config(config_path):
    json_file = open(config_path, 'rb')
    data = json.load(json_file)
    json_file.close()
    return data

def main():
    define("port", default=8321, help="Server port", type=int)
    define("verbose", default=False, help="Prints debugging statements")
    define("tabix_executable", default=".", help="Path to tabix executable")
    define("vcf_config", default=".", help="Path to VCF configuration JSON")

    tornado.options.parse_command_line()

    logging.info("Starting Tornado web server on http://localhost:%s" % options.port)
    logging.info("--tabix_executable=%s" % options.tabix_executable)
    logging.info("--vcf_config=%s" % options.vcf_config)
    logging.info("--verbose=%s" % options.verbose)

    # Try loading the VCF mapping
    logging.info("Loading VCF configuration file \'" + options.vcf_config + "\'...")
    try:
        global DATA_MAP
        DATA_MAP = parse_vcf_config(options.vcf_config)
    except Exception as e:
        logging.error("Failed to load VCF configuration file:")
        logging.error(e)
        sys.exit(1)

    application = tornado.web.Application([
        (r"/(\w+)/(X|Y|M|\d{1,2})/(\d+)", TabixLookup)
    ], **settings)
    application.listen(options.port, **server_settings)
    tornado.ioloop.IOLoop.instance().start()


if __name__ == "__main__":
    main()






