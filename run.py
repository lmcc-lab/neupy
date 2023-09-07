from databases.load_databases import load_nubase, load_fyu235thermal, logging
from pprint import pprint

nubase2020, config_file =load_nubase(log_level=logging.ERROR)
fyu235 = load_fyu235thermal()

print(nubase2020)