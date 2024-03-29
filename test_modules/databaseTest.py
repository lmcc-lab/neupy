from databases.load_databases import *


def nubase_selector_test():
    nubase = load_nubase('./databases/')
    print(rename_and_simplify_nubase_columns(nubase))


def make_db_test():
    nubase, fission, config = load_data()
    nubase = rename_and_simplify_nubase_columns(nubase)
    make_database(nubase, fission, debug=True)


def load_config_test():
    pprint(load_config('databases/'))