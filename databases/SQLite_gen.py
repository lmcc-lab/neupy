import sqlite3
import pandas as pd
from yaml import safe_load
from tqdm import tqdm
from pprint import pprint
import math
from typing import Tuple

"""
Transform Parsed_nubase2016.xlsx file into a SQL database. Data will be
normalised into SI units, 
"""


def load_nubase(path: str) -> pd.DataFrame:
    """
    Load nubase file

    @params
    path: str, path where Parsed_nubase2016.xlsx file is contained.
    """
    nubase = pd.read_excel(path + 'Parsed_nubase2016.xlsx', header = 0)
    return nubase


def load_config(config_path: str) -> dict:
    """
    Load config.yaml file

    @params
    path: str, path where config.yaml file is contained.

    @returns
    fiss_data: dict
    """
    with open(config_path + 'config.yaml', 'r') as file:
        config = safe_load(file)
    return config


def load_fyu235thermal(path: str) -> pd.DataFrame:
    """
    Load fyu235thermal.txt file, and sorted by atomic number "A".

    @params
    path: str, path where fyu235thermal.txt file is contained.

    @returns
    fiss_data: pd.DataFrame
    """
    fiss_data = pd.read_csv(path + 'fyu235thermal.txt', sep="	", header=0)
    fiss_data.columns = ["Z", "A", "Level", "YI", "YI uncert"]
    fiss_data = fiss_data.sort_values('A')
    return fiss_data

def isNaN(num):
    """
    Checks if any input is nan. Checking val against itself will return True if it isn't a nan, and visa versa.
    
    @return
    None if input is NaN
    input otherwise
    """
    if num!=num:
        return None
    return num


def convert_half_life_to_SI(half_life, unit, uncertainty, config):
    try:
        # First check that 
        half_life = float(half_life)
        uncertainty = float(uncertainty)
        prefix = None
        if len(unit) > 1:
            prefix, unit = unit
        conversion = config['prefix'].get(prefix, 1) * config['unit'].get(unit, 1)
        half_life = half_life * conversion
        uncertainty = uncertainty * conversion
        unit = 's'
        
    except ValueError:
        pass
    return half_life, unit, uncertainty


def make_float(variable):
    try:
        return float(variable)
    except ValueError:
        return variable


def load_data(root_path: str='./databases/') -> Tuple[pd.DataFrame, pd.DataFrame, dict]:
    r"""
    Load the 3 files for creating our SQL database: Parsed_nubase2016.xlsx, fyu235thermal.txt and config.yaml

    @params:
    root_path: str, root path where all files are kept. Default is './databases/

    @returns
    nubase: pd.DataFrame, nubase dataframe.
    fission: pd.DataFrame, fission dataframe.
    config: dict, config file

    """
    # load files. This will fail if root_path is incorrect
    print("loading nubase")
    nubase = load_nubase(root_path)
    print("loading done")
    print("loading config")
    config = load_config(root_path)
    print("loading done")
    print("loading fission data")
    fission = load_fyu235thermal(root_path)
    print("loading done")

    return nubase, fission, config


def show_unique_types(column: list, debug: bool=False) -> set:
    """
    Given a column of values, get the set of the unique values. If there are floats then only output
    float as one of the values. This can be used to identify the different datatypes/symbols used in
    a pandas dataframe

    @params
    column: list, The column of values. For example, a dataframe with column names A, B, C ..., then column is
            df.A.values
    debug: bool, default False. If True it prints out the unique values in the function itself.
    
    @returns
    set of unique values, and float class name if there are floats.
    """
    unique_values = set([float if isinstance(make_float(val), float) else val for val in column])
    if debug:
        print(f"Unique nubase entries in column")
        pprint(unique_values)
    return unique_values


def rename_and_simplify_nubase_columns(nubase: pd.DataFrame) -> pd.DataFrame:
    """
    This function just changes the column names to make them clearer. 

    The new column names are at times abreviated, but this abbreviation is
    explicitly written in one of the columns. For example, half_life_hl has the
    abreviation at the end of the column name `hl'. `hl' is then used in other columns
    such as hl_uncert and hl_unit so it can be easily read as half_life_uncert and
    half_life_unit repectively.

    All other columns are deleted. 

    @params
    nubase: pd.DataFrame, nubase database

    @returns:
    pd.DataFrame, simplified with clearer headers.
    """
    # Go row by row through nubase database. This is for readablity
    rename_columns = {'A': 'A',
                      'Z': 'Z',
                      'element': 'element',
                      'id': 'id',
                      'Excess energy': 'excess_energy_en', 
                      'uncertainity': 'en_uncert', 
                      'extrapolated': 'en_extrapolated',
                      'Half life': 'half_life_hl', 
                      'Half life uncert': 'hl_uncert', 
                      'unit': 'hl_unit', 
                      'relation': 'hl_relation', 
                      'extrapolated2': 'hl_extrapolated',
                      'J pi': 'j_pi', 
                      'Decay mode': 'decay_mode_dm', 
                      'relation7': 'dm_relation',
                      'extrapolated5': 'dm_extrapolated', 
                      'Excitation energy': 'isomer_excitation_energy_iee', 
                      'Excitation uncert': 'iee_uncert', 
                      'extrapolated10': 'iee_extrapolated', 
                      'Half life 2': 'isomer_half_life_ihl', 
                      'Half life uncert 2': 'ihl_uncert',
                      'unit 2': 'ihl_unit',
                      'relation13': 'ihl_relation',
                      'extrapolated12': 'ihl_extrapolated', 
                      'Decay mode 2': 'isomer_decay_mode_idm',
                      'Branch ratio 2': 'idm_branch_ratio', 
                      'Branch ratio uncert 2': 'idm_branch_ratio_uncert',
                      'relation20': 'idm_relation', 
                      'extrapolated17': 'idm_extrapolated', 
                      'J pi 2': 'isomer_j_pi', 
                      'T': 'isomer_T'}

    #columns that aren't in rename_columns
    unused_columns = [column for column in nubase.columns if column not in rename_columns.keys()]
    
    # rename columns of nubase df
    nubase.rename(columns=rename_columns, inplace=True)

    # Create a decay mode mask to check if value is isotope abundance or not (True if it is)
    decay_mode_check_mask = nubase['decay_mode_dm'] != 'is'
    # Creat it's anti-mask
    anti_mask = nubase['decay_mode_dm'] == 'is'
    # Use mask to define branch ratio column
    nubase['dm_branch_ratio'] = nubase['Branch ratio/Isoptope abund'].where(decay_mode_check_mask, None)
    # and it's uncertainty
    nubase['dm_branch_ratio_uncert'] = nubase['Branch ratio/Isoptope abund uncert'].where(decay_mode_check_mask, None)
    # Likewise, use anti-mask to define the isotope abundance column
    nubase['isotope_abundance_ia'] = nubase['Branch ratio/Isoptope abund'].where(anti_mask, None)
    # and it's uncertainty
    nubase['ia_uncert'] = nubase['Branch ratio/Isoptope abund uncert'].where(anti_mask, None)

    # remove unused/unchaged columns
    nubase.drop(unused_columns, axis=1, inplace=True)

    return nubase


def make_new_table(df: pd.DataFrame, table_name:str, root_path: str = './databases/'):
    # start connection to sql database, or create it if it doesn't exist
    con = sqlite3.connect(root_path + '/nuclide_data.db')
    # start a cursor
    cur = con.cursor()

    try:
        # columns of table
        columns = ', '.join([val.replace(' ', '_') for val in list(df.columns)])
        # begin table, giving it columns
        cur.execute("""CREATE TABLE IF NOT EXISTS {0}({1})""".format(table_name, columns))

        # make data ready for inserting into table
        data = [tuple(row) for row in df.values]

        # make placeholders
        placeholders = ', '.join(['?']*len(data[0]))

        # insert data into nuclide
        cur.executemany(f"INSERT INTO {table_name} VALUES({placeholders})", data)

        # commit changes
        con.commit()
    
        # close cursor and connection
        cur.close()
        con.close()
    
    except sqlite3.OperationalError as e:
        cur.close()
        con.close()
        raise(e)


def database_debugging(root_path):
    con = sqlite3.connect(root_path + '/nuclide_data.db')
    # start a cursor
    cur = con.cursor()
    res = cur.execute("SELECT name FROM sqlite_master")
    tables = res.fetchall()
    print("Tables: ", tables)
    print("")
    for table in tables:
        if isinstance(table, tuple):
            table = table[0]
        print(table)
        cur.execute(f"SELECT * FROM {table}")
        table_cols = list(map(lambda x: x[0], cur.description))
        print(f"{table} column names:")
        pprint(table_cols)
        print("")
    
        res = cur.execute(f"SELECT {', '.join(table_cols)} FROM {table} ORDER BY A")
        first_entry = res.fetchmany(10)
        print("First 10 entries")
        for entry in first_entry:
            print(entry)


def make_database(nubase: pd.DataFrame, fission: pd.DataFrame, debug=False, **kwargs) -> None:
    """
    Take the nubase and fission dataFrames and save them as tables of a SQLite database. This will allow for faster
    loading in the future.

    @params
    nubase: pd.DataFrame, the nubase data
    fission: pd.DataFrame, the fission yield data
    debug: bool, prints debugging information.
    
    kwargs: passed onto make_new_table, or used if debug is True
            ---------------------------
            root_path: str, the path to the databases
    """
    # Make a new table for each of the databases
    make_new_table(nubase, 'nubase', **kwargs)
    make_new_table(fission, 'fission', **kwargs)


    # debugging output, print table names, table column names, and first 10 rows.
    if debug:
        database_debugging(kwargs.get('root_path', './databases/'))
        

if __name__ == "__main__":
    make_database(print_set_nubase_col='Branch ratio/Isoptope abund')