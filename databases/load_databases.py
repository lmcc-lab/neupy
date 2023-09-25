import pandas as pd
from yaml import safe_load
from tqdm import tqdm
from pprint import pprint
import math
from typing import Tuple
import logging
import os
import io
import numpy as np

"""
Transform Parsed_nubase2016.xlsx file into a SQL database. Data will be
normalised into SI units, 
"""

def load_nubase(path: str='databases/', filename: str='nubase2020.txt', log_level=logging.WARNING) -> pd.DataFrame:
    """
    We load the nubase2020.txt file, and we use the information in the config.yaml file
    to get the columns, and convert the data into the necessary format (float or string).
    The config file is assumed to be in the same path as the database itself.
    """
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT)
    log = logging.getLogger('load_nubase')
    log.setLevel(log_level)
    with open(f"{path}{filename}", 'r') as f:
        # read lines of nubase2020
        nubase_lines = f.readlines()
    
    # load the config.yaml file
    config = load_config(path)
    # get the column descriptions for nubase2020
    column_desc = config[filename.replace('.txt','')]['columns']

    # prepare a dictionary with the appropriate headers(keys)
    entry = {meta_data['quantity']: [] for meta_data in column_desc.values()}
    
    # line by line, ignoring any commented lines starting with #. 
    for line in nubase_lines:
        if line[0] == '#':
            continue
        # go through the column range and other information stored in the config file
        for column_range, meta_data in column_desc.items():
            # convert start and end column numbers to ints
            first, second = [int(num) for num in column_range.split(':')]
            
            heading = meta_data['quantity']
            # get the value of the column, removing any spaces and new line markers
            value = line[first-1: second].replace(' ', '') if heading != 'BR' else line[first-1: second].replace('\n','')
            # get format and quantity information for this column
            form = meta_data['format']
            # convert value to float if value is not ''.
            systematics_flag = False
            equality_flag = False
            equality = ''
            if ('f' in form and value != '') and not value.isspace():

                if '#' in value:
                    value = value.replace('#', '')
                    systematics_flag = True

                if value[0] in ['>', '<', '~']:
                    equality_flag = True
                    equality = value[0]
                    value = value[1:]
                try:
                    value = float(value)
                except ValueError:
                    log.warning(f'WARNING: Value {value} was attempted to be converted to a float. Please check if this is incorrect')
            if systematics_flag:
                value = (value, '#')
            if equality_flag:
                value = (value, equality)
            entry[heading].append(value)
    # Convert entries to df
    df = pd.DataFrame(entry)
    df["AZI"] = df["AAA"].astype(str) + df["ZZZi"].astype(str)
    df["A"] = df["AAA"].astype(int)
    df["Z"] = df["ZZZi"].apply(lambda ZZZi: ZZZi[:3]).astype(int)
    df["i"] = df["ZZZi"].apply(lambda ZZZi: ZZZi[-1]).astype(int)
    df.index = df['AZI']

    return df, config


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


def load_fyu235thermal(path: str='databases/') -> pd.DataFrame:
    """
    Load fyu235thermal.txt file, and create a new column in AZI format.

    @params
    path: str, path where fyu235thermal.txt file is contained.

    @returns
    fiss_data: pd.DataFrame
    """
    fiss_data = pd.read_csv(path + 'fyu235thermal.txt', sep="	", header=0)
    fiss_data.columns = ["Z", "A", "Level", "YI", "YI uncert"]
    fiss_data['AZI'] = fiss_data['A'].apply(lambda A: ''.join(["0"]*(3-len(f"{A}")))+f"{A}") + fiss_data['Z'].apply(lambda Z: ''.join(["0"]*(3-len(f"{Z}"))) + f"{Z}") + fiss_data['Level'].apply(lambda level: f"{level}")
    return fiss_data


def load_all_fy_databases(path: str='databases/') -> dict:
    """
    Load any txt file of the form fy{}.txt file, add them to a dictionary with key {}
    and add a column of AZI, assuming all fy files are in the same format (sep = "   ", header=0 with
    columns "Z", "A", "Level", "YI", "YI uncert")

    @params
    path: str, path where fy{}.txt files is contained.

    @returns
    fiss_data: Dict[str: pd.DataFrame]
    """
    files = os.listdir(path)
    db = {}
    for f in files:
        if '.' not in f:
            continue
        f_split = f.split('.', 1)
        if f_split[1] != 'txt':
            continue
        if f[:2] != 'fy':
            continue
        key = f_split[0][2:]
        fiss_data = break_ENSDF_db(path+f, save_file=False)
        db[key] = fiss_data
    return db


def unique_decay_modes(nubase: pd.DataFrame):
    BR = nubase['BR']
    unique_dm = set()
    split_symbols = ['=','>','<', '~', '?']
    for row in BR:
        dm = row.split(';')
        for d in dm:
            for ss in split_symbols:
                if ss not in d:
                    continue
                d = d.split(ss)[0]
                if d == 'B':
                    print(row)
            unique_dm.add(d)

def break_ENSDF_db(filename: str, save_file=True):
    with open(filename, 'r') as f:
        lines = f.readlines()

    neutron_energies_index = [i for i, line in enumerate(lines) if 'Neutron energy' in line]
    if len(neutron_energies_index) == 0:
        lines.remove('----------- ------- ------- ----------- -----------\n')
        df = pd.read_csv(io.StringIO(''.join(lines)), delimiter=r"\s{2,}", engine='python')
        df = df.rename(columns={'FPS': 'Level'})
        df['Level'] = df['Level'].astype(int)
        df['ZAFP'] = df['ZAFP'].astype(int)
        df['A'] = (df['ZAFP'] % 1000).astype(int)
        df['Z'] = ((df['ZAFP'] - df['ZAFP']%1000)/1000).astype(int)
        df['AZI'] = df['A'].apply(lambda A: ''.join(["0"]*(3-len(f"{A}")))+f"{A}") + df['Z'].apply(lambda Z: ''.join(["0"]*(3-len(f"{Z}"))) + f"{Z}") + df['Level'].apply(lambda level: f"{level}")
        df = df.set_index('AZI', drop=True)
        df = df.drop(columns=['ZAFP', 'PRODUCT'], axis=1)
        tables = {'independent': df}
        return tables
    
    tables = {float(lines[row_index].split(' = ', 1)[1].replace('\n', "").replace(' ', '')):
              lines[row_index+2: neutron_energies_index[i+1]-1]
              if i+1 < len(neutron_energies_index) else lines[row_index+2:]
              for i, row_index in enumerate(neutron_energies_index)}
    
    for ne, table in tables.items():
        table.remove('----------- ------- ------- ----------- -----------\n')
        if save_file:
            with open(filename.replace('.txt', '')+f'_ne={ne}.txt', 'w+') as f:
                f.writelines(table)
        df = pd.read_csv(io.StringIO(''.join(table)), delimiter=r"\s{2,}", engine='python')
        df = df.rename(columns={'FPS': 'Level'})
        df['Level'] = df['Level'].astype(int)
        df['ZAFP'] = df['ZAFP'].astype(int)
        df['A'] = (df['ZAFP'] % 1000).astype(int)
        df['Z'] = ((df['ZAFP'] - df['ZAFP']%1000)/1000).astype(int)
        df['AZI'] = df['A'].apply(lambda A: ''.join(["0"]*(3-len(f"{A}")))+f"{A}") + df['Z'].apply(lambda Z: ''.join(["0"]*(3-len(f"{Z}"))) + f"{Z}") + df['Level'].apply(lambda level: f"{level}")
        df = df.set_index('AZI', drop=True)
        df = df.drop(columns=['ZAFP', 'PRODUCT'], axis=1)
        tables[ne] = df
    return tables




if __name__ == "__main__":
    tables = break_ENSDF_db('databases/fyU235.txt', save_file=False)