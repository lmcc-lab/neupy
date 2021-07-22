'''
Naming module for converting NU isotope format to readable and seperable formats.

Author: Liam McClelland

Last Edited: 31/05/21

'''


import pandas as pd
import dependencies

nubase = pd.read_excel('./databases/Parsed_nubase2016.xlsx')
def readable(iso):
    '''
    @param iso, int in NU format

    returns readable isotope, eg 2H-1
    '''
    A = nubase['A'].tolist()
    Aindex = A.index((iso//10000)%1000)
    Z = nubase['Z'].loc[Aindex:].tolist()
    Zindex = Z.index(iso//10000000)
    return str(nubase['id'].loc[Aindex+Zindex])+'-'+str(int(iso//10000000))

def seperate(iso):
    '''
    @param iso, int in NU format

    returns list [Z, A, Level]
    '''
    return [iso//10000000, (iso//10000)%1000, iso%10000]

def short_readable(iso):
    '''
    @param iso, int in NU format

    returns shortened readable isotope, eg 2H
    '''
    A = nubase['A'].tolist()
    Aindex = A.index((iso//10000)%1000)
    Z = nubase['Z'].loc[Aindex:].tolist()
    Zindex = Z.index(iso//10000000)
    return str(nubase['id'].loc[Aindex+Zindex])

def NU(Z,A,level):
    '''
    @param Z, int - proton number
    @param A, int - Atomic number
    @param level, int - energy level

    returns NU, formatted in one int.
    '''
    return int(Z*10000000+A*10000+level)

def readableChain(decayChain):
    '''
    @param decayChain, list with embedded daughter lists for each generation. In NU format

    '''
    read = []
    for i,gen in enumerate(decayChain):
        if type(gen) == list:
            read.append([])
            for _,daughter in enumerate(gen):
                readableDaughter = readable(daughter)
                read[i].append(readableDaughter)
        else:
            read.append(readable(gen))
    return read

def readableToSeperate(readableName):
    '''
    Convert a normal readable element to a seperated list of [Z, A, level]
    for example, 135I-53 converts to [53, 135, 0]
    
    parameters
    ------------------------
    readableName, str - readable name in format Z+ID+"-"+A
    '''
    elements = dependencies.split(readableName)
    for i, val in enumerate(elements):
        try:
            int(val)
        except ValueError:
            Z_index = i
            break
    Z = int(''.join(elements[:Z_index]))
    A = int(readableName.split('-',1)[1])
    return [A, Z, 0]

