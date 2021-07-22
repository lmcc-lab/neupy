'''
This script requires the following dependencies and functions.
'''
import sympy
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import pandas as pd
import ast

def nMax(n, array):
    '''
    Get the nth maximum value from array
    
    @param - n, int. The nth value from array
    @param - array, array. The array in question

    returns nth maximum value

    '''
    for i in range(n-1):
        del array[array.index(max(array))]
    return max(array)

def split(word): 
    '''
    Split a word into characters

    @param - word, str.

    returns list of characters
    '''
    return [char for char in word]

def remove_elem(List, element):
    '''
    Removes element from list

    @param - List, list
    @param - element, any type

    returns list without element
    '''
    try:
        if element in List:
            List.remove(element)
    except TypeError:
        pass
    return List

def conv_str_to_list(source):
    '''
    Convert a string representation of a list into a list.
    First tries a literal evaluation, if that doesn't work then
    it does it manually.

    @param - source, str
    
    returns list as translated from string.
    '''
    try:
        source = ast.literal_eval(source)
    except:
        source = [[i.replace(']','')] for i in source.split('[[',1)[1].split(']]')[0].split(', [')]
        for i, sub in enumerate(source):
            sub = sub[0].split(', ')
            try:
                sub = [float(j) for j in sub]
            except ValueError:
                pass
            source[i] = sub
    return source

def counter(val, listobj):
    '''
    counts number of times value is in list

    @param - val, any type
    @param - listobj, list

    reutrns count as int
    '''
    count = 0
    for i, elem in enumerate(listobj):
        if elem == val:
            count += 1
    return count


def derivative(f,a):
    '''Compute the difference formula for f'(a) with step size h.

    Parameters
    ----------
    f : function
        Vectorized function of one variable
    a : number
        Compute derivative at x = a

    Returns
    -------
    float
        Difference formula:
            central: f(a+h) - f(a-h))/2h
            forward: f(a+h) - f(a))/h
            backward: f(a) - f(a-h))/h            
    '''
    dy = np.zeros(len(f))
    for i,_ in enumerate(f):
        if i > 0 and i < len(a)-1:
            h = a[i+1]-a[i-1]
            dy[i] = (f[i+1]-f[i-1])/h
        elif i==0:
            h = a[i+1]-a[i]
            dy[i] = (f[i+1]-f[i])/h
        else:
            h = a[i-1]-a[i]
            dy[i] = (f[i-1]-f[i])/h
    return dy

def find_index(array, value):
    '''
    Find the index of a value in an array

    parameters 
    ----------
    array : array
    value : type of value that is 
            consistent with elements of string

    Returns
    -------
    index

    '''
    array = np.array(array)
    idx = (np.abs(array-value)).argmin()
    return idx

def integral(x, y):
    '''
    Numerically integrate
    '''
    yd = np.zeros(len(x))
    for i, xval in enumerate(x):
        if i+1<len(x):
            miny = min(y[i],y[i+1])
            maxy = max(y[i],y[i+1])
            base = x[i+1]-xval
            yd[i] = base*miny+0.5*(maxy-miny)*base
        else:
            miny = min(y[i],y[i-1])
            maxy = max(y[i],y[i-1])
            base = x[i-1]-xval
            yd[i] = base*miny+0.5*(maxy-miny)*base
    return yd
