# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 09:00:45 2015

@author: xz804
"""
import numpy as np

def ReadInp(InpName="PhaFra.inp"):
    """
    read Abaqus-generated inp file

    Call ReadNode
    Call ReadElement

    Node_info       node information
    Element_info    element information
    """



    inp_info = []
    fp = open(InpName, 'r+')
    inp_info = fp.readlines()
    
    Node_info = ReadNode(inp_info)
    Element_info = ReadElement(inp_info)

    fp.close()
    
    return Node_info, Element_info
    
def IsNumber(line):
    try:
        int(line[0])
        return True
    except ValueError:
        return False
        
def ReadNode(inp_info):
    # read the nodes information
    node_info = []
    
    for line in inp_info:
        if line.startswith("*Node"):
            index = inp_info.index(line)
            break
        
    for line in inp_info[index+1:]:
        itemID = inp_info.index(line)
        line = line.split(',')
        if IsNumber(line):
            del inp_info[itemID]
            node_info.append([float(i) for i in line])
        else:
            break
        
    node_info = np.array(node_info)
    return node_info

def ReadElement(inp_info):
    # read the elements information

    element_info = []
    for line in inp_info:
        if line.startswith("*Element"):
            index = inp_info.index(line)
            break
    for line in inp_info[index+1:]:
        itemID = inp_info.index(line)
        line = line.split(',')
        if IsNumber(line):
            del inp_info[itemID]
            element_info.append([int(i) for i in line])
        else:
            break
    element_info = np.array(element_info)
    return element_info