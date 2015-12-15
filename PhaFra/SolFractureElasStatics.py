# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 09:00:45 2015

@author: xz804
"""

import numpy as np
import SolElasStatics as SolFEM
import SolPhaseField  as SolPha
#from DataStruct import *


def SolFractureElasStatics(FemInf,MatInf,MacroInf):
    """
    read Abaqus-generated inp file

    Call ReadNode
    Call ReadElement

    Node_info       node information
    Element_info    element information
    """
    
    #-----  solve Elastic Static Problem in 2D
    FemInf = SolFEM.SolElasStatics(FemInf,MatInf,MacroInf)
    
    #------ solve Phase Field
    FemInf = SolPha.SolPhaseField(FemInf,MatInf,MacroInf)

    return FemInf
    
    
    
    
