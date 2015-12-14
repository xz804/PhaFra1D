# -*- coding: utf-8 -*-
"""
Data Structure

@author: Xue Zhang
"""
import numpy as np

class DataStructFem:

    Nmat    = 1

    UBound  = []    # Boundary condition for displacement
    FBound  = []    # Boundary condition for forces

    def __init__(self,EleType,Node_xy,Elem):

        if (EleType=='Tri6'):
            self.Npoint = 6
            self.Ndim   = 2
            self.Ngauss = 3
            self.Nsigma = 4
            self.Neps   = 4

        self.EleType = EleType
        self.Node_xy = Node_xy
        self.Elem    = Elem

        self.Nnode    = int(Node_xy.shape[0])
        self.Nelem    = int(Elem.shape[0])

        self.U_tot    = np.zeros((self.Ndim*self.Nnode,),dtype=np.float)
        self.Phase    = np.zeros((1*self.Nnode,),dtype=np.float)

        self.Sigma    = np.zeros((self.Nelem,self.Ngauss*self.Nsigma),dtype=np.float)
        self.Eps      = np.zeros((self.Nelem,self.Ngauss*self.Neps),dtype=np.float)

        self.Hmax    = np.zeros((self.Nelem,self.Ngauss),dtype=np.float)





class DataStructMaterial:

    PhaseEta = 1e-20

    def __init__(self,E,Gamma,Density,GcI,GcII,Lc,PhaseEta,ProbType='PlaneStrain',PhaseType='PhaseIso'):

        self.E        = E
        self.Gamma    = Gamma
        self.Density  = Density

        self.GcI      = GcI
        self.GcII     = GcII
        self.Lc       = Lc
        self.PhaseEta = PhaseEta

        self.ProbType = ProbType
        self.PhaseType = PhaseType     #'PhaseIso', 'PhaseAmor', 'PhaseMiehe'


class DataStructMacro:

    TotalAnalysis   = True


    def __init__(self):
        pass
