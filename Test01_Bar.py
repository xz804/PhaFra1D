# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 08:58:06 2015

@author: xz804
"""

import os
import PhaFra
import numpy as np

def main():
    print "Current working dir: %s \n" % os.getcwd()

    #------  Mesh node and Mesh element information stored in DataStructFem
    Node_xy, Elem = PhaFra.ReadInp("CrackBar.inp")

    FemInf   = PhaFra.DataStruct.DataStructFem('Tri6',Node_xy,Elem)

    #-----   Material information stored in DataStructMaterial
    E        = 210.
    Gamma    = 0.0
    Density  = 0.
    ProbType = "PlaneStrain"
    PhaseType = "PhaseIso"
    GcI      = 2.7e-3
    GcII     = GcI
    Lc       = 0.01
    PhaseEta = 0.001
    MatInf = PhaFra.DataStruct.DataStructMaterial(E,Gamma,Density,GcI,GcII,Lc,PhaseEta,ProbType,PhaseType)

    #-----   Boundary Condition information stored in DataStructFem
    UBound = []
    for inode in range(int(Node_xy.shape[0])):
        if  Node_xy[inode,1]<-0.999999:
            UBound.append([inode,"Ux",0.])
            UBound.append([inode,"Uy",0.])
        elif Node_xy[inode,1]>0.999999:
            UBound.append([inode,"Ux",0.005])   # incremental Ux
            UBound.append([inode,"Uy",0.])
    FemInf.UBound = UBound
    #-----   Macro control information stored in DataStructMacro
    MacroInf = PhaFra.DataStruct.DataStructMacro()

    ###########################################################
    print ("Initialization  ----  completed \n")
    ###########################################################

    nstep = 60

    for istep in range(nstep):

        NUBound = len(UBound)
        for i in range(NUBound):
            if FemInf.UBound[i][2] >0.00:
                FemInf.UBound[i][2] = 0.005*(istep+1)


        FemInf = PhaFra.SolFractureElasStatics(FemInf,MatInf,MacroInf)

        print ("========istep = %d completed=========", istep)



if __name__ == "__main__":
    main()
