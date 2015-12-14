# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 09:00:45 2015

@author: xz804
"""
import numpy as np

from scipy import sparse

from scipy.sparse.linalg import dsolve

import Ele01 as EleTask


def SolElasStatics(FemInf,MatInf,MacroInf):

    #  PhaseCoef_Gauss = np.zeros((Nelem,NGP),dtype=np.float)
    PhaseCoef_Gauss = EleTask.TranPhase(FemInf.Elem,FemInf.Nelem,FemInf.Phase,MatInf.PhaseEta)

    #------- form global stiffness matrix K_g
    K_g = EleTask.GlobalStiffness(FemInf.Nnode,FemInf.Node_xy,FemInf.Nelem,FemInf.Elem,MatInf,PhaseCoef_Gauss)

    K_g = K_g.tocsr()

    print ("Revise Stiffness -----begin")

    #------- impose Essential boundary condition
    # index for the free node  -- initialization

    Nessb = len(FemInf.UBound)
    Ndof = 2*FemInf.Nnode - Nessb
    T_index = EleTask.CalDofIndex(FemInf.Nnode, Nessb,FemInf.UBound, Ndof)

    # reformulation K*U=Fext
    #---- U = T*Ur + Uf   Ur:(2*Ndof,1)    Uf: (2*Nnode,1)
    #---- K(T*Ur + Uf) = F
    #---- Transpose(T)*K*T * Ur = Transpose(T) * (F - K*Uf)

    K_gR = K_g.dot(T_index)
    K_gR = (T_index.T).dot(K_gR)

    EssBound = FemInf.UBound
    dUf = np.zeros((2*FemInf.Nnode,),dtype=np.float)
    for i in range(Nessb):
        if EssBound[i][1] == 'Ux':
            dUf[2*EssBound[i][0]]   = EssBound[i][2]
        elif EssBound[i][1] == 'Uy':
            dUf[2*EssBound[i][0]+1] = EssBound[i][2]

#    RHS = T*K_g.dot(-1*dUf)
    RHS = K_g.dot(-1*dUf)
    RHS = (T_index.T).dot(RHS)

    print ("UMFSolve -----begin")

    dUr = dsolve.spsolve(K_gR,RHS,use_umfpack=True)
    print ("UMFSolve -----end")
    dU = T_index.dot(dUr) + dUf

    dSigma, dEps = EleTask.CaldSigmadEps(dU, FemInf.Node_xy, FemInf.Elem, FemInf.Nnode, FemInf.Nelem, MatInf, PhaseCoef_Gauss)

    FemInf.U_tot = dU
    FemInf.Sigma = dSigma
    FemInf.Eps   = dEps

    FemInf.Hmax = EleTask.CalHmax(FemInf,MatInf,PhaseCoef_Gauss)


    Fn = EleTask.CalFn(FemInf.Node_xy, FemInf.Elem, FemInf.Nnode, FemInf.Nelem, FemInf.Sigma)
    F_output=0.
    for i in range(FemInf.Nnode):
        if(FemInf.Node_xy[i,1]>0.9999):
            F_output = F_output + Fn[2*i]

    print ("Dis = %f, Reaction Force Fy= %f", FemInf.U_tot[2401*2], F_output)
    return FemInf