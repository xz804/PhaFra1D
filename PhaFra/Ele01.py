## 6-node Triangle element  Tri6
## 3 gauss point

import numpy as np
import math
from DataStruct import *
from scipy import sparse


def TranPhase(Elem,Nelem,Phase,PhaseEta):

    NGP = 3
    NPE = 6

    w,r,s,t = QuadraturePoint()

    Phase_Gauss = np.zeros((Nelem,NGP),dtype=np.float)
    for ielement in range(Nelem):

        Phase_e = np.zeros((NPE,),dtype=np.float)
        for inode in range(NPE):
            Phase_e[inode]     =   Phase[Elem[ielement,inode+1]-1]

        for igauss in range(NGP):

            N_Phase = np.array([(2.*r[igauss]-1.)*r[igauss], (2.*s[igauss]-1.)*s[igauss], (2.*t[igauss]-1.)*t[igauss],
                                4.*r[igauss]*s[igauss],      4.*s[igauss]*t[igauss],      4.*t[igauss]*r[igauss]])

            Phase_Gauss[ielement,igauss] = N_Phase.dot(Phase_e)

    Phase_tem = 1.- Phase_Gauss
    PhaseCoef_Gauss = Phase_tem*Phase_tem+PhaseEta
    return PhaseCoef_Gauss

def GlobalStiffness(Nnode,Node_xy,Nelem,Elem,MatInf,PhaseCoef_Gauss):

    NPE = 6

    Kg_irow  = np.zeros((12*12*Nelem,),dtype = np.int)
    Kg_icol  = np.zeros((12*12*Nelem,),dtype = np.int)
    Kg_V     = np.zeros((12*12*Nelem,),dtype = np.float)
    K_g_index = 0
    for ielement in range(Nelem):

        if ielement%20000 == 0:
            print ("ielement = %d", ielement)

        Eldof = np.zeros((NPE*2,),dtype=np.int)

        for inode in range(NPE):
            Eldof[inode*2]   = (Elem[ielement,inode+1]-1)*2
            Eldof[inode*2+1] = (Elem[ielement,inode+1]-1)*2+1

        Node_xy_e = np.zeros((6,2),dtype=np.float)

        for inode in range(NPE):
            Node_xy_e[inode,:] = Node_xy[Elem[ielement,inode+1]-1,1:3]

        # calculate element stiffness matrix K_e(NPE*2,NPE*2)

        K_e = CalKe_Elas(ielement, Node_xy_e, MatInf, PhaseCoef_Gauss[ielement,:])

        # Generate irow icol and V for Constructing Global stiffness matrix

        for i in range(12):
            for j in range(12):
                Kg_irow[K_g_index] = Eldof[i]
                Kg_icol[K_g_index] = Eldof[j]
                Kg_V   [K_g_index] = K_e[i,j]
                K_g_index = K_g_index + 1

    #========== loop over K_e ----finished


    #------- form global stiffness matrix K_g
    K_g = sparse.coo_matrix((Kg_V,(Kg_irow,Kg_icol)),shape=(Nnode*2,Nnode*2))

    return K_g

def CalKe_Elas(ielement, Node_xy_e, MatInf, Phase_e):
    """
    Calculation of element stiffness matrix
    """

    #------Gauss point
    NGP = 3
    w,r,s,t = QuadraturePoint()

    if ielement == 0 or ielement == 1 or ielement == 2 or ielement == 3:
        E     = MatInf.E/10.
    else:
        E     = MatInf.E

    Gamma     = MatInf.Gamma
    ProbType  = MatInf.ProbType
    PhaseType = MatInf.PhaseType

    if ProbType == "PlaneStrain":
        D_e = np.array([[1.-Gamma,    Gamma,              0.],
                        [   Gamma, 1.-Gamma,              0.],
                        [      0.,       0., (1-2.*Gamma)/2.]],dtype=np.float) * E/(1+Gamma)/(1-2.*Gamma)
    elif ProbType == "PlaneStress":
        D_e = np.array([[       1,    Gamma,              0.],
                        [   Gamma,       1.,              0.],
                        [      0.,       0.,    (1-Gamma)/2.]],dtype=np.float) * E/(1-Gamma*Gamma)

    K_e = np.zeros((12,12),dtype=np.float)

    # ----- loop over Gauss point
    for igauss in range(NGP):

        dNdr = np.array([[4.*r[igauss]-1.,                               0.,   -3.+4.*(r[igauss]+s[igauss]),
                             4.*s[igauss],                   -4.*s[igauss],  4.-8.*r[igauss]- 4.*s[igauss]],
                         [              0.,                 4.*s[igauss]-1.,    -3.+4.*(r[igauss]+s[igauss]),
                             4.*r[igauss],  4.-4.*r[igauss]- 8.*s[igauss],                    -4.*r[igauss]]],dtype=np.float)

        J    = np.dot(dNdr,Node_xy_e)
        detJ = J[0,0]*J[1,1]-J[1,0]*J[0,1]
        invJ = np.array([[    J[1,1], -1.*J[0,1]],
                         [-1.*J[1,0],     J[0,0]]])/detJ

        dNdx = np.dot(invJ,dNdr)
        Bu_e = np.zeros((3,12),dtype=np.float)

        for i in range(6):
            Bu_e[0,2*i]   = dNdx[0,i]
            Bu_e[1,2*i+1] = dNdx[1,i]
            Bu_e[2,2*i]   = dNdx[1,i]
            Bu_e[2,2*i+1] = dNdx[0,i]

        BtD  =  np.dot(Bu_e.T,D_e)
        BtDB =  np.dot(BtD,Bu_e)

        K_e  =  K_e + w[igauss]*detJ*BtDB*Phase_e[igauss]


    return K_e

def CalFn(Node_xy, Elem, Nnode, Nelem, Sigma):
    """
    Node_xy  Node information for the Element
    NGP        Gauss Integration point
    """

    #------Gauss point
    NGP = 3
    NPE = 6
    w,r,s,t = QuadraturePoint()

    Fn = np.zeros((2*Nnode,),dtype=np.float)
    for ielement in range(Nelem):

        Eldof = np.zeros((NPE*2,),dtype=np.int)
        for inode in range(NPE):
            Eldof[inode*2]   = (Elem[ielement,inode+1]-1)*2
            Eldof[inode*2+1] = (Elem[ielement,inode+1]-1)*2+1

        Node_xy_e = np.zeros((6,2),dtype=np.float)

        for inode in range(NPE):
            Node_xy_e[inode,:] = Node_xy[Elem[ielement,inode+1]-1,1:3]

        Fn_e = np.zeros((2*NPE,),dtype=np.float)
        for iGpoint in range(NGP):
            dNdr = np.array([[4.*r[iGpoint]-1.,                               0.,   -3.+4.*(r[iGpoint]+s[iGpoint]),
                                 4.*s[iGpoint],                   -4.*s[iGpoint],  4.-8.*r[iGpoint]- 4.*s[iGpoint]],
                             [              0.,                 4.*s[iGpoint]-1.,    -3.+4.*(r[iGpoint]+s[iGpoint]),
                                 4.*r[iGpoint],  4.-4.*r[iGpoint]- 8.*s[iGpoint],                    -4.*r[iGpoint]]],dtype=np.float)

            J    = np.dot(dNdr,Node_xy_e)
            detJ = J[0,0]*J[1,1]-J[1,0]*J[0,1]
            invJ = np.array([[    J[1,1], -1.*J[0,1]],
                             [-1.*J[1,0],     J[0,0]]])/detJ

            dNdx = np.dot(invJ,dNdr)
            Bu_e = np.zeros((3,12),dtype=np.float)

            for i in range(6):
                Bu_e[0,2*i]   = dNdx[0,i]
                Bu_e[1,2*i+1] = dNdx[1,i]
                Bu_e[2,2*i]   = dNdx[1,i]
                Bu_e[2,2*i+1] = dNdx[0,i]

            Sigma_local = np.zeros((3,),dtype=np.float)
            Sigma_local[0] = Sigma[ielement,iGpoint*4+0]
            Sigma_local[1] = Sigma[ielement,iGpoint*4+1]
            Sigma_local[2] = Sigma[ielement,iGpoint*4+3]

            BtSigma = (Bu_e.T).dot(Sigma_local)
            Fn_e  =  Fn_e + w[iGpoint]*detJ*BtSigma

        for i in range(12):
            Fn[Eldof[i]] = Fn[Eldof[i]] + Fn_e[i]

    return Fn

def CalDofIndex(Nnode, Nessb, EssBound, Ndof):

    FreeDof_index = np.zeros((2*Nnode,),dtype=np.int)  #  -1: ux is described;  -2: uy is described

    index = 0
    for i in range(Nessb):
        if EssBound[i][1] == 'Ux':
            index = index - 1
            FreeDof_index[EssBound[i][0]*2] = index
        elif EssBound[i][1] == 'Uy':
            index = index - 1
            FreeDof_index[EssBound[i][0]*2+1] = index
        else:
            print ("Error in FreeDof_index[i][1]")

    index = 0
    for i in range(2*Nnode):
        if FreeDof_index[i] == 0:
            FreeDof_index[i] = index
            index = index + 1

    #T_index for the recover of x from xr:  x = T_index*xr
    T_irow = np.zeros((Ndof,),dtype=np.int)
    T_icol = np.zeros((Ndof,),dtype=np.int)
    T_V    = np.ones((Ndof,),dtype=np.int)

    index = 0
    for i in range(2*Nnode):
        if FreeDof_index[i] >= 0:
            T_irow[index] = i
            T_icol[index] = FreeDof_index[i]
            index = index + 1


    T_index = sparse.coo_matrix((T_V,(T_irow,T_icol)),shape=(2*Nnode,Ndof))


    return T_index


def CaldSigmadEps(dU, Node_xy, Elem, Nnode, Nelem, MatInf, PhaseCoef_Gauss):
    """
    Calculation of Stress increment and Strain increment
    """

    ProbType  = MatInf.ProbType
    PhaseType = MatInf.PhaseType
    NGP = 3
    NPE = 6

    # Gauss Point
    w,r,s,t = QuadraturePoint()

    dSigma = np.zeros((Nelem,12),dtype=np.float)
    dEps   = np.zeros((Nelem,12),dtype=np.float)

    for ielement in range(Nelem):

        if ielement == 0 or ielement == 1 or ielement == 2 or ielement == 3:
            E     = MatInf.E/10.
        else:
            E     = MatInf.E

        Gamma = MatInf.Gamma

        if ProbType == "PlaneStrain":
            D_e = np.array([[1.-Gamma,    Gamma,              0.],
                            [   Gamma, 1.-Gamma,              0.],
                            [      0.,       0., (1-2.*Gamma)/2.]],dtype=np.float) * E/(1+Gamma)/(1-2.*Gamma)
        elif ProbType == "PlaneStress":
            D_e = np.array([[       1,    Gamma,              0.],
                            [   Gamma,       1.,              0.],
                            [      0.,       0.,    (1-Gamma)/2.]],dtype=np.float) * E/(1-Gamma*Gamma)

        Eldof = np.zeros((NPE*2,),dtype=np.int)
        for inode in range(NPE):
            Eldof[inode*2]   = (Elem[ielement,inode+1]-1)*2
            Eldof[inode*2+1] = (Elem[ielement,inode+1]-1)*2+1

        Node_xy_e = np.zeros((6,2),dtype=np.float)
        for inode in range(NPE):
            Node_xy_e[inode,:]   = Node_xy[Elem[ielement,inode+1]-1,1:3]

        dU_local  = np.zeros((12,),dtype=np.float)
        for i in range(12):
            dU_local[i] = dU[Eldof[i]]

        for igauss in range(NGP):

            dNdr = np.array([[4.*r[igauss]-1.,                               0.,   -3.+4.*(r[igauss]+s[igauss]),
                                 4.*s[igauss],                   -4.*s[igauss],  4.-8.*r[igauss]- 4.*s[igauss]],
                             [              0.,                 4.*s[igauss]-1.,    -3.+4.*(r[igauss]+s[igauss]),
                                 4.*r[igauss],  4.-4.*r[igauss]- 8.*s[igauss],                    -4.*r[igauss]]],dtype=np.float)

            J    = np.dot(dNdr,Node_xy_e)
            detJ = J[0,0]*J[1,1]-J[1,0]*J[0,1]
            invJ = np.array([[    J[1,1], -1.*J[0,1]],
                             [-1.*J[1,0],     J[0,0]]])/detJ

            dNdx = np.dot(invJ,dNdr)
            Bu_e = np.zeros((3,12),dtype=np.float)

            for i in range(6):
                Bu_e[0,2*i]   = dNdx[0,i]
                Bu_e[1,2*i+1] = dNdx[1,i]
                Bu_e[2,2*i]   = dNdx[1,i]
                Bu_e[2,2*i+1] = dNdx[0,i]

            dEps_local   = Bu_e.dot(dU_local)
            dSigma_local = D_e.dot(dEps_local)*PhaseCoef_Gauss[ielement,igauss]


            index = 4*igauss
            dEps[ielement,index]     = dEps_local[0]
            dEps[ielement,index+1]   = dEps_local[1]
            dEps[ielement,index+3]   = dEps_local[2]
            dSigma[ielement,index+0] = dSigma_local[0]
            dSigma[ielement,index+1] = dSigma_local[1]
            dSigma[ielement,index+3] = dSigma_local[2]

            if ProbType == "PlaneStrain":
                dEps[ielement,index+2]   = 0.
                dSigma[ielement,index+2] = Gamma*(dSigma[ielement,index+0]+dSigma[ielement,index+1])
            elif ProbType == "PlaneStress":
                dEps[ielement,index+2]   = -1.*Gamma*(dSigma[ielement,index+0]+dSigma[ielement,index+1])/E
                dSigma[ielement,index+2] = 0.


    return dSigma,dEps

def CalHmax(FemInf,MatInf,PhaseCoef_Gauss):

    if MatInf.PhaseType == 'PhaseAmor':
        H_tem = CalHmaxAmor (FemInf,MatInf)
    elif MatInf.PhaseType == 'PhaseIso':
        H_tem = CalHmaxIso  (FemInf,MatInf,PhaseCoef_Gauss)
    elif MatInf.PhaseType == 'PhaseMiehe':
        H_tem = CalHmaxMiehe(FemInf,MatInf)
    else:
        print("wrong phase option in CalHmax \n")

    Hmax = FemInf.Hmax
    NGP = 3
    for ielement in range(FemInf.Nelem):
        for igauss in range(NGP):
            if H_tem[ielement,igauss] > Hmax[ielement,igauss]:
                Hmax[ielement,igauss] = H_tem[ielement,igauss]

    return Hmax

def CalHmaxIso(FemInf,MatInf,PhaseCoef_Gauss):

    NGP = 3

    Nelem = FemInf.Nelem

    Eps   = FemInf.Eps
    Sigma = FemInf.Sigma

    Hmax = FemInf.Hmax

    for ielement in range(Nelem):
        for igauss in range(NGP):
            Eps_local = np.array([[Eps[ielement,4*igauss]],
                                  [Eps[ielement,4*igauss+1]],
                                  [Eps[ielement,4*igauss+3]]],dtype=np.float)

            Sigma_local = np.array([[Sigma[ielement,4*igauss]],
                                    [Sigma[ielement,4*igauss+1]],
                                    [Sigma[ielement,4*igauss+3]]],dtype=np.float)

            H_tem = (Sigma_local.T).dot(Eps_local)*0.5/PhaseCoef_Gauss[ielement,igauss]

            Hmax[ielement,igauss] = H_tem

    return Hmax


def CalHmaxAmor(FemInf,MatInf):

    NGP = 3

    E     = MatInf.E
    Gamma = MatInf.Gamma
    K     = E/3.0/(1-2.0*Gamma)
    mu    = E/2.0/(1+Gamma)

    Nelem = FemInf.Nelem
    Eps   = FemInf.Eps
    Hmax  = FemInf.Hmax

    Gfactor = MatInf.GcI/MatInf.GcII

    for ielement in range(Nelem):
        for igauss in range(NGP):
            Eps_local = np.array([[Eps[ielement,4*igauss  ]],
                                  [Eps[ielement,4*igauss+1]],
                                  [Eps[ielement,4*igauss+2]],
                                  [Eps[ielement,4*igauss+3]/2.0]],dtype=np.float)
            Eps_tr  = Eps_local[0]+Eps_local[1]+Eps_local[2]
            Eps_dev = np.array([[Eps[ielement,4*igauss  ]-Eps_tr/3.0],
                                [Eps[ielement,4*igauss+1]-Eps_tr/3.0],
                                [Eps[ielement,4*igauss+2]-Eps_tr/3.0],
                                [Eps[ielement,4*igauss+3]/2.0]],dtype=np.float)

            H_tem = 0.0
            if Eps_tr>=0:
                H_tem = 0.5*K*Eps_tr*Eps_tr

            H_tem = H_tem + mu*(Eps_dev[0]*Eps_dev[0] + Eps_dev[1]*Eps_dev[1] + Eps_dev[2]*Eps_dev[2] + 2.0*Eps_dev[3]*Eps_dev[3])*Gfactor

            Hmax[ielement,igauss] = H_tem

    return Hmax


def CalHmaxMiehe(FemInf,MatInf):

    NGP = 3

    E     = MatInf.E
    Gamma = MatInf.Gamma                      # Poisson ratio
#    K     = E/3.0/(1-2.0*Gamma)               # Bulk modulus
    Lame  = E*Gamma/(1.+Gamma)/(1.-2.*Gamma)  # Lame constant
    mu    = E/2.0/(1+Gamma)                   # Shear modulus

    Nelem = FemInf.Nelem
    Eps   = FemInf.Eps
    Hmax  = FemInf.Hmax

    Gfactor = MatInf.GcI/MatInf.GcII

    for ielement in range(Nelem):
        for igauss in range(NGP):
            Eps_local = np.array([[Eps[ielement,4*igauss  ]],
                                  [Eps[ielement,4*igauss+1]],
                                  [Eps[ielement,4*igauss+2]],
                                  [Eps[ielement,4*igauss+3]/2.0]],dtype=np.float)

            tem1 = Eps_local[0] + Eps_local[1]
            tem2 = 4.0*Eps_local[3]*Eps_local[3]+(Eps_local[0] - Eps_local[1])*(Eps_local[0] - Eps_local[1])

            EigenV1 = 0.5*(tem1+math.sqrt(tem2))
            EigenV2 = 0.5*(tem1-math.sqrt(tem2))

            H_tem = 0.0

            Eps_tr = EigenV1+EigenV2
            if Eps_tr>0:
                H_tem = H_tem + 0.5*Lame*Eps_tr*Eps_tr

            if EigenV1>0:
                H_tem = H_tem + mu*EigenV1*EigenV1*Gfactor

            if EigenV2>0:
                H_tem = H_tem + mu*EigenV2*EigenV2*Gfactor


            Hmax[ielement,igauss] = H_tem


    return Hmax

def QuadraturePoint():

    w = np.array([0.1666666667,0.1666666667,0.1666666667])   # weight
    r = np.array([0.6666666667,0.1666666667,0.1666666667])   # x
    s = np.array([0.1666666667,0.6666666667,0.1666666667])   # y
    t = np.array([0.1666666667,0.1666666667,0.6666666667])   # 1-x-y

    return w, r, s, t