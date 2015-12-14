import numpy as np

from scipy import sparse

from scipy.sparse.linalg import dsolve

def SolPhaseField(FemInf,MatInf,MacroInf):
    Node_xy = FemInf.Node_xy
    Elem    = FemInf.Elem
    Hmax   = FemInf.Hmax
    GcI     = MatInf.GcI
    Lc      = MatInf.Lc

    Nnode = FemInf.Nnode
    Nelem = FemInf.Nelem

    NPE  = FemInf.Npoint  # Node Number per element
    NGP  = FemInf.Ngauss  # Number of Gauss Point per element

    w,r,s,t = QuadraturePoint()

    Kg_irow  = np.zeros((6*6*Nelem,),dtype = np.int)
    Kg_icol  = np.zeros((6*6*Nelem,),dtype = np.int)
    Kg_V     = np.zeros((6*6*Nelem,),dtype = np.float)
    K_g_index = 0

    Fn = np.zeros((Nnode),dtype=np.float)

    for ielement in range(Nelem):

        Node_xy_e = np.zeros((6,2),dtype=np.float)
        for inode in range(NPE):
            Node_xy_e[inode,:]   = Node_xy[Elem[ielement,inode+1]-1,1:3]

        Eldof = np.zeros((NPE,),dtype=np.int)

        for inode in range(NPE):
            Eldof[inode]   = Elem[ielement,inode+1]-1

        K1_e = np.zeros((6,6),dtype=np.float)
        K2_e = np.zeros((6,6),dtype=np.float)
        Fn_e = np.zeros((6,),dtype=np.float)
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
            Bs_e = dNdx


            Ns_e = np.matrix([(2.*r[igauss]-1.)*r[igauss], (2.*s[igauss]-1.)*s[igauss], (2.*t[igauss]-1.)*t[igauss],
                                4.*r[igauss]*s[igauss],      4.*s[igauss]*t[igauss],      4.*t[igauss]*r[igauss]])



            K1_e = K1_e + w[igauss]*detJ*(Ns_e.T*Ns_e)*(Hmax[ielement,igauss]+GcI/4./Lc)

            Ndx=np.matrix(Bs_e[0,:])
            Ndy=np.matrix(Bs_e[1,:])

            K2_e = K2_e + w[igauss]*detJ*(Ndx.T*Ndx+Ndy.T*Ndy)*GcI*Lc

            Fn_e = Fn_e + w[igauss]*detJ*1.*Hmax[ielement,igauss]*Ns_e

        for i in range(6):
            Fn[Eldof[i]] = Fn[Eldof[i]] + Fn_e[0,i]
            for j in range(6):
                Kg_irow[K_g_index] = Eldof[i]
                Kg_icol[K_g_index] = Eldof[j]
                Kg_V   [K_g_index] = K1_e[i,j] + K2_e[i,j]
                K_g_index = K_g_index + 1

    K_g = sparse.coo_matrix((Kg_V,(Kg_irow,Kg_icol)),shape=(Nnode,Nnode))

    K_g = K_g.tocsc()
    Phase = dsolve.spsolve(K_g,Fn,use_umfpack=True)

    FemInf.Phase = Phase
    return FemInf

def QuadraturePoint():

    w = np.array([0.1666666667,0.1666666667,0.1666666667])   # weight
    r = np.array([0.6666666667,0.1666666667,0.1666666667])   # x
    s = np.array([0.1666666667,0.6666666667,0.1666666667])   # y
    t = np.array([0.1666666667,0.1666666667,0.6666666667])   # 1-x-y

    return w, r, s, t