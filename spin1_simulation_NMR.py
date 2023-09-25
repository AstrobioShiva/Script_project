from libraries import *
from functions import Rabc, Quad, ChemShift, AntiShift
#*******************************Input Parameters********************************************
Nptx = 100
# w0 = float(input('Enter Larmor Frequency for nuclei in MHz:')) #Larmor fre in MHz #for N14
B0 = 14.1 #magnetic field
w0 = 43.4 #input for 14N
# print(w0)
wX = w0*10**6 #actual freq in Hz

ntheta  = 500
nphi = 500

dangle  = 2*np.pi/Nptx

#Spin Quantum Number
Ispin = 1.

#Quadrupolar Coupling Tensor
#coupling values for NAV (taken from paper https://pubmed.ncbi.nlm.nih.gov/22027340/)

CQ_M = -3.03 #CQ in MHz
# CQ = 0
Qeta = 0.45 #eta of Q

# Symmetric 2nd-rank chemical shift anisotropy (CSA)   tensor
# Siso_ppm = 0.     #isotropic chemical shift + offset(ppm)
delta_ppm = 110.0  #chemical shift anisotropy (CSA) (ppm)
eta = 0.415   #eta of CSA


Sxy = 50; Sxz = 50; Syz = 50

#************************************************************************************************

freq1D = np.zeros((Nptx+1,1))
freq2D = np.zeros((Nptx+1,1))
freqSUM = np.zeros((Nptx+1,1))
XX = np.zeros((Nptx+1,1))
freqDIFF = np.zeros((Nptx+1,1))
qacs = np.zeros((Nptx+1,1))
qcsa = np.zeros((Nptx+1,1))
qacs_anti = np.zeros((Nptx+1,1))


CQ = CQ_M*10**6/(2*Ispin*(2*Ispin-1))
Sxy = Sxy*w0; Sxz = Sxz*w0; Syz = Syz*w0;

delta = delta_ppm*w0


#*****Relative Tensor Orientations
#input parameters
#       {a,b,c)       {zeta,lamda,nu}         {alpha,beta,gama}           {phi,theta, 0}
# CSA===========>Quad==================>X-tal=======================>Gon=================>Lab

a, b, c = 0*np.pi/180, 0*np.pi/180, 0*np.pi/180
zeta, lamda, nu = 0*np.pi/180, 0*np.pi/180, 0*np.pi/180
alpha, beta, gama = 0*np.pi/180, 0*np.pi/180, 0*np.pi/180

# tensor parameter at PAS
QPAS = np.zeros((3, 3))
Csa = np.zeros((3, 3))
Acs = np.zeros((3, 3))

#Quadrupolar interactions
QPAS[0, 0] = (Qeta - 1) * CQ / 2
QPAS[1, 1] = -(1 + Qeta) * CQ / 2
QPAS[2, 2] = CQ

# 2nd-rank chemical shift anisotropy (CSA)

Csa[0, 0] = (eta - 1) * delta / 2
Csa[1, 1] = -(1 + eta) * delta / 2
Csa[2, 2] = delta
##############################################

Siso_ppm = 0
Siso = Siso_ppm*w0
##############################################

# 1st-rank antisymmetric chemical shift (ACS)

Acs[0, 1] = Sxy;  Acs[0, 2] = Sxz
Acs[1, 0] = -Sxy; Acs[1, 2] = Syz
Acs[2, 0] = -Sxz; Acs[2, 1] = -Syz

# Combining CSA and ACS
CSPAS = Csa + Acs


#Transformation between different frames
U = Rabc(a, b, c)
CsaQ = np.matmul(np.matmul(U, Csa), np.linalg.inv(U))
AcsQ = np.matmul(np.matmul(U, Acs), np.linalg.inv(U))
CSQ = np.matmul(np.matmul(U, CSPAS), np.linalg.inv(U))

U = Rabc(zeta, lamda, nu)
QX = np.matmul(np.matmul(U, QPAS), np.linalg.inv(U))
CsaQX = np.matmul(np.matmul(U, CsaQ), np.linalg.inv(U))
AcsQX = np.matmul(np.matmul(U, AcsQ), np.linalg.inv(U))
CSQX = np.matmul(np.matmul(U, CSQ), np.linalg.inv(U))

U = Rabc(alpha, beta, gama) 
QXG = np.matmul(np.matmul(U, QX), np.linalg.inv(U)); 
CsaQXG = np.matmul(np.matmul(U, CsaQX), np.linalg.inv(U));
AcsQXG = np.matmul(np.matmul(U, AcsQX), np.linalg.inv(U)); 
CSQXG = np.matmul(np.matmul(U, CSQX), np.linalg.inv(U));

def data_file():
    
    text = ['z rotation', 'y rotation', 'x rotation']
    df = [0]*3
    for i in (range(len(text))):
        ang = 0 #starting angle
        for j in range(0, Nptx+1):
            

            atheta = [np.pi/2, ang, -ang] #-z, y, -x rotation
            aphi = [-ang, 0, np.pi/2]
            
            theta = atheta[i] 
            phi = aphi[i]
            
            ct = np.cos(theta); st = np.sin(theta); s2t = 2*ct*st; c2t = ct*ct-st*st;
            cp = np.cos(phi); sp = np.sin(phi); s2p = 2*cp*sp; c2p = cp*cp-sp*sp;
            
            R2m2Q, R2m1Q, R20Q, R2p1Q, R2p2Q = Quad(QXG,ct, st, s2t, c2t, cp ,sp ,c2p ,s2p)
            R2m1cs, R20cs, R2p1cs = ChemShift(CsaQXG,ct, st, s2t, c2t, cp, sp, c2p, s2p)
            R1m1acs, R1p1acs = AntiShift(AcsQXG,ct, st, cp, sp);

        #********************change made from original code according to theory in next line**************************
            HCSA1 = (1/3*Siso + R20cs)
            HQCSA = -0.5/(2*Ispin*(2*Ispin - 1))*(R2m1Q*R2p1cs + R2p1Q*R2m1cs)/wX;
            HQACS = 0.5/(2*Ispin*(2*Ispin - 1))*(R2m1Q*R1p1acs - R2p1Q*R1m1acs)/wX;
            HQ1 = 0.25/(2*Ispin*(2*Ispin - 1))*R20Q
            HQ2 = 0.5/(2*Ispin*(2*Ispin - 1))**2*(R2m2Q*R2p2Q - R2m1Q*R2p1Q)/wX


        #********************change made from original code according to theory in next line**************************
            freq1D[j] = 3*np.real(HQ1) + np.real(HQ2) + np.real(HCSA1) + 3*np.real(HQCSA) + 3*np.real(HQACS) # 1 <-> 0
            freq2D[j] = -3*np.real(HQ1 )+ np.real(HQ2) + np.real(HCSA1) + -3*np.real(HQCSA) + -3*np.real(HQACS); # 0 <-> -1

            freqSUM[j] = freq1D[j]+freq2D[j]           #1 <-> -1 transition 
            freqDIFF[j] = freq1D[j]-freq2D[j]          #1 <-> -1 transition 
            

            qcsa[j] = 3*np.real(HQCSA)
            qacs[j] = 3*np.real(HQACS)
            qacs_anti[j] = 3*np.real(-HQACS)
        
            
            XX[j]=ang*180/np.pi
            ang = ang + dangle
        
        data = np.column_stack((XX, np.real(freqSUM), np.real(freqDIFF), np.real(freq1D), np.real(freq2D), np.real(qacs), np.real(qcsa), np.real(qacs_anti)))

        df[i]= pd.DataFrame(data, columns=['angle', 'freq_sum', 'freq_diff', 'freq_1D', 'freq_2D', 'acs', 'csa', 'acs_anti'])



    return df
    


# def data_file2():
#     df = dict()
    
#     # df = [0]*3
#     for i in range(0,3):
#         ang = 0 #starting angle
#         for j in range(0, Nptx+1):
#             text = ['z', 'y', 'x']

#             atheta = [np.pi/2, ang, -ang] #-z, y, -x rotation
#             aphi = [-ang, 0, np.pi/2]
            
#             theta = atheta[i] 
#             phi = aphi[i]
            
#             ct = np.cos(theta); st = np.sin(theta); s2t = 2*ct*st; c2t = ct*ct-st*st;
#             cp = np.cos(phi); sp = np.sin(phi); s2p = 2*cp*sp; c2p = cp*cp-sp*sp;
            
#             R2m2Q, R2m1Q, R20Q, R2p1Q, R2p2Q = Quad(QXG,ct, st, s2t, c2t, cp ,sp ,c2p ,s2p)
#             R2m1cs, R20cs, R2p1cs = ChemShift(CsaQXG,ct, st, s2t, c2t, cp, sp, c2p, s2p)
#             R1m1acs, R1p1acs = AntiShift(AcsQXG,ct, st, cp, sp);

#             HCSA1 = (1/3*Siso + R20cs)
#             HQCSA = -0.5/(2*Ispin*(2*Ispin-1))*(R2m1Q*R2p1cs+R2p1Q*R2m1cs)/wX;
#             HQACS = 0.5/(2*Ispin*(2*Ispin-1))*(R2m1Q*R1p1acs-R2p1Q*R1m1acs)/wX;
#             HQ1 = 0.25/(2*Ispin*(2*Ispin-1))*R20Q
#             HQ2 = 0.5/(2*Ispin*(2*Ispin-1))**2*(R2m2Q*R2p2Q-R2m1Q*R2p1Q)/wX

#         #change made from original code according to theory in next line
#             freq1D[j] = 3*np.real(HQ1) + np.real(HQ2) + np.real(HCSA1) + 3*np.real(HQCSA) + 3*np.real(HQACS) # 1 <-> 0
#             freq2D[j] = -3*np.real(HQ1 )+ np.real(HQ2) + np.real(HCSA1) + -3*np.real(HQCSA) + -3*np.real(HQACS); # 0 <-> -1

#             freqSUM[j] = freq1D[j]+freq2D[j]           #1 <-> -1 transition 
#             freqDIFF[j] = freq1D[j]-freq2D[j]          #1 <-> -1 transition 
            

#             qcsa[j] = 3*np.real(HQCSA)
#             qacs[j] = 3*np.real(HQACS)
            
#             XX[j]=ang*180/np.pi
#             ang = ang + dangle
        
#         data = np.column_stack((XX, np.real(freqSUM), np.real(freqDIFF), np.real(freq1D), np.real(freq2D), np.real(qacs), np.real(qcsa)))
#         df[f'{i}']= pd.DataFrame(data, columns=['angle', 'freq_sum', 'freq_diff', 'freq_1D', 'freq_2D', 'ACS', 'CSA'])
        
#     return df
    

# if __name__ =='__main__':
#     df=data_file()
#     # print(df['z'])
#     x=df[2]['freq_sum']
#     print(x)
#     # for key,value in df.items():
#     #     print(key, df)



# # len(text)
# # df=data_file()