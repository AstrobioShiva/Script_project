from libraries import *
from functions import Rabc, Quad, ChemShift, AntiShift
#*******************************Input Parameters********************************************
Nptx = 64
B0 = 14.1 #magnetic field
w0 = 192.5 #input for 11B
# print(w0)
wX = w0*10**6 #actual freq in Hz

dangle  = np.pi/Nptx

#Spin Quantum Number
Ispin = 3/2

########### Quadrupolar coupling tensor #########################
#coupling values provided to match the magnitude of simulated with experimental data
# CQ_M = 0.01 #CQ in MHz
# Qeta = 0.9979 #eta of Q

# # Symmetric 2nd-rank chemical shift anisotropy (CSA)   tensor
# # Siso_ppm = 0 #isotropic chemical shift + offset(ppm) #value for Siso_ppm to match the magnitude of simulated with experimental data
# Siso_ppm = 0
   
# delta_ppm = -10  #chemical shift anisotropy (CSA) (ppm)
# eta = -0.5  #eta of CSA

CQ_M = 2.182
Qeta = 0.915 
Siso_ppm = 0
csa_ppm = -0
eta = 0.3  

######### Antisymmetric 1st-rank chemical shift tensor ##########
Sxy_set = [500]; Sxz_set = [500]; Syz_set = [500]                #Use ACS values for different sites

for k in (range(len(Sxy_set))):
    Sxy = Sxy_set[k]; Sxz = Sxz_set[k]; Syz = Syz_set[k]
#************************************************************************************************

    freq3212 = np.zeros((Nptx+1,1))
    freq1212 = np.zeros((Nptx+1,1))
    freq1232 = np.zeros((Nptx+1,1))
    freqSUM = np.zeros((Nptx+1,1))
    XX = np.zeros((Nptx+1,1))
    freqDIFF = np.zeros((Nptx+1,1))
    qacs = np.zeros((Nptx+1,1))
    qcsa = np.zeros((Nptx+1,1))
    qacs_anti = np.zeros((Nptx+1,1))


    CQ = CQ_M*10**6/(2*Ispin*(2*Ispin-1))
    Sxy = Sxy*w0; Sxz = Sxz*w0; Syz = Syz*w0;

    csa = csa_ppm*w0
    Siso = Siso_ppm*w0 

    #*****Relative Tensor Orientations
    #input parameters
    #       {a,b,c)       {zeta,lamda,nu}         {alpha,beta,gama}           {phi,theta, 0}
    # CSA===========>Quad==================>X-tal=======================>Gon=================>Rot. Frame

    a, b, c = 362.3*np.pi/180, 87.8*np.pi/180, 75.2*np.pi/180
    zeta, lamda, nu = 82.1*np.pi/180, 37*np.pi/180, 335*np.pi/180
    alpha, beta, gamma = 0*np.pi/180, 0*np.pi/180, 0*np.pi/180

    # tensor parameter at PAS
    QPAS = np.zeros((3, 3))
    Csa = np.zeros((3, 3))
    Acs = np.zeros((3, 3))

    #Quadrupolar interactions
    QPAS[0, 0] = (Qeta - 1) * CQ / 2
    QPAS[1, 1] = -(1 + Qeta) * CQ / 2
    QPAS[2, 2] = CQ

    # 2nd-rank chemical shift anisotropy (CSA)

    Csa[0, 0] = (eta - 1) * csa / 2  + Siso
    Csa[1, 1] = -(1 + eta) * csa / 2 + Siso
    Csa[2, 2] = csa + Siso


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

    U = Rabc(alpha, beta, gamma) 
    QXG = np.matmul(np.matmul(U, QX), np.linalg.inv(U)); 
    CsaQXG = np.matmul(np.matmul(U, CsaQX), np.linalg.inv(U));
    AcsQXG = np.matmul(np.matmul(U, AcsQX), np.linalg.inv(U)); 
    CSQXG = np.matmul(np.matmul(U, CSQX), np.linalg.inv(U));


            
    text = ['z rotation', 'y rotation', 'x rotation']
    df = [0]*len(text)
    file_path = ('/home/shiva/WMU/PhD/Scripts/phd_project/Python/NMR/Data_simulation/')
    for i in (range(len(text))):
        ang = 0 #starting angle
        for j in range(0, Nptx+1):
            
            aphi = [-(ang), 0, np.pi/2]   #changes to be made from original code as initial position for x, z orientation should match the diagram in Vosegaard et. al [-x rot: (phi = -pi/2, theta = 0, 0)] [-z rot: (phi = pi/2, theta = -pi/2, 0)]?
            atheta = [np.pi/2, ang, -ang] #-z, y, -x rotation
            
            theta = atheta[i] 
            phi = aphi[i]
            
            ct = np.cos(theta); st = np.sin(theta); s2t = 2*ct*st; c2t = ct*ct-st*st;
            cp = np.cos(phi); sp = np.sin(phi); s2p = 2*cp*sp; c2p = cp*cp-sp*sp;
            
            R2m2Q, R2m1Q, R20Q, R2p1Q, R2p2Q = Quad(QXG,ct, st, s2t, c2t, cp ,sp ,c2p ,s2p)
            R2m1cs, R20cs, R2p1cs = ChemShift(CsaQXG,ct, st, s2t, c2t, cp, sp, c2p, s2p)
            R1m1acs, R1p1acs = AntiShift(AcsQXG,ct, st, cp, sp);

            #********************change made from original code according to theor415y in next line**************************
            HCSA1 = (R20cs + 1/3*(CsaQXG[2, 2] + CsaQXG[0,0] + CsaQXG[1,1]))                                                # Siso = 1/3*(CsaQXG[2, 2] + CsaQXG[0,0] + CsaQXG[1,1]); Siso added to CSA tensor
            HQCSA = -(0.5/(2*Ispin*(2*Ispin-1)))*(R2m1Q*R2p1cs+R2p1Q*R2m1cs)/wX; #change made based on equations in doc *****multiplied factor of -0.5/2I(2I-1)
            HQACS = -(0.5/(2*Ispin*(2*Ispin-1)))*(R2p1Q*R1m1acs - R2m1Q*R1p1acs)/wX; #change made based on equations in doc *****multiplied factor of -0.5/2I(2I-1)

            
            HQ1 = 1/(2*Ispin*(2*Ispin-1))*R20Q
            HQ2a= (0.5/(2*Ispin*(2*Ispin-1))**2)*(R2m2Q*R2p2Q)/wX                         
            HQ2b = (0.5/(2*Ispin*(2*Ispin-1))**2)*(R2m1Q*R2p1Q)/wX

            # HQ2a= 0.5*(R2m2Q*R2p2Q)/wX                         
            # HQ2b = 0.5*(R2m1Q*R2p1Q)/wX 

        #********************change made from original code according to theory in next line**************************
            #****************coefficients are from code nmr_eq_coefficients*********************************************
            freq3212[j] = 6*(np.real(HQ1) + np.real(HQCSA) + np.real(HQACS)) + 1*np.real(HCSA1) + 0*np.real(HQ2a) + -12*np.real(HQ2b)      # 3/2 <-> 1/2
            freq1212[j] = 0*(np.real(HQ1) + np.real(HQCSA) + np.real(HQACS)) + 1*np.real(HCSA1) + 6*np.real(HQ2a) + 12*np.real(HQ2b) # 1/2 <-> -1/2
            freq1232[j] = -6*(np.real(HQ1) + np.real(HQCSA) + np.real(HQACS)) + 1*np.real(HCSA1) + 0*np.real(HQ2a) + -12*np.real(HQ2b)  #-1/2   <-> -3/2                                                                                    #-1/2<->-3/2

            freqSUM[j] = freq3212[j]+freq1212[j] + freq1232[j]          #3/2 <-> 1/2 + 1/2 <-> -1/2 + -1/2   <-> -3/2  transition 
            freqDIFF[j] = freq3212[j]-freq1232[j]                     #3/2 <-> 1/2 - -1/2   <-> -3/2  transition 
            
            #the coeffient are taken from coefficients in freq equation
            qcsa[j] = 6*np.real(HQCSA)
            qacs[j] = 6*np.real(HQACS)
            qacs_anti[j] = 6*np.real(-HQACS)
        
            
            XX[j]=(ang*180/np.pi)
            ang = ang + dangle
        
        data = np.column_stack((XX, np.real(freqSUM), np.real(freqDIFF), np.real(freq3212), np.real(freq1212), np.real(freq1232) ,np.real(qacs), np.real(qcsa), np.real(qacs_anti)))

        df[i]= pd.DataFrame(data, columns=['angle', 'freq_sum', 'freq_diff', 'freq_3212', 'freq_1212', 'freq_1232','acs', 'csa', 'acs_anti'])

        df[i].to_csv(f'{file_path}{text[i]}_set{k+1}_B.csv', index=True)