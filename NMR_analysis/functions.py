import numpy as np

def Rabc(alfa1, beta1, gama1): #Euler Rotation Matrix
    U = np.zeros((3, 3))
    U[0, 0] = np.cos(alfa1)*np.cos(beta1)*np.cos(gama1) - np.sin(alfa1)*np.sin(gama1)
    U[0, 1] = np.sin(alfa1)*np.cos(beta1)*np.cos(gama1) + np.cos(alfa1)*np.sin(gama1)
    U[0, 2] = -np.sin(beta1)*np.cos(gama1)
    U[1, 0] = -np.cos(alfa1)*np.cos(beta1)*np.sin(gama1) - np.sin(alfa1)*np.cos(gama1)
    U[1, 1] = -np.sin(alfa1)*np.cos(beta1)*np.sin(gama1) + np.cos(alfa1)*np.cos(gama1)
    U[1, 2] = np.sin(beta1)*np.sin(gama1)
    U[2, 0] = np.cos(alfa1)*np.sin(beta1)
    U[2, 1] = np.sin(alfa1)*np.sin(beta1)
    U[2, 2] = np.cos(beta1)
    return U

def Quad(QXG, ct, st, s2t, c2t, cp, sp, c2p, s2p): #Quadrupolar Interaction Function
    sq15 = np.sqrt(1.5)
    sq6 = np.sqrt(6)
    #**********************change made from original code according to theory in next line***************************
    R20Q = 0.25*(QXG[2,2]*(3*ct*ct-1)+2*QXG[0,2]*s2t*cp+2*QXG[1,2]*s2t*sp+(QXG[0,0]-QXG[1,1])*st*st*c2p+2*QXG[0,1]*st*st*s2p) #factor of 1/4 multiplied
    R2q1 = 0.75*QXG[2,2]*s2t
    R2q2r = QXG[0,2]*c2t*cp
    R2q2i = QXG[0,2]*ct*sp
    R2q3r = QXG[1,2]*c2t*sp
    R2q3i = QXG[1,2]*ct*cp
    R2q4r = 0.25*(QXG[0,0]-QXG[1,1])*s2t*c2p
    R2q4i = 0.5*(QXG[0,0]-QXG[1,1])*st*s2p
    R2q5r = 0.5*QXG[0,1]*s2t*s2p
    R2q5i = QXG[0,1]*st*c2p
    R2p1Q = R2q1 - R2q2r + 1j*R2q2i - R2q3r - 1j*R2q3i - R2q4r + 1j*R2q4i - R2q5r - 1j*R2q5i
    R2m1Q = -R2q1 + R2q2r + 1j*R2q2i + R2q3r - 1j*R2q3i + R2q4r + 1j*R2q4i + R2q5r - 1j*R2q5i;
  
    R4q1=0.75*QXG[2,2]*st*st;
    R4q2r=0.5*QXG[0,2]*s2t*cp; 
    R4q2i=QXG[0,2]*st*sp;
    R4q3r=0.5*QXG[1,2]*s2t*sp; 
    R4q3i=QXG[1,2]*st*cp;
    R4q4r=0.25*(QXG[0,0]-QXG[1,1])*(1+ct*ct)*c2p; 
    R4q4i=0.5*(QXG[0,0]-QXG[1,1])*ct*s2p;
    R4q5r=0.5*QXG[0,1]*(1+ct*ct)*s2p; 
    R4q5i=QXG[0,1]*ct*c2p;
    R2p2Q=R4q1-R4q2r+1j*R4q2i-R4q3r-1j*R4q3i+R4q4r-1j*R4q4i+R4q5r+1j*R4q5i;
    R2m2Q=R4q1-R4q2r-1j*R4q2i-R4q3r+1j*R4q3i+R4q4r+1j*R4q4i+R4q5r-1j*R4q5i;
    
    return [R2m2Q, R2m1Q, R20Q, R2p1Q, R2p2Q]
    

def ChemShift(CsaQXG, ct, st, s2t, c2t, cp, sp, c2p, s2p):
    sq15 = np.sqrt(1.5)
    #***************************change made from original code according to theory in next line********************************
    R20cs = 0.5 *(1/3*(2 * CsaQXG[2, 2] - CsaQXG[0, 0] - CsaQXG[1, 1]) * (3*ct ** 2 - 1) +2 * CsaQXG[0, 2] * s2t * cp + 2 * CsaQXG[1, 2] * s2t * sp +(CsaQXG[0, 0] - CsaQXG[1, 1]) * st ** 2 * c2p + 2 * CsaQXG[0, 1] * st ** 2 * s2p)
    
    R2cs0 = 0.25 * (2 * CsaQXG[2, 2] - CsaQXG[0, 0] - CsaQXG[1, 1]) * s2t
    R2cs1r = CsaQXG[0, 2] * c2t * cp
    R2cs1i = CsaQXG[0, 2] * ct * sp
    R2cs2r = CsaQXG[1, 2] * c2t * sp
    R2cs2i = CsaQXG[1, 2] * ct * cp
    R2cs3r = 0.25 * (CsaQXG[0, 0] - CsaQXG[1, 1]) * s2t * c2p
    R2cs3i = 0.5 * (CsaQXG[0, 0] - CsaQXG[1, 1]) * st * s2p
    R2cs4r = 0.5 * CsaQXG[0, 1] * s2t * s2p
    R2cs4i = CsaQXG[0, 1] * st * c2p
    R2p1cs = R2cs0 - R2cs1r + 1j * R2cs1i - R2cs2r - 1j * R2cs2i - R2cs3r + 1j * R2cs3i - R2cs4r - 1j * R2cs4i
    R2m1cs = -R2cs0 + R2cs1r + 1j * R2cs1i + R2cs2r - 1j * R2cs2i + R2cs3r + 1j * R2cs3i + R2cs4r - 1j * R2cs4i
    
    return [R2m1cs, R20cs, R2p1cs]

#**********************************Changes made to ACS due to modified equations***************************************
def AntiShift(AcsQXG, ct, st, cp, sp):
    R1p1acs = -1j*AcsQXG[0,1]*st + AcsQXG[0,2]*cp - 1j*AcsQXG[0,2]*ct*sp + AcsQXG[1,2]*sp + 1j*AcsQXG[1,2]*ct*cp
    R1m1acs = 1j*AcsQXG[0,1]*st + AcsQXG[0,2]*cp + 1j*AcsQXG[0,2]*ct*sp + AcsQXG[1,2]*sp - 1j*AcsQXG[1,2]*ct*cp
    
    return[R1m1acs, R1p1acs]

def fourier3(x, a, b, c):
    return a + b*np.cos(2*x*np.pi/180.) + c*np.sin(2*x*np.pi/180.)
def fourier5(x, a, b, c, d, e):
    return a + b*np.cos(2*x*np.pi/180.) + c*np.sin(2*x*np.pi/180.) + d*np.cos(4*x*np.pi/180.) + e*np.sin(4*x*np.pi/180)