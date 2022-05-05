###########################################################
#                    Hollenbach.py
###########################################################
# Alle input er i SI enheter
# Forklaring til noen utvalgte parametre:
# Vvec: Velocity vector [m/s]. Something like this:[9.8,10,11,12,13,13.5]
# L = Length between perpendiculars [m]
# Lwl = Length of waterline [m]
# Los = Length over Surface [m] (se kompendiet)
# B = Beam [m]
# TF = Dypgang ved FP [m]
# TA = Dypgang ved AP [m]
# CB = Block coefficient
# S = Wetted surface of hull
# Dp = Propelldiameter [m]
# Nrud = Antall ror [-]
# NBrac = Antall braketter [-]
# NBoss = Antall propellboss [-]
# NThr = Antall tunnelthrustere [-]

# Output er en matrise der første kollonne er hastighetsverdiene man puttet
# inn, og de neste kolonnene gir motstand og effektbehov for hver
# hastighet.

import numpy as np
import matplotlib.pyplot as plt

def hollenbach(Vsvec, L, Lwl, Los, B, TF, TA, CB, S, Dp, NRud, NBrac, NBoss, NThr):

    T = (TF + TA)/2
    Fi = (CB/L)*((B/2)*(TF + TA))**0.5
    k = 0.6*Fi + 145*Fi**3.5

    rho = 1025          # Density, seawater
    gravk = 9.81        # Gravity
    nu = 1.1395E-6      # Viscosity

    # Calculation of 'Froude length', Lfn:
    if Los/L < 1:
        Lfn = Los
    elif Los/L >= 1 and Los/L < 1.1:
        Lfn = L + 2/3*(Los - L)
    elif Los/L >= 1.1:
        Lfn = 1.0667*L
    
    ###########################################################
    # Constants from Hollenbachs paper:
    ###########################################################

    # 'Mean' resistance coefficients
    a = np.array([-0.3382, 0.8086, -6.0258, -3.5632, 9.4405, 0.0146, 0, 0, 0, 0])
    # a1 means a[0] and so on (Python arrays starts at 0)
    b = np.array([[-0.57424, 13.3893, 90.5960],
        [4.6614, -39.721, -351.483],
        [-1.14215, -12.3296, 459.254]])    # b12 means b[0,1]
    d = np.array([0.854, -1.228, 0.497])
    e = np.array([2.1701, -0.1602])
    f = np.array([0.17, 0.20, 0.60])
    g = np.array([0.642, -0.635, 0.150])

    # 'Minimum' resistance coefficients
    a_min = np.array([-0.3382, 0.8086, -6.0258, -3.5632, 0, 0, 0, 0, 0, 0])
    b_min = np.array([[-0.91424, 13.3893, 90.5960],
        [4.6614, -39.721, -351.483],
        [-1.14215, -12.3296, 459.254]])
    d_min = np.array([0, 0, 0])
    e_min = np.array([1, 0])
    f_min = np.array([0.17, 0.2, 0.6])
    g_min = np.array([0.614, -0.717, 0.261])
    
    # Preallocating vectors to later store results for plotting
    CFsvec = np.zeros(Vsvec.size)
    CRvec = np.zeros(Vsvec.size)
    C_Tsvec = np.zeros(Vsvec.size)
    R_T_meanvec = np.zeros(Vsvec.size)
    CR_minvec = np.zeros(Vsvec.size)
    C_Ts_minvec = np.zeros(Vsvec.size)
    R_T_minvec = np.zeros(Vsvec.size)
    P_E_meanvec = np.zeros(Vsvec.size)
    P_E_minvec = np.zeros(Vsvec.size)
        

    cc = 0
    # Loop over velocities
    for Vs in Vsvec:

        # Froude's number
        Fn = Vs/np.sqrt(gravk*Lfn)

        Fnkrit = np.dot(d, np.array([1, CB, CB**2]))
        c1 = Fn/Fnkrit
        c1_min = Fn/Fnkrit

        Rns = Vs*L/nu                       # Reynold's number for ship
        CFs = 0.075/(np.log10(Rns) - 2)**2  # ITTC friction line for ship

        # Calculation of C_R for given ship
        # Mean value

        CRFnkrit = np.max(np.array([[1.0],[(Fn/Fnkrit)**c1]]))

        kL = e[0]*L**e[1]

        # There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which
	      # is corrected in this formula by dividing by 10
        CRstandard = np.dot(np.array([1, CB, CB**2]), np.dot(b, np.array([1, Fn, Fn**2])/10))
        
        CR_hollenbach = CRstandard*CRFnkrit*kL*np.prod(np.array([T/B, B/L, Los/Lwl, Lwl/L, (1+(TA-TF)/L),
            Dp/TA, (1+NRud), (1+NBrac), (1+NBoss), (1+NThr)]**a))
        
        CR = CR_hollenbach*B*T/S    			# Resistance coefficient, scaled for wetted surface
        C_Ts = CFs + CR                         # Total resistance coeff. ship 
        R_T_mean = C_Ts*rho/2*Vs**2*S           # Total resistance to the ship

        ###########################################################
        # Minimum values

        # There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which
	      # is corrected in this formula by dividing by 10
        CRstandard_min = np.dot(np.array([1, CB, CB**2]), np.dot(b_min, np.array([1, Fn, Fn**2]))/10)

        CR_hollenbach_min = CRstandard_min*np.prod(np.array([T/B, B/L, Los/Lwl, Lwl/L, (1+(TA-TF)/L),
			  Dp/TA, (1+NRud), (1+NBrac), (1+NBoss), (1+NThr)]**a_min))
        
        CR_min = CR_hollenbach_min*B*T/S

        # Total resistance coefficient of the ship 
        C_Ts_min = CFs + CR_min			
        # Total resistance	
        R_T_min = C_Ts_min*rho/2*Vs**2*S

        # Propulsion power
        P_E_mean = R_T_mean*Vs      # [W]
        P_E_min = R_T_min*Vs        # [W]

        print('***********************************************')
        print(f'Vs = {Vs/0.5144:3.1f} knots')
        print('***********************************************')
        print('Mean values')
        print(f'CRh: {CR:0.3E}')
        print(f'CF: {CFs:0.3E}')
        print(f'CT: {C_Ts:0.3E}')
        print(f'RT: {R_T_mean:0.3E}')
        print('***********************************************')
        print('Minimum values')
        print(f'CRh: {CR_min:0.3E}')
        print(f'CF: {CFs:0.3E}')
        print(f'CT: {C_Ts_min:0.3E}')
        print(f'RT: {R_T_min:0.3E}')
        print('***********************************************')

        # Store results for plotting
        CFsvec[cc] = CFs
        CRvec[cc] = CR
        C_Tsvec[cc] = C_Ts
        R_T_meanvec[cc] = R_T_mean
        CR_minvec[cc] = CR_min
        C_Ts_minvec[cc] = C_Ts_min
        R_T_minvec[cc] = R_T_min
        P_E_meanvec[cc] = P_E_mean
        P_E_minvec[cc] = P_E_min

        cc += 1
 

    return np.array([Vsvec, R_T_meanvec, R_T_minvec, P_E_meanvec, P_E_minvec])
  


