import numpy as np
def geo(z,alpha_t,beta,rb,ra,rf,x,m,mt,n,al,i):
    '''Function to generate the tooth profile of involute gears according to \n
    LITVIN'''
    rhoF = 0.38*m
    ## VARIABLES ASSIGNMENT ===================================================
    z = z[i]
    rb = rb[i]
    ra = ra[i]
    rf = rf[i]
    x = x[i]
    ## INVOLUTE FUNCTION ======================================================
    def inv(alpha):
        return np.tan(alpha) - alpha
    ## TOOTH ROTATION =========================================================
    # alpha = np.pi/9
    r = mt*z/2
    # aPRIM = 0.5*m*(np.pi/2 + 2*x*np.tan(alpha))/r
    # aBASE = aPRIM + inv(alpha_t)
    aROOT = np.pi/z
    ## POINT A ================================================================
    # aA = ((np.pi - 5*np.tan(alpha_t))*mt/4 - rhoF*(1-np.sin(alpha_t))/np.cos(alpha_t))
    aA = np.pi*mt/4 - m*np.tan(alpha_t)-rhoF*np.cos(alpha_t) 
    bA = 1.25*m - x*m - rhoF  
    alphaG = np.arctan(np.tan(alpha_t) - 4*(m-x*m)/(mt*z*np.sin(2*alpha_t)))
    rG = r*np.cos(alpha_t)/(np.cos(alphaG))
    ## ROOT FILLET ============================================================
    thetaR = np.linspace(0,np.pi/2-alpha_t,n)
    phi = (aA-bA*np.tan(thetaR))/r
    xf = rhoF*np.sin(thetaR-phi) + aA*np.cos(phi) - bA*np.sin(phi) + r*(np.sin(phi) - phi*np.cos(phi))
    yf = -rhoF*np.cos(thetaR-phi) - aA*np.sin(phi) - bA*np.cos(phi) + r*(np.cos(phi) + phi*np.sin(phi)) 
    xF = xf*np.cos(aROOT) - yf*np.sin(aROOT)
    yF = xf*np.sin(aROOT) + yf*np.cos(aROOT)
    rG = np.sqrt(xF[-1]**2+yF[-1]**2)
    ## INVOLUTE ===============================================================
    thI = np.arccos(rb/rG)
    thF = np.arccos(rb/ra)
    theta = np.linspace(thI+inv(thI),thF+inv(thF),n)
    xe = rb*(np.sin(theta) - theta*np.cos(theta))
    ye = rb*(np.cos(theta) + theta*np.sin(theta))
    aBASE = -np.arctan(xF[-1]/yF[-1])+inv(thI)
    xI = xe*np.cos(aBASE) - ye*np.sin(aBASE)
    yI = xe*np.sin(aBASE) + ye*np.cos(aBASE)
    ## DEDDENDUM CIRCLE =======================================================
    xd = np.linspace(-rf*np.sin(aROOT), xF[0],n)
    yd = np.sqrt(rf**2-xd**2)
    ## ADDENDUM CIRCLE ========================================================
    xa = np.linspace(xI[-1], -xI[-1],n//5)
    ya = np.sqrt(ra**2-xa**2)
    ## CONCATENATE
    xGEO = np.concatenate((xd,xF,xI,xa,-xI[::-1],-xF[::-1],-xd[::-1]), axis=0)
    yGEO = np.concatenate((yd,yF,yI,ya, yI[::-1], yF[::-1], yd[::-1]), axis=0)
    xRF = np.concatenate((xd,xF), axis=0)
    yRF = np.concatenate((yd,yF), axis=0)
    return xGEO, yGEO, xRF, yRF, xI, yI

