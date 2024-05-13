# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 14:07:06 2018
Spectral analysis tools.

By no stretch of the imagination is it a port of the matlab tools (yet)

For now it does have a collection of sea state spectra and some manipulation tools

It also contains some nice examples of passing function handles that i liked while programming this

Feel free to add to this toolbox


@author: Thijs van Krimpen

Functions
*********
"""
from numpy import sqrt, exp, log, pi, array, cos, arange, logical_and, sign
from scipy import trapz
import matplotlib.pyplot as plt
import copy
#spectrum function module comments to follow
#
# UPDATE:
# 1.0 author Thijs van Krimpen


def raise_value(error):
    #wrapper function for raising a value error
    raise ValueError(error)
    
def Tz2Tp(spectrumname):
    """
    Parameters
    ----------
    spectrumname: str
        allowable values PM, JONSWAP, CSIR, TOSETHAUGEN
    
    Returns
    -------
    
    Tz2Tp: lambda
    
    """
    #returns a function that can convert the Tz to Tp, requires a spectrum type
    Tz2Tp={}
    Tz2Tp['PM']=lambda Tz:Tz*1.407
    Tz2Tp['JONSWAP']=lambda Tz,gamma_val:Tz/(0.6673 + 0.05037 * gamma_val - 0.006230 * gamma_val**2 + 0.0003341 * gamma_val**3)
    Tz2Tp['CSIR']=lambda Tz: raise_value('CISR must be Tp') 
    Tz2Tp['TOSETHAUGEN']=lambda Tz: raise_value('TOSETHAUGEN must be Tp') 
    return Tz2Tp[spectrumname]
    
def gamma_value(Period,PeriodType,Hs) :   
    """
    
    Parameters
    ----------
        Period: float
            [s]
        PeriodType: str
            Tz or Tp
    
    Returns
    -------
    gamma: float
    returns an appropriate as per dnv-rp-c205 gamma value
    # input period, period type (tp/Tz) and Hs
    """
    gamma_valstart = 0
    gamma_val   = 9999
    
    # Loop to calculate Tz/Tp ratio and gamma_val
    # Loop is repeated until gamma_val does not change anymore
    
    while round(gamma_valstart*10000)/10000 != round(gamma_val*10000)/10000:
        gamma_valstart = gamma_val # gamma_valstart is results of previous loop
        if PeriodType == 'Tz':
            Tp=Tz2Tp('JONSWAP')(Period,gamma_valstart)
        else:
            Tp=Period
        
        if Tp/sqrt(Hs)<3.6:
            gamma_val = 5
        elif Tp/sqrt(Hs)>=3.6 and Tp/sqrt(Hs)<= 5:
            gamma_val = exp(5.75-1.15*Tp/sqrt(Hs))
        elif Tp/sqrt(Hs)>5:
            gamma_val = 1
        else:
            gamma_val = 3.3
        
    return gamma_val

#ToDo pep8
def JONSWAP(gamma_val,Period, Periodtype,w,Spec,Hs):  
    """
    creates a Jonswap spectrum definition
     Parameters:
     ------------
     gamma_val:
         peakedness factor, 
     Period:
         WavePeriod,
     Periodtype:
         Tz or Tp
     w:
         frequencies at which the density spectrum is assessed,
     spec:
         Spectrum dictionairy 
     Hs:
         significant wave height
    
    Returns
    -------
    spectrum dictionairy contianing fields:
    'density'= spectrum density assessed at frequencies
    """
    
    if Periodtype=='Tz':
        Tp=Tz2Tp('JONSWAP')(Period)
    else:
        Tp=Period
    if gamma_val==0:
        gamma_val=gamma_value(Period,Periodtype,Hs)
    
    omega=w
    
    omegap          = 2*pi/Tp # Definition of peak frequency
    sigma=copy.deepcopy(omega)
    sigma[omega<=omegap]= 0.07*omega[omega<=omegap]
    sigma[omega>omegap]=0.09*omega[omega>omegap]
    A = exp(-((omega/omegap-1)/(sigma*sqrt(2)))**2)
    Spec['density']    = ((1-0.287*log(gamma_val)) * 5/16 * Hs**2 * omegap**4 * omega**(-5)* exp(-5/4 * (omega/omegap)**(-4))* gamma_val**A)
    Spec['gamma']=gamma_val
    
def PM(gamma_val,Period, Periodtype,w,Spec,Hs):  
    """
    #PM spectrum
     Parameters:
     ------------
     gamma_val:
         peakedness factor, 
     Period:
         WavePeriod,
     Periodtype:
         Tz or Tp
     w:
         frequencies at which the density spectrum is assessed,
     spec:
         Spectrum dictionairy 
     Hs:
         significant wave height
    
    Returns
    -------
    spectrum dictionairy contianing fields:
    'density'= spectrum density assessed at frequencies
    """
    if Periodtype=='Tz':
        Tp=Tz2Tp('PM')(Period)
    else:
        Tp=Period
    omega=w
    T1           = 0.772*Tp #mean centroid wave periode
    Spec['density'] = ((173*Hs**2/T1**4)*omega**(-5)*exp(-692/T1**4*omega**(-4)))
#    gamma_val        = 0

def CSIR(gamma_val,Period, Periodtype,w,Spec,Hs):  
    """
    CSIR spectrum
     Parameters:
     ------------
     gamma_val:
         peakedness factor, 
     Period:
         WavePeriod,
     Periodtype:
         Tz or Tp
     w:
         frequencies at which the density spectrum is assessed,
     spec:
         Spectrum dictionairy 
     Hs:
         significant wave height
    
    Returns
    -------
    spectrum dictionairy contianing fields:
    'density'= spectrum density assessed at frequencies
    """
    
    if gamma_val > 5 or gamma_val < 3:
        print('gamma_val outside range for CSIR spectrum')
    if Periodtype=='Tz':
        Tp=Tz2Tp('CSIR')(Period)
    else:
        Tp=Period
    omega=w

    m0           = (Hs/4 )**2;
    omegap       = 2*pi/Tp;
    sigma=copy.deepcopy(omega)
    sigma[omega<=omegap]= 0.10*omega[omega<=omegap]
    sigma[omega>omegap]=0.14*omega[omega>omegap]
    q            = exp(-((omega/omegap-1)**2)/(2*sigma**2));
    density      = 13.5/(2*pi) * ((omega/omegap)**(-9/4))* exp(-0.3*(omega/omegap)**(-7))* gamma_val**q;
    
    # Force that m0=Hs**2/16 = int S(w)dw / normalize spectrum
    m0_actual    = trapz(density,omega);

    Spec['density'] = density* m0/m0_actual;
    

def TORSETHAUGEN(gamma_val,Period, Periodtype,w,Spec,Hs): 
    """
    #torsethaugen spectrum
     Parameters:
     ------------
     gamma_val:
         peakedness factor, 
     Period:
         WavePeriod,
     Periodtype:
         Tz or Tp
     w:
         frequencies at which the density spectrum is assessed,
     spec:
         Spectrum dictionairy 
     Hs:
         significant wave height
    
    Returns
    -------
    spectrum dictionairy contianing fields:
    'density'= spectrum density assessed at frequencies
    """
    if Periodtype=='Tz':
        Tp=Tz2Tp('TORSETHAUGEN')(Period)
    else:
        Tp=Period
    if Hs>11 or Hs>max((Tp/3.6)**2, (Tp-2)*12/11):
        print('**************************************************')
        print('Warning: Hs is outside the valid range')
        print('The validity of the spectral density is questionable')
        print('Hs: ' + str(Hs) +' Tp:' +str(Tp))
    
    if Tp>20 or Tp<3:
        print('**************************************************')
        print('Warning: Tp is outside the valid range')
        print('The validity of the spectral density is questionable')
        print('Hs: '+ str(Hs)+ ' Tp:'+str(Tp))
    
    omega=w
    
    omega_p=2*pi/Tp
    
    Af=6.6
    Ae=2.0
    Au=25
    A10=0.7
    A1=0.5
    Kg=35.0
    A20=0.6
    A2=0.3
    B1=2.0
    G0=3.26
    Tpf=Af*Hs**(1/3)
    Tl=Ae*Hs**0.5# lower limit
    Tu=Au# upper limit
    #Non-dimensional scales for the spectral peak period
    Epsi_l = min(max((Tpf-Tp)/(Tpf-Tl),0),1) #wind sea
    Epsi_u = min(max((Tp-Tpf)/(Tu-Tpf),0),1) #Swell
    
    if Tp<=Tpf:
#        se 1 # Wind dominated seas
        Rw=min(1,(1-A10)*exp(-(Epsi_l/A1)**2)+A10) # formulae(31)
        Hw1=Rw*Hs #significant waveheight wind
        
        Tpw1=Tp  # primary peak period
        Hw2=(1-Rw**2)**0.5*Hs # significant waveheight swell
        Tpw2=Tpf+B1 # Secondary peak period
        H1=Hw1
        H2=Hw2
        Tp1=Tpw1
        Tp2=Tpw2
        gamma_val=max(Kg*((2*pi/9.81)*Hw1/Tpw1**2)**(6/7),1)
    else:
#        case 0 # Swell dominated sea
        Rs=(1-A20)*exp(-(Epsi_u/A2)**2)+A20
        Hs1=Rs*Hs#significant waveheight swell
        Tps1=Tp # primary peak period
        Hs2= (1-Rs**2)**0.5*Hs
        Tps2=Af*Hs2**(1/3)
        H1=Hs1
        H2=Hs2
        Tp1=Tps1 
        Tp2=Tps2
        gamma_val=max((Kg*((2*pi/9.81)*Hs/(Tpf**2))**(6/7))*(1+6*Epsi_u),1)
    
    
    E1 = (1/16)*H1**2*Tp1
    E2=(1/16)*H2**2*Tp2
    f1n = w/2/pi*Tp1
    f2n=w/2/pi*Tp2
    Ar = (1+1.1*(log(gamma_val))**1.19)/gamma_val
    # if omega<omega_p  Sigma=0.07 else Sigma=0.09; 
    ind1 = w<omega_p;
    ind2 = not( ind1);
    Sigma=copy.deepcopy(omega)
    Sigma[ind1]=0.07 
    Sigma[ind2]=0.09;
    S_1n=G0*Ar*(f1n**(-4))*(exp(-(f1n)**(-4)))*gamma_val**(exp(-(1/2/Sigma**2)*(f1n-1)**2));
    S_2n=G0*(f2n**(-4))*exp(-(f2n**(-4)));
    
    Spectrum_f = E1*S_1n+E2*S_2n;
    Spec['density'] = Spectrum_f/2/pi;    


def Spectrum(type, Hs,  Periodtype, Period, w, gamma_val=0):
    """
    wrapper function for 4 different spectrum definitions
    Parameters:
    -----------
    type:
        spectrum name, 
    Hs:
        wave height, 
    Periodtype:
        period type (Tp/Tz), 
    period:
        [s]
    w:
        frequencies
    gamma_val:
        peakedness factor default = 0
        
    returns
    ---------
        a spectrum dictionary 
    """
    spec={'type':type}
    types=['JONSWAP','PM','CSIR','TORSETHAUGEN']
    
    if type not in types:
        raise ValueError('unrecognized spectrum type, type may be:'.join(type+',' for type in types))
    
    spec['frequencies']=w
    
    spectra={
            'JONSWAP':JONSWAP,
             'PM':PM,
            'CSIR':CSIR,
            'TORSETHAUGEN':TORSETHAUGEN
             }

    spectra[type](gamma_val,Period, Periodtype,w,spec,Hs) 
    return spec


def crappy_forward(fwd_v,Spectrum,direction):
    """
    #the author acknowledges that this is a bad way to treat encounter frequency
    #However it does treat the amount of energy that the vessel encounters at a 
    #particular frequency correctly. 
    #Note that there is not a particularly good way to do this short of running AQWA with forward speed maybe?
    """
    direction=direction*pi/180
    g=9.81
    omega=Spectrum['frequencies']
    r=Spectrum['density']
    dwdt0=g/(fwd_v*cos(direction))
    
    
    dwdt=1-2*omega*fwd_v*cos(direction)/g
    we=omega-omega**2/g*fwd_v*cos(direction)
    
    ind1=omega<dwdt0/2
    ind2=logical_and(omega>dwdt0/2, omega<dwdt0)
    ind3=omega>dwdt0
    w1=abs(we[ind1])
    w2=abs(we[ind2])
    w3=abs(we[ind3])
    r1=abs(r[ind1]/dwdt[ind1])
    r2=abs(r[ind2]/dwdt[ind2])
    r3=abs(r[ind3]/dwdt[ind3])
    
#    #include a cap on the assymptote
    r1[r1>5*max(Spectrum['density'])]=5*max(Spectrum['density'])
    r2[r2>5*max(Spectrum['density'])]=5*max(Spectrum['density'])
    r3[r3>5*max(Spectrum['density'])]=5*max(Spectrum['density'])
    
    
    spec={}
    spec['type']=Spectrum['type']+'_fwd'
    spec['fwd']=True
    spec['frequencies']=[w1,w2,w3]
    spec['density']=[r1,r2,r3]


    return spec
    

def less_crappy_forward(fwd_v,base_spec,direction):
    """
    #the author acknowledges that this is a bad way to treat encounter frequency
    #However it does treat the amount of energy that the vessel encounters at a 
    #particular frequency correctly. 
    #Note that there is not a particularly good way to do this short of running AQWA with forward speed maybe?
    """
    direction0=direction
    direction = direction*pi/180
    g = 9.81
    w = base_spec['frequencies']
    Sz = base_spec['density']
    
    #limit where waves are moving more slowly than the ship
    fwd_lim = g / (fwd_v*cos(direction))
    ind1 = w < fwd_lim/2
    ind2 = logical_and(w >= fwd_lim/2, w < fwd_lim)
    ind3 = w >= fwd_lim

    #encounter spectrum
    dwedw = 1 - 2*w*fwd_v * cos(direction) / g
    we = w - w**2/g * fwd_v * cos(direction)
    Sze = Sz / dwedw
    
    #include a cap on the assymptote
    cap = 5 * max(base_spec['density'])
    Sze[abs(Sze)>cap] = sign(Sze[abs(Sze)>cap])*cap
    
    #create wave spectrum components for OrcaFlex
    spec = {}
    spec['type'] = base_spec['type'] + '_fwd'
    spec['fwd'] = True
    spec['frequencies'] = [we[ind1], we[ind2], -we[ind3]]
    spec['density'] = [Sze[ind1], -Sze[ind2], -Sze[ind3]]
    spec['direction'] = [direction0, direction0, (180+direction0*180/pi) % 360]

    return spec
     
    
# #### test script to generate and plot spectra
#tp=9.0
#w=arange(0.5*2*pi/tp,10*2*pi/tp,(10*2*pi/tp-0.5*2*pi/tp)/1000)    
#spec=Spectrum(type='JONSWAP',Hs=1,Period=9,Periodtype='Tp',w=w,gamma_val=2.0)
#spec_v=crappy_forward(3,spec,0)
#
#
#plt.plot(spec['frequencies'],spec['density'])
#plt.xlim(0.5*2*pi/tp,10*2*pi/tp)
#
#plt.plot(spec_v['frequencies'][0],spec_v['density'][0])
#plt.plot(spec_v['frequencies'][1],spec_v['density'][1])
#plt.plot(spec_v['frequencies'][2],spec_v['density'][2])

#Spectrum(type='PM',Hs=1,Period=9,Periodtype='Tp',w=array([0.4,0.5]),gamma_val=2.0)
#Spectrum(type='CSIR',Hs=1,Period=9,Periodtype='Tp',w=array([0.4,0.5]),gamma_val=4.0)
#Spectrum(type='TORSETHAUGEN',Hs=1,Period=9,Periodtype='Tp',w=array([0.4,0.5]),gamma_val=2.0)
