# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:15:43 2022

@author: Sander
"""

#### Experimental validation functions and classes 

import collections
try:
    from collections import abc
    collections.MutableMapping = abc.MutableMapping
except:
    pass
import glob 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import OrcFxAPI as ofx
import pickle
import os
from scipy import signal
#### The most recent input class needs to be copied into the same main file 
import InputClasses

import h5py
import matplotlib.pyplot as plt


def read_marin_data(setup, run):
        
    # # fill in setup and run from excel sheet
    # setup = '002'
    # run = '009_01'
    location=r'C:\Users\Sander\PW2\Engineering - ENG\01_Hydromechanics\Project Management\07_ResearchPartners\01_MARIN\02_TKIWindWaveInteration\03_FTPDataMARIN\01_Raw Data'
    
    # test number
    test = '34110_01OB_01_' + setup + '_' + run
    
    # measurement file
    fileName =os.path.join(location, test[14:17], test + '.h5m')
    
    # inlezen van h5m file
    fileData = h5py.File(fileName, 'r')
    Data100Hz = fileData.get(list(fileData.keys())[0])
    Data200Hz = fileData.get(list(fileData.keys())[1])
    
    #### The 100Hz and 200Hz have differnt data sets 
    #### 100Hz Platform motions
    #### 200Hz Waves and mooring forces
    
    return Data100Hz, Data200Hz

def print_data_MARINOF(MARINBody,OFbody, iDOF):
    print('Ca: MARIN, OF')
    print(MARINBody.Ca[iDOF,iDOF].round(3),OFbody.Ca[iDOF,iDOF].round(3))
    print('LinearDamping: MARIN, OF')
    print(MARINBody.Bl[iDOF].round(3),OFbody.Bl[iDOF].round(3))
    print('QuadraticDamping: MARIN, OF')
    print(MARINBody.Bq[iDOF].round(3),OFbody.Bq[iDOF].round(3))
    print('Tn: MARIN, OF')
    print(MARINBody.TnDecay[iDOF].round(3),OFbody.TnDecay[iDOF].round(3))
    print('Total Mass: MARIN, OF')
    print(MARINBody.TotalMass[iDOF][iDOF].round(3),OFbody.TotalMass[iDOF][iDOF].round(3))
    print('q: MARIN, OF')
    print(MARINBody.Q[iDOF].round(3),OFbody.Q[iDOF].round(3))
    print('p: MARIN, OF')
    print(MARINBody.P[iDOF].round(3),OFbody.P[iDOF].round(3))
    
def print_data_2_comp(Body1,Body2, Tag1, Tag2, iDOF1, iDOF2):
    print('Ca: ' + Tag1 +', ' +Tag2)
    print(Body1.Ca[iDOF1,iDOF1].round(3),Body2.Ca[iDOF2,iDOF2].round(3))
    print('LinearDamping: ' + Tag1 +', ' +Tag2)
    print(Body1.Bl[iDOF1].round(3),Body2.Bl[iDOF2].round(3))
    print('QuadraticDamping: ' + Tag1 +', ' +Tag2)
    print(Body1.Bq[iDOF1].round(3),Body2.Bq[iDOF2].round(3))
    print('Tn: ' + Tag1 +', ' +Tag2)
    print(Body1.TnDecay[iDOF1].round(3),Body2.TnDecay[iDOF2].round(3))
    print('Total Mass: ' + Tag1 +', ' +Tag2)
    print(Body1.TotalMass[iDOF1][iDOF1].round(3),Body2.TotalMass[iDOF2][iDOF2].round(3))
    # print('q: ' + Tag1 +', ' +Tag2)
    # print(Body1.Q[iDOF1].round(3),Body2.Q[iDOF2].round(3))
    # print('p: ' + Tag1 +', ' +Tag2)
    # print(Body1.P[iDOF1].round(3),Body2.P[iDOF2].round(3))


def select_time_trace_from_MARIN_data(StTime,EndTime,scale, M100Hz):
    #### returns a new time array, motions 2d array (6dof, time) for the selected time and at the correct scale 
    #### Time is at the scale time
    signals = ['X TARGET','Y TARGET','Z TARGET','ROLL','PITCH','YAW']
    indst=(np.abs(M100Hz['Time'][:]*scale**0.5-StTime)).argmin()
    indend=(np.abs(M100Hz['Time'][:]*scale**0.5-EndTime)).argmin()
    Motions=np.zeros([indend-indst,6])
    for i in range(6):
        if i<3:
            Motions[:,i]=M100Hz[signals[i]][indst:indend]*scale
        else:
            Motions[:,i]=M100Hz[signals[i]][indst:indend]
    Time=M100Hz['Time'][indst:indend]*scale**0.5-StTime
    return Motions, Time 

def select_time_trace_from_OF(StTime, EndTime, Motions, RawData):
    indstOF=(np.abs(RawData['Environment']['time']-StTime)).argmin()
    indendOF=(np.abs(RawData['Environment']['time']-EndTime)).argmin()  
    MotionsOF=Motions[indstOF:indendOF,:]
    Time=RawData['Environment']['time'][indstOF:indendOF]-StTime
    return MotionsOF, Time

def phase_sync(ListOfMotions, ListOfTimes, SyncDOF, Distance=2):
    #### This function will phase sync to the first peak occuring of the SyncDOF selected to the first motion of the list
    #### Motions shouldbe in a 2d array standard form 
    #### Returns new timehistories in which motions are synced
    indPeaks=[]
    for iM, Mot in enumerate(ListOfMotions):
        dt=ListOfTimes[iM][1]-ListOfTimes[iM][0]
        ind, peakinfo=signal.find_peaks(Mot[:,SyncDOF], distance=Distance/dt)    
        indPeaks.append(ind)

    TimesOut=[]
    for iM, Time in enumerate(ListOfTimes):
        dTime=ListOfTimes[0][indPeaks[0][0]]-Time[indPeaks[iM][0]]
        TimesOut.append(Time+dTime)
    
    return TimesOut
        
def POI_translation(POI, TimeTrace):
    ## this can be used to calculate the velocity, acceleration or force at a point based on the 6dof local input
    ## velocity acceleration or force will be output at the POI in the local axis system
    ## This can also be used for a linearized POI translation best for small angles
    TranslationMatrix=np.matrix([[1,0,0,      0, POI[2], -POI[1]],
                                 [0,1,0,-POI[2],       0, POI[0]],
                                 [0,0,1, POI[1], -POI[0],     0]])
    TimeTrace[:,3:6]=TimeTrace[:,3:6]*np.pi/180
    POIOut=np.matrix(TimeTrace)*TranslationMatrix.transpose()
    POIOut=np.hstack([np.array(POIOut), TimeTrace[:,3:6]*180/np.pi])
    return POIOut   

    
class HydBody:
    
    def __init__(self):
        self.Mass=np.zeros([6,6])
        self.StiffnessHydStat=np.zeros([6,6])
        self.StiffnessMooring=np.zeros([6,6])
        # self.Stiffness=np.zeros([6,6])
        self.DOF=['Surge', 'Sway', 'Heave', 'Roll','Pitch','Yaw']
        self.Unit=['m','m','m','deg','deg','deg']
        self.TnDecay=np.zeros([6])
        self.WnDecay=np.zeros([6])
        self.TotalMass=np.zeros([6,6])
        self.AddedMass=np.zeros([6,6])
        self.Ca=np.zeros([6,6])
        self.Q=np.zeros([6])
        self.P=np.zeros([6])
        self.Bl=np.zeros([6])
        self.Bq=np.zeros([6])
        
    def set_mass(self,Mass):
        for i in range(3):
            self.Mass[i,i]=Mass
        # self.WnDecay=np.zeros([6])
        
    
    def analyze_decay(self,Time,Signal,Tag,Moored=False,axes=False, iDOF=2, Distance=3):
        #### Note DOF is in python indexing
        #### function makes a plot of the signal, peak finder and a scatter of the peak periods
        #### distance is the time inbetween periods used for filtering
        
        dt=Time[1]-Time[0]  
        indpeaks, peakinfo=signal.find_peaks(Signal, distance=Distance/dt)
        # indpeaksmin, peakinfomin=signal.find_peaks(Signal*-1, distance=Distance/dt)
        
        TimePeaks=Time[indpeaks]
        Peaks=Signal[indpeaks]
        dTimePeaks=[]
        dxperxm=[]
        xm=[]
        for i in range(len(TimePeaks)-1):
            dTimePeaks.append(TimePeaks[i+1]-TimePeaks[i])
            xm.append((Peaks[i]+Peaks[i+1])/2)
            dxperxm.append((Peaks[i]-Peaks[i+1])/xm[i])
            
        #### plot signal and natural period
        if axes:
            fig=axes[0]
            ax1=axes[1]
            ax1a=axes[2]
            ax2=axes[3]
        else:
            fig, (ax1,ax2)=plt.subplots(1,2)   
            ax1a=ax1.twinx() 
            
        ax1.plot(Time-Time[0],Signal, label=Tag)
        ax1.scatter(Time[indpeaks]-Time[0],Signal[indpeaks], marker='+')

        ax1.set_xlabel('Time[s]')
        ax1.set_ylabel('d'+self.DOF[iDOF]+ ' [' + self.Unit[iDOF] +']')   
        ax1.legend()
        #    
        ax1a.scatter(Time[indpeaks[:-1]]-Time[0 ],np.array(dTimePeaks),marker='x')
        ax1a.set_ylabel('Period between peaks [s]')
        
        #### plot linear regression
        ax2.scatter(xm,dxperxm, label='peaks')
        ax2.set_ylabel('delta peaks / mean,(n-(n+1))/mean(n,n+1) [m]')
        ax2.set_xlabel('mean peaks, (n+(n+1))/2) [m]')
        
        fig.set_size_inches([13,5])
        # plt.subplot_tool()
        # plt.show()
        plt.subplots_adjust(left=0.1,right=0.95,wspace=0.35)
        
        self.TnDecay[iDOF]=np.mean(dTimePeaks)
        self.WnDecay[iDOF]=2*np.pi/self.TnDecay[iDOF]
        if Moored:
            self.TotalMass[iDOF,iDOF]=(self.StiffnessHydStat[iDOF,iDOF]+self.StiffnessMooring[iDOF,iDOF])/self.WnDecay[iDOF]**2
        else:
            self.TotalMass[iDOF,iDOF]=self.StiffnessHydStat[iDOF,iDOF]/self.WnDecay[iDOF]**2
        self.AddedMass[iDOF,iDOF]=self.TotalMass[iDOF,iDOF]-self.Mass[iDOF,iDOF]
        self.Ca[iDOF,iDOF]=self.AddedMass[iDOF,iDOF]/self.Mass[iDOF,iDOF]
           
        Q, P=np.polyfit(xm,dxperxm,1)
        ax2.plot(xm, Q*np.array(xm)+P, label=str(P.round(3))+' + '+str(Q.round(3))+'*x')
        ax2.legend()
        self.Q[iDOF]=Q
        self.P[iDOF]=P
        
        self.Bl[iDOF]=2*P*(self.TotalMass[iDOF,iDOF])/self.TnDecay[iDOF]
        self.Bq[iDOF]=3/8*Q*(self.TotalMass[iDOF,iDOF])
        
                      
        return plt
        
    def calc_mass_and_stiffness(self, Sim, PlatformName):
        #### This function should calc the hydrostatic stiffness and the bodies inertia 
        #### This is only for simple ridig body models
        
        #### Make mass objects of the floater and platform
        Platform=self.model[PlatformName]
        
        #### Add mass objects to get a CoG and intertia properties around CoG
        
        #### Get floater diameter locations
        
        #### Calc hydrostatic stiffness
        
        #### Create a that includes the floater positions and seperate masses of the bodies
        

    
    
    # def POI_translation_NLEstimate(self, POI, TimeTrace):
    #     ## this can be used to calculate the velocity, acceleration or force at a point based on the 6dof local input
    #     ## velocity acceleration or force will be output at the POI in the local axis system
    #     ## This can also be used for a linearized POI translation
    #     TimeTrace1=np.zero()
    #     TimeTraxc
        
        
    #     TranslationMatrix=np.matrix([[1,0,0,      0, POI[2], -POI[1]],
    #                                   [0,1,0,-POI[2],       0, POI[0]],
    #                                   [0,0,1, POI[1], -POI[0],     0]])
    #     POIOut=np.matrix(TimeTrace)*TranslationMatrix.transpose()
    #     POIOut=np.array(POIOut)
    #     return POIOut   
        
