# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 09:11:33 2022

@author: Sander
"""
import sys
sys.path.insert(1, 'Classes')
import OrcFxAPI as ofx
import numpy as np
import numpy.matlib
import os
import shutil
import pandas as pd
from openpyxl import Workbook, load_workbook
import pickle
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


class Plot:
    def __init__(self, model):
        self.model=model
        self.WavePeriods=model.periods
        self.BodyNames=model.BodyName
        self.Headings=model.headings
        self.CMap=plt.get_cmap("tab10")
        self.LineStyle='-'
        self.Color=self.CMap(0)
        self.Xlim=False
        
    def dispRAOs(self, BodySelection, DoF_ind, Dir_ind):
        DoF=['Surge', 'Sway', 'Heave','Roll', 'Pitch', 'Yaw']
        Unit=['m/m','m/m','m/m','deg/m','deg/m','deg/m']
        
        Cnt=DoF_ind
        plt.figure('Disp RAO '+ DoF[DoF_ind] +' ' + 'WvDr ' + str(self.headings[Dir_ind]))
        for i, Body in enumerate(BodySelection):
            if DoF_ind<3:
                plt.plot(self.WavePeriods, abs(self.model.displacementRAOs[Dir_ind,:,Cnt]), label=Body)
            else:
                plt.plot(self.WavePeriods, abs(self.model.displacementRAOs[Dir_ind,:,Cnt])*180/np.pi, label=Body)
            Cnt=Cnt+6
        
        plt.title('Disp RAO '+ DoF[DoF_ind] +' ' + 'WvDr ' + str(self.Headings[Dir_ind]))
        plt.xlabel('WavePeriod[s]')
        plt.ylabel(Unit[DoF_ind])
        plt.legend()
        plt.grid()
    
    def A_B_WavePeriod(self, Data, Dof1, Dof2, YLabel='AddedMass [t/m]', Label=''):
        plt.figure(str(Dof1)+str(Dof2) +' '+YLabel )
        plt.plot(self.WavePeriods, Data[:,Dof1,Dof2], label=Label, color=self.Color,linestyle=self.LineStyle)
        plt.xlabel('Motion Period[s]')
        plt.ylabel(YLabel)
        plt.title(str(Dof1)+str(Dof2) +' '+YLabel)
        if self.Xlim:
            plt.xlim(self.Xlim)
    
    def field_points(self, Dir_ind, Freq_ind):

        Title='FieldPoint RAO amp ' + str(self.model.periods[Freq_ind]) +'s ' + str(self.model.headings[Dir_ind]) + 'deg'
        plt.figure(Title)
        X=np.unique(self.model.FieldPointX)
        Y=np.unique(self.model.FieldPointY)
        Z_2D=np.zeros((len(Y),len(X)))
        c=0
        for ix, x in enumerate(X):
            for iy, y in enumerate(Y):
                Z_2D[iy,ix]=abs(self.model.fieldPointRAO[Dir_ind,Freq_ind,c])
                c=c+1 
                
        plt.contourf(X,Y,Z_2D, levels=np.arange(0.6,1.45,0.05))
        plt.colorbar()
        plt.title(Title)
        plt.axis('equal')
        plt.grid()
        
        