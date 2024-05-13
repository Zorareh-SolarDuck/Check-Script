# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:03:28 2021

@author: Sander
"""
import collections
try:
    from collections import abc
    collections.MutableMapping = abc.MutableMapping
except:
    pass
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
plt.rcParams['figure.max_open_warning'] = 1000                                              
from scipy.spatial.transform import Rotation as R
import re
import time
import copy

class GeneralFunctions:
    def filterdf(self, df, Headers, Selection, Type='Contains'):
        ## filter df based on a selection of headers, makes str contains or equals filter from   
        for i, GenS in enumerate(Headers): 
            if Type=='Contains':
                try:
                    df=df[df[GenS].str.contains(Selection[i])]
                except:
                    df=df[df[GenS]==Selection[i]]
            elif Type=='Exact':
                df=df[df[GenS]==Selection[i]]
        return df
    

class GetData:
    def __init__(self, model,Compose, SpecifiedTimePeriod=1):
        self.model=model
        self.SpecifiedTimePeriod=SpecifiedTimePeriod
        self.RawData={}
        self.translate_coordinates=Compose.translate_coordinates
        self.rotate_coordinates_around_z=Compose.rotate_coordinates_around_z
        
    def move_incomplete(self):
        dirname=os.path.dirname(self.model.latestFileName)
        filename=os.path.basename(self.model.latestFileName)
        incomplete=os.path.join(dirname,'Incomplete')
        
        if not os.path.exists(incomplete):
            os.makedirs(incomplete)
        
        shutil.move(self.model.latestFileName,os.path.join(incomplete,filename))
        shutil.move(self.model.latestFileName[0:-4]+'.yml', os.path.join(incomplete,filename[0:-4]+ '.yml'))
        print('moved ' + filename + ' to incomplete folder')
        
    def environment(self, N, LastTwoWaves=True):
            env=self.model.environment
            self.RawData['Environment']={}
            self.RawData['Environment']['Hex']=N
            self.RawData['Environment']['WaterDepth']=env.WaterDepth
            self.RawData['Environment']['WaveDirection']=env.WaveDirection
            self.RawData['Environment']['CurrentSpeed'] =env.CurrentSpeedAtSurface
            self.RawData['Environment']['CurrentDirection']=env.RefCurrentDirection
            self.RawData['Environment']['time']=self.model.general.TimeHistory('Time',self.SpecifiedTimePeriod)
            self.RawData['Environment']['WindSpeed']=env.WindSpeed
            self.RawData['Environment']['WindDirection']=env.WindDirection
            
            self.RawData['Environment']['WaveType']=env.WaveType
            if env.WaveType=='JONSWAP':
                self.RawData['Environment']['Tp']=env.WaveTp
                self.RawData['Environment']['Hs']=env.WaveHs
                if env.WaveNumberOfSpectralDirections>1:
                    self.RawData['Environment']['WaveSpreadingCoefficent']=env.WaveDirectionSpreadingExponent
                elif env.WaveNumberOfSpectralDirections==1:
                    self.RawData['Environment']['WaveSpreadingCoefficent']=None
            if env.WaveType=='Cnoidal' or env.WaveType=="Stokes' 5th" or env.WaveType=='Dean stream' or env.WaveType=='Airy':
                self.RawData['Environment']['WavePeriod']=env.WavePeriod
                self.RawData['Environment']['HMax']=round(env.WaveHeight,3)
                TEnd=self.RawData['Environment']['time'][-1]
                if LastTwoWaves:
                    self.SpecifiedTimePeriod=ofx.SpecifiedPeriod(TEnd-2*env.WavePeriod, TEnd)
                    self.RawData['Environment']['time']=self.model.general.TimeHistory('Time',self.SpecifiedTimePeriod)
   
    def get_timehistory_and_position(self, OFBuoy, Variable, POI, PlatPos, Phi):
                Dict={}
                OFPOI= ofx.oeBuoy(POI)
                Dict['TimeHistory']=OFBuoy.TimeHistory(Variable, self.SpecifiedTimePeriod, OFPOI)
                Dict['General']={}
                # Dict['General']['StaticZ']=OFBuoy.StaticResult('Z', OFPOI)
                Dict['General']['Local Pos']=POI
                Dict['General']['Global Pos']=PlatPos+self.rotate_coordinates_around_z(POI, Phi)
                return Dict
            
    def endA_endB_coordinates(self, Object):
        EndA=[Object.EndAX,Object.EndAY,Object.EndAZ]
        EndB=[Object.EndBX,Object.EndBY,Object.EndBZ]
        return EndA, EndB
    
    def initialXYZ(self, Object):
        XYZ=np.array([Object.InitialX,Object.InitialY,Object.InitialZ])
        return XYZ
    
    def global_to_local(self,i, GlobalPosPOI):
        ## Can only be run once platform positions in RawData are known
        ## i is the platform number
        Local1=GlobalPosPOI-self.RawData['Platform'][i]['General']['Initial Pos']
        return self.rotate_coordinates_around_z(Local1, self.RawData['Platform'][i]['General']['Initial Orientation'])
            
    def local_to_global(self, i, LocalPosPOI):
        ## Can only be run once platform positions in RawData are known
        ## i is the platform number
        Local=self.rotate_coordinates_around_z(LocalPosPOI, self.RawData['Platform'][i]['General']['Initial Orientation'])
        return Local+self.RawData['Platform'][i]['General']['Initial Pos']
    
    
    def platform_raw(self, Layout, Platform,MotionsFromStaticPosPOI=False,SkidCoG=False):
        print('AirGap Probes selected on edge length Input class')
        
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0:
                self.RawData['Platform']={}
                self.Split=True
                try:
                    OFPlatform=self.model['Platform_'+str(i)]
                    self.Split=False
                except:
                    pass
        if self.Split:
            self.RawData['Environment']['PlatformSplit']=True
        else:
            self.RawData['Environment']['PlatformSplit']=False
                
        
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0:  
                if self.Split:
                    for j in [0,120,240]:
                        OFPlatform=self.model['Platform_'+str(i)+'_'+str(j)]
                        self.get_platform_raw_main(OFPlatform, Layout, Platform, i,j,MotionsFromStaticPosPOI=MotionsFromStaticPosPOI,SkidCoG=SkidCoG)
                        # self.get_platform_raw_directional(OFPlatform, Platform, i,j)
                        self.RawData['Platform'][str(i)+'_'+str(j)]['General']['Split']=True
                else:
                    OFPlatform=self.model['Platform_'+str(i)]
                    self.get_platform_raw_main(OFPlatform, Layout, Platform, i,MotionsFromStaticPosPOI=MotionsFromStaticPosPOI,SkidCoG=SkidCoG)
                    # self.get_platform_raw_directional(OFPlatform, Platform, i)
                    self.RawData['Platform'][str(i)]['General']['Split']=False          
    
    def platform_raw_stability(self, Layout,pfnum, Platform,MotionsFromStaticPosPOI=False):
        print('AirGap Probes selected on edge length Input class')
        
        for i in range(len(Layout.Types)):
            i = pfnum
            self.RawData['Platform']={}
            OFPlatform=self.model['Platform_'+str(i)]
            self.Split=False
            
        if self.Split:
            self.RawData['Environment']['PlatformSplit']=True
        else:
            self.RawData['Environment']['PlatformSplit']=False
                
        
        for i in range(len(Layout.Types)):
            i = pfnum
            if Layout.Types[i] != 0:  
                if self.Split:
                    for j in [0,120,240]:
                        OFPlatform=self.model['Platform_'+str(i)+'_'+str(j)]
                        self.get_platform_raw_main(OFPlatform, Layout, Platform, i,j,MotionsFromStaticPosPOI=MotionsFromStaticPosPOI)
                        # self.get_platform_raw_directional(OFPlatform, Platform, i,j)
                        self.RawData['Platform'][str(i)+'_'+str(j)]['General']['Split']=True
                else:
                    OFPlatform=self.model['Platform_'+str(i)]
                    self.get_platform_raw_main(OFPlatform, Layout, Platform, i,MotionsFromStaticPosPOI=MotionsFromStaticPosPOI)
                    # self.get_platform_raw_directional(OFPlatform, Platform, i)
                    self.RawData['Platform'][str(i)]['General']['Split']=False   
                    
                    
    def get_platform_raw_main(self,OFPlatform, Layout, Platform, i, j=300,MotionsFromStaticPosPOI=False,SkidCoG=False):
        if j<300:
            Name=str(i)+'_'+str(j)
        else:
            Name=str(i)
        
        self.RawData['Platform'][Name]={}
        ## Extract General platform info
        self.RawData['Platform'][Name]['General']={}
        PlatInitialPos=[OFPlatform.InitialX, OFPlatform.InitialY, OFPlatform.InitialZ]
        self.RawData['Platform'][Name]['General']['Initial Pos']=PlatInitialPos
        Phi=OFPlatform.InitialRotation3
        self.RawData['Platform'][Name]['General']['Initial Orientation']=Phi
        self.RawData['Platform'][Name]['General']['PlatformType']=Layout.Types[i]
        self.RawData['Platform'][Name]['General']['Drawing']=np.asarray([OFPlatform.VertexX, OFPlatform.VertexY,OFPlatform.VertexZ])
        self.RawData['Platform'][Name]['General']['CoG']=[OFPlatform.CenterOfMassX, OFPlatform.CenterOfMassY, OFPlatform.CenterOfMassZ]
        
        GenDict={}
        # GenDict['Platform']=Name
        GenDict['Local Pos']=[0,0,0] ## Local position used as reference for the motions and accelrations
        GenDict['Global Pos']=PlatInitialPos
        GenDict['PlatformType']=Layout.Types[i]
        ## Extract Motions at local reference origin these need to at this reference to calculate realtive pitch
        # self.RawData['Platform'][Name]['General']['POI Motions']=[0,0,0]
        Variables='X', 'Y', 'Z', 'Rotation 1', 'Rotation 2', 'Rotation 3'
        self.RawData['Platform'][Name]['Positions']={}
        self.RawData['Platform'][Name]['Positions']['Origin']={}
        self.RawData['Platform'][Name]['Positions']['Origin']['TimeHistory']=OFPlatform.TimeHistory(Variables,self.SpecifiedTimePeriod)
        self.RawData['Platform'][Name]['Positions']['Origin']['General']=GenDict.copy()
        ## Extract Accelerations at local reference origin these acceleration do not contain the gravity component, this can extract also if desired.
        # self.RawData['Platform'][Name]['General']['POI Accelerations']=[0,0,0]
        Variables=['x acceleration rel. g', 'y acceleration rel. g', 'z acceleration rel. g', 'x angular acceleration', 'y angular acceleration','z angular acceleration']
        self.RawData['Platform'][Name]['Accelerations']={}
        self.RawData['Platform'][Name]['Accelerations']['Origin']={}
        self.RawData['Platform'][Name]['Accelerations']['Origin']['TimeHistory']=OFPlatform.TimeHistory(Variables,self.SpecifiedTimePeriod)
        self.RawData['Platform'][Name]['Accelerations']['Origin']['General']=GenDict.copy()
        if SkidCoG:
            self.RawData['Platform'][Name]['Accelerations']['SkidCoG']={}
            self.RawData['Platform'][Name]['Accelerations']['SkidCoG']['TimeHistory']=OFPlatform.TimeHistory(Variables,self.SpecifiedTimePeriod,ofx.oeBuoy(np.array(SkidCoG)))
            
            GenDictSkidCoG=GenDict.copy()
            GenDictSkidCoG['Local Pos']=SkidCoG
            GenDictSkidCoG['Global Pos']=np.array(PlatInitialPos)+self.rotate_coordinates_around_z(np.array(SkidCoG),np.cos(self.RawData['Platform'][Name]['General']['Initial Orientation']*np.pi/180))
            self.RawData['Platform'][Name]['Accelerations']['SkidCoG']['General']=GenDictSkidCoG    
            
        ## Extract Motions, these are platform motions from static position
        Variables=['Dynamic x','Dynamic y','Dynamic z', 'Dynamic Rx', 'Dynamic Ry','Dynamic Rz']
        self.RawData['Platform'][Name]['MotionsFromStaticPos']={}
        self.RawData['Platform'][Name]['MotionsFromStaticPos']['Origin']={}
        self.RawData['Platform'][Name]['MotionsFromStaticPos']['Origin']['TimeHistory']=OFPlatform.TimeHistory(Variables,self.SpecifiedTimePeriod)
        self.RawData['Platform'][Name]['MotionsFromStaticPos']['Origin']['General']=GenDict.copy()
        if MotionsFromStaticPosPOI:
            self.RawData['Platform'][Name]['MotionsFromStaticPos']['POI']={}
            self.RawData['Platform'][Name]['MotionsFromStaticPos']['POI']['TimeHistory']=OFPlatform.TimeHistory(Variables,self.SpecifiedTimePeriod,ofx.oeBuoy(np.array(MotionsFromStaticPosPOI)))
            GenDictPOI=GenDict.copy()
            GenDictPOI['Local Pos']=MotionsFromStaticPosPOI
            GenDictPOI['Global Pos']=np.array(PlatInitialPos)+self.rotate_coordinates_around_z(np.array(MotionsFromStaticPosPOI),np.cos(self.RawData['Platform'][Name]['General']['Initial Orientation']*np.pi/180))
            self.RawData['Platform'][Name]['MotionsFromStaticPos']['POI']['General']=GenDictPOI
            
        
        
    # def get_platform_raw_directional(self, OFPlatform, Platform, i,j):
        ### This function was set-up to split platforms
        Phi=self.RawData['Platform'][Name]['General']['Initial Orientation']
        PlatInitialPos=self.RawData['Platform'][Name]['General']['Initial Pos']
        
        ## Extract Normals for relative pitch motions
        for k in [0, 120, 240]:
            Normal=self.rotate_coordinates_around_z([1,0,0], k)   
            self.RawData['Platform'][Name]['Positions']['POI Normal '+str(k)]={}
            self.RawData['Platform'][Name]['Positions']['POI Normal '+str(k)]['General']={}
            self.RawData['Platform'][Name]['Positions']['POI Normal '+str(k)]['General']['Local Pos']=Normal
            POI=ofx.oeBuoy(Normal)
            Variables='X','Y','Z'
            self.RawData['Platform'][Name]['Positions']['POI Normal '+str(k)]['TimeHistory']=OFPlatform.TimeHistory(Variables,self.SpecifiedTimePeriod, POI)
            
        ## Extract Platform AirGaps
        self.RawData['Platform'][Name]['AirGap']={}
        
        for iAGP in range(len(Platform.AirGapProbesX)):   
            if j<300:
                POI=Platform.AirGapProbesX[iAGP], Platform.AirGapProbesY[iAGP], Platform.AirGapProbesZ[iAGP]
                POI0=self.translate_coordinates(Platform, POI, -j)
                if POI0[0]-0.1<=POI0[1]*np.tan(30*np.pi/180)+Platform.Length/3 and POI0[1]+0.1>=0:
                    self.RawData['Platform'][Name]['AirGap'][str(iAGP)]=self.get_timehistory_and_position(OFPlatform,'Sea surface clearance', POI, PlatInitialPos, Phi)           
            else:    
                POI=Platform.AirGapProbesX[iAGP], Platform.AirGapProbesY[iAGP], Platform.AirGapProbesZ[iAGP]
                self.RawData['Platform'][Name]['AirGap'][str(iAGP)]=self.get_timehistory_and_position(OFPlatform,'Sea surface clearance', POI, PlatInitialPos, Phi)

    def floater_raw(self, FloaterExtra=[''], MarineGrowthT=0):
        self.RawData['Floater']={}
        for i in self.RawData['Platform'].keys():
            PlatPos=self.RawData['Platform'][i]['General']['Initial Pos']
            PlatPhi=self.RawData['Platform'][i]['General']['Initial Orientation']
            for jj, ex in enumerate(FloaterExtra): ## This line has been added to account for extra floaters once used in a senstivity
                for ii,k  in enumerate([0, 120, 240]):
                    if self.RawData['Platform'][i]['General']['Split']:
                        Name='Floater'+str(i)+ ex
                    else:
                        Name='Floater'+str(i)+'_'+ str(k)+ ex
                    self.RawData['Floater'][Name]={}
                    OFFloater=self.model[Name]
                    POI=[OFFloater.InitialX,OFFloater.InitialY,OFFloater.InitialZ]
                    Dict={}
                    Dict['Local Pos']=POI  ## This is the connection location to the platform 
                    Dict['Global Pos']=PlatPos+self.rotate_coordinates_around_z(POI, PlatPhi) #### This is the connection location
                    # Dict['DiameterAboveWater']=OFFloater.CylinderOuterDiameter[0]
                    # Dict['DiameterBelowWater']=OFFloater.CylinderOuterDiameter[-1]
                    Dict['TopDiameter']=OFFloater.CylinderOuterDiameter[0]
                    Dict['ClashDiameter']=OFFloater.CylinderOuterDiameter[-2]
                    Dict['Local Orientation']=OFFloater.InitialRotation3
                    Dict['Global Orientation']=OFFloater.InitialRotation3+PlatPhi
                    Dict['MarineGrowthThicknes']=MarineGrowthT
                    Dict['Platform']=i
                    Dict['PlatformType']=self.RawData['Platform'][i]['General']['PlatformType']
                    self.RawData['Floater'][Name]['General']=Dict
                    ## Extract Connection forces
                    # ShortName=['Fx', 'Fy,' 'Fz', 'Mrx, ']
                    # for iV, V in enumerate(['Connection x force','Connection y force','Connection z force', 'Connection x moment', 'Connection y moment', 'Connection z moment']):
                    Variables=['Connection x force','Connection y force','Connection z force', 'Connection x moment', 'Connection y moment', 'Connection z moment']
                    ConnectionForces=OFFloater.TimeHistory(Variables,self.SpecifiedTimePeriod)
                    HorzForce=np.sqrt(ConnectionForces[:,0]**2+ConnectionForces[:,1]**2)
                    BendingMoment=np.sqrt(ConnectionForces[:,3]**2+ConnectionForces[:,4]**2)
                    self.RawData['Floater'][Name]['ConnectionForce']=np.hstack((ConnectionForces,HorzForce[:,None],BendingMoment[:,None]))
                    #### Extra Floater Airgap and clashing positions
                    # MarineGrowth=(Dict['DiameterBelowWater']-Dict['DiameterAboveWater'])/2
                    #### The section below has been coded for the taper without marine growth this can be manually adjusted the input
                    POI=[0,0,OFFloater.StackBaseCentreZ+MarineGrowthT]
                    self.RawData['Floater'][Name]['AirGap']=self.get_timehistory_and_position(OFFloater, 'Sea surface clearance',  POI, Dict['Global Pos'], 0)
                    Variables='X', 'Y', 'Z'
                    POI=[0,0,OFFloater.StackBaseCentreZ+OFFloater.CylinderLength[-1]+MarineGrowthT]               
                    self.RawData['Floater'][Name]['Position']=self.get_timehistory_and_position(OFFloater, Variables, POI, Dict['Global Pos'], 0)
                    Variables='Declination'
                    self.RawData['Floater'][Name]['Declination']=self.get_timehistory_and_position(OFFloater, Variables, POI, Dict['Global Pos'], 0)   
                    
    def constraints_raw(self, Coupling):
        self.RawData['Cons']={}
        for i in self.RawData['Platform'].keys():
            i0=i.split('_')[0]
            
            #### Check if the ref constraint is found and anaylze it. For future versions it will be better to loop trough the OrcaFlex objects, this is a general note for all the post processing. 
            try:
                Exists=True
                Con=self.model['Con_Ref_'+i]
            except:
                Exists=False
            if Exists:
                self.constraint_raw_details(i0, i, '0' , '1', 'Ref', Con,'Con_Ref_'+i)
            
            for Dir in ['0', '120', '240']:
                for Loc in ['1','-1']:
                    for Type in ['In', 'Out', 'Mid', 'Cen', 'Frld']:
                        Name='Con_' + Type +'_' + i0 +'_' + Dir+'_'+Loc                        
                        try:
                            Con=self.model[Name]
                            Exists=True
                        except:
                            Exists=False
                      
                        if Exists:
                            self.constraint_raw_details(i,i0,Dir,Loc, Type, Con, Name)
                            
    def constraints_raw_stability(self, Coupling,pfnum):
        self.RawData['Cons']={}
        for i in self.RawData['Platform'].keys():
            i0=i.split('_')[0]
            
            #### Check if the ref constraint is found and anaylze it. For future versions it will be better to loop trough the OrcaFlex objects, this is a general note for all the post processing. 
            try:
                Exists=True
                # Con=self.model['Con_Ref_'+i]
                Con = self.model['Constraint1']
            except:
                Exists=False
            if Exists:
                self.constraint_raw_details_stability(i0, i, None, Con,'Constraint1',pfnum)
            
            # for Dir in ['0', '120', '240']:
            #     for Loc in ['1','-1']:
            #         for Type in ['In', 'Out', 'Mid', 'Cen', 'Frld']:
            #             Name='Con_' + Type +'_' + i0 +'_' + Dir+'_'+Loc                        
            #             try:
            #                 Con=self.model[Name]
            #                 Exists=True
            #             except:
            #                 Exists=False
                      
            #             if Exists:
            #                 self.constraint_raw_details(i,i0,Dir,Loc, Type, Con, Name)
    def constraint_raw_details_stability(self,i,i0,Type, Con, Name,pfnum):
        #### This function has been added to be able to also process the reference constraint witht the same code
        self.RawData['Cons'][Name]={}
        General={}

        LocalP=[Con.InitialX, Con.InitialY, Con.InitialZ]
        if Type=='Ref':
            Connection=i
        else:
            # Connection=Con.Connection[9:]
            Connection = str(pfnum)
        
        GlobalP=self.rotate_coordinates_around_z(LocalP,self.RawData['Platform'][Connection]['General']['Initial Orientation'])+ self.RawData['Platform'][Connection]['General']['Initial Pos']        
        General['Global Pos']= GlobalP
        General['Local Pos']=LocalP
        General['Platform']=Connection
        # General['Orientation']=self.RawData['Platform'][Connection]['General']['Initial Orientation']##(np.mod(float(Connection),2)-1)*180
        General['Type']=Type
        General['PlatformType']=self.RawData['Platform'][i]['General']['PlatformType']
        self.RawData['Cons'][Name]['General']=General
        Variables=['x', 'y','z','Rx', 'Ry', 'Rz','In-frame connection Lx moment', 'In-frame connection Ly moment', 'In-frame connection Lz moment']
        self.RawData['Cons'][Name]['Motions']=Con.TimeHistory(Variables,self.SpecifiedTimePeriod)
                        
    def constraint_raw_details(self,i,i0, Dir,Loc,Type, Con, Name):
        #### This function has been added to be able to also process the reference constraint witht the same code
        self.RawData['Cons'][Name]={}
        General={}

        LocalP=[Con.InitialX, Con.InitialY, Con.InitialZ]
        if Type=='Ref':
            Connection=i
        else:
            Connection=Con.Connection[9:]
        
        GlobalP=self.rotate_coordinates_around_z(LocalP,self.RawData['Platform'][Connection]['General']['Initial Orientation'])+ self.RawData['Platform'][Connection]['General']['Initial Pos']        
        General['Global Pos']= GlobalP
        General['Local Pos']=LocalP
        General['Platform']=Connection
        # General['Orientation']=self.RawData['Platform'][Connection]['General']['Initial Orientation']##(np.mod(float(Connection),2)-1)*180
        General['Type']=Type
        General['PlatformType']=self.RawData['Platform'][i]['General']['PlatformType']
        self.RawData['Cons'][Name]['General']=General
        if Type=='Ref':
            Variables=['x', 'y','z','Rx', 'Ry', 'Rz']
            self.RawData['Cons'][Name]['Motions']=Con.TimeHistory(Variables,self.SpecifiedTimePeriod)
        else:    
            Variables=['In-frame connection Lx force', 'In-frame connection Ly Force', 'In-frame connection Lz Force', 'In-frame connection Lx moment', 'In-frame connection Ly moment', 'In-frame connection Lz moment']
            self.RawData['Cons'][Name]['ConnectionForce']=Con.TimeHistory(Variables,self.SpecifiedTimePeriod)
        
        ## Calculate the equivalent cen load in the local axis system of the platform for type 2 platforms
        if Type == 'Mid':
            Name1='Con_Cen_' + i0 +'_' + Dir+'_'+Loc
            self.RawData['Cons'][Name1]={}
            self.RawData['Cons'][Name1]['General']=General.copy()
            self.RawData['Cons'][Name1]['General']['Type']='Cen'
            ## Local position is establised at a later point in the coupling raw processing
            Forces1=self.RawData['Cons']['Con_In_' + i0 +'_' + Dir+'_'+Loc ]['ConnectionForce']
            Forces2=self.RawData['Cons']['Con_Mid_' + i0 +'_' + Dir+'_'+Loc ]['ConnectionForce']
            Forces3=self.RawData['Cons']['Con_Out_' + i0 +'_' + Dir+'_'+Loc ]['ConnectionForce']
            self.RawData['Cons'][Name1]['ConnectionForce']=Forces1+Forces2+Forces3
            
    def coupling_raw(self, Coupling):
        if self.Split:
            Coupling.split_connections()
            
        self.RawData['Coupling']={}
        self.RawData['Gap']={}
        self.RawData['CouplingClash']={}
        for i in  range(len(Coupling.Names)):
            Name=Coupling.Names[i]
            ## This checks if link is present and not left out for example for alt link confgurations
            try:
                Link=self.model[Name]
                Exists=True
            except:
                try:
                    ind=Name.find('_')
                    Name=Name[:ind] + Name[ind+1:]
                    Link=self.model[Name]
                    Exists=True
                except:
                    Exists=False
                    
            
            if Exists:
                
                self.RawData['Coupling'][Name]={}
                self.RawData['Coupling'][Name]['General']={}    
                self.RawData['Coupling'][Name]['General']['Type']= Coupling.Types[i]
                self.RawData['Coupling'][Name]['General']['Coupling Direction']=Coupling.Directions[i]
                                
                if Link.typeName=='Link':
                    ## Extract dynamic results            
                    Variables='Tension', 'Length', 'Velocity'            
                    for V in Variables:
                        self.RawData['Coupling'][Name][V]=Link.TimeHistory(V,self.SpecifiedTimePeriod)
                    self.RawData['Coupling'][Name]['General']['LinDamping']=Link.DampingConstant
                    self.RawData['Coupling'][Name]['General']['Stiffness']=Link.Stiffness
                    self.RawData['Coupling'][Name]['General']['UnstretchedLength']=Link.UnstretchedLength
                    if Link.UnstretchedLength==None:
                        Tension=np.array(Link.SpringTension)
                        Length=np.array(Link.SpringLength)
                        ind=np.argmin(Tension**2)
                        self.RawData['Coupling'][Name]['General']['UnstretchedLength']=Length[ind]
                        self.RawData['Coupling'][Name]['General']['Stiffness']='NonLin ~' +str(round((Tension[-1]-Tension[0])/(Length[-1]-Length[0])))
                elif Link.typeName=='Line':
                    #### This is only for the bar if it is a line at the moment
                    Variables='Effective Tension',  'y bend moment', 'Shear force'
                    for V in Variables:
                        # self.RawData['Beams'][NameA][V]=Beam.TimeHistory(V,self.SpecifiedTimePeriod, ofx.oeEndA)
                        self.RawData['Coupling'][Name][V]=Link.TimeHistory(V,self.SpecifiedTimePeriod, ofx.oeEndB)
                    self.RawData['Coupling'][Name]['General']['LinDamping']=0
                    self.RawData['Coupling'][Name]['General']['UnstretchedLength']=Link.Length[0]
                    LineType=Link.LineType[0]
                    self.RawData['Coupling'][Name]['General']['Stiffness']=self.model[LineType].EA/Link.Length[0]
                    
 
                LocalA=[Link.EndAX,Link.EndAY,Link.EndAZ]
                LocalB=[Link.EndBX,Link.EndBY,Link.EndBZ]
                
                if (LocalA - np.array([0,0,0])).all() < 1e-10:
                    LocalA = [0,0,0]
                    
                if LocalA==[0,0,0]:
                    EndAConec=Link.EndAConnection
                    EndBConec=Link.EndBConnection
                    EndBConec = EndBConec.replace('_Link_','_Cen_')
                    EndBConec = EndBConec.replace('_Bar_','_Cen_')
                    self.RawData['Coupling'][Name]['General']['Global Pos']= self.RawData['Cons'][EndAConec]['General']['Global Pos']
                    self.RawData['Coupling'][Name]['General']['Local Pos']=self.RawData['Cons'][EndAConec]['General']['Local Pos']
                    self.RawData['Coupling'][Name]['General']['Local Pos End B']=self.RawData['Cons'][EndBConec]['General']['Local Pos']
                    self.RawData['Coupling'][Name]['General']['Connection A']=EndAConec
                    self.RawData['Coupling'][Name]['General']['Connection B']=EndBConec
                    self.RawData['Coupling'][Name]['General']['Platform']= self.RawData['Cons'][EndAConec]['General']['Platform']
                    self.RawData['Coupling'][Name]['General']['Platform End B']= self.RawData['Cons'][EndBConec]['General']['Platform']
                    # self.RawData['Coupling'][Name]['General']['PlatformType']= self.RawData['Cons'][EndAConec]['General']['PlatformType']
                    GlobalB=self.RawData['Cons'][EndBConec]['General']['Global Pos']
                    ## This is to set the correct positions for the non-physical center constraint 
                    if  self.RawData['Coupling'][Name]['General']['Type']=='Mid':
                        # Name= EndAConn
                        seperator='_'
                        val=EndAConec.split('_')
                        val[1]='Cen'
                        Name1=seperator.join(val)
                        self.RawData['Cons'][Name1]['General']['Global Pos']=self.RawData['Cons'][EndBConec]['General']['Global Pos']
                        deltapos=self.RawData['Cons'][EndBConec]['General']['Global Pos']-self.RawData['Cons'][EndAConec]['General']['Global Pos']
                        self.RawData['Cons'][Name1]['General']['Local Pos']=self.RawData['Cons'][EndAConec]['General']['Local Pos']+deltapos*[-1,-1,0]
                    # self.RawData['Cons'][]
                else:         
                    EndAConec=Link.EndAConnection[9:]
                    EndBConec=Link.EndBConnection[9:]
                    GlobalA=self.rotate_coordinates_around_z(LocalA,self.RawData['Platform'][EndAConec]['General']['Initial Orientation'])+ self.RawData['Platform'][EndAConec]['General']['Initial Pos']        
                    GlobalB=self.rotate_coordinates_around_z(LocalB,self.RawData['Platform'][EndBConec]['General']['Initial Orientation'])+ self.RawData['Platform'][EndBConec]['General']['Initial Pos']        
                    self.RawData['Coupling'][Name]['General']['Global Pos']= GlobalA
                    self.RawData['Coupling'][Name]['General']['Local Pos']=LocalA
                    self.RawData['Coupling'][Name]['General']['Local Pos End B']=LocalB
                    self.RawData['Coupling'][Name]['General']['Platform']=EndAConec
                    self.RawData['Coupling'][Name]['General']['Platform End B']=EndBConec
                    self.RawData['Coupling'][Name]['General']['Connection A']=Link.EndAConnection
                    self.RawData['Coupling'][Name]['General']['Connection B']=Link.EndBConnection
                    # self.RawData['Coupling'][Name]['General']['Platform']=EndAConec                    
                
                #### This whole section has been commented to speed up the pp process
                if Name.find('Mid')==0:
                    self.RawData['Coupling'][Name]['General']['Global Pos']= GlobalB
                    # Extract data for clashing check of coupling
                    OFPlatformA=self.model[Coupling.EndAConnections[i]]
                    OFPlatformB=self.model[Coupling.EndBConnections[i]]   
                    self.RawData['CouplingClash'][Name[3:]]={}
                    self.RawData['CouplingClash'][Name[3:]]['General']=self.RawData['Coupling'][Name]['General'].copy()
                    Variables=['X', 'Y', 'Z']
                    # self.RawData['CouplingClash'][Name[3:]]['EndA']=OFPlatformA.TimeHistory(Variables,self.SpecifiedTimePeriod, ofx.oeBuoy(self.RawData['Coupling'][Name]['General']['Local Pos']))
                    # self.RawData['CouplingClash'][Name[3:]]['EndB']=OFPlatformB.TimeHistory(Variables,self.SpecifiedTimePeriod, ofx.oeBuoy(self.RawData['Coupling'][Name]['General']['Local Pos End B']))
                    
                    # Extract gap motions for the relative motions
                    # ga motion is now still based on position defined in the input, this can be adjusted
                    POIA=ofx.oeBuoy(Coupling.PosProbeEndAs[i])
                    POIB=ofx.oeBuoy(Coupling.PosProbeEndBs[i])
                    self.RawData['Gap'][Name[3:]]={}
                    self.RawData['Gap'][Name[3:]]['EndA']=OFPlatformA.TimeHistory(Variables,self.SpecifiedTimePeriod, POIA)
                    self.RawData['Gap'][Name[3:]]['EndB']=OFPlatformB.TimeHistory(Variables,self.SpecifiedTimePeriod, POIB)
                    self.RawData['Gap'][Name[3:]]['General']={}
                    self.RawData['Gap'][Name[3:]]['General']['Direction']=Coupling.Directions[i]
                    LocalA=Coupling.PosProbeEndAs[i]
                    LocalB=Coupling.PosProbeEndBs[i]
                    EndAConec=str(Coupling.EndAConnections[i])[9:]
                    EndBConec=str(Coupling.EndBConnections[i])[9:]
                    GlobalA=self.rotate_coordinates_around_z(LocalA,self.RawData['Platform'][EndAConec]['General']['Initial Orientation'])+ self.RawData['Platform'][EndAConec]['General']['Initial Pos']         
                    GlobalB=self.rotate_coordinates_around_z(LocalB,self.RawData['Platform'][EndBConec]['General']['Initial Orientation'])+ self.RawData['Platform'][EndBConec]['General']['Initial Pos']        
                    self.RawData['Gap'][Name[3:]]['General']['EndA Reference']=EndAConec
                    self.RawData['Gap'][Name[3:]]['General']['Local Pos End A']=LocalA
                    self.RawData['Gap'][Name[3:]]['General']['EndB Reference']=EndBConec
                    self.RawData['Gap'][Name[3:]]['General']['Local Pos End B']=LocalB
                    self.RawData['Gap'][Name[3:]]['General']['Global Plot Pos']=(GlobalA+GlobalB)/2  
    
   
    def mooring_raw(self, Mooring):
        ## Note that is only for a tether mooring system
        self.RawData['Mooring']={}
        for i, Plat in enumerate(Mooring.ConnectedPlatformsTotal):
            for il, l in enumerate(['-1','1']):
                for m in ['1','2']:
                    for di in [0,120,240]: 
                        try:
                            name = 'M'+ m +'_' + str(Plat) +'_' + str(di)
                            Mooring=self.model[name]
                            Exists1=True
                            try: 
                                a=Mooring.LineType
                                Type='Line'
                            except: 
                                Type='Tether'
                        except:
                            try:
                                name='M'+ m +'_' + str(Plat) +'_' + str(di)+ '_' + l
                                Mooring=self.model[name]
                                Exists1=True
                                try: 
                                    a=Mooring.LineType
                                    Type='Line'
                                except: 
                                    Type='Tether'
                            except:
                                Exists1=False
                        # try: 
                        #     name='M'+ m +'_' + str(Plat) +'_' + str(di)+ '_' + l
                        #     Mooring=self.model[name]
                        #     Exists1=True
                        #     try: 
                        #         a=Mooring.LineType
                        #         Type='Line'
                        #     except: 
                        #         Type='Tether'
                        # except:
                        #     Exists1=False
                        

                        if Exists1: 
                            if Type=='Tether':                 
                                # elif m=='2':
                                #     name='M'+ m +'_' + str(Plat) + '_' 
                                self.RawData['Mooring'][name]={}
                                self.RawData['Mooring'][name]['Tension']=Mooring.TimeHistory('Tension',self.SpecifiedTimePeriod)
                                General={}
                                General['Stiffness']=Mooring.Stiffness
                                Connec=Mooring.EndAConnection
                                if Connec.find('3DBuoy')== 0:
                                    Shackel=self.model[Connec]
                                    General['Global Pos']=[Shackel.InitialX, Shackel.InitialY, Shackel.InitialZ]
                                    General['Vlink']='Yes'
                                elif Connec.find('Platform')== 0:
                                    ConnecLocal=[Mooring.EndAX, Mooring.EndAY, Mooring.EndAZ]
                                    General['Global Pos']=self.rotate_coordinates_around_z(ConnecLocal,self.RawData['Platform'][Connec[9:]]['General']['Initial Orientation'])+ self.RawData['Platform'][Connec[9:]]['General']['Initial Pos'] 
                                    General['Vlink']='No'
                                elif Connec.find('Con_Frld')== 0:
                                    #### Sets the mooring location, please note that this based on the FrldLocation so the raw_cons functions should be run first
                                    General['Local Pos']=self.RawData['Cons'][Connec]['General']['Local Pos'].copy()
                                    General['Global Pos']=self.RawData['Cons'][Connec]['General']['Global Pos'].copy()
                                    General['Vlink']='No'
                                General['UnstrechtedLength']=Mooring.UnstretchedLength
                                General['Type']='Tether'
                                self.RawData['Mooring'][name]['General']=General
                                
                            elif Type=='Line':               
                                # elif m=='2':
                                #     name='M'+ m +'_' + str(Plat) + '_' 
                                self.RawData['Mooring'][name]={}
                                # TensionRange=Mooring.RangeGraph("Effective Tension", self.SpecifiedTimePeriod)
                                # TensionRange=TensionRange.Max
                                # self.RawData['Mooring'][name]['Tension']=np.max(TensionRange)
                                self.RawData['Mooring'][name]['Tension']=Mooring.TimeHistory('Effective Tension',self.SpecifiedTimePeriod, ofx.oeEndA)
                                
                                # print(Mooring.Representation, Mooring.CumulativeLength[-1])
                                if  Mooring.Representation != 'Analytic catenary':
                                    NormTensionRange=Mooring.RangeGraph("Normalised Tension", self.SpecifiedTimePeriod) #, ofx.arSpecifiedArclengths(0, Mooring.CumulativeLength[-1]))                                    
                                    NormTensionRange = pd.DataFrame(NormTensionRange.Max)
                                    if m == '1':
                                        NormTensionRange.drop(NormTensionRange.tail(2).index, inplace=True) # to drop anchor UC, currently problem with number of nodes in Orcaflex
                                    elif m == '2': # some problem in number of nodes of M2 lines connected to buoy
                                        NormTensionRange.drop(NormTensionRange.tail(1).index, inplace=True) 
                                    tmp_max_node = NormTensionRange[0].idxmax()+1   
                                    # print(name, tmp_max_node)
                                    self.RawData['Mooring'][name]['Normalized_Tension']=Mooring.TimeHistory('Normalised Tension',self.SpecifiedTimePeriod, ofx.oeNodeNum(tmp_max_node)) #ofx.oeEndA) 

                                #### Section below has been added for specific line outputs in the mooring
                                if m == '1':
                                    self.RawData['Mooring'][name]['TensionBottom']=Mooring.TimeHistory('Effective Tension',self.SpecifiedTimePeriod, ofx.oeEndB)
                                    self.RawData['Mooring'][name]['UpliftAngle']=Mooring.TimeHistory('Declination',self.SpecifiedTimePeriod, ofx.oeEndB)-90
                                    self.RawData['Mooring'][name]['UpliftForce']=Mooring.TimeHistory('End GZ force',self.SpecifiedTimePeriod, ofx.oeEndB)*-1
                                    self.RawData['Mooring'][name]['TouchDownPoint']=np.sum(Mooring.Length)-Mooring.TimeHistory('Arc length',self.SpecifiedTimePeriod, ofx.oeTouchdown)
                                    #### Watch out the location of the line end has been based on the fact there is only a single end chain
                                    # LineLocation=ofx.oeArcLength(np.sum(Mooring.Length[0:2]))
                                    # self.RawData['Mooring'][name]['SeaBedClearanceSynLineEnd']=Mooring.TimeHistory('Vertical seabed clearance',self.SpecifiedTimePeriod, LineLocation)
                                    if  Mooring.Representation != 'Analytic catenary':
                                        self.RawData['Mooring'][name]['Anchor_Normalized_Tension']=Mooring.TimeHistory('Normalised Tension',self.SpecifiedTimePeriod, ofx.oeEndB)                                  
                                       
                                        MooringLineType = {} 
                                        MooringSeaBedClearance = {}
                                        for i,linetype in enumerate(Mooring.LineType):    
                                            if linetype in MooringLineType:
                                                linetype += '_1'
                                            MooringLineType[linetype] = Mooring.Length[i]
                                            if i == 0:
                                                arcLengthRange = [0,Mooring.Length[i]]
                                            else:
                                                arcLengthRange = [np.sum(Mooring.Length[:i]),np.sum(Mooring.Length[:i+1])]
                                            Data = Mooring.RangeGraph('Vertical seabed clearance',self.SpecifiedTimePeriod, arclengthRange=ofx.arSpecifiedArclengths(arcLengthRange[0], arcLengthRange[-1]))
                                            MooringSeaBedClearance[linetype] = pd.concat([pd.DataFrame(Data.X),pd.DataFrame(Data.Max),pd.DataFrame(Data.Mean),pd.DataFrame(Data.Min),pd.DataFrame(Data.StdDev)],axis=1)
                                            MooringSeaBedClearance[linetype].columns = ['X','Max','Mean','Min','StdDev']
                                        self.RawData['Mooring'][name]['SeaBedClearance']=MooringSeaBedClearance   
                                General={}
                                # MoorLine={}
                                #### Calc lines equivalent stiffness
                                k_inv=[]
                                try: 
                                    for i,LineType in enumerate(Mooring.LineType):
                                        LineType=self.model[LineType]      
                                        k_inv.append(1/(LineType.EA/Mooring.Length[i]))
                                        k=1/(np.sum(k_inv))
                                        General['Stiffness']= k*np.sum(Mooring.Length)

                                except:
                                    General['Stiffness']='NonLinear'
                                #### End of calc mooring stiffness
                                # Connec=Mooring.EndAConnection
                                if m == '1':
                                    Connec=Mooring.EndAConnection
                                elif m == '2':
                                    Connec=Mooring.EndBConnection
                                # print(Connec)
                                if Connec.find('Con_Frld') == 0:
                                    #### Sets the mooring location, please note that this based on the FrldLocation so the raw_cons functions should be run first
                                    General['Local Pos']=self.RawData['Cons'][Connec]['General']['Local Pos'].copy()
                                    General['Global Pos']=self.RawData['Cons'][Connec]['General']['Global Pos'].copy()
                                    General['Vlink']='No'
                                
                                elif Connec.find('Platform') == 0:
                                    ConnecLocal=[Mooring.EndAX, Mooring.EndAY, Mooring.EndAZ]
                                    General['Global Pos']=self.rotate_coordinates_around_z(ConnecLocal,self.RawData['Platform'][Connec[9:]]['General']['Initial Orientation'])+ self.RawData['Platform'][Connec[9:]]['General']['Initial Pos'] 
                                    General['Vlink']='No'
                                    
                                elif Connec.find('VLinkBuoy') == 0: ### VLinkBuoy exists
                                    Shackel=self.model[Connec]
                                    General['Global Pos']=[Shackel.InitialX, Shackel.InitialY, Shackel.InitialZ]
                                    General['Vlink']='Yes'

                                General['UnstrechtedLength']=np.sum(Mooring.Length)
                                General['Type']='Line'
                                self.RawData['Mooring'][name]['General']=General
                                # #### End of calc mooring stiffness
                                # Connec=Mooring.EndAConnection
                                # # if Connec.find('3DBuoy')== 0:
                                # if Connec.find('VLinkBouy') != -1:
                                #     Shackel=self.model[Connec]
                                #     General['Global Pos']=[Shackel.InitialX, Shackel.InitialY, Shackel.InitialZ]
                                #     General['Vlink']='Yes'
                                # elif Connec.find('Platform')== 0:
                                #     ConnecLocal=[Mooring.EndAX, Mooring.EndAY, Mooring.EndAZ]
                                #     General['Global Pos']=self.rotate_coordinates_around_z(ConnecLocal,self.RawData['Platform'][Connec[9:]]['General']['Initial Orientation'])+ self.RawData['Platform'][Connec[9:]]['General']['Initial Pos'] 
                                #     General['Vlink']='No'
                                # elif Connec.find('Con_Frld')== 0:
                                #     #### Sets the mooring location, please note that this based on the FrldLocation so the raw_cons functions should be run first
                                #     General['Local Pos']=self.RawData['Cons'][Connec]['General']['Local Pos'].copy()
                                #     General['Global Pos']=self.RawData['Cons'][Connec]['General']['Global Pos'].copy()
                                #     General['Vlink']='No'
                                # elif m in ['1']: #Connec.find('Buoy')== 0 or Connec.find('ConAft')== 0 :
                                #     #### set for the towig line location (just for later data extract with the same headers)
                                #     General['Local Pos']=[0,0,0]
                                #     General['Global Pos']=[0,0,0]
                                #     General['Vlink']='No'
                                # General['UnstrechtedLength']=np.sum(Mooring.Length)
                                # General['Type']='Line'
                                # self.RawData['Mooring'][name]['General']=General
                                
               
    def beams_raw(self,Platform=None):
        self.RawData['Beams']={}
        for i in self.RawData['Platform'].keys():
            i0=i.split('_')[0]
            for Dir in ['0', '120', '240']:
                for n in ['0', '1', '2']:
                    try:
                        Name='Beam_' + i0 +'_' + Dir +'_'+ n 
                        NameA= Name+'_EndA'                       
                        NameB= Name +'_EndB' 
                        NameCen= Name+'_Center'
                        Beam=self.model[Name]
                        NameMember = ['Member1','Member2','Member3','Member4']
                        # NameMemberA = [NameA + mem for mem in NameMember]
                        # NameMemberB = [NameB + mem for mem in NameMember]
                    except:
                        pass
                    GenA={}
                    GenA['LineType']=Beam.LineType[0]
                    GenA['Platform']=str(i)
                    GenB=GenA.copy()
                    GenC=GenA.copy()
                    GenA['Connection']=Beam.EndAConnection
                    GenB['Connection']=Beam.EndBConnection
                    GenC['Connection']=Name+'_Center'
                    [LocalA, LocalB] =self.endA_endB_coordinates(Beam)
                    if Beam.EndAConnection.split('_')[0]=='Node':       
                        GlobalA0=self.initialXYZ(self.model[GenA['Connection']])
                        ## Local pos is the actual position of the connection
                        GenA['Local Pos']=self.global_to_local(i,GlobalA0)
                    else: 
                        GlobalA0=self.local_to_global(i, LocalA)
                        GenA['Local Pos']=LocalA
                    if Beam.EndBConnection.split('_')[0]=='Node': 
                        GlobalB0=self.initialXYZ(self.model[GenB['Connection']])
                        GenB['Local Pos']=self.global_to_local(i,GlobalB0)                      
                    else: 
                        GlobalB0=self.local_to_global(i, LocalB)
                        GenB['Local Pos']=LocalB
                    ## Global pos is altered for plotting purpose
                    GenA['Global Pos']=np.mean([np.mean([GlobalA0,GlobalB0],axis=0),GlobalA0], axis=0)
                    GenB['Global Pos']=np.mean([np.mean([GlobalA0,GlobalB0], axis=0),GlobalB0], axis=0)
                    GenC['Local Pos']=np.mean([GenA['Local Pos'],GenB['Local Pos']],axis=0)
                    GenC['Global Pos']=np.mean([GenA['Global Pos'],GenB['Global Pos']], axis=0)
                                       
                    self.RawData['Beams'][NameA]={}
                    self.RawData['Beams'][NameB]={}
                    self.RawData['Beams'][NameCen]={}
                    self.RawData['Beams'][NameA]['General']=GenA
                    self.RawData['Beams'][NameB]['General']=GenB
                    self.RawData['Beams'][NameCen]['General']=GenC
                    
                    Variables='Effective Tension', 'x bend moment', 'y bend moment', 'Shear force', 'Declination'

                    for V in Variables:
                        self.RawData['Beams'][NameA][V]=Beam.TimeHistory(V,self.SpecifiedTimePeriod, ofx.oeEndA)
                        self.RawData['Beams'][NameB][V]=Beam.TimeHistory(V,self.SpecifiedTimePeriod, ofx.oeEndB)
                        
                    # for im, memb in enumerate(NameMember):
                    self.RawData['Beams'][NameA][NameMember[0] + ' Force'] = self.RawData['Beams'][NameA]['Effective Tension']/4 + self.RawData['Beams'][NameA]['x bend moment']/Platform.TrussWidth + self.RawData['Beams'][NameA]['y bend moment']/Platform.TrussHeight
                    self.RawData['Beams'][NameA][NameMember[1] + ' Force'] = self.RawData['Beams'][NameA]['Effective Tension']/4 - self.RawData['Beams'][NameA]['x bend moment']/Platform.TrussWidth + self.RawData['Beams'][NameA]['y bend moment']/Platform.TrussHeight
                    self.RawData['Beams'][NameA][NameMember[2] + ' Force'] = self.RawData['Beams'][NameA]['Effective Tension']/4 - self.RawData['Beams'][NameA]['x bend moment']/Platform.TrussWidth - self.RawData['Beams'][NameA]['y bend moment']/Platform.TrussHeight
                    self.RawData['Beams'][NameA][NameMember[3] + ' Force'] = self.RawData['Beams'][NameA]['Effective Tension']/4 + self.RawData['Beams'][NameA]['x bend moment']/Platform.TrussWidth - self.RawData['Beams'][NameA]['y bend moment']/Platform.TrussHeight
                    
                    
                    
                    Variables= ['Sea surface clearance'] #, 'x relative velocity', 'y relative velocity'
                    #### should add node 2 somewhere
                    for V in Variables:
                        self.RawData['Beams'][NameCen][V]=Beam.TimeHistory(V,self.SpecifiedTimePeriod, ofx.oeNodeNum(2))
                    
                    # tmp=self.RawData['Beams'][NameCen]['x relative velocity'].copy()
                    # tmp[self.RawData['Beams'][NameCen]['Sea surface clearance']<0]=0
                    # self.RawData['Beams'][NameCen]['x relative velocity air']=tmp
                    
                    # tmp=self.RawData['Beams'][NameCen]['y relative velocity'].copy()
                    # tmp[self.RawData['Beams'][NameCen]['Sea surface clearance']<0]=0
                    # self.RawData['Beams'][NameCen]['y relative velocity air']=tmp
                    
                    # tmp=self.RawData['Beams'][NameCen]['x relative velocity'].copy()
                    # tmp[self.RawData['Beams'][NameCen]['Sea surface clearance']>0]=0
                    # self.RawData['Beams'][NameCen]['x relative velocity water']=tmp
                    
                    # tmp=self.RawData['Beams'][NameCen]['y relative velocity'].copy()
                    # tmp[self.RawData['Beams'][NameCen]['Sea surface clearance']>0]=0
                    # self.RawData['Beams'][NameCen]['y relative velocity water']=tmp
                                            
class ProcessData:
    def __init__(self, RawData, Compose, ExtractionInstance=None):
        self.ProData={}## just as preliminary self property
        self.StatisticsDuration=3 #hours
        self.ExtractionInstance=ExtractionInstance ## Seconds
        self.RawData=RawData
        self.rotate_coordinates_around_z=Compose.rotate_coordinates_around_z
        self.translate_coordinates=Compose.translate_coordinates
        # self.DurationSimulation=RawData['Environment']['time'][-1]
        self.IncFEMLabel=True
        
    def environment(self):
        self.ProData['Environment']=self.RawData['Environment'].copy()
        self.Time=self.ProData['Environment']['time']
        if self.ExtractionInstance!= None:    
            self.ExtractionIndex=np.where(np.array(self.RawData['Environment']['time'])==self.ExtractionInstance)
            self.ProData['Environment']['ExtractionTime']=self.ExtractionInstance
            # self.ProData['Environment']['ExtractionIndex']=self.ExtractionIndex
            
        del self.ProData['Environment']['time']
        
        try:                                        
            self.ProData['Environment']['Name']= ("Wdep%.2f" % (self.ProData['Environment']['WaterDepth']) +
                                                  "m Hs%.2f" % (self.ProData['Environment']['Hs']) +
                                                  "m Tp%.2f" % (self.ProData['Environment']['Tp']) + 
                                                  "s Wdir%.2f" % (self.ProData['Environment']['WaveDirection']) +
                                                  "deg Cv%.2f" % (self.ProData['Environment']['CurrentSpeed']) +
                                                  "m/s Cdir%.2f" % (self.ProData['Environment']['CurrentDirection'])+
                                                  "deg")
                                                    # "_Wndir%.2f" % (0))
        except:                                     
            self.ProData['Environment']['Name']= ("Wdep%.2f" % (self.ProData['Environment']['WaterDepth']) +
                                                  "m HMax%.2f" % (self.ProData['Environment']['HMax']) +
                                                  "m T%.2f" % (self.ProData['Environment']['WavePeriod']) + 
                                                  "s Wdir%.2f" % (self.ProData['Environment']['WaveDirection']) +
                                                  "deg Cv%.2f" % (self.ProData['Environment']['CurrentSpeed']) +
                                                  "m/s Cdir%.2f" % (self.ProData['Environment']['CurrentDirection'])+
                                                  "deg")                               
            
    def find_down_crossings(self, TimeTrace, Mean):
        data=TimeTrace-Mean
        pos = data > 0
        return (pos[:-1] & ~pos[1:]).nonzero()[0]
    
    def POI_translation(self, POI, TimeTrace):
        ## this can be used to calculate the velocity, acceleration or force at a point based on the 6dof local input
        ## velocity acceleration or force will be output at the POI in the local axis system
        ## This function is called up so many times that it significnatly slows the process 
        # a=time.perf_counter() 
        TranslationMatrix=np.array([[1,0,0,      0, POI[2], -POI[1]],
                                     [0,1,0,-POI[2],       0, POI[0]],
                                     [0,0,1, POI[1], -POI[0],     0]])
        # b=time.perf_counter()
        # POIOut=np.array(TimeTrace)*TranslationMatrix.transpose()
        # POIOut=np.matmul(np.array(TimeTrace),TranslationMatrix)
        # c=time.perf_counter()
        # # # POIOut=np.matmul(np.array(TimeTrace),TranslationMatrix)
        # print(b-a)
        # print(c-b)
        return TimeTrace@TranslationMatrix.transpose()
    
        
    # def calc_rayleigh_max(self, TimeTrace):
    #     Mean=np.mean(TimeTrace,0)
    #     CrossInd, _ = self.find_down_crossings(TimeTrace, Mean)
    #     CrossingPeriod=self.DurationSimulation/len(CrossInd)
    #     n=self.StatisticsDuration/CrossingPeriod
    #     Max=Dict['Mean']+Dict['Sig Amplitude']*np.sqrt(np.log(n)/2)
    #     return Max
        
    def POI_translation_position(self, i, POI, Tag):
        ## This can be used to recalculate position on a vessel based on orcaflex conventions and the motions of the origin
        ## In this manner it is possible to calc positions soley on recalculation and no need to exract everything directly from the simulation
        self.RawData['Platform'][i]['Positions'][Tag]={}
        self.RawData['Platform'][i]['Positions'][Tag]['General']=self.RawData['Platform'][i]['Positions']['Origin']['General'].copy()
        self.RawData['Platform'][i]['Positions'][Tag]['General']['Local Pos']=POI
        self.RawData['Platform'][i]['Positions'][Tag]['General']['Global Pos']=np.array(self.RawData['Platform'][i]['Positions']['Origin']['General']['Global Pos'])+self.rotate_coordinates_around_z(POI, self.RawData['Platform'][i]['General']['Initial Orientation'])
        TimetraceO=self.RawData['Platform'][i]['Positions']['Origin']['TimeHistory']
        r = R.from_euler('zyx', np.flip(TimetraceO[:,3:6], axis=1), degrees=True)
        self.RawData['Platform'][i]['Positions'][Tag]['TimeHistory']=r.apply(POI)+TimetraceO[:,:3]
    
    
    def statistics(self, TimeTrace, Distance=False):
        Dict={}
        Dict['Max']=np.max(TimeTrace,0)
        Dict['Min']=np.min(TimeTrace,0)
        Dict['StDev']=np.std(TimeTrace,0)
        Dict['Mean']=np.mean(TimeTrace,0)
        Dict['Sig Amplitude']= 2*Dict['StDev']
        Dict['AbsMax']=np.max([Dict['Max'] ,np.abs(Dict['Min'])],0)
        if Distance:
            Dict['Distance']=np.sum(abs(TimeTrace[1:]-TimeTrace[:-1]))
        if self.ExtractionInstance!=None:
            Dict['InstantVal']=TimeTrace[self.ExtractionIndex][0]
        # Dict['RayleighMax'+str(self.StatisticsDuration)+'h']=                
        return Dict     
        
    def constraints(self):
        self.ProData['Cons']={}
        for i in self.RawData['Cons'].keys():
            self.ProData['Cons'][i]={}
            self.ProData['Cons'][i]['General']=self.RawData['Cons'][i]['General']
            if i.find('Ref')==-1:
                # Empty=np.zeros(np.shape(self.RawData['Cons'][i]['ConnectionForce'][:,0:3]))
                self.RawData['Cons'][i]['ConnectionForce']= np.hstack((self.RawData['Cons'][i]['ConnectionForce'][:,0:6], np.linalg.norm(self.RawData['Cons'][i]['ConnectionForce'][:,0:3],axis=1)[:,None]))
                #### Caclualte angles of the loads wrt to platform, this does not work well yet
                ConForce=self.RawData['Cons'][i]['ConnectionForce'][:,0:3]
                Angles=np.zeros(np.shape(ConForce))
                Angles[:,0]=np.arctan(ConForce[:,2]/ConForce[:,1])*180/np.pi
                Angles[:,1]=np.arctan(ConForce[:,2]/ConForce[:,0])*180/np.pi
                Angles[:,2]=np.arctan(ConForce[:,1]/ConForce[:,0])*180/np.pi
                Angles[np.isnan(Angles)]=0
                self.RawData['Cons'][i]['Angles']=Angles
            ## Get statistics
            Variables=self.RawData['Cons'][i].keys()
            for V in Variables:
                if V!='General':
                    self.ProData['Cons'][i][V]={}
                    self.ProData['Cons'][i][V]['Statistics']=self.statistics(self.RawData['Cons'][i][V])           
            if self.IncFEMLabel:
                NameSplt=i.split('_')
                Name='_'
                self.ProData['Cons'][i]['General']['FEM Label']=Name.join(NameSplt[:2])+'_'+Name.join(NameSplt[3:])
                
                
    def coupling(self, RawData):
        self.ProData['Coupling']={}
        for i in RawData['Coupling'].keys():
            self.ProData['Coupling'][i]={}
            self.ProData['Coupling'][i]['General']=RawData['Coupling'][i]['General']
            self.ProData['Coupling'][i]['General']['PlatformType']=RawData['Platform'][RawData['Coupling'][i]['General']['Platform']]['General']['PlatformType']      
            # Calc and add power
            if self.RawData['Coupling'][i]['General']['LinDamping']:
                self.ProData['Coupling'][i]['Power']={}
                self.ProData['Coupling'][i]['Power']['TimeTrace']=RawData['Coupling'][i]['General']['LinDamping']*RawData['Coupling'][i]['Velocity']**2
                self.ProData['Coupling'][i]['Power']['Statistics']=self.statistics(self.ProData['Coupling'][i]['Power']['TimeTrace'], Distance=True)
            # Extract other data
            # Variables='Tension', 'Length', 'Velocity'
            for V in self.RawData['Coupling'][i].keys():
                if V!='General':
                    self.ProData['Coupling'][i][V]={}
                    self.ProData['Coupling'][i][V]['Statistics']=self.statistics(RawData['Coupling'][i][V], Distance=True)     


    def get_excursions(self, i, POI, Tag):    
        ## Gets excursions, this only works if the tag is the same as the tag for a certain position 
        self.RawData['Platform'][i]['Excursions'][Tag]={}
        General=self.RawData['Platform'][i]['Positions'][Tag]['General'].copy()
        self.RawData['Platform'][i]['Excursions'][Tag]['General']=General
        Excursions=self.RawData['Platform'][i]['Positions'][Tag]['TimeHistory'][:,:3]-General['Global Pos']
        Horz=(Excursions[:,0]**2+Excursions[:,1]**2)**0.5    
        self.RawData['Platform'][i]['Excursions'][Tag]['TimeHistory']=np.hstack((Excursions,np.zeros(Excursions.shape),Horz[:,None]))    
       
    def platform(self, Platform, StructuralOutputPOI_Acc=False, StructuralOutputPOI_Name=False, Include_Accelerations=False, PositionPOIs=False, PositionPOIs_Name=False):
        self.ProData['Platform']={}
        for i in self.RawData['Platform'].keys():
            self.ProData['Platform'][i]={}
            self.ProData['Platform'][i]['General']=self.RawData['Platform'][i]['General']
            
            if Include_Accelerations or StructuralOutputPOI_Name: 
                #### This process takes a lot of time, therefore it is not normally inlcuded
                # First calc acceleration timetraces and add them to the RawData 
                for iP, XP in enumerate(Platform.AirGapProbesX):   
                    POI=[Platform.AirGapProbesX[iP], Platform.AirGapProbesY[iP], Platform.AirGapProbesZ[iP]] #### This lines is not safe!! needs to be adjusted the platform height can be read from the drawing output 
                    ### Test if the POI is in the correct subplatform
                    if self.RawData['Platform'][i]['General']['Split']:
                        POI0=self.translate_coordinates(Platform, POI, -float(i.split('_')[1]))
                        if POI0[0]<=POI0[1]*np.tan(30*np.pi/180)+Platform.Length/3 and POI0[1]>=0:                
                            Gen={}
                            Gen['Local Pos']=POI
                            Gen['Global Pos']=self.RawData['Platform'][i]['General']['Initial Pos']+self.rotate_coordinates_around_z(POI, self.RawData['Platform'][i]['General']['Initial Orientation'])
                            Gen['PlatformType']=self.RawData['Platform'][i]['General']['PlatformType']
                            self.RawData['Platform'][i]['Accelerations']['Probe'+str(iP)]={}
                            self.RawData['Platform'][i]['Accelerations']['Probe'+str(iP)]['General']=Gen 
                            self.RawData['Platform'][i]['Accelerations']['Probe'+str(iP)]['TimeHistory']=self.POI_translation(POI, self.RawData['Platform'][i]['Accelerations']['Origin']['TimeHistory'])
                    else:
                        Gen={}
                        Gen['Local Pos']=POI
                        Gen['Global Pos']=self.RawData['Platform'][i]['General']['Initial Pos']+self.rotate_coordinates_around_z(POI, self.RawData['Platform'][i]['General']['Initial Orientation'])
                        Gen['PlatformType']=self.RawData['Platform'][i]['General']['PlatformType']
                        self.RawData['Platform'][i]['Accelerations']['Probe'+str(iP)]={}
                        self.RawData['Platform'][i]['Accelerations']['Probe'+str(iP)]['General']=Gen 
                        self.RawData['Platform'][i]['Accelerations']['Probe'+str(iP)]['TimeHistory']=self.POI_translation(POI, self.RawData['Platform'][i]['Accelerations']['Origin']['TimeHistory'])
              
                if StructuralOutputPOI_Acc:
                    for iSOPOI, SOPOI in enumerate(StructuralOutputPOI_Acc):
                        Gen={}
                        Gen['Local Pos']=SOPOI
                        Gen['Global Pos']=self.RawData['Platform'][i]['General']['Initial Pos']+self.rotate_coordinates_around_z(SOPOI, self.RawData['Platform'][i]['General']['Initial Orientation'])
                        Gen['PlatformType']=self.RawData['Platform'][i]['General']['PlatformType']
                        self.RawData['Platform'][i]['Accelerations'][StructuralOutputPOI_Name[iSOPOI]]={}
                        self.RawData['Platform'][i]['Accelerations'][StructuralOutputPOI_Name[iSOPOI]]['General']=Gen 
                        self.RawData['Platform'][i]['Accelerations'][StructuralOutputPOI_Name[iSOPOI]]['TimeHistory']=self.POI_translation(SOPOI, self.RawData['Platform'][i]['Accelerations']['Origin']['TimeHistory'])
            
            ## Extract positions at a certain point (CoG)
            POI=self.RawData['Platform'][i]['General']['CoG']
            self.POI_translation_position(i, POI, 'CoG')

            ## Extract exursions from pre-set netural point (GlobalPosition)
            self.RawData['Platform'][i]['Excursions']={}
            self.get_excursions( i, POI, 'CoG')
            self.get_excursions(i, [0,0,0], 'Origin')      
            
            if PositionPOIs:
                for iPOI, POI in enumerate(PositionPOIs):
                    self.POI_translation_position(i, POI, PositionPOIs_Name[iPOI])
                    self.get_excursions( i, POI, PositionPOIs_Name[iPOI])
            
            for V in ['Positions', 'Accelerations', 'MotionsFromStaticPos', 'Excursions']:
                self.ProData['Platform'][i][V]={}
                for k in self.RawData['Platform'][i][V].keys():
                    if V!='AirGap' and k!='POI Normal 0' and k!='POI Normal 120' and k!='POI Normal 240':
                        self.ProData['Platform'][i][V][k]={}
                        self.ProData['Platform'][i][V][k]['General']=self.RawData['Platform'][i][V][k]['General']
                        self.ProData['Platform'][i][V][k]['Statistics']=self.statistics(self.RawData['Platform'][i][V][k]['TimeHistory'])            
                        
    def floater_beams(self, RawData, Mooring=None, Type='Floater'):
        self.ProData[Type]={}
        for i in RawData[Type].keys():
            self.ProData[Type][i]={}
            self.ProData[Type][i]['General']=RawData[Type][i]['General']
            self.ProData[Type][i]['General']['PlatformType']=RawData['Platform'][RawData[Type][i]['General']['Platform']]['General']['PlatformType']   
            if Type == 'Floater' and Mooring:
                Split=i.split('Floater')
                NameSplt=Split[1].split('_')
                #### This if statement recognizes floaters that are at the edge
                if any(((np.array(Mooring.ConnectedPlatformsTotal)== float(NameSplt[0])) & (np.array(Mooring.OrientationLocalTotal)==float(NameSplt[1]))) |
                    (((np.array(Mooring.ConnectedPlatformsTotal)== float(NameSplt[0])) & (np.array(Mooring.OrientationLocalTotal)==np.mod(float(NameSplt[1])-120,360))))):
                    self.ProData[Type][i]['General']['OuterFloater']='Yes'
                else:
                    self.ProData[Type][i]['General']['OuterFloater']='No'
            # elif Type == 'Beams':
            #     #### get member force
            #     Member1_Force = RawData[Type][i]['Effective Tension']/4 + RawData[Type][i]['x bend moment']/
                
            # for V in ['ConnectionForce']:
            for V in RawData[Type][i].keys():
                if V!='General' and V!='AirGap' and V!='Position' and V!='Declination':
                    self.ProData[Type][i][V]={}
                    self.ProData[Type][i][V]['Statistics']=self.statistics(RawData[Type][i][V])    
                elif V=='Position':
                    self.ProData[Type][i][V]={}
                    self.ProData[Type][i][V]['Statistics']=self.statistics(RawData[Type][i][V]['TimeHistory'])
                elif V == 'Declination' and Type == 'Beams':
                    RawData[Type][i][V] = RawData[Type][i][V] - 90 #### Declination = 90 means point in the XY plane
                    self.ProData[Type][i][V]={}
                    self.ProData[Type][i][V]['Statistics']=self.statistics(RawData[Type][i][V])
            if self.IncFEMLabel:
                if Type=='Floater':
                    Split=i.split('Floater')
                    NameSplt=Split[1].split('_')
                    Name='_'
                    self.ProData[Type][i]['General']['FEM Label']='Floater_'+Name.join(NameSplt[1:])
                                      
    def airgap(self, RawData):
        self.ProData['AirGap']={}
        for i in RawData['Platform'].keys():    
            self.ProData['AirGap'][i]={}
            for j in RawData['Platform'][i]['AirGap'].keys(): 
                self.ProData['AirGap'][i][j]={}
                self.ProData['AirGap'][i][j]['General']= RawData['Platform'][i]['AirGap'][j]['General']
                self.ProData['AirGap'][i][j]['General']['Reference']='Platform' 
                self.ProData['AirGap'][i][j]['General']['PlatformType']= RawData['Platform'][i]['General']['PlatformType']
                self.ProData['AirGap'][i][j]['General']['UniqueName']=i+'_'+j
                TimeTrace=RawData['Platform'][i]['AirGap'][j]['TimeHistory'] 
                self.ProData['AirGap'][i][j]['Statistics']=self.statistics(TimeTrace)
                CrossInd =self.find_down_crossings(TimeTrace,0)
                self.ProData['AirGap'][i][j]['Statistics']['Num Contacts']=len(CrossInd)
        for j in RawData['Floater'].keys(): 
            i=RawData['Floater'][j]['General']['Platform']
            self.ProData['AirGap'][i][j]={}
            self.ProData['AirGap'][i][j]['General']=RawData['Floater'][j]['AirGap']['General']
            self.ProData['AirGap'][i][j]['General']['Reference']='Floater'
            self.ProData['AirGap'][i][j]['General']['PlatformType']= RawData['Platform'][i]['General']['PlatformType']
            self.ProData['AirGap'][i][j]['General']['UniqueName']=i+'_'+j
            TimeTrace=RawData['Floater'][j]['AirGap']['TimeHistory']
            self.ProData['AirGap'][i][j]['Statistics']=self.statistics(TimeTrace)
            CrossInd =self.find_down_crossings(TimeTrace,0)
            self.ProData['AirGap'][i][j]['Statistics']['Num Contacts']=len(CrossInd)
                  
    def gap_details(self, RawData):
        self.ProData['Gap']={}
        for i in RawData['Gap'].keys():
            MotionsA=RawData['Gap'][i]['EndA']
            MotionsB=RawData['Gap'][i]['EndB']
            DistanceVector=MotionsB-MotionsA
            NormalVectors=[[1,0,0], [0,1,0], [0,0,1]] 
            Distances=[]
            for N in NormalVectors:
                N=self.rotate_coordinates_around_z(N, RawData['Gap'][i]['General']['Direction'])
                Distance=np.matmul(N,DistanceVector.transpose())
                Distances.append(Distance)
            Distances=np.asarray(Distances)
            self.ProData['Gap'][i]={} 
            self.ProData['Gap'][i]['TimeTrace']=Distances.transpose() ## These detail are all describing delta width, height and shear translation 
            self.ProData['Gap'][i]['Statistics']=self.statistics(Distances.transpose())  
            self.ProData['Gap'][i]['General']=RawData['Gap'][i]['General']  
            self.ProData['Gap'][i]['General']['Variable']='Gap'

    def clash_details(self, Layout):
        self.ProData['Clash']={}
 
        #### A lot has been commmented here to remove the coupling clash calcualtions
        #### This angle represent the distance of the roation point to the platform
        # Angle2=0.09103477803741532
        
        for i in self.RawData['CouplingClash'].keys():
            Plat=self.RawData['CouplingClash'][i]['General']['Platform End B']
            Orient=self.RawData['CouplingClash'][i]['General']['Coupling Direction']            
            # #### First calculate the different between coupling end and plane normal of the traiangle with a reference point
            # PlaneNormal=self.RawData['Platform'][Plat]['Positions']['POI Normal '+str(Orient)]['TimeHistory']-self.RawData['Platform'][Plat]['Positions']['Origin']['TimeHistory'][:,:3]
            # POI=self.RawData['CouplingClash'][i]['EndA']
            # POR=self.RawData['CouplingClash'][i]['EndB'] ## This is the reference that is on the plane 
            # DiffPOIPOR= np.subtract(POI,POR)
            # ClashLength=np.sum(-1*PlaneNormal*DiffPOIPOR, 1)
            # #### Calculate the angle and add andlge due to connection poitions on T1 platform use platform height to calculate clahsing distance
            # Angle1=np.arcsin(ClashLength/np.linalg.norm(DiffPOIPOR, axis=1))
            # self.RawData['CouplingClash'][i]['Clearance']=self.RawData['CouplingClash'][i]['General']['Local Pos End B'][2]*np.sin(Angle1+Angle2)
            # self.ProData['Clash']['Coupling'+i]={}
            # self.ProData['Clash']['Coupling'+i]['General']=self.RawData['CouplingClash'][i]['General'].copy()
            # self.ProData['Clash']['Coupling'+i]['General']['Type']='Coupling'
            # self.ProData['Clash']['Coupling'+i]['Clearance']={}           
            # self.ProData['Clash']['Coupling'+i]['Clearance']['Statistics']=self.statistics(ClashLength)  
            # del self.ProData['Clash']['Coupling'+i]['General']['LinDamping']
            # del self.ProData['Clash']['Coupling'+i]['General']['Stiffness']
            # del self.ProData['Clash']['Coupling'+i]['General']['UnstretchedLength']
            
            ## Get the floater clash -> Horizontal tip distance - diameter
            PlatA=self.RawData['CouplingClash'][i]['General']['Platform'].split('_')[0]
            PlatB=Plat.split('_')[0]
            if i.find('-1')==-1:
                OrientA=np.mod(float(Orient)+120,360)
                OrientB=Orient
            else:
                OrientA=Orient
                OrientB=np.mod(float(Orient)+120,360)
            NameA='Floater'+PlatA+'_'+str(round(OrientA))
            NameB='Floater'+PlatB+'_'+str(round(OrientB))
            
            ## Suggest changeing NameA and NameB such that clashing is plotted better
            ## Added np.mod as NameA and NameB are fliped 
            self.calc_floater_clash(NameB, NameA, Orient , i, Extra='T12')
            
            #### Commented for efficieny, this is not needed for this desing 
            # ## Get the floater clash for T2 floaters
            # if OrientA==120:
            #     NameB='Floater'+ str(int(float(PlatA)+2))+'_0'
            #     self.calc_floater_clash(NameA, NameB, 270, i, Extra='T22')
            # if OrientA==240: 
            #     NameB='Floater' + str(int(float(PlatA)-Layout.Matrix.shape[1]-1))+'_120'
            #     self.calc_floater_clash(NameA, NameB, 30, i,Extra='T22')
            # if OrientA==0: 
            #     NameB='Floater' + str(int(float(PlatA)+Layout.Matrix.shape[1]-1))+'_240'
            #     self.calc_floater_clash(NameA, NameB, 150, i,Extra='T22')
            
    def calc_bottom_clash(self, ConeHeight):
        #### Calculate the floaters bottom clearance
        for FloatName in self.RawData['Floater'].keys():
            Name='Ground_'+ FloatName.split('Floater')[1]
            Gen={}
            Gen['Type']='Ground'
            Gen['Coupling Direction']=FloatName.split('_')[1]
            Gen['Global Pos']=self.RawData['Floater'][FloatName]['Position']['General']['Global Pos']
            Gen['Local Pos']=self.RawData['Floater'][FloatName]['Position']['General']['Local Pos']
            Gen['Local Pos End B']=np.array([0,0,0])
            Gen['Connection A']=FloatName
            Gen['Conenction B']='Ground'            
            # Gen['LinDamping']='-'
            Gen['Platform']=self.RawData['Floater'][FloatName]['General']['Platform']
            Gen['Platform End B']='-'
            # Gen['Stiffness']='-'
            
            # Gen['UnstretchedLenght']=0
            self.ProData['Clash'][Name]={}
            self.ProData['Clash'][Name]['General']=Gen
            self.ProData['Clash'][Name]['Clearance']={}
            Clearance=self.RawData['Floater'][FloatName]['Position']['TimeHistory'][:,2]-ConeHeight+self.ProData['Environment']['WaterDepth']
            self.ProData['Clash'][Name]['Clearance']['Statistics']=self.statistics(Clearance)
            
    def calc_floater_clash(self, NameA, NameB, Orient, i, Extra=''):
        check=False
        Orient=np.mod(Orient+180,360)
        try:
            DXY=self.RawData['Floater'][NameA]['Position']['TimeHistory'][:,:2]-self.RawData['Floater'][NameB]['Position']['TimeHistory'][:,:2]
            check=True
        except:
            pass
        if check:
            HorRadA=np.cos(self.RawData['Floater'][NameA]['Declination']['TimeHistory']*np.pi/180)*self.RawData['Floater'][NameA]['General']['ClashDiameter']/2
            HorRadB=np.cos(self.RawData['Floater'][NameB]['Declination']['TimeHistory']*np.pi/180)*self.RawData['Floater'][NameB]['General']['ClashDiameter']/2
            ClashLength=np.linalg.norm(DXY, axis=1)-HorRadA-HorRadB
            
            self.RawData['FloaterClash']={}
            self.RawData['FloaterClash'][NameA+'_'+str(Orient)]={}
            # Gen=self.ProData['Clash']['Coupling'+i]['General'].copy()
            
            Gen=self.RawData['CouplingClash'][i]['General'].copy()
            Gen['Type']='Coupling'
            del Gen['LinDamping']
            del Gen['Stiffness']
            del Gen['UnstretchedLength']
            
            Gen['Connection A']=NameA
            Gen['Connection B']=NameB
            Gen['Local Pos']=self.RawData['Floater'][NameA]['General']['Local Pos']
            Gen['Local Pos End B']=self.RawData['Floater'][NameB]['General']['Local Pos']
            NormalVector=self.rotate_coordinates_around_z([1,0,0],float(Orient))
            Gen['Global Pos']=  self.RawData['Floater'][NameA]['Position']['General']['Global Pos']+ self.RawData['Floater'][NameA]['General']['ClashDiameter']*0.5*NormalVector
            Gen['Type']='Floater'+'_'+Extra
            self.RawData['FloaterClash'][NameA+'_'+str(Orient)]['General']=Gen
            self.RawData['FloaterClash'][NameA+'_'+str(Orient)]['Clearance']=ClashLength
            
            self.ProData['Clash'][NameA+'_'+str(Orient)]={}
            self.ProData['Clash'][NameA+'_'+str(Orient)]['General']=Gen
            self.ProData['Clash'][NameA+'_'+str(Orient)]['Clearance']={}
            self.ProData['Clash'][NameA+'_'+str(Orient)]['Clearance']['Statistics']=self.statistics(ClashLength)
        
        
           
    def relative_angles(self):
        ## First calculate pitch angles in each plane based on normal vectors
        for i in self.RawData['Platform'].keys():
            for j in [0, 120, 240]:
                #### Find phi y for rotated axis system ( pitch)
                DeltaZ=self.RawData['Platform'][i]['Positions']['POI Normal '+str(j)]['TimeHistory'][:,2]-self.RawData['Platform'][i]['Positions']['Origin']['TimeHistory'][:,2]
                PhiY=np.arcsin(-1*DeltaZ)*180/np.pi 
                #### Find phi z for rotated axis system (relative yaw) ### Yaw code not working, don't understand why??
                # DeltaY=self.RawData['Platform'][i]['Positions']['POI Normal '+str(j)]['TimeHistory'][:,1]-self.RawData['Platform'][i]['Positions']['Origin']['TimeHistory'][:,1]
                # DeltaX=self.RawData['Platform'][i]['Positions']['POI Normal '+str(j)]['TimeHistory'][:,0]-self.RawData['Platform'][i]['Positions']['Origin']['TimeHistory'][:,0]
                # PhiZ=np.arctan2(DeltaY,DeltaX)*180/np.pi 
                #### Find phi x for rotated axis system (relative roll)
                POI=self.rotate_coordinates_around_z([0,1,0],j)
                self.POI_translation_position( i, POI, 'POI Normal y' + str(j))
                DeltaZ=self.RawData['Platform'][i]['Positions']['POI Normal y'+str(j)]['TimeHistory'][:,2]-self.RawData['Platform'][i]['Positions']['Origin']['TimeHistory'][:,2]
                PhiX=np.arcsin(DeltaZ)*180/np.pi 
                
                self.RawData['Platform'][i]['Pitch'+str(j)]={}
                self.RawData['Platform'][i]['Roll'+str(j)]={}
                # self.RawData['Platform'][i]['Yaw'+str(j)]={}
                self.RawData['Platform'][i]['Pitch'+str(j)]['TimeHistory']=PhiY
                self.RawData['Platform'][i]['Roll'+str(j)]['TimeHistory']=PhiX
                # self.RawData['Platform'][i]['Yaw'+str(j)]['TimeHistory']=PhiZ
                # self.ProData['Platform'][i]['Pitch'+str(j)]['Statistics']=self.statistics(Phi)
                
        ## Calculate realtive pitch based on sequence also used for gap calculations
        self.ProData['RelativeAngle']={}
        self.RawData['RelativeAngle']={}
        for i in self.RawData['Gap'].keys():
            if i.find('-')==-1:
                self.ProData['RelativeAngle'][i]={}
                for t in ['Pitch', 'Roll']:#, 'Yaw']:
                    PlatA=self.RawData['Gap'][i]['General']['EndA Reference']
                    PlatB=self.RawData['Gap'][i]['General']['EndB Reference']
                    Direction=self.RawData['Gap'][i]['General']['Direction']
                    
                    if t=='Yaw':
                        #### Thiss approach does not seem to work well...
                        AngleA=self.RawData['Platform'][PlatA]['MotionsFromStaticPos']['Origin']['TimeHistory'][:,5]
                        AngleB=self.RawData['Platform'][PlatB]['MotionsFromStaticPos']['Origin']['TimeHistory'][:,5]
                        AngleRelative=AngleA-AngleB
                    else:
                        AngleA=self.RawData['Platform'][PlatA][t+str(Direction)]['TimeHistory']
                        AngleB=self.RawData['Platform'][PlatB][t+str(Direction)]['TimeHistory']
                        AngleRelative=AngleA+AngleB
                    ### raw data section currently commented could be uncommented for instant value calcs
                    # self.RawData['RelativePitch'][i][t]={}
                    Gen={}
                    Gen['Connection A']=PlatA
                    Gen['Connection B']=PlatB
                    Gen['RelativeAngleDirection']=Direction
                    Gen['Global Pos']=(self.RawData['Gap'][i]['General']['Global Plot Pos']+ self.RawData['Gap'][i[:-1]+'-1']['General']['Global Plot Pos'])/2
                    # self.RawData['RelativePitch'][i][t]['General']=Gen
                    # self.RawData['RelativePitch'][i][t]['TimeHistory']=PitchRelative
                    
                    self.ProData['RelativeAngle'][i][t]={}                                
                    self.ProData['RelativeAngle'][i][t]['Statistics']=self.statistics(AngleRelative)
                    self.ProData['RelativeAngle'][i]['General']=Gen
                
    def sytheticLine_statistics(self,df,Slice=False):
        Dict={}
        if Slice:
            df = df.iloc[:-1]
        Dict['Max']=np.max(df['Max'],0)
        Dict['Min']=np.min(df['Min'],0)
        Dict['StDev']=np.std(df['StdDev'],0)
        Dict['Mean']=np.mean(df['Mean'],0)
        Dict['Sig Amplitude']= 2*Dict['StDev']
        Dict['AbsMax']=np.max([Dict['Max'] ,np.abs(Dict['Min'])],0)
        return Dict
    
    def mooring(self):
        self.ProData['Mooring']={}
        for i in self.RawData['Mooring'].keys():
            self.ProData['Mooring'][i]={}
            self.ProData['Mooring'][i]['General']=self.RawData['Mooring'][i]['General']
            for j in self.RawData['Mooring'][i].keys():
                # print(j)
                if j!='General' and j!='SeaBedClearance': ### @min also for the line section otherwise 
                    self.ProData['Mooring'][i][j]={}
                    self.ProData['Mooring'][i][j]['Statistics']=self.statistics(self.RawData['Mooring'][i][j])
                if j == 'SeaBedClearance':
                    for k in self.RawData['Mooring'][i][j].keys():
                        if k.find('KappaPolyester') != -1:
                            if k.find('_1') != -1:
                                name = 'SynLineShort'
                                Slice = True
                            else:
                                name = 'SynLineLong'
                                Slice = False
                            self.ProData['Mooring'][i][j+name]={}
                            self.ProData['Mooring'][i][j+name]['Statistics']=self.sytheticLine_statistics(self.RawData['Mooring'][i][j][k],Slice=Slice)
                    
class ToExcel(GeneralFunctions):
    ## This class has been made to output data to excel such that pivot charts can be made
    ## The idea is that each functions creates its own sheet in the output Excel
    
    def __init__(self, DataFiles=False):
        self.Sims=DataFiles
        
    def get_environment(self, ProData):
        Out=[]
        Headers=ProData['Environment'].keys()
        for key in Headers:
            Out.append(ProData['Environment'][key]) 
        return Headers, Out
    
        
    def get_data(self, Dict, Length=False, Type=False):
        Out=[]
        Headers=Dict.keys()
        Header2=[]
        ## HorzForce has been added for the floater connection was the easiest way to implement it
        Extra=(' X',' Y', ' Z', ' rX', ' rY', ' rZ', 'Horz', 'BendMoment')
        if Type=='Cons':
            Extra=( ' X',' Y', ' Z', ' rX', ' rY', ' rZ', ' Total')
            Length=7
        for key in Headers:
            A=np.array(Dict[key])
            if not Length:
                if A.size>1:
                    for i in range(len(A)):
                        Header2.append(key+Extra[i])
                        Out.append(Dict[key][i])
                else: 
                    Out.append(Dict[key])
                    Header2.append(key)    
            else: 
                for i in range(Length):
                    Header2.append(key+Extra[i])
                    if i<len(A):
                        Out.append(Dict[key][i])
                    else:
                        Out.append(None)
        return Header2, Out
            
    def get_statistics_gap(self, Dict):
        Out=[]
        Headers=Dict.keys()
        Headers2=[]
        for key in Headers:
            Out.extend(Dict[key])
            for j in ['PERH', 'PAR', 'PERV']:
                Headers2.append(key + ' ' + j)          
        return Headers2, Out
    
    def save_excel(self, DataFrame, ExcelFile, SaveExcel=True, SheetName='Sheet1'):
        if SaveExcel:
            writer = pd.ExcelWriter(ExcelFile)
            DataFrame.to_excel(writer,SheetName)
            
            #### cell backgroud setting only applicable to single directionary
            workbook=writer.book
            worksheet=writer.sheets[SheetName]
            ##define background color
            # formatred=workbook.add_format({'bg_color':'red'})
            # # col = DataFrame.columns.get_loc("Pass [Yes / No]")
            # col = np.where(DataFrame.columns == "Pass [Yes / No]")[0] ## work for columns with same name
            # for ii in col:
                # for loc, val in enumerate(DataFrame.iloc[:,ii]):
                    # if val == 'No':
                        # worksheet.write(loc+1, ii+1, DataFrame.iloc[loc,ii], formatred)
            writer.close()
            
    def append_excel_sheet(self, DataFrame, ExcelFile, SheetName='Sheet1'):
        book = load_workbook(ExcelFile)
        writer = pd.ExcelWriter(ExcelFile,  engine = 'openpyxl')
        writer.book = book
           
        # writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
        
        DataFrame.to_excel(writer, SheetName)
        writer.close()
        
    def get_extras(self, Sim):
        SimBase=os.path.basename(Sim)
        # print(Sim)
        array=SimBase.split('_')
        ## This is still hard coded such that it can be ajusted if file naming is incorrect
        x=2
        HeadOut=[]
        Out=[]
        while array[x].find('Wdp')!= 0 and array[x].find('Hs')!=0:
            # print(x)
            HeadOut.append(array[x])
            x=x+1
            Out.append(array[x])  
            x=x+1
        return HeadOut, Out
    
    def set_headers(self,Headers,EnvHeader,StatHeader,GenHeader,HeadOutEx):
        Headers.extend(EnvHeader)
        Headers.extend(StatHeader)
        Headers.extend(GenHeader)
        Headers.extend(HeadOutEx)
        
    def air_gap(self, ExcelFile, SaveExcel) :
        Out3=[]
        for Sim in self.Sims:
            # print(Sim)
            ProData=pickle.load(open(Sim,'rb'))
            [EnvHeader, EnvOut]= self.get_environment(ProData)
            Out2=[]
            for i in ProData['AirGap'].keys():
                Out1=[] 
                for j in ProData['AirGap'][i].keys(): 
                    ## Ad basis info
                    Out=[os.path.basename(Sim)[:-8], i, j, 'AirGap']
                    # Add Environemnt 
                    Out.extend(EnvOut)
                    # Statistical in 
                    StatHeader, StatOut= self.get_data(ProData['AirGap'][i][j]['Statistics'])
                    Out.extend(StatOut)
                    # Add ExtraData
                    GenHeader, GenOut = self.get_data(ProData['AirGap'][i][j]['General'])
                    Out.extend(GenOut)
                    HeadOutEx, OutEx= self.get_extras(Sim)
                    Out.extend(OutEx)
                    Out1.append(Out)
                Out2.extend(Out1)
            Out3.extend(Out2)
        Headers=['SimName', 'Platform', 'ProbeName', 'Variable']
        self.set_headers(Headers,EnvHeader,StatHeader,GenHeader,HeadOutEx)
        DataFrame=pd.DataFrame(Out3, columns=Headers)
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl') 
        self.save_excel(DataFrame, ExcelFile, SaveExcel)
        return DataFrame
        
    def coupling_gap(self, ExcelFile, SaveExcel) :
        Out3=[]
        for Sim in self.Sims:
            ProData=pickle.load(open(Sim,'rb'))
            [EnvHeader, EnvOut]= self.get_environment(ProData)
            Out1=[]
            for i in ProData['Gap'].keys():
                ## Ad basis info
                Out=[os.path.basename(Sim)[:-8], i]
                # Add Environemnt 
                Out.extend(EnvOut)
                # Statistical in 
                StatHeader, StatOut= self.get_statistics_gap(ProData['Gap'][i]['Statistics'])
                Out.extend(StatOut)
                # Add extra data 
                GenHeader, GenOut = self.get_data(ProData['Gap'][i]['General'])
                Out.extend(GenOut)
                HeadOutEx, OutEx= self.get_extras(Sim)
                Out.extend(OutEx)
                Out1.append(Out)
            Out3.extend(Out1)
        Headers=['SimName', 'ProbeName']
        self.set_headers(Headers,EnvHeader,StatHeader,GenHeader,HeadOutEx)
        DataFrame=pd.DataFrame(Out3, columns=Headers)
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl') 
        self.save_excel(DataFrame, ExcelFile, SaveExcel)
        return DataFrame

    def mooring(self, ExcelFile, SaveExcel):
        Out3=[]
        for Sim in self.Sims:
           ProData=pickle.load(open(Sim,'rb'))
           [EnvHeader, EnvOut]= self.get_environment(ProData)
           Out2=[]       
           for i in ProData['Mooring'].keys():
               ## Ad basis info
               Out=[os.path.basename(Sim)[:-8],i, 'Tension']
               # Add Environemnt 
               Out.extend(EnvOut)
               # Statistical in 
               StatHeader, StatOut= self.get_data(ProData['Mooring'][i]['Statistics'])
               Out.extend(StatOut)
               # Add Extra Data
               GenHeader, GenOut = self.get_data(ProData['Mooring'][i]['General'])
               Out.extend(GenOut)
               HeadOutEx, OutEx= self.get_extras(Sim)
               Out.extend(OutEx)
               Out2.append(Out)
           Out3.extend(Out2)       
        Headers=['SimName', 'UniqueName', 'Variable']
        self.set_headers(Headers,EnvHeader,StatHeader,GenHeader,HeadOutEx)
        DataFrame=pd.DataFrame(Out3, columns=Headers)
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl') 
        self.save_excel(DataFrame, ExcelFile, SaveExcel)
        return DataFrame
    
    def stress_raos(self, ExcelFile, Nodes, PlatformType, Type='PlatformStress', SaveExcel=False):
        ## This for a floater or platform type pro data has to do with the way data is saved
        Out1=[]
        for Sim in self.Sims:
            ProData=pickle.load(open(Sim,'rb'))
            [EnvHeader, EnvOut]= self.get_environment(ProData)
            for i in ProData[Type].keys():
                if ProData['Platform'][i+'_0']['General']['PlatformType']==PlatformType:
                    for var in ProData[Type][i].keys():
                        if var !='General' and var!='Pitch120' and var!='Pitch0' and var!='Pitch240':
                            # for loc in ProData[Type][i][var].keys():
                                # for i in enumerate(ProData[Buoy][i][j]['Statistics']['AbsMax']):
                            Out2=[os.path.basename(Sim)[:-8], 'Platform'+i, var, i]
                            Out2.extend(EnvOut)
                            
                            #### Section to get the right ampltidue data
                            ind1=ProData[Type][i][var]['P1']['Statistics']['AbsMax']>ProData[Type][i][var]['P3']['Statistics']['AbsMax']
                            ind2=ProData[Type][i][var]['P3']['Statistics']['AbsMax']>ProData[Type][i][var]['P1']['Statistics']['AbsMax']
                            #### This has been updated and needs to be checked in something like a presentation
                            ind3=(ProData[Type][i][var]['P1']['Statistics']['AbsMax']>ProData[Type][i][var]['P3']['Statistics']['AbsMax']) & (ProData[Type][i][var]['P1']['Statistics']['Min']<ProData[Type][i][var]['P3']['Statistics']['Min']*-1)
                            ind4=(ProData[Type][i][var]['P3']['Statistics']['AbsMax']>ProData[Type][i][var]['P1']['Statistics']['AbsMax']) & (ProData[Type][i][var]['P3']['Statistics']['Max']>ProData[Type][i][var]['P1']['Statistics']['Max']*-1)
                            
                            StressRange=np.zeros(np.shape(ProData[Type][i][var]['P1']['Statistics']['AbsMax']))
                            StressRange[ind1]=(ProData[Type][i][var]['P1']['Statistics']['Max']-ProData[Type][i][var]['P1']['Statistics']['Min'])[ind1]
                            StressRange[ind2]=(ProData[Type][i][var]['P3']['Statistics']['Max']-ProData[Type][i][var]['P3']['Statistics']['Min'])[ind2]
                            StressRange[ind3]=(ProData[Type][i][var]['P1']['Statistics']['Max']-ProData[Type][i][var]['P3']['Statistics']['Min'])[ind3]
                            StressRange[ind4]=(ProData[Type][i][var]['P1']['Statistics']['Max']-ProData[Type][i][var]['P3']['Statistics']['Min'])[ind4]
                            
                            StessRAO=StressRange/ProData['Environment']['HMax']
                            # StatHead=Nodes 
                            Out2.extend(StessRAO)
                            # StatHeader, StatOut= self.get_data(ProData[Type][i][var][loc]['Statistics'], Length=7)
                            # Out2.extend(StatOut)
                            # GenHeader, GenOut = self.get_data(ProData[Type][i][var][loc]['General'])
                            # Out2.extend(GenOut)
                            HeadOutEx, OutEx= self.get_extras(Sim)
                            Out2.extend(OutEx)
                            Out1.append(Out2)    
                        
        Headers=['SimName', 'UniqueName','Variable', 'Platform']
        self.set_headers(Headers,EnvHeader,Nodes,[],HeadOutEx)
        DataFrame=pd.DataFrame(Out1, columns=Headers)   
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl')                       
        self.save_excel(DataFrame, ExcelFile, SaveExcel)   
        return DataFrame   
    
    def floater_cons_coupling_clash(self, ExcelFile, Type, SaveExcel=False):
        ## This for a floater or platform type pro data has to do with the way data is saved
        Out3=[]
        for Sim in self.Sims:
            ProData=pickle.load(open(Sim,'rb'))
            if not (Type=='Beams' and not ProData['Environment']['PlatformSplit']):
                [EnvHeader, EnvOut]= self.get_environment(ProData)
                for i in ProData[Type].keys():
                    for j in ProData[Type][i].keys():
                        if j !='General' and j!='HalfCycles':
                            # and j!='General1' and j!='Pitch120' and j!='Pitch0' and j!='Pitch240':
                            # for i in enumerate(ProData[Buoy][i][j]['Statistics']['AbsMax']):
                            Out=[os.path.basename(Sim)[:-8],i,j]
                            Out.extend(EnvOut)
                            if Type=='Floater':    
                                StatHeader, StatOut= self.get_data(ProData[Type][i][j]['Statistics'], Length=8)
                            else:
                                StatHeader, StatOut= self.get_data(ProData[Type][i][j]['Statistics'], Type=Type)
                            Out.extend(StatOut)
                            GenHeader, GenOut = self.get_data(ProData[Type][i]['General'])
                            Out.extend(GenOut)
                            HeadOutEx, OutEx= self.get_extras(Sim)
                            Out.extend(OutEx)
                            Out3.append(Out)
        Headers=['SimName', 'UniqueName','Variable']
        self.set_headers(Headers,EnvHeader,StatHeader,GenHeader,HeadOutEx)
        DataFrame=pd.DataFrame(Out3, columns=Headers)
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl')
        self.save_excel(DataFrame, ExcelFile, SaveExcel)   
        return DataFrame
    
        
    def platform(self, ExcelFile, Type, SaveExcel=False):
        ## This for a floater or platform type pro data has to do with the way data is saved
        Out1=[]
        for Sim in self.Sims:
            ProData=pickle.load(open(Sim,'rb'))
            [EnvHeader, EnvOut]= self.get_environment(ProData)
            for i in ProData[Type].keys():
                for var in ProData[Type][i].keys():
                    if var !='General' :
                        for loc in ProData[Type][i][var].keys():
                            # for i in enumerate(ProData[Buoy][i][j]['Statistics']['AbsMax']):
                            Out2=[os.path.basename(Sim)[:-8], 'Platform'+i+'_'+loc , var, i]
                            Out2.extend(EnvOut)
                            StatHeader, StatOut= self.get_data(ProData[Type][i][var][loc]['Statistics'], Length=7)
                            Out2.extend(StatOut)
                            GenHeader, GenOut = self.get_data(ProData[Type][i][var][loc]['General'])
                            Out2.extend(GenOut)
                            HeadOutEx, OutEx= self.get_extras(Sim)
                            Out2.extend(OutEx)
                            Out1.append(Out2)                    
        Headers=['SimName', 'UniqueName','Variable', 'Platform']
        self.set_headers(Headers,EnvHeader,StatHeader,GenHeader,HeadOutEx)
        DataFrame=pd.DataFrame(Out1, columns=Headers)   
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl')                       
        self.save_excel(DataFrame, ExcelFile, SaveExcel)   
        return DataFrame   
    
    def node_stress(self, ExcelFile,  SaveExcel=False):
        Out1=[]
        for Sim in self.Sims:
            StressData=pickle.load(open(Sim,'rb'))
            [EnvHeader, EnvOut]= self.get_environment(StressData)
            for i in StressData.keys():
                if i != 'Environment':
                    for j in StressData[i].keys(): 
                        if j!= 'General':
                            Out2=[os.path.basename(Sim)[:-8], i,j]
                            Headers=['SimName',  'FATName','Platform']
                            Out2.extend(EnvOut)
                            # Max=[max(arr) for arr in StressData[i][j]['HalfCyclesY']]
                            Out2.extend([max(max(arr) for arr in StressData[i][j]['HalfCyclesY'])])
                            # Out2.extend([max(Max), Max])
                            # Out2.extend([min(min(arr) for arr in StressData[i][j]['HalfCyclesY'])])
                            StatHeader=['Max']
                            # NodeHeader=StressData[i]['General']['Nodes']
                            HeadOutEx, OutEx= self.get_extras(Sim)
                            Out2.extend(OutEx)
                            Out1.append(Out2) 
            print(os.path.basename(Sim))
        self.set_headers(Headers,EnvHeader,StatHeader,[], HeadOutEx)
        DataFrame=pd.DataFrame(Out1, columns=Headers)   
        DataFrame.to_pickle(ExcelFile[:-4]+'pkl')                       
        self.save_excel(DataFrame, ExcelFile, SaveExcel)   
        
    def extract_scatters(self, df, Variable, Statistic, Tag='', GeneralSelector=None, GeneralSelection=None, dfIreg=None, Axes=None):
        df['Hs'] = df['Hs'].round(2)
        df['Tp'] = df['Tp'].round(2)
        WaveHeights = np.sort(df['Hs'].unique())
        WavePeriods = np.sort(df['Tp'].unique())
        if Variable == '':
            df = df
        else:
            df = df[df['Variable']==Variable]
        if GeneralSelector:    
            df=self.filterdf(df, GeneralSelector, GeneralSelection)
            # if 'SkidLocation' in GeneralSelection or 'POI' in GeneralSelection or 'Origin' or 'CouplingPos' or 'CouplingAcc' in GeneralSelection:
            #     Type = 'Contains'
            # else:
            #     Type = 'Exact'
            # for i, GenS in enumerate(GeneralSelector): 
            #     if Type=='Contains':
            #         try:
            #             df=df[df[GenS].str.contains(GeneralSelection[i])]
            #         except:
            #             df=df[df[GenS]==GeneralSelection[i]]
            #     elif Type=='Exact':
            #         df=df[df[GenS]==GeneralSelection[i]]
        Means = [idx for idx, val in enumerate(Statistic) if val.find('Mean')!=-1]
        for im in Means:
            Statistic.insert(im+1,'m-'+Statistic[im]) ### to incude min mean
        dfTable = {}
        for extvar in Statistic:
            dfTable[extvar] = pd.DataFrame(index=WaveHeights,columns=WavePeriods)
            for wp in df['Tp'].unique():
                df1 = df[df['Tp'] == wp]
                col = dfTable[extvar].columns.get_loc(wp)
                for wh in np.round(df1['Hs'].unique(),2):
                    df2 = df1[df1['Hs']== wh]
                    # if CouplingPosOrACC:
                    #     dfnew = pd.DataFrame()
                    #     for pf in df2['Platform'].unique():
                    #         # print(pf)
                    #         df3 = (df2[df2['Platform']==pf])
                    #         for var in df3['UniqueName']:
                    #             if pf.find('_0') != -1 and (var.find('_240_-1')!=-1 or var.find('_0_1')!=-1):
                    #                 dfnew = pd.concat([dfnew,df3[df3['UniqueName']==var]])
                    #             elif pf.find('_120') != -1 and (var.find('_0_-1')!=-1 or var.find('_120_1')!=-1):
                    #                 dfnew = pd.concat([dfnew,df3[df3['UniqueName']==var]])
                    #             elif pf.find('_240') != -1 and (var.find('_120_1')!=-1 or var.find('_240_1')!=-1):
                    #                 dfnew = pd.concat([dfnew,df3[df3['UniqueName']==var]])
                            
                            
                    #         if pf.find('_0') != -1:
                    #             for var in df3['UniqueName']:
                    #                 if var.find('_240_-1')!=-1 or var.find('_0_1')!=-1:
                    #                     dfnew = pd.concat([dfnew,df3[df3['UniqueName']==var]])
                    #                 elif var.find('_0-1') !=-1 or var.find('_120-1') != -1:
                    #                     dfnew = pd.concat([dfnew,df3[df3['UniqueName']==var]])
                                
                    #             # Newdf.append(df2[df2['UniqueName']==[j for j in df2['UniqueName'] if j.find('_240_-1')!=-1 or j.find('_0-1')!=-1]])
                    if extvar.find('Min')!=-1:
                        df3 = df2[extvar].min()
                    elif extvar.find('m-Mean')!=-1:
                        df3 = df2[extvar.split('-')[1]].min()
                    else:
                        df3 = df2[extvar].max()
                    ind = dfTable[extvar].index.get_loc(wh)
                    dfTable[extvar].iloc[ind,col] = round(df3,2)
            # dfTable[extvar] = dfTable[extvar].fillna(0) ### Replcae NAN by 0
        if len(Statistic[0].split(' ')) == 1:
            # Pos = Statistic[0].split(' ')[0]
            Pos = ''
        elif len(Statistic[0].split(' ')) == 2:
            Pos = Statistic[0].split(' ')[1]
            Variable = Variable+' '
        writer = pd.ExcelWriter(os.path.join(self.ExcelPath,'Scatters'+'_'+Tag+Variable+Pos+'.xlsx'), engine='xlsxwriter')
        for key in dfTable.keys():
            df =  dfTable[key]
            if key.find('Min')!=-1:
                sheetname = 'Min ' +key
            elif key.find('m-Mean')!=-1:
                sheetname = 'Min '+key.split('-')[1]
            else:
                sheetname = 'Max '+ key
            df.to_excel(writer, sheet_name=sheetname, index=True)
        writer.close()  
        # df=self.filterdf(df, GeneralSelector, GeneralSelection)
        # a = 1
        # Plot.plot_scatter(dfPlatform, Variable='MotionsFromStaticPos', Statistic='Max rY', function=np.max, Tag='Platform T1', Criteria='Max Platform T1 rY', GeneralSelector=['PlatformType'], GeneralSelection=[1])
        # Plot.plot_scatter(dfPlatform, Variable='Accelerations', Statistic='AbsMax Z', function=np.max, Tag='Skid' ,GeneralSelector=['Platform', 'UniqueName'], GeneralSelection=['2_240', 'SkidLocation'] )
        # WaveHeightBins = np.round(np.arange(0,df['Hs'].max()+1,1),decimals=2)
        # WavePeriodBins = np.round(np.arange(1,df['Tp'].max()+1.5,1.5),decimals=2)
        # WpCol = []
        # for iwp,wp in enumerate(WavePeriodBins):
        #     if iwp>0:
        #         wpcol = str(WavePeriodBins[iwp-1]) + '-' + str(wp)
        #         WpCol.append(wpcol)
        # WhInd = []
        # for iwh,wh in enumerate(WaveHeightBins):
        #     if iwh > 0:
        #         whind = str(WaveHeightBins[iwh-1])+'-'+str(wh)
        #         WhInd.append(whind)
        # WaveHeights = np.round(np.sort(df['Hs'].unique()),2)
        # WavePeriods = np.round(np.sort(df['Tp'].unique()),2)
        
        # dfTable = {}
        # df1 = df[df['Variable'] == Variable]
        # for extvar in ExtrVars:
        #     dfTable[extvar] = pd.DataFrame(index=WaveHeights,columns=WavePeriods)
        #     for wp in np.round(df1['Tp'].unique(),2):
        #         df2 = df1[df1['Tp'] == wp]
        #         col = dfTable[extvar].columns.get_loc(wp)
        #         for wh in df2['Hs'].unique():
        #             df3 = df2[df2['Hs']== wh]
        #             if extvar == 'Min':
        #                 df4 = df3[extvar].min()
        #             else:
        #                 df4 = df3[extvar].max()
        #             ind = dfTable[extvar].index.get_loc(wh)
        #             dfTable[extvar].iloc[ind,col] = round(df4,2)
        #     dfTable[extvar] = dfTable[extvar].fillna(0)
        # writer = pd.ExcelWriter(os.path.join(self.ExcelPath,'Scatters'+'_'+Tag+'.xlsx'), engine='xlsxwriter')
        # for key in dfTable.keys():
        #     df =  dfTable[key]
        #     if key  == 'Min':
        #         sheetname = key + '_Min'
        #     else:
        #         sheetname = key + '_Max'
        #     df.to_excel(writer, sheet_name=sheetname, index=True)
        # writer.close()  
            
class Plot(GeneralFunctions):
    
    def __init__(self, Compose, LocSims, ProDataReference=None):
            self.rotate_around_z=Compose.rotate_coordinates_around_z
            self.SizeMarker=75
            self.LineWidth=1
            self.LocSims=LocSims
            self.ScatXaxis='WavePeriod'
            self.ScatLegend='HMax'
            self.HsOrHMax='HMax'
            self.DirectionMarkers=True
            self.LegendIndex=range(20)
            self.LegendLabel=False
            self.OveralSelector=None
            self.OveralSelection=None
            self.HScatterIndex=0
            self.OutputSummary=[]
            self.OutputHeaders=[]
            self.PlotRAO=False
            self.LineStyle='-'
            self.CMap=plt.get_cmap("tab10")
            self.ProDataReference=ProDataReference #### This is a prodata file that is used to make the background for localization plots
            self.PlotLocalMax=False
            self.PreviousULSCrit= False
            self.PreviousULSCritTag=None
            self.Hex2Rows = False
            self.Shape = 'Hex'
    
    # def save_excel(self, DataFrame, ExcelFile, SaveExcel=True, SheetName='Sheet1'):
    #     if SaveExcel:
    #         writer = pd.ExcelWriter(ExcelFile)
    #         DataFrame.to_excel(writer,SheetName)
    #         writer.close()
    
    def set_unit(self, Variable):
        if Variable.find('Length')!=-1 or Variable.find('dLength')!=-1 or Variable.find('Clearance')!=-1  or Variable.find('Position')!=-1  or Variable.find('AirGap')!=-1  or Variable.find('HMax')!=-1 or Variable.find('Excursions')!=-1 or Variable.find('TouchDown')!=-1:
            Unit='[m]'
        elif Variable.find('Normalized')!=-1:
            Unit='[-]'                                          
        elif Variable.find('Tension')!=-1 or Variable.find('force')!=-1 or Variable.find('Force')!=-1 :
            Unit='[kN]'
        elif Variable.find('oment')!=-1:
            Unit='[kNm]'
        elif Variable=='ConnectionForce' :
            Unit='[kN] or [kNm]'
        elif Variable=='Power':
            Unit='[kW]'
        elif Variable=='Accelerations': 
            Unit='[m/s^2]'
        elif Variable=='RelativePitch' or Variable=='WaveDirection' or Variable.find('Angle')!=-1 or Variable=='Pitch':
            Unit='[deg]'
        elif Variable=='WavePeriod'or Variable=='Tp':
            Unit='[s]'
        elif Variable=='MotionsFromStaticPos' or Variable=='Motions':
            Unit='[m] or [deg]'
        elif Variable.find('Velocity')!=1:
            Unit='[m/s]'
        else:
            Unit='Undefined'
            
        return Unit
           
    def plot_floater(self, PosGlobal, Diameter):
        rad=np.arange(0,2*np.pi,0.1)
        X=Diameter/2*np.sin(rad)+PosGlobal[0]
        Y=Diameter/2*np.cos(rad)+PosGlobal[1]
        plt.plot(X,Y, color='black',  linewidth=self.LineWidth)
    
    def plot_background(self, ProData):
        # plt.close()
        plt.figure('1')
        for i in ProData['Platform'].keys():
            InitialPos=ProData['Platform'][i]['General']['Initial Pos']
            Phi=ProData['Platform'][i]['General']['Initial Orientation']
            Coord=ProData['Platform'][i]['General']['Drawing']    
            Coord=self.rotate_around_z(Coord.transpose(), Phi).transpose()
            X=Coord[0,:][Coord[2,:]==0]+InitialPos[0]
            Y=Coord[1,:][Coord[2,:]==0]+InitialPos[1]    
            # X=Coord[0,:][Coord[2,:]==0]+InitialPos[0]
            # Y=Coord[1,:][Coord[2,:]==0]+InitialPos[1]  
            X=np.append(X,X[0])
            Y=np.append(Y,Y[0])
            plt.plot(X,Y, color='black', linewidth=self.LineWidth)
        for j in ProData['Floater'].keys():
            Diameter=ProData['Floater'][j]['General']['TopDiameter']
            PosGlobal=ProData['Floater'][j]['General']['Global Pos']
            self.plot_floater(PosGlobal,Diameter)            
            
    def set_axis(self, X, Y, ProData,dfEnvDir=False):
        plt.axis('equal')
        plt.grid()
        plt.colorbar()
        plt.xlabel('Global X [m]')
        plt.ylabel('Global Y [m]')
        plt.legend()
        fig = plt.gcf()
        # fig.set_size_inches(16.84,  9.82)
        fig.set_size_inches(16.84,  9.82)
        # plt.ylim((np.min(Y)-5, np.max(Y)+10))
        
        if isinstance(dfEnvDir, pd.Series):
            Hmax=dfEnvDir.max()
            for i, H in enumerate(dfEnvDir):
                if H == 0 or Hmax == 0:
                    pass
                else:
                    ArrowAmp=np.max(Y)/5*H/Hmax
                    ArrowAmp0=np.max(Y)/5
                    plt.arrow(-ArrowAmp0-ArrowAmp0/2,0, ArrowAmp*np.cos(dfEnvDir.index[i]*np.pi/180),ArrowAmp*np.sin(dfEnvDir.index[i]*np.pi/180), head_width=ArrowAmp0/10, head_length=ArrowAmp0/5)
                    if self.Shape == 'Row':
                        plt.arrow(-ArrowAmp0-ArrowAmp0/2,0,15*np.sin(np.radians(60)),-15*np.cos(np.radians(60)), width=0.1,head_width=2, head_length=3,color='red')
                        plt.text(-ArrowAmp0-ArrowAmp0/2 + 15*np.sin(np.radians(60))+4, -14, 'N', fontsize=12,color='red')
                    elif self.Shape == 'Hex':
                        plt.arrow(-ArrowAmp0-ArrowAmp0/2,0,0,-15, width=0.1,head_width=2, head_length=3,color='red')
                        plt.text(-ArrowAmp0-ArrowAmp0/2, -25, 'N', fontsize=12,color='red')
        else: 
            #### This selected if we make only a plot of the pro-data
            plt.text(np.min(X),np.max(Y)+5, ProData['Environment']['Name'], fontsize=12)                
            ## Set wave arrow
            ArrowAmp=np.max(Y)/5
            plt.arrow(-ArrowAmp-ArrowAmp/2,0, ArrowAmp*np.cos(ProData['Environment']['WaveDirection']*np.pi/180),ArrowAmp*np.sin(ProData['Environment']['WaveDirection']*np.pi/180), head_width=ArrowAmp/10, head_length=ArrowAmp/5)

        
    def save_fig(self, Folder,FigTitle):
        FigPath=os.path.join(self.LocSims, Folder)        
        if not os.path.exists(FigPath):
            os.makedirs(FigPath)
        plt.savefig(os.path.join(FigPath, FigTitle +'.png'),bbox_inches='tight')#,transparent=True)    
                
    def air_gap(self,ProData, Sim, Statistic, Type='Platform'):
        self.plot_background(ProData)
        X=[]
        Y=[]
        Z=[]
        Out=[]
        for i in ProData['AirGap'].keys():
            for j in ProData['AirGap'][i].keys():
                Data=ProData['AirGap'][i][j]
                X.append(Data['General']['Global Pos'][0])
                Y.append(Data['General']['Global Pos'][1])
                Z.append(Data['General']['Global Pos'][2])
                Out.append(Data['Statistics'][Statistic])  
        X=np.asarray(X)
        Y=np.asarray(Y)
        Z=np.asarray(Z) 
        Out=np.asarray(Out)
        
        Unit='[m]'
        if Statistic=='Num Contacts':
            Unit=''
        label= Statistic + ' Air Gap ' + Unit 
        
        if Type=='Platform':
            plt.scatter(X[Z>0], Y[Z>0],s=self.SizeMarker, c=Out[Z>0], label=label)
        elif Type=='Floater':
            plt.scatter(X[Z<0], Y[Z<0],s=self.SizeMarker, c=Out[Z<0], label=label)
        
        FigTitle=  Statistic + 'ArGp'+ Type + '' +   os.path.basename(Sim)[12:-32] 
        plt.title(FigTitle)
        
        self.set_axis(X,Y, ProData)
        self.save_fig('AirGap', Sim, FigTitle)
        
    def coupling_clash_platform(self, ProData, Sim, Variable, Statistic, Type='Coupling', index=None):
        self.plot_background(ProData)
        Unit=self.set_unit(Variable)       
        
        X=[]
        Y=[]
        Out=[]
        for i in ProData[Type].keys():
            if Type=='Platform':
                for j in ProData[Type][i][Variable].keys():
                    X.append(ProData[Type][i][Variable][j]['General']['Global Pos'][0])
                    Y.append(ProData[Type][i][Variable][j]['General']['Global Pos'][1])
                    Out.append(ProData[Type][i][Variable][j]['Statistics'][Statistic][index])
            else:
                try:
                    if index==None:
                        Out.append(ProData[Type][i][Variable]['Statistics'][Statistic])    
                    else:
                        Out.append(ProData[Type][i][Variable]['Statistics'][Statistic][index]) 
                    if Type=='Floater': #### Quick fix need to adjust the post processing to global pos
                        X.append(ProData[Type][i]['General']['Global Pos Connection'][0])
                        Y.append(ProData[Type][i]['General']['Global Pos Connection'][1])
                    else:
                        X.append(ProData[Type][i]['General']['Global Pos'][0])
                        Y.append(ProData[Type][i]['General']['Global Pos'][1])
                        
                    
                except:
                    pass
        
        Axis=['X','Y', 'Z', 'rx','ry','rz', 'horz','bend moment']
        
        if Statistic.find('Instant')!=-1:
            Statistic=Statistic+'_'+ str(round(ProData['Environment']['ExtractionTime'],1))+'s'
        
        if index==None:
            label=Type+ ' ' + Variable + ' ' + Statistic +' ' + Unit
            FigName=Type + '_'+ Variable + '_' + Statistic + '_' + os.path.basename(Sim)[:-16]
        else:
            label=Type+ ' ' + Variable + ' ' +Axis[index]+' '  + Statistic +' ' + Unit
            FigName=Type +'_'+ Variable + '_'+ Axis[index]+'_'+ Statistic + '_'+ os.path.basename(Sim)[:-16]
                    
        plt.scatter(X,Y,s=self.SizeMarker,c=Out, label=label)
        FigTitle=FigName.split('_Wdp')[0]
        plt.title(FigTitle)
        
        self.set_axis(X,Y, ProData)
        self.save_fig(Type, Sim, FigName)
                       
    def gap(self, ProData, Sim, Type, Statistic):
        ind=['PERH', 'PAR', 'PERV'].index(Type)
        self.plot_background(ProData)            
        X=[]
        Y=[]
        Out=[]
        for i in ProData['Gap'].keys():
            X.append(ProData['Gap'][i]['General']['Global Plot Pos'][0])
            Y.append(ProData['Gap'][i]['General']['Global Plot Pos'][1])
            Out.append(ProData['Gap'][i]['Statistics'][Statistic][ind])
            
        label='Gap ' + Type + ' ' + Statistic +' [m]'
        plt.scatter(X,Y,s=self.SizeMarker,c=Out, label=label)
       
        FigTitle='Gap '+ Type + ' ' + Statistic + ' ' + os.path.basename(Sim)[:-14]
        plt.title( FigTitle)
        
        self.set_axis(X,Y, ProData)
        self.save_fig('Gap', Sim, FigTitle)
        
    def mooring(self, ProData, Sim, Statistic):
        self.plot_background(ProData)            
        X=[]
        Y=[]
        Out=[]
        for i in ProData['Mooring'].keys():
            X.append(ProData['Mooring'][i]['General']['Global Plot Pos'][0])
            Y.append(ProData['Mooring'][i]['General']['Global Plot Pos'][1])
            Out.append(ProData['Mooring'][i]['Statistics'][Statistic])
            
        label='Mooring ' + Statistic +' [kN]'
        plt.scatter(X,Y,s=self.SizeMarker,c=Out, label=label)
       
        FigTitle='Mooring ' + Statistic + ' ' + os.path.basename(Sim)[:-14]
        plt.title( FigTitle)
        
        self.set_axis(X,Y, ProData)
        self.save_fig('Mooring', Sim, FigTitle)
        
    def plot_scatter(self, df, Variable, Statistic, function, Criteria=None,  Tag='', GeneralSelector=None, GeneralSelection=None, dfIreg=None, Axes=None):
        # plt.close()
        
        if self.PlotRAO:
            function=np.max
            Tag='RAO '+Tag
        
        if self.ScatXaxis=='WaveFrequency':
            df['WaveFrequency']=2*np.pi/df['WavePeriod']
        
        OutSummary=[Tag, Variable, Statistic,str(function)[11:15],  str(GeneralSelector), str(GeneralSelection), self.OveralSelector, self.OveralSelection]
        OutHeaders=['Description','Variable', 'Statistic', 'PivotTable Function',  'GeneralSelector','GeneralSelection', 'OveralSelector','OveralSelection']
        
        if self.OveralSelector and GeneralSelection:        
            GeneralSelector.extend(self.OveralSelector)
            GeneralSelection.extend(self.OveralSelection)
        elif self.OveralSelector and not GeneralSelection :
            GeneralSelector=self.OveralSelector
            GeneralSelection=self.OveralSelection
        else:
            pass
        
        df=df[df['Variable']==Variable]     
        if GeneralSelector:                         
            df=self.filterdf(df, GeneralSelector, GeneralSelection)
                
        # Depth=df['WaterDepth'].unique()
        UnitY=self.set_unit(Variable)
        UnitX=self.set_unit(self.ScatXaxis)
        
        if Axes:
            fig=Axes[0]
            ax=Axes[1]
        else:
            fig,ax=plt.subplots()
            
        if self.HOverScatter:
            ax2=ax.twinx()
        
        if Criteria: 
            Criteria=self.dfULSCriteria[self.dfULSCriteria.ULSCriteria==Criteria]
            # print(Criteria)
            if Criteria.Value.values == 'None':
                CritValue = None
            else:
                if Criteria.Value.values == 'Previous ULS':
                    if self.PreviousULSCrit:
                        Criteria = self.dfPreviousULS
                        dfc1 = Criteria[Criteria['Description']==Tag]
                        dfc2 = dfc1[dfc1['Variable']==Variable]
                        dfc3 = dfc2[dfc2['Statistic']==Statistic]
                        if function==np.max:
                            CritValue=max(dfc3[self.PreviousULSCritTag+'-MG'].values,dfc3[self.PreviousULSCritTag].values)
                            # CritValue=max(dfc3['41t-T1Inv-MG'].values,dfc3['41t-T1Inv'].values)
                        else:
                            CritValue=min(dfc3[self.PreviousULSCritTag+'-MG'].values,dfc3[self.PreviousULSCritTag].values)
                            # CritValue=min(dfc3['41t-T1Inv-MG'].values,dfc3['41t-T1Inv'].values)     
                    else: 
                        CritValue=None
                else:
                    CritValue = Criteria.Value
                    CritValue = float(CritValue)
            # Criteria = self.select_criteria(Tag,Variable,Statistic,dfCriteria)
        else: 
            CritValue=None
        
        if CritValue:
            ax.axhline(y = CritValue, color = 'r', linestyle = '--')
            
        LegendU=df[self.ScatLegend].unique()      
        
        c=0
        WD=0
        
        for i, Leg in enumerate(LegendU): 
            if i in self.LegendIndex:  
                df1=df[df[self.ScatLegend]==Leg]  
                # if self.Hex2Rows:
                #     if Leg == 'SD230t-MG-FltMG0.1-Gap-1.0':
                #         Leg = 'Hex2'
                #     elif Leg == 'SD230t-MG-ResetMooring':
                #         Leg = 'Layout1-Row2'
                #     elif Leg == 'SD230t-MG-ReOrientModel':
                #         Leg = 'Layout2-Row3'
                        
                if self.PlotLocalMax:
                    self.ProDataReference=pickle.load(open(os.path.join(self.LocSims, df1['SimName'].iloc[0]+'_Pro.pkl'),"rb"))
                    self.plot_max_overal(df1,Variable,Statistic,function, Tag, GeneralSelector=GeneralSelector, GeneralSelection=GeneralSelection, Directory=os.path.join(self.ScatLegend+self.ScatXaxis +'Plot',Leg))            
                
                if self.PlotRAO:
                    #### Plot RAO is always about the amplitude extreme(max-min)/2
                    #### RAO plots therefore should therefore be added to the mean value
                    # df1[Statistic]=2*df1[Statistic]/df1['HMax']  
                    
                    if Statistic.find('Min')==0:
                        StatSplit=Statistic.split('Min')                        
                    elif Statistic.find('Max')==0:
                        StatSplit=Statistic.split('Max')
                    elif Statistic.find('AbsMax')==0:
                        StatSplit=Statistic.split('AbsMax')
                    StatMax='Max'+StatSplit[1]
                    StatMin='Min'+StatSplit[1]
                    df1[Statistic]=2*((df1[StatMax]-df1[StatMin])/2)/df1['HMax']  
                    
                df2=pd.pivot_table(df1, values=[Statistic], index=self.ScatXaxis, aggfunc=function)
                    
                # dfSimName=pd.pivot_table(df1, values=[Statistic,'SimName'], index=self.ScatLegend, aggfunc=function)
                if function==np.max:
                    ind=df1[Statistic].idxmax()
                elif function==np.min:
                    ind=df1[Statistic].idxmin()
   
                if self.LegendLabel:
                    p=ax.plot(df2[Statistic],linestyle=self.LineStyle,color=self.CMap(c),label=str(i)+' ' + self.ScatLegend + ' '+ str(self.LegendLabel[i]))    
                else:
                    p=ax.plot(df2[Statistic],linestyle=self.LineStyle, color=self.CMap(c), label=str(i)+' ' + self.ScatLegend + ' '+ str(Leg))
                    OutSummary.append(float(function(df2, axis=0)))
                    OutHeaders.append(Leg)
                OutSummary.append(df1['SimName'][ind])
                OutHeaders.append('ExtremeSimName')
                c=c+1
                #### check criteria
                # OutSummary.append(CritValue)
                if CritValue is not None:
                    OutSummary.append(CritValue)
                    if function==np.max:
                        if float(function(df2, axis=0)) <= CritValue:
                            OutSummary.append('Yes')
                        else:
                            OutSummary.append('No')
                    else:
                        if float(function(df2, axis=0)) >= CritValue:
                            OutSummary.append('Yes')
                        else:
                            OutSummary.append('No')
                else:
                    OutSummary.append('None')
                    OutSummary.append('-')
                OutHeaders.append('Criteria Value')
                OutHeaders.append('Pass [Yes / No]')
                
                if self.DirectionMarkers:
                    for j in df1['WaveDirection'].unique():
                        df2=df1[df1['WaveDirection']==j]
                        df3=pd.pivot_table(df2, values=Statistic, index=self.ScatXaxis, aggfunc=function)
                        ax.scatter(df3.index,df3[Statistic], marker=(2,2,270+j), color=p[0].get_color())
                
                #OutSummary.append(df1['SimName'][ind])
                #OutHeaders.append('ExtremeSimName')
                #c=c+1
                
                # if i==self.HScatterIndex:'
                WD1=df2['WaterDepth'].unique()[0]
                
                if self.HOverScatter and WD!=WD1:   
                    df3=pd.pivot_table(df1, values=self.HsOrHMax, index=self.ScatXaxis, aggfunc=np.max)        
                    ax2.plot(df3,'^--', label=self.HsOrHMax +str(df2['WaterDepth'].unique())+'m', color='gray')
                    WD=df2['WaterDepth'].unique()[0]
        
        ax.set_xlabel(self.ScatXaxis + ' '+ UnitX)
        ax.set_ylabel(Tag+' '+ Statistic +' ' + Variable+ ' '+ UnitY)
        ax2.set_ylabel(self.HsOrHMax +'[m]')
        ax2.legend()
        
        OutSummary.append(UnitY)
        OutHeaders.append('Unit')                            
        Depth=df1['WaterDepth'].unique()
        
        if self.PlotRAO:
            title= Tag +  ' '+  Statistic + 'Variable'  #+' Dpth'+ str(Depth)+ 'm' 
        else:
            title= Tag +  ' '+  Variable + ' ' + Statistic + ' '+ str(function)[10:15] #+ ' Dpth'+ str(Depth)+ 'm' 
        ax.set_title(title)
        ax.grid()
                                
        if dfIreg is not None:
            df=dfIreg[dfIreg['Variable']==Variable] 
            if GeneralSelector: 
                df=self.filterdf(df, GeneralSelector, GeneralSelection)
                
            Depth1=df['WaterDepth'].unique()
            if Depth1!=Depth:
                print('Depths regular and irregular waves to not match')
            if len(df[self.ScatLegend].unique())!=1:
                print('Two sets of irregular waves plotted, thus plot incorrect')
                
            for i, Leg in enumerate(LegendU): 
                if i in self.LegendIndex:
                    df1=df[df[self.ScatLegend]==Leg]  
                    if self.PlotLocalMax:
                        self.plot_max_overal(df1,Variable,Statistic,function, Tag, GeneralSelector=GeneralSelector, GeneralSelection=GeneralSelection, Directory=os.path.join(self.ScatLegend+self.ScatXaxis +'Plot',Leg))
            
                    for j in df1['WaveDirection'].unique():
                        df2=df1[df1['WaveDirection']==j]
                        df3=pd.pivot_table(df2, values=Statistic, index=['Tp', 'Hs'], columns='sd', aggfunc=function)   
                        for it in range(len(df3.index)):
                            self.Env.Tp=df3.index[it][0]
                            self.Env.Hs=df3.index[it][1]
                            self.Env.calc_Hmax_rayleigh()
                            # self.Env.calc_Tass()
                            self.Env.Tass=self.Env.Tp*0.9
                            ax.scatter(np.repeat(self.Env.Tass, len(df3.iloc[it])), df3.iloc[it], marker=(2,2,270+j), color='black')
                            ax.scatter(self.Env.Tass, df3.iloc[it].mean(), marker=(2,2,270+j), color='red')
                    ax.plot(self.Env.Tass,df3.iloc[it].mean(),'k',label='IrregSeeds '+str(Leg))
                    ax.plot(self.Env.Tass,df3.iloc[it].mean(),'red',label='IrregMean '+ str(Leg))
                
        ax.legend()       
        if self.PlotRAO:
            self.save_fig(self.ScatLegend+self.ScatXaxis +'_RAO_Plot',  title)
        if dfIreg is not None:
            self.save_fig(os.path.join(self.LocIrregular, self.ScatLegend+self.ScatXaxis +'Plot'),  title)
        else:
            self.save_fig(self.ScatLegend+self.ScatXaxis +'Plot',  title)
        self.OutputSummary.append(OutSummary)
        self.OutputHeaders.append(OutHeaders)
        

        return [fig, ax]
    
    def plot_max_overal(self, df, Variable, Statistic, aggfunc,  Tag='',GeneralSelector=None, GeneralSelection=None, Directory='MaxOveralPlots'):
        #### This function makes a localized plot a given varible for a given dataframe
        #### It plots the maximum in the data frame
        #### Filter DataFrame
        if GeneralSelector:
            df=self.filterdf(df, GeneralSelector, GeneralSelection)
        
        #### Plot background
        # fig1=plt.figure('1')
        self.plot_background(self.ProDataReference)  
        ####
        if aggfunc==np.min:
            groupidx=df.groupby('UniqueName')[Statistic].idxmin()
        elif aggfunc==np.max:
            groupidx=df.groupby('UniqueName')[Statistic].idxmax()
    
        dfPlot=df.loc[groupidx] 
        X=dfPlot['Global Pos X']           
        Y=dfPlot['Global Pos Y']           
        Z=dfPlot[Statistic]
        
        Unit=self.set_unit(Variable)
        Title=Tag+ ' ' + Variable + ' ' + Statistic + ' '+ str(aggfunc)[11:15]+' ' + Unit
        FigTitle= Tag+ ' ' + Variable + ' ' + Statistic + ' '+ str(aggfunc)[11:15]
        
        if GeneralSelector: 
            plt.scatter(X,Y,s=self.SizeMarker,c=Z, label='Overal max or min ' + Unit)
        else: 
            plt.scatter(X,Y,s=self.SizeMarker,c=Z, label='Overal max or min ' + Unit + ' ' + str(GeneralSelector) + str(GeneralSelection))
        
        plt.title(Title)
        
        dfEnvDir=df.groupby('WaveDirection')[self.HsOrHMax].max()
        
        #### flip axis to match the overview to the design axis
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        self.set_axis(X,Y, self.ProDataReference, dfEnvDir=dfEnvDir)
        self.save_fig(Directory, '0LctnPlt '+ FigTitle)
        dfPlot.to_excel(os.path.join(self.LocSims, Directory, '01Lctn'+FigTitle+'.xlsx'))
        # plt.close(fig1)
        plt.close()
        
class StucturalLoadSelection(GeneralFunctions): 
        
    def __init__(self):
        self.a=1
        self.LCCounter=0 #### This has been set as an easy way to give loadcase numbers 

    def rotate_coordinates_around_z(self, Pos, Phi,):
        #Input phi is in degrees
        Phi=Phi*np.pi/180
        Pos=np.matrix(Pos)
        RotMat=np.matrix([[np.cos(Phi), -np.sin(Phi), 0], [np.sin(Phi), np.cos(Phi),0],[0,0,1]])
        NewPos=np.squeeze(np.asarray(RotMat*Pos.transpose()))
        NewPos=NewPos.transpose()
        return NewPos      
    
    # This first function locates platform and attribute with highst load from a series of simulations      
    # def select_max_loads_T1_T2(self, df, Statistic, LoadCaseName, MaxOrMin='Max'):
    #     dfOut=pd.DataFrame()
    #     c=0
    #     for iPT, PlatType in enumerate([1,2]):
    #         try:
    #             df1=df[df['PlatformType']==PlatType]
    #             if MaxOrMin=='Max':
    #                 ind=df1[Statistic].idxmax()
    #             elif MaxOrMin=='Min':
    #                 ind=df1[Statistic].idxmin()
    #             dfOut=pd.concat([dfOut,pd.DataFrame(df1.loc[ind,:]).T],ignore_index=True)
    #             dfOut.loc[c,'LoadCase']=LoadCaseName
    #             c=c+1
    #         except:
    #             print(LoadCaseName+ ' not found for T' + str(PlatType))
    #     dfOut['LC no']=self.LCCounter
    #     self.LCCounter=self.LCCounter+1
    #     return dfOut
    
    def select_max_loads_T1_T2(self, df, Statistic, LoadCaseName):
        dfOut=pd.DataFrame()
        for iPT, PlatType in enumerate([1,2]):
            try:
                df1=df[df['PlatformType']==PlatType]
                if Statistic.find('AbsMax')!=-1:
                    try: 
                        IndexExtra=Statistic[6:]
                    except:
                        IndexExtra=''
                    dfOut1=self.select_df_absmax(df1, index_extra=IndexExtra, LC=LoadCaseName)
                elif Statistic.find('Max')!=-1:
                    try: 
                        IndexExtra=Statistic[3:]
                    except:
                        IndexExtra=''
                    dfOut1=self.select_df_max(df1, index_extra=IndexExtra, LC=LoadCaseName)
                elif Statistic.find('Min')!=-1:
                    try: 
                        IndexExtra=Statistic[3:]
                    except:
                        IndexExtra=''
                    dfOut1=self.select_df_min(df1, index_extra=IndexExtra, LC=LoadCaseName)
                dfOut=pd.concat([dfOut,dfOut1])
            except:
                print(LoadCaseName+ ' not found for T' + str(PlatType)) 
        dfOut['LC no']=self.LCCounter
        self.LCCounter=self.LCCounter+1    
        return dfOut
    
    def select_max_floater_draft(self, df, Statistic, LoadCaseName):
        dfOut=pd.DataFrame()
        for iPT, PlatType in enumerate([1,2]):
            try:
                df1=df[df['PlatformType']==PlatType]
                if Statistic.find('AbsMax')!=-1:
                    try: 
                        IndexExtra=Statistic[6:]
                    except:
                        IndexExtra=''
                    dfOut1=self.select_df_absmax(df1, index_extra=IndexExtra, LC=LoadCaseName)
                elif Statistic.find('Max')!=-1:
                    try: 
                        IndexExtra=Statistic[3:]
                    except:
                        IndexExtra=''
                    dfOut1=self.select_df_max(df1, index_extra=IndexExtra, LC=LoadCaseName)
                elif Statistic.find('Min')!=-1:
                    try: 
                        IndexExtra=Statistic[3:]
                    except:
                        IndexExtra=''
                    dfOut1=self.select_df_min(df1, index_extra=IndexExtra, LC=LoadCaseName)
                dfOut=pd.concat([dfOut,dfOut1])
            except:
                print(LoadCaseName+ ' not found for T' + str(PlatType)) 
        dfOut['LC no']=self.LCCounter
        self.LCCounter=self.LCCounter+1    
        return dfOut
    
    def set_lc_No(self, df):
        df.reset_index( drop=True, inplace=True)
        df.reset_index( inplace=True)
        df = df.rename(columns={'index': 'LC no'})
        return df
    
    def select_df_max(self, df1, index_extra='', LC=None):
        #### selec max loads of data frame from T1 and T2 
        dfOut=pd.DataFrame()
        ind=df1['Max'+index_extra].idxmax()
        dfOut=pd.concat([dfOut,pd.DataFrame(df1.loc[ind,:]).T],ignore_index=True)
        dfOut['DOF Ind']=self.find_loadcase_index(index_extra)
        dfOut['MaxMin']='Max'
        if LC:
            dfOut['LoadCase']=LC
        return dfOut
    
    def select_df_min(self, df1, index_extra='', LC=None):
        #### selec max loads of data frame from T1 and T2 
        dfOut=pd.DataFrame()
        ind=df1['Min'+index_extra].idxmin()      
        dfOut=pd.concat([dfOut,pd.DataFrame(df1.loc[ind,:]).T],ignore_index=True)
        dfOut['DOF Ind']=self.find_loadcase_index(index_extra)
        dfOut['MaxMin']='Min'
        if LC:
            dfOut['LoadCase']=LC  
        return dfOut
    
    def select_df_absmax(self, df1, index_extra='', LC=None):
        #### selec max loads of data frame from T1 and T2 
        dfOut=pd.DataFrame()
        ind=df1['AbsMax'+index_extra].idxmax()
        dfOut=pd.concat([dfOut,pd.DataFrame(df1.loc[ind,:]).T],ignore_index=True)
        dfOut['DOF Ind']=self.find_loadcase_index(index_extra)
        dfOut['MaxMin']='AbsMax'
        if LC:
            dfOut['LoadCase']=LC
        return dfOut
      
    def select_sim_maxima(self, df1, index_extra=''):
        #### slect the maximum df line from a results df
        #### df selection should be done beforehand
        
        Sims=[]
        FailComponent1=[]
        # FailComponent2=[]
        dfFailing=pd.DataFrame()
        
        for Sim in df1['SimName'].unique():
            df2=self.filterdf(df1, ['SimName'], [Sim])
            Sims.append(Sim)
            ind=df2['AbsMax'+index_extra].idxmax()
            # ind2=df2['AbsMax'+index_extra].nlargest(2).index[-1]
            
            FailComponent1.append(df2['UniqueName'][ind])
            # FailComponent2.append(df2['UniqueName'][ind])
            dfFailing=pd.concat([dfFailing, pd.DataFrame(df1.loc[ind,:]).T])
            # dfFailing=pd.concat([dfFailing, pd.DataFrame(df1.loc[ind2,:]).T])
        
        return dfFailing
      
    def find_loadcase_index(self,LoadCase):
        ind_mapping = {
            'X': 0,
            'Y': 1,
            'Z': 2,
            'rX':3,
            'rY':4,
            'rZ':5,
            'Horz':6,
            'Total':6,
            'BendMoment':7,
            '':0,}
        
        for key, value in ind_mapping.items():
            if key in LoadCase:
                return value
            
    def get_time_instance_with_dof(self, df, iS, RawData, TimeTrace):
        if df.loc[iS,'MaxMin'].find('AbsMax')!=-1: 
            ind=(np.argmax(np.abs(TimeTrace[:,df.loc[iS,'DOF Ind']]))) 
            # ind=(np.argmax(np.abs(TimeTrace[df.loc[iS,'DOF Ind']]))) 
        elif df.loc[iS,'MaxMin'].find('Max') != -1:
            ind=(np.argmax(TimeTrace[:,df.loc[iS,'DOF Ind']])) 
        elif df.loc[iS,'MaxMin'].find('Min') != -1:
            ind=(np.argmin(TimeTrace[:,df.loc[iS,'DOF Ind']])) 
        return RawData['Environment']['time'][ind]

    ## timestamp functions locate timestamps of max loads 
    def timestamp_platformvar(self, df, iS, RawData):
        Probe=df.loc[iS,'UniqueName'].split('_')[-1]
        TimeTrace=RawData['Platform'][df.loc[iS,'Platform']][df.loc[iS,'Variable']][Probe]['TimeHistory']
        return self.get_time_instance_with_dof(df, iS, RawData, TimeTrace)           

    def timestamp_beam_loads(self, df, iS, RawData):
        TimeTrace=RawData['Beams'][df.loc[iS,'UniqueName']][df.loc[iS,'Variable']]
        return self.get_time_instance_with_dof(df, iS, RawData, TimeTrace)
    
    def timestamp_con_loads(self, df, iS, RawData):
        TimeTrace=RawData['Cons'][df.loc[iS,'UniqueName']][df.loc[iS,'Variable']]
        return self.get_time_instance_with_dof(df, iS, RawData, TimeTrace)
    
    def timestamp_floater_loads(self, df, iS, RawData):
        TimeTrace=RawData['Floater'][df.loc[iS,'UniqueName']]['ConnectionForce']
        return self.get_time_instance_with_dof(df, iS, RawData, TimeTrace)   
    
    def timestamp_floater_draft(self,df,iS,RawData):
        TimeTrace=RawData['Floater'][df.loc[iS,'UniqueName'][6:]]['AirGap']['TimeHistory']
        return RawData['Environment']['time'][np.argmin(TimeTrace)]
    
    ## extract fucntions extract desired data for the desired platform
    
    def extract_headers_and_lc(self, dfSims, dfInstant, iS, Headers):
        ind=[]
        for i in dfInstant['Platform']:
            if i.split('_')[0]==dfSims.loc[iS,'Platform'].split('_')[0]:
                ind.append(True)
            else: 
                ind.append(False)
        # for i in dfInstant['Platform']:
            
        # df1=dfInstant[dfInstant['Platform']==dfSims.loc[iS,'Platform']]
        df1=dfInstant[ind]
        df2=df1.loc[:,Headers]
        df2['LoadCase']= dfSims.loc[iS,'LoadCase']
        df2['LC no']= dfSims.loc[iS,'LC no']
        df2['MaxMin']= dfSims.loc[iS,'MaxMin']
        df2['DOF Ind']= dfSims.loc[iS,'DOF Ind']
        return df2
            
    def extract_con_data(self, df, dfCons, iS):
        ## Extract the desired data at the desired platform
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime', 'PlatformType','Max X', 'Max Y',   'Max rX', 'Max rY', 'Max rZ','Max Z', 'Max Total',
                 'Min X', 'Min Y', 'Min Z', 'Min rX', 'Min rY', 'Min rZ', 'Min Total',  'InstantVal X','InstantVal Y',
                 'InstantVal Z', 'InstantVal rX','InstantVal rY','InstantVal rZ','InstantVal Total', 'Local Pos X','Local Pos Y','Local Pos Z','FEM Label']
        dfOut=self.extract_headers_and_lc(df, dfCons, iS, Headers)
        return dfOut
    
    def extract_floater_data(self, df, dfFloater, iS):
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime','PlatformType', 'InstantVal X','InstantVal Y','InstantVal Z',
                  'InstantVal rX','InstantVal rY','InstantVal rZ','InstantValHorz', 'InstantValBendMoment', 'Max X', 'Max Y', 'Max Z',
                   'Max rX', 'Max rY', 'Max rZ', 'MaxHorz','MaxBendMoment', 'Min X', 'Min Y', 'Min Z',
                   'Min rX', 'Min rY', 'Min rZ', 'MinHorz', 'MinBendMoment','Local Pos X','Local Pos Y','Local Pos Z', 'FEM Label']           
        dfOut=self.extract_headers_and_lc(df, dfFloater, iS, Headers)
        return dfOut
    
    def extract_platform_data(self, df, dfPlatform, iS):
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime', 'PlatformType', 'InstantVal X','InstantVal Y','InstantVal Z',
                  'InstantVal rX','InstantVal rY','InstantVal rZ', 'Max X', 'Max Y', 'Max Z',
                   'Max rX', 'Max rY', 'Max rZ', 'Min X', 'Min Y', 'Min Z',
                   'Min rX', 'Min rY', 'Min rZ', 'Local Pos X','Local Pos Y','Local Pos Z']
        dfOut=self.extract_headers_and_lc(df, dfPlatform, iS, Headers)
        return dfOut
    
    def extract_coupling_data(self, df, dfPlatform, iS):
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime','Type', 'PlatformType', 'InstantVal', 
                   'Max', 'Min',  'UnstretchedLength', 'Stiffness', 'LinDamping', 'Local Pos X',	
                   'Local Pos Y', 'Local Pos Z', 'Connection A', 'Connection B']
        dfOut=self.extract_headers_and_lc(df, dfPlatform, iS, Headers)
        return dfOut
    
    def extract_coupling_data_2(self, df, dfInstant, iS):
        #### This is to extraxt data from the adjacted platform of interest
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime','Type', 'PlatformType', 'InstantVal', 
                   'Max', 'Min',  'UnstretchedLength', 'Stiffness', 'LinDamping', 'Local Pos X',	
                   'Local Pos Y', 'Local Pos Z', 'Connection A', 'Connection B']
        ind=[]
        for i in dfInstant['Connection B']:
            if i==df.loc[iS,'UniqueName']:
                ind.append(True)
            else: 
                ind.append(False)
        # for i in dfInstant['Platform']:
            
        # df1=dfInstant[dfInstant['Platform']==dfSims.loc[iS,'Platform']]
        df1=dfInstant[ind]
        df2=df1.loc[:,Headers]
        df2['LoadCase']= df.loc[iS,'LoadCase']
        df2['LC no']= df.loc[iS,'LC no']
        return df2
    
    def extract_beam_data(self, df, dfPlatform, iS):
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime', 'PlatformType', 'InstantVal', 
                   'Max', 'Min', 'Local Pos X',	'Local Pos Y', 'Local Pos Z']
        dfOut=self.extract_headers_and_lc(df, dfPlatform, iS, Headers)
        return dfOut
    
    
    def extract_airgap_data(self, df, dfPlatform, iS):
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime','PlatformType', 'ProbeName', 'InstantVal', 
                    'Max', 'Min',  'Local Pos X',	
                    'Local Pos Y', 'Local Pos Z']
        dfOut=self.extract_headers_and_lc(df, dfPlatform, iS, Headers)
        return dfOut

    def extract_relangle_data(self, df, dfInstant, iS,dfOut):
        #### This function only works after extracting the cylinder loads from extract_coupling_data_2
        Headers=['SimName', 'UniqueName', 'Variable', 'ExtractionTime',  'InstantVal', 
                    'Max', 'Min',  'Connection A', 'Connection B']
        ind=[]
        PlatAOut=df.loc[iS, 'Platform'].split('_')[0]
        PlatBOut=list(dfOut['UniqueName'])[-1].split('_')[1]
        PlatOut=sorted([PlatAOut,PlatBOut])
        for i, ConA in enumerate(dfInstant['Connection A']):
            PlatA=list(dfInstant['Connection A'])[i].split('_')[0]
            PlatB=list(dfInstant['Connection B'])[i].split('_')[0]
            Plat=sorted([PlatA,PlatB])
            if Plat==PlatOut:
                ind.append(True)
            else: 
                ind.append(False)
        df1=dfInstant[ind]
        df2=df1.loc[:,Headers]
        df2['LoadCase']= df.loc[iS,'LoadCase']
        df2['LC no']= df.loc[iS,'LC no']
        
        return df2

class Frequency_domain_anaylsis:
    #### This class is specifically for FD analysis
    #### This is for creating RAOs, reading in scatter, and doing FD calulations based on the dict structure of OF objects
    #### This varies from the fatigue analysis functions as they work with  a different stucture that accomodates stress nodes
    
    def __init__(self, LocSims, SLS):
        self.LocSims=LocSims
        #### set of wave frequencies that shall be interpolated
        self.WavePeriods=np.concatenate([np.arange(2.5,14.5,0.25),np.arange(15,20,0.5), np.arange(20,30,3)])#, np.arange(20,30,3)
        self.WaveFrequencies=2*np.pi/self.WavePeriods
        self.filterdf=SLS.filterdf
                
    def read_scatter(self, ScatterPath):
        #### Reads a scatter text file, colum and row indeces can be used for optimization
        fileObject = open(ScatterPath, "r")
        data = pd.read_csv(fileObject,sep='\t',skiprows=(0,1,2), header=(1))
        
        cind1=2
        cind2=-11
        rind1=0
        rind2=-2
        
        Hss=data.iloc[rind1:rind2,0]
        Hss=[re.findall(r"[+]?\d*\.\d+|\d+",i) for i in Hss]
        Hss=np.asarray(Hss, dtype=float)
        Hss=Hss.mean(axis=1)
        # self.Hss=Hss
        
        Tps=np.array(data.columns[cind1:cind2])
        Tps=[re.findall(r"[+]?\d*\.\d+|\d+",i) for i in Tps]
        Tps=np.asarray(Tps, dtype=float)
        Tps=Tps.mean(axis=1)
        # self.Tps=Tps
        
        Scatter=np.array(data.iloc[rind1:rind2,cind1:cind2])
       
        Scatter0=Scatter.copy()
        Scatter0[Scatter=='-']=0
        Scatter0=np.asarray(Scatter0,dtype=float)
        
        Scatter025=Scatter.copy()
        Scatter025[Scatter025==0]=0.0025
        Scatter025[Scatter025=='-']=0
        Scatter025=np.asarray(Scatter025,dtype=float)
        
        ScatterDict={}
        ScatterDict['Hss']=Hss
        ScatterDict['Tps']=Tps
        ScatterDict['Scatter']=Scatter
        ScatterDict['Scatter0']=Scatter0
        ScatterDict['Scatter0025']=Scatter025
        ScatterDict['ScatterRaw']=data
        
        return ScatterDict
    
    def gen_spect_dict(self,ScatterDict):
        import spectral
        SpecDict={}
        for iTp, Tp in enumerate(ScatterDict['Tps']):
            for iHs, Hs in enumerate(ScatterDict['Hss']):
                #### below is an import statement, this takes up all sea-sates that could happen
                #### in the fatigue functions scatter 0 should as there calculation time is more critical
                if ScatterDict['Scatter'][iHs,iTp]!='-':
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]=spectral.Spectrum('JONSWAP', Hs, 'Tp', Tp, self.WaveFrequencies, gamma_val=3.3)
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['Hs']=Hs
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['Tp']=Tp
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['TpInd']=iTp
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['HsInd']=iHs
                    # plt.plot(SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['frequencies'], SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['density'])
        print('SpecDict Generated')
        SpecDict['General']={}
        for key in ScatterDict.keys():
            SpecDict['General'][key]=ScatterDict[key].copy()
        SpecDict['General']['WaveFrequencies']=self.WaveFrequencies
        return SpecDict    
    
    def extract_raos(self, df, Variable, StatisticAddon, dfName, OverWrite=True):
        #### This function extracts and save ROAs for a set of compoenent for a single variable and statistic    
        FileN=dfName+'_'+ Variable+'_'+StatisticAddon
        
        if os.path.isfile(os.path.join(self.LocSims,'RAOS_'+FileN+'.pkl')) and not OverWrite:
            print('RAOs Loaded')
            DictRAO=pickle.load(open(os.path.join(self.LocSims,'RAOS_'+FileN+'.pkl'),'rb'))
          
        else:
            if StatisticAddon.find('Total')!=-1 or StatisticAddon.find('Horz')!=-1 or StatisticAddon.find('BendMoment')!=-1:
                df.loc[:,'RAO']=(df['Max'+StatisticAddon].copy()*2)/df['HMax'].copy()
            else:
                df.loc[:,'RAO']=(df['Max'+StatisticAddon].copy()-df['Min'+StatisticAddon].copy())/df['HMax'].copy()
            DictRAO={}
            DictRAO['General']={}
            for Lim in [ '_1m','_BrLm']:#'_1m',
                if Lim=='_1m':
                    df1=df[df['HMax']==1] 
                elif Lim=='_BrLm':
                    df1=df[df['HMax']>1]
                print(Lim)
                               
                
                for UN in df['UniqueName'].unique():   
                    if Lim=='_1m':
                        DictRAO[UN]={}
        
                    for V in [Variable]:#df1['Variable'].unique():
                        if Lim=='_1m':
                            DictRAO[UN][V]={}
                        dfSpec=self.filterdf(df1,['UniqueName','Variable'],[UN,V])
                        RAO1m=[]
                        for WD in sorted(dfSpec['WaveDirection'].unique()):
                            Pivot=pd.pivot_table(dfSpec[dfSpec['WaveDirection']==WD], values=['RAO','HMax'], index=['WavePeriod'], aggfunc=np.mean)
                            ## average value of with and without current is taken in the pivot above
                            # RAO=np.interp(WavePeriods, np.array(Pivot.index),Pivot['RAO'])
                            RAO1m.append(np.interp(self.WavePeriods, np.array(Pivot.index),Pivot['RAO']))
                            #### test interpolation
                            # plt.plot(WavePeriods, RAO)
                            # plt.plot(Pivot['RAO'])
                        DictRAO[UN][V]['RAO'+ Lim]=np.array(RAO1m)  
                            
                    DictRAO['General']['WaveDirections'+Lim]=sorted(dfSpec['WaveDirection'].unique())
                    DictRAO['General']['WavePeriods_OrcaFlex'+Lim]=np.array(Pivot.index)
                    DictRAO['General']['WavePeriods'+Lim]=self.WavePeriods
                    DictRAO['General']['WaveFrequencies_OrcaFlex'+Lim]=2*np.pi/np.array(Pivot.index)
                    DictRAO['General']['WaveFrequencies'+Lim]=self.WaveFrequencies
                    DictRAO['General']['WaveHeights_OrcaFlex'+Lim]=np.array(Pivot['HMax'])
                    DictRAO['General']['WaveHeights'+Lim]=np.interp(self.WavePeriods, np.array(Pivot.index),Pivot['HMax'])
                    DictRAO['General']['FileName']=FileN
              
            pickle.dump(DictRAO,open(os.path.join(self.LocSims,'RAOS_'+FileN+'.pkl'),'wb'))
        return DictRAO
        
    def calc_responses(self, DictResp, SpecDict, SavePKLandExcel=True, Threshold=False):
        #### Dict RAO becomes dict response
        WaveDirections=DictResp['General']['WaveDirections_1m']
        Tps=SpecDict['General']['Tps']
        Hss=SpecDict['General']['Hss']
        WaveFrequencies=SpecDict['General']['WaveFrequencies']
        if Threshold:
            DictResp['General']['Threshold']=Threshold
        
        WH1m=DictResp['General']['WaveHeights_1m']
        WHBrLm=DictResp['General']['WaveHeights_BrLm']

        WHBrLmTpss=np.interp(0.9*np.array(Tps), self.WavePeriods, WHBrLm)

        SigTotal=np.zeros([len(Hss),len(Tps),len(WaveDirections),len(DictResp.keys())-1])
        TzTotal=np.zeros([len(Hss),len(Tps),len(WaveDirections),len(DictResp.keys())-1])
        NumPoTperHTotal=np.zeros([len(Hss),len(Tps),len(WaveDirections),len(DictResp.keys())-1])
        m0Total=np.zeros([len(Hss),len(Tps),len(WaveDirections),len(DictResp.keys())-1])

        for iUN, UN in enumerate(DictResp.keys()):
            if UN!='General':
                for V in DictResp[UN].keys():
                    DictResp[UN][V]['HsTp']={}
                    DictResp[UN][V]['SigSummary']=np.zeros([len(Hss),len(Tps), len(WaveDirections)])
                    DictResp[UN][V]['TzSummary']=np.zeros([len(Hss),len(Tps),len(WaveDirections)])
                    DictResp[UN][V]['m0Summary']=np.zeros([len(Hss),len(Tps),len(WaveDirections)])
                    DictResp[UN][V]['CyclesPerHourSumary']=np.zeros([len(Hss),len(Tps),len(WaveDirections)])
                    DictResp[UN][V]['PoTSummary']=np.zeros([len(Hss),len(Tps),len(WaveDirections)])
                    DictResp[UN][V]['NumPoTperHrSummary']=np.zeros([len(Hss),len(Tps),len(WaveDirections)])
                    for HsTp in SpecDict.keys():
                        if HsTp!='General':
                            ## Interp Hs is based on the ratio 2(Hs-H1m)/(HBrLm-H1m) H1m is hardcoded as 1
                            tmp=SpecDict[HsTp]
                            Ratio=(SpecDict[HsTp]['Hs']*2-1)/(WHBrLmTpss[SpecDict[HsTp]['TpInd']]-1)
                            if Ratio>1:
                                a=1
                            Ratio=np.max([Ratio,0])
                            Ratio=np.min([Ratio,1])
                            # Ratio=0
                            tmp['RAOintp']=(DictResp[UN][V]['RAO_BrLm']-DictResp[UN][V]['RAO_1m'])*Ratio+DictResp[UN][V]['RAO_1m']
                                                
                            tmp['RespDens']=tmp['RAOintp']**2*SpecDict[HsTp]['density']
                            
                            # if SpecDict[HsTp]['Tp']==10.5:
                            #     plt.figure(1)
                            #     plt.plot(WaveFrequencies,tmp['RAOintp'].transpose())
                            #     plt.plot(WavePeriods,tmp['RAOintp'].transpose())
                            #     plt.figure(2)
                            #     plt.plot(WaveFrequencies,SpecDict[HsTp]['density'])
                            #     plt.figure(3)
                            #     plt.plot(WaveFrequencies,tmp['RespDens'].transpose())
                            #     a=1
                            
                            tmp['m2']=np.trapz(tmp['RespDens']*WaveFrequencies**2, WaveFrequencies)*-1
                            tmp['m0']=np.trapz(tmp['RespDens'],WaveFrequencies)*-1
                            # tmp['m0']=np.trapz(np.flip(tmp['RespDens'],axis=1),np.flip(WaveFrequencies))
                            tmp['Tz']=2*np.pi*(tmp['m0']/tmp['m2'])**0.5
                            # tmp['CyclesPerHour']=3600/tmp['Tz']
                            tmp['Sig']=2*(tmp['m0']**0.5)
                            DictResp[UN][V]['HsTp'][HsTp]=tmp.copy()
                            if Threshold:
                                tmp['PoT']= np.exp(-Threshold**2/(2*tmp['m0'])) #### Probability over threshold
                                tmp['NumPoTperHr']=3600/tmp['Tz']*tmp['PoT']
                                DictResp[UN][V]['PoTSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]= tmp['PoT']
                                DictResp[UN][V]['NumPoTperHrSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['NumPoTperHr']
                    
                            DictResp[UN][V]['SigSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['Sig']
                            DictResp[UN][V]['TzSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['Tz']
                            DictResp[UN][V]['m0Summary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['m0']
                            # DictResp[UN][V]['CyclesPerHourSumary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['CyclesPerHour']
                                        
        #       #### This only works with one variable
                SigTotal[:,:,:,iUN-1]=DictResp[UN][V]['SigSummary']
                TzTotal[:,:,:,iUN-1]=DictResp[UN][V]['TzSummary']
                m0Total[:,:,:,iUN-1]=DictResp[UN][V]['m0Summary']
                
                if Threshold:
                    NumPoTperHTotal[:,:,:,iUN-1]=DictResp[UN][V]['NumPoTperHrSummary']
                    DictResp['General']['NumPoTperHTotal']=NumPoTperHTotal                 
           
                    
        DictResp['General']['SigTotal']=SigTotal
        DictResp['General']['TzTotal']=TzTotal
        DictResp['General']['m0Total']=m0Total
        
        if SavePKLandExcel == True: 
            self.save_signifcant_response_pickle_and_excel(DictResp, SpecDict, Threshold=Threshold)
        
        return DictResp
        
    def save_signifcant_response_pickle_and_excel(self, DictResp, SpecDict, SavePKL=True, SaveScatterXLS=True, Threshold=False):
        #### Max component max direction
        SigMax=np.max(np.max( DictResp['General']['SigTotal'],axis=3),axis=2)
        #### Mean component mean direction
        TzMean=np.mean(np.mean( DictResp['General']['TzTotal'],axis=3), axis=2)
        SigMean=np.mean(np.mean(DictResp['General']['SigTotal'],axis=3),axis=2)
                
        DictResp['General']['SigMax']=SigMax
        DictResp['General']['TzMean']=TzMean
        DictResp['General']['SigMean']=SigMean
        
        Tps=SpecDict['General']['Tps']
        Hss=SpecDict['General']['Hss']
        dfSigMax=pd.DataFrame(SigMax, index=Hss, columns=Tps)
        dfSigMean=pd.DataFrame(SigMean, index=Hss, columns=Tps)
        dfTz=pd.DataFrame(TzMean, index=Hss, columns=Tps)
        
        
        if SavePKL:
            pickle.dump(DictResp,open(os.path.join(self.LocSims,'Resp_'+DictResp['General']['FileName']+'.pkl'),'wb'))
            
        if SaveScatterXLS:
            TE=ToExcel()
            FileN='Scatter_'+DictResp['General']['FileName']+'.xlsx'
            TE.save_excel(SpecDict['General']['ScatterRaw'], os.path.join(self.LocSims, FileN), SheetName='InputScatter')
            TE.append_excel_sheet(dfTz, os.path.join(self.LocSims,FileN), SheetName='Tz-MeanCmpnnt-MeanDrc')
            TE.append_excel_sheet(dfSigMax, os.path.join(self.LocSims,FileN), SheetName='SigVal-MaxCmpnnt-MaxDrc')                   
            TE.append_excel_sheet(dfSigMean, os.path.join(self.LocSims,FileN), SheetName='SigVal-MeanCmpnnt-MeanDrc')    
            if Threshold:   
                NPoTMax=np.max(np.max(DictResp['General']['NumPoTperHTotal'],axis=3),axis=2)
                dfNPoTMax=pd.DataFrame(NPoTMax, index=Hss, columns=Tps) 
                TE.append_excel_sheet(dfNPoTMax, os.path.join(self.LocSims,FileN), SheetName='NPeakOvThresh-MaxComp-MaxDir')   
    
        return DictResp
        
    def calc_seastate_cycle_hist(self, DictResp, Duration):
        #### Calculates the load spectrum of a given ResponseDict
        from matplotlib.ticker import LogLocator
        
        #### Select the m0, Tz and create bins
        m0=DictResp['General']['m0Total']
        Tz=DictResp['General']['TzTotal']
        Sig=DictResp['General']['SigTotal']
        
        MaxLoad=np.max(DictResp['General']['SigTotal'])*2
        Bins=np.linspace(0,MaxLoad,100)
        dBin=Bins[1]-Bins[0]
        BinsCenter=Bins[:-1]+(dBin)/2
        
        
        Duration=Duration
        Cycles=Duration/Tz
        
        #### Calculate histogram for amount of cycles in each bin for each sea-state, component, direction
        pdf=BinsCenter[:,None,None,None,None]/m0*np.exp(-(BinsCenter[:,None,None,None,None]**2)/(2*m0))
        CycleHist=pdf*Cycles[None,:]/np.sum(pdf, axis=0)[None,:]
        CycleHist[np.isnan(CycleHist)] = 0
        
        return CycleHist, Bins, BinsCenter
        
    def calc_yearly_cycle_hist(self, DictResp, SpecDict):
        #### This function now assumes a mean direction        
        Duration=3600
        
        CycleHist, Bins, BinsCenter =self.calc_seastate_cycle_hist(DictResp,  Duration)
        
        # #### Start test with single m0 code ####
        # def rayleigh_pdf(x,m0):
        #     return x/m0*np.exp(-x**2/(2*m0))

        # P_n=rayleigh_pdf(BinsCenter,m0[7,10,3,3])
        # Cycles1D=Cycles[7,10,3,3]       
        # Hist=P_n/np.sum(P_n)*Cycles1D
        
        # Hs=2*(m0[7,10,3,3])**0.5
        
        # plt.plot(BinsCenter,Hist)
        # plt.plot(BinsCenter,CycleHist[:,7,10,3,3])
        # #### End test with single m0 code ####
        
        #### Calculate the number of cylces in each load bin by summing over the seastate scatter 
        Scatter=SpecDict['General']['Scatter0025'] #### This is the scatter that contains an occurance of 0.005 procent for all unknowns, this is to stay conservative
        Scatter=Scatter*1/np.sum(Scatter)
        
        CycleHistYear=CycleHist*Scatter[None,:,:,None,None]*365*24*3600/Duration  
        #### Now take the average cycles for the 5 directions, assumed each direction happens 1/5th of the time
        CycleHistYear1=np.sum(CycleHistYear,axis=3)/CycleHistYear.shape[3]
        #### Now add over all the sea-states to get the total load cylces spectrum per component
        CycleHistYearSum=np.sum(np.sum(CycleHistYear1, axis=1),axis=1)
        plt.semilogx(CycleHistYearSum,BinsCenter,'-x')
        plt.xlim(0.001,10**7)
        plt.legend(list(DictResp.keys())[1:])
        plt.grid()
        # plt.grid(True, which="minor", axis="x")
        # ax = plt.gca()
        # ax.grid(which='major')
        

        plt.ylabel('Load[kN]')
        plt.xlabel('Cycles[-]')
        plt.title('LoadCyclesPerYear')
        plt.ylabel('Load [kN]')
        plt.xlabel('Cycles per year [-]')
        plt.savefig(os.path.join(self.LocSims,'CyclePlot_'+DictResp['General']['FileName']+'png'))
        
        # plt.ylim(0,1000)
        # plt.savefig(os.path.join(self.LocSims,'CyclePlot1000Cycle_'+DictResp['General']['FileName']+'png'))
        
        df_hist = pd.DataFrame(CycleHistYearSum, columns=list(DictResp.keys())[1:] , index=[f'{Bins[i]:.2f} - {Bins[i+1]:.2f}' for i in range(len(Bins)-1)])
        df_hist.to_excel(os.path.join(self.LocSims,'LoadCycles_'+DictResp['General']['FileName']+'.xlsx'))
    
    
    
    
    
    
