# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 13:53:58 2022

@author: Sander
"""

import OrcFxAPI as ofx
import numpy as np

class Gen: 
    
    def __init__(self, model):
        self.model=model
        
    def environment(self, Env):
        self.model.WaterDepth=Env.Depth
        self.model.WavesReferredToBy='period (s)' 
        self.model.NumberOfPeriodsOrFrequencies=len(Env.OrcaWavePeriods)
        self.model.PeriodOrFrequency=Env.OrcaWavePeriods
        self.model.NumberOfWaveHeadings=len(Env.OrcaWaveHeadings)
        self.model.WaveHeading=Env.OrcaWaveHeadings
        
    def calculation_and_output(self):
        self.model.ValidatePanelArrangement='Yes'
        self.model.DivideNonPlanarPanels ='Yes'
        self.model.SolveType='Potential and source formulations'
        self.model.HasResonanceDampingLid ='No'
        
    def translate_inertia_matrix(self, Mass, CoG, I):
        #### Translates the interitia matrix to the origin
        # I0=np.matrix([[I[0],0,0],[0,I[1],0],[0,0,I[2]]])
        I0=np.matrix(I)
        Tmat=np.matrix([[CoG[1]**2+CoG[2]**2, -1*CoG[0]*CoG[1]   , -1*CoG[0]*CoG[2]],
                        [-1*CoG[0]*CoG[1]   , CoG[0]**2+CoG[2]**2, -1*CoG[1]*CoG[2]],
                        [-1*CoG[0]*CoG[2]   ,  -1*CoG[1]*CoG[2]  , CoG[0]**2+CoG[1]**2]])
        I1=I0+Mass*Tmat
        return I1
        
    
    def bodies_and_inertia(self, Layout, Floater, Platform, MeshFile, Constraints=[False,False,False,False,False,False]):
        PlatNames=[]
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0:       
                PlatNames.append('Platform'+str(i))
        self.model.NumberOfBodies=len(PlatNames)
        self.model.BodyName=PlatNames
        self.model.BodyMeshPositionX=Layout.PlatformsX[Layout.Types!=0]
        self.model.BodyMeshPositionY=Layout.PlatformsY[Layout.Types!=0]
        self.model.BodyMeshHeading=Layout.Orientation[Layout.Types!=0]
        
        DraftList=[]
        
        for PlNm in PlatNames:
            self.model.SelectedBody=PlNm
            self.model.BodyOrcaFlexImportSymmetry='None'
            self.model.BodyMeshFileName=MeshFile
            self.model.BodyOrcaFlexImportLength=Platform.Length
            self.model.BodyMeshFormat='Wamit gdf'
            self.model.BodyAddInteriorSurfacePanels='Yes'
            self.model.BodyInteriorSurfacePanelMethod='Triangulation method'
            
            self.model.BodyInertiaSpecifiedBy='Matrix'
            self.model.BodyCentreOfGravityx=Platform.CoG[0]
            self.model.BodyCentreOfGravityy=Platform.CoG[1]    
            
            i=1
            if np.mod(float(PlNm[8:]),2)==0: 
                i=-1
                
            self.model.BodyMass=Platform.Mass#+Platform.CouplingMassContribution*i
            # Draft=self.model.BodyMass/(1.025*Floater.AreaAxial*Floater.NumberOfFloaters)
            # DraftList.append(Draft)
            # self.model.BodyCentreOfGravityz=Platform.CoG[2]+Floater.Height-Draft
            
            self.model.BodyCentreOfGravityz=Platform.CoG[2]+Floater.Height-Floater.Draft
            CoG=np.concatenate([Platform.CoG[0:2], [ self.model.BodyCentreOfGravityz]])
            
            Ixx=Platform.Ixx#+Platform.CouplingInertiaContribution[0]*i
            Iyy=Platform.Iyy#+Platform.CouplingInertiaContribution[1]*i
            Izz=Platform.Izz#+Platform.CouplingInertiaContribution[2]*i
            I1 =self.translate_inertia_matrix(self.model.BodyMass,CoG, [Ixx, Iyy, Izz])
            self.model.BodyInertiaTensorRx[0]=I1[0,0]
            self.model.BodyInertiaTensorRy[:2]=I1[0:2,1]
            self.model.BodyInertiaTensorRz=I1[:,2]
            
            ### Setup constraints
            self.model.BodyFixedDOFx=Constraints[0]
            self.model.BodyFixedDOFy=Constraints[1]
            self.model.BodyFixedDOFz=Constraints[2]
            self.model.BodyFixedDOFRx=Constraints[3]
            self.model.BodyFixedDOFRy=Constraints[4]
            self.model.BodyFixedDOFRz=Constraints[5]
            
            
        # self.model.BodyMeshPositionZ=np.array(DraftList)*-1
        self.model.BodyMeshPositionZ=np.ones(len(PlatNames))*Floater.Draft*-1
        
    def field_points(self, xlist, ylist):
        self.model.NumberOfFieldPoints=len(xlist)
        self.model.FieldPointX=xlist
        self.model.FieldPointY=ylist

                
        
        