
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 14:22:51 2021

@author: Sander Morshuiss
"""
import sys
sys.path.insert(1, 'Classes')
import collections
try:
    from collections import abc
    collections.MutableMapping = abc.MutableMapping
except:
    pass
from openpyxl import load_workbook
import OrcFxAPI as ofx
import numpy as np
import numpy.matlib
import InputClasses as Input
import OrcaFlexConstructionClasses
import DataExtractionClasses
import os
# import copy
import pandas as pd
import copy
# import PostProcessingClasses as PPC
import multiprocessing
import glob
import matplotlib.pyplot as plt
import winsound
import re

def read_env(Path,ExcelName,SelectWaves=False):
    dfEnv=pd.read_csv(os.path.join(Path, ExcelName))
    if SelectWaves:
        dfEnv1 = pd.DataFrame()
        for i in dfEnv['WDir']:
            if i == 0 or i == 30 or i >= 210 and i<= 330: #### select high waves
                dfEnv1 = pd.concat([dfEnv1,dfEnv[dfEnv['WDir']==i]])
        dfEnv = dfEnv1.reset_index(drop=True)
    WDirs=dfEnv['WDir'][:]
    WDirs=np.array(WDirs.astype(float))
    return dfEnv,WDirs

if __name__ == '__main__': 
    #### the if statement is because of how multiprocessing works
    #### it is here to protect any code block from re-executing in a child process.
    #### before go model running, it is better to check everyhting right
    RunSims = False
    LargerFloterOnEdgePlf = False#This turn on option to have larger floater diameter of floaters on edge platforms
    EdgeFltD = 3 ## This only works when LargerFloterOnEdgePlf=True
    
    BridalMooring = True
    AdjustBridalCornerAnchorPos = True #### this turns on adjustment of bridal mooring outskirts

    EstiMass = False ### this turns on option to apply estimate mass to all platforms
    HeavierEdgePlatforms = False ### apply esti mass on the edged platforms
    BoxTruss = True
    
    
    FloaterDia = 3
    FltHPlus = 1 ### allow extra 2.5m airgap, bottom floater can be lower
    GapChange = 0 ### Platform gap change only applied on coupling center x position
    MoorLengthChange = 6 # Length adjustment for aged/MG line
    
    Threads = 15
    
    ##### below is the paramatric setting for models

    LocSims=r'Z:\01_Hydro\56_240228_BoxTruss'
    LayoutShape = [9,3]
    Update = 'ULS'
    FolderName = r'04_Moor_Rev1\Test' 
    
    #### Set up environment ####
    dfEnv,WDirs = read_env(r'Z:\01_Hydro\56_240228_BoxTruss','DirectionalHMax1.csv',SelectWaves=False)
    Hmax=np.max(dfEnv['50 yr rtp DHI'])
    # WDirsOF = (np.mod((360-WDirs)+90,360))
    WDirsOF = [0,15,30,45,60]
    # WDirsOF = np.arange(120,181,15)
    # Ts = [6.5,7,7.5,8,8.5,10]
    Ts = [6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.4]
    # Ts=[6.5,7,7.5,8.5 ,9.5,10.5, 11.7]
    #### Hs and Tp given for Agnes
    # Hs = 6.74
    # Tp = 9.73
    # Ts0 = Tp*0.9 ### 8.76s
    
    # Ts=[6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10]
    # Ts = np.arange(3,9.5,0.5)
    if not RunSims:
        WDirsOF = [WDirsOF[0]]
        Ts = [Ts[0]]
        
    LocCylinders=r"X:\01_Hydromechanics\Reference Info\05_CouplingSpringsEwoud\20221108_GasSpringConceptCalculationsV3_CylinderLeak.xlsx"
    wbCylinders=load_workbook(LocCylinders, data_only=True)
    
    #### Activate all the classes
    General, Environment, Platform, Floater, DragPlate, Coupling, Mooring, Layout = OrcaFlexConstructionClasses.Activate.all_input(Input)
    # Wing=Input.Wing(Platform)

    Environment.Depth=29.3
    IrregularEnv=False
    
    Alt=70
    
    LocSims = os.path.join(LocSims,FolderName)
    if not os.path.exists(LocSims):
        os.makedirs(LocSims) 
    # Hmax=12.42
    OutNum = 0 # for Model Check Summary
    for iMG, MG in enumerate([False]):
    
        # MG=False    
        print("MG = ",MG)
        if MG:
            MGEx='-MG'
            LengthChange = MoorLengthChange
        else:
            MGEx='-NMG'
            LengthChange = 0

        FloaterMG = MG
        #### Activate gen and compose class by allocating base model with line types    
        model=ofx.Model(Mooring.LocationLineTypes)
        Gen=OrcaFlexConstructionClasses.Gen(model)
        Compose=OrcaFlexConstructionClasses.Compose(model)  
            
        General.ImplicitConstantTimeStep=0.025 ### Made smaller to 0.025 for mooring anaylsis with lines 
        General.TargetLogSampleInterval=0.05
        
        #### First loop is over the platform type
        for iN, N in enumerate([LayoutShape]):
            Layout.calc_triangle(N[0],N[1]) 
            #### Rest some parameters for the loop
            Platform=Input.Platform()
            Platform.EdgeLength=24.784
            Platform.EdgeLengthOuter=36.9
            Platform.TrussZOffset=0
            Platform.TrussHeight = 3.4
            Platform.Height=Platform.TrussHeight+Platform.TrussZOffset
            Platform.calc_details_basic()
            Platform.SkidPlatform = 0 ###0: no skid considered
            Platform.BoxTruss = BoxTruss
            #### Set platform masses for mass checks
            Platform.MassT1=56.776
            Platform.MassT2=61.889
            Platform.BeamMass = 5.6+2.38+1.232 #### sum of L1L, L2, L3
            
            Platform.TrussLayer1=Input.MassObject(Platform.BeamMass*1,[Platform.GeoCoG[0],0, 1.7])
            Platform.TrussLayer2=Input.MassObject(Platform.BeamMass*0,[Platform.GeoCoG[0],0, 1.7])
            
            Platform.MassCornersTotal=1.6
            Platform.MassTopFloaterTotal=18.20*0.36#### This mass will be added to the corner mass, Ratio based on Merganser
            Platform.ZTopFloater=1.36
            Platform.MassCouplingsTotal_T2=17.155 + 1.5 #### bridges mass is added on coupling 
            Platform.MassCouplingsTotal_T1=12.042 + 1.5 #### *3/2 was done as Henk only included 2 set of coupling, in parametric code 6 are taken into account
            Platform.MassPanelsAndGrating=14.222 
            Platform.PanelsAndGrating=Input.MassObject(Platform.MassPanelsAndGrating,[Platform.GeoCoG[0],0, 4.177])
            # Platform.PanelsAndGrating.set_filled_triangle_interia(Platform.EdgeLengthOuter)
            Platform.MassNodesTotal=0
            Platform.AdjustMassOfOuterPlatforms=True
            Platform.MassFairlead=(0.611)/2 ### +200 is an estimate for the chains

            #### Section used define masses for use of ridid bodies, masses are copied from the compound properties in OrcaFlex
            # dfW = pd.read_excel(os.path.join(LocSims, "OST1-SYS-003-REV0 Mass Calculation 20230210.xlsx"), sheet_name='RealizedMass', header=None)
        
            # Platform.MassObj_Layer4=Input.MassObject(0, [0,0,0])
        
            # #### T1 represents the Platform without coupling, T2 the Platfom with coupling    
            # Platform.MassObj_T2_AllMinL4Flt=Input.MassObject(dfW.iloc[5,2], np.array(dfW.iloc[8,2:5]))
            # Platform.MassObj_T2_AllMinL4Flt.I=np.array( dfW.iloc[9:12,2:5])
             
            # Platform.MassObj_T1_AllMinL4Flt=Input.MassObject(dfW.iloc[5,7], np.array(dfW.iloc[8,7:10]))
            # Platform.MassObj_T1_AllMinL4Flt.I=np.array( dfW.iloc[9:12,7:10])
            
            # Platform.MassObj_T1=Platform.add_mass_objects([Platform.MassObj_T1_AllMinL4Flt,Platform.MassObj_Layer4])
            # Platform.MassObj_T2=Platform.add_mass_objects([Platform.MassObj_T2_AllMinL4Flt,Platform.MassObj_Layer4])
            
            # Platform.MassNoFloaters_T1=Platform.MassObj_T1.Mass
            # Platform.MassNoFloaters_T2=Platform.MassObj_T2.Mass
            
            # #### End section rigid body
           
            
            
            #### Split platforms for structural model
            
            Platform.split_platform_in_three_initiate()
            Platform.BeamAlpha_T1=3
            Platform.BeamAlpha_T2=3
            Platform.TranslationT2=np.array([0,0,0])
            Platform.split_platform_in_three_calc_details()
            Platform.find_platforms(Layout)
            
            Coupling=Input.Coupling()            
            Coupling.calc_details()
            
            Coupling.Position['Center']=np.array([-5.3,10.9, 2.65]) # wrt T1 due connect to T1
            Coupling.Len['In']=4.282 ## also out mid        
            Coupling.Len['Out']=Coupling.Len['In']
            Coupling.Len['Mid']=4.667
            
            Coupling.Position['In']=np.array([-0.3, 9.3, 3.4])#np.array([-0.270,	9.392,	3.692])
            Coupling.Position['Out']=np.array([-0.3, 12.5, 3.4])#np.array([-0.270,	11.992,	3.692])
            Coupling.Position['Mid']=np.array([-0.379,10.9,0.03])#np.array([-0.492,	10.692,	0.181])
            
            Coupling.Gap= 9.5 + Coupling.Position['Center'][0] + Coupling.Position['Mid'][0]# -Coupling.Position['Center'][0] - (-Coupling.Position['Mid'][0])#3.3# x distance betweeen Postion 'center' and Position 'mid'  #### 3.5 set for better statics stability
            Coupling.DetailedModel=False ### this true to include tripods in coupling modelling

            Platform.set_coupling_inertia_contribution(Coupling, Compose, CouplingZPos=3.1)
            Platform.AdjustMassOfEdgePlatforms = HeavierEdgePlatforms
            PlatformEdge = copy.deepcopy(Platform)
            
            Mooring100t=Input.Mooring(Platform, Coupling, Environment)
            Mooring100t.Representation = 'Finite element' # 'Finite element' or 'Analytic catenary' method for mooring lines
            Mooring100t.LocationLineTypes=r"X:\01_Hydromechanics\Calculations\OrcaFlex\01_LineTypes\20231127_BoxTrussUpd.dat"
        
            # StiffType='KappaPolyester88mm-SD244t' ### StaticDynamic or QuasiStatic
            StiffTypeMain='KappaPolyester143mm-SD700t'#'KappaPolyester120mm-SD454t' #'KappaPolyester88mm-SD244t'#'KappaPolyester120mm-SD454t' 
            StiffTypeMain2='Chain100mm'
            StiffTypeMain3='Swivel100mm'            
            StiffTypeVLink='Chain80mm'#'KappaPolyester88mm-SD244t' 
            StiffTypeVLink2='KappaPolyester120mm-SD454t'#'KappaPolyester88mm-SD244t' 
            # Mooring100t.M1LineTypes=[StiffType,StiffType,  StiffType, StiffType, StiffType]   
            Mooring100t.M1LineTypes=[StiffTypeMain,StiffTypeMain2,StiffTypeMain3]### no clumpweight and chains, synthetic line split at mooring buoy
            Mooring100t.VLinkLineTypes = [StiffTypeVLink,StiffTypeVLink2] ### synthetic line type
            Mooring100t.VLink = BridalMooring
            if Mooring100t.VLink:
                Mooring100t.AdjustBridalCornerAnchorPos = AdjustBridalCornerAnchorPos
            # VLinkMidLength = 45#45
            VlinkLength = 46
            
            # MooringPosition=np.array([-7.572, 9.922, 0.223])#### This position is based on information from Niels
            # MooringPosition[0:2]=MooringPosition[0:2]+Platform.GeoCoG[0:2]
            ### adjust the mooring position to the coupling center pos
            # MooringPosition[2] = Coupling.Position['Center'][2]
            Mooring100t.Position=np.array([0, 11.583,0])
            Mooring100t.Angle=0
            Mooring100t.ConnectedType=[1,2]
            # Mooring100t.ConnectedDirection=[240]
            Mooring100t.MooringBuoy='Yes'
            # Mooring100t.Depth=Environment.Depth-0.1
            # Mooring100t.VLinkLength = (VLinkMidLength**2+Mooring100t.Position[1]**2)**0.5
            VLinkMidLength = (VlinkLength**2 - Mooring100t.Position[1]**2)**0.5
            
            if MG:
                Mooring100t.M1LineTypes=[StiffTypeMain+'-MG',StiffTypeMain2+'-MG',StiffTypeMain3+'-MG']#,StiffType+'-MG',  StiffType+'-MG', StiffType+'-MG', StiffType+'-MG']       
                Mooring100t.VLinkLineTypes = [StiffTypeVLink+'-NMG-Corrod',StiffTypeVLink2+'-NMG-Aged'] ### synthetic line type
                # Mooring100t.VLinkLineTypes = ['OFChain38mmGrade3','OFChain38mmGrade3'] ### synthetic line type
            Mooring100t.Type='Line'
            Mooring100t.M1SegementLengths=[2,1,0.2]#,2.5,2.5,2.5,2.5] #### This is only needed if a line is used, not needed for tether mooring
            Mooring100t.M2SegementLengths=[1,2]
            if Mooring100t.VLink:
                Mooring100t.M1Lengths=[280-VLinkMidLength,125,0.2]#[10,225-VLinkMidLength,25,6.0,27-3]
                Mooring100t.M2Lengths = [10-LengthChange, VlinkLength-10]
                Mooring100t.RadiusAnchor=np.sum(Mooring100t.M1Lengths)+VLinkMidLength-4 # for design V3 -1m
            else:
                Mooring100t.M1Lengths=[280,125,0.2]#[10,225-VLinkMidLength,25,6.0,27-3]
                Mooring100t.RadiusAnchor=np.sum(Mooring100t.M1Lengths)-4
            
            
            # MooringPosition=np.array([-7.572, 9.922, 0.223])#### This position is based on information from Niels
            # MooringPosition[0:2]=MooringPosition[0:2]+Platform.GeoCoG[0:2]
            # ### adjust the mooring position to the coupling center pos
            # # MooringPosition[2] = Coupling.Position['Center'][2]
            # Mooring100t.Position=MooringPosition 
            # Mooring100t.Angle=0
            # Mooring100t.ConnectedType=[1,2]
            # # Mooring100t.ConnectedDirection=[240]
            # Mooring100t.MooringBuoy='Yes'
            # Mooring100t.Depth=Environment.Depth-0.1
            # Mooring100t.VLinkLength = (VLinkMidLength**2+Mooring100t.Position[1]**2)**0.5
            
            Layout.calc_details(Platform, Coupling) 
            Mooring100t.calc_connected_platforms(Layout)#, Directions=[0])  
            # StiffType='KappaPolyester54mm-SD80t'
            # Mooring80t=copy.deepcopy(Mooring100t)
            # Mooring80t.M1Lengths=[10,175,25,4.5,27-1.5]
            # Mooring80t.RadiusAnchor=np.sum(Mooring80t.M1Lengths)-2
            # Mooring80t.M1LineTypes=['OFChain38mmGrade3',StiffType,  StiffType, 'ConcreteChainLarge', 'OFChain38mmGrade3']   
            # if MG:
            #     Mooring80t.M1LineTypes=['OFChain38mmGrade3',StiffType+'-MG',  StiffType+'-MG', 'ConcreteChainLarge', 'OFChain38mmGrade3MG-Dcm']        
            
            # StiffType='KappaPolyester46mm-SD55t'
            # Mooring55t=copy.deepcopy(Mooring100t)
            # Mooring55t.M1Lengths=[10,125,25,3,27]
            # Mooring55t.RadiusAnchor=np.sum(Mooring55t.M1Lengths)-2
            # Mooring55t.M1LineTypes=['OFChain38mmGrade3',StiffType,  StiffType, 'ConcreteChainLarge', 'OFChain38mmGrade3']   
            # if MG:
            #     Mooring55t.M1LineTypes=['OFChain38mmGrade3',StiffType+'-MG',  StiffType+'-MG', 'ConcreteChainLarge', 'OFChain38mmGrade3MG-Dcm']        
                
        
            MooringBuoy1=Input.MooringBuoy(Mooring100t)
            MooringDiaRatio = float(re.compile(r'\d+').findall(StiffTypeMain)[0])/60 # 60m was synthetic mooring diameter used in merganser
            # float(StiffTypeMain.split('-')[0][-5:-2])/60 ### 76mm is the current using mooring line diameter, 60m was synthetic mooring diameter used in merganser
            MooringBuoy1.NetBuoyancy= 7.5 #12.15 #15.845 #2.3*MooringDiaRatio**2
            # MooringBuoy1.MassRatio=628/550 ###1.14
            MooringBuoy1.Height = 1.5 #1.00*MooringDiaRatio
            MooringBuoy1.MassRatio= 2.5/7.5 #3170/12150 #3685/15845
            MooringBuoy1.BuoyLocation=[np.sum(Mooring100t.M1Lengths[-2:])]#([np.sum(Mooring100t.M1Lengths[-3:])])#### This should be in brackets
            # MooringBuoy1.BuoyLocation=np.array([Mooring100t.M1Lengths[-1]])#([np.sum(Mooring100t.M1Lengths[-3:])])#### This should be in brackets
            MooringBuoy1.calc_details_cylinder_buoyant()
            if MG:
                MooringBuoy1.add_marine_growth_cylindershape(0.05,1.3) 
            MooringBuoy1.Offset=MooringBuoy1.Height/2 
            MooringBuoy1.set_name()          
        
            # MooringBuoy2=Input.MooringBuoy(Mooring80t)
            # MooringBuoy2.NetBuoyancy=2.3
            # MooringBuoy2.MassRatio=628/550
            # MooringBuoy2.BuoyLocation=np.array([np.sum(Mooring100t.M1Lengths[-3:])])#### This should be in brackets
            # MooringBuoy2.calc_details_sphere_buoyant()
            # if MG:
            #     MooringBuoy2.add_marine_growth(0.1,1.3)  
            # MooringBuoy2.Offset=MooringBuoy2.Height/2         
            # MooringBuoy2.set_name()
            
            # MooringBuoy3=Input.MooringBuoy(Mooring55t)
            # MooringBuoy3.NetBuoyancy=1.3
            # MooringBuoy3.MassRatio=628/550
            # MooringBuoy3.BuoyLocation=np.array([np.sum(Mooring100t.M1Lengths[-3:])])#### This should be in brackets
            # MooringBuoy3.calc_details_sphere_buoyant()
            # if MG:
            #     MooringBuoy3.add_marine_growth(0.1,1.3)  
            # MooringBuoy3.Offset=MooringBuoy2.Height/2         
            # MooringBuoy3.set_name()
               
            # PTTool=Input.MooringBuoy(Mooring100t)
            # PTTool.SpecificGravity=7.8
            # PTTool.Mass=0.15
            # PTTool.calc_details_sphere_specificgravity()
            # #### This should be in brackets 
            # PTTool.ReferenceEnd='End A'
            # PTTool.Name='PTTool'
              
            
            for k in [1.25]: #### stiffness factor of 1.25 is applied on the gas spring to account for leakage effects
                    
                #### Reset floater properties for in the loop
                Floater=Input.Floater()
                Floater.Height=10.74+FltHPlus
                Floater.MassTop=0
                Floater.MassTempex=0
                Floater.HeightTempex=6
                Floater.ZTempex=10.74-7.162###7.162 #### Adjusted for floater reference system
                
                Floater.MassBottom=18.20*0.64/3
                Floater.Diameter=FloaterDia
                Floater.Mass=Floater.MassTop+Floater.MassBottom+Floater.MassTempex
                # Floater.CdAxial = 1.2

                MassCheck=(Floater.Mass*3+Platform.MassNoFloaters_T2-Platform.MassT2)**2<0.01
                print('T2MassCheck='+str(MassCheck))
                print('T2 Mass Input = ', Platform.MassT2)
                print('T2 Mass Calc = ',round(Floater.Mass*3+Platform.MassNoFloaters_T2,3))
                MassCheck=(Floater.Mass*3+Platform.MassNoFloaters_T1-Platform.MassT1)**2<0.01
                print('T1MassCheck='+str(MassCheck))
                print('T1 Mass Input = ', Platform.MassT1)
                print('T1 Mass Calc = ',round(Floater.Mass*3+Platform.MassNoFloaters_T1,3))
                                
                Platform.MassNoFloatersMean=(Platform.MassNoFloaters_T1+Platform.MassNoFloaters_T2)/2
                # Floater.CoG = np.array([0,0,6])
                Floater.calc_details(Platform.MassNoFloatersMean)
                
                
                for DampingRatio in [0.25]:#1,2]:
                ## Automatically change the draft of the floaters
                    ## Due to recent update apply marine growth should always be used
                    if FloaterMG:
                        MGThickness = 0.05
                    else:
                        MGThickness=0 ##0.05
                    
                    if MGThickness==0:
                        CdBelowWL=0.8
                    else:   
                        CdBelowWL=1.05
                    
                    Dia=Floater.Diameter
                    Rio = Dia/2.25
                    
                    Floater.Spacing = 20.8
                    
                    Floater_Alpha=(Platform.EdgeLength-Floater.Spacing)/2
                    
                    Floater.HeightBottomCone = 0.75
                    Floater.VBottomCone = 3.5
                    Floater.apply_double_taper(Dia, Dia, Dia ,Platform.Height, Floater.Height-Floater.HeightBottomCone, 3,8) #### Diameters, Heigths, Segements (top to bottom)
                    # Floater.VBottomCone=1.489
                    Floater.calc_taper_draft_and_inertia(Platform.MassNoFloatersMean, Platform)                                                    
                    Floater.apply_taper_marine_growth(Platform,DensityGrowth=1.325, Thickness=MGThickness, Cd=CdBelowWL)
                   
                    Floater.PosBottom[0:2]=[1.15,Floater.Spacing/2]#[Floater_Alpha*np.tan(30*np.pi/180),Platform.EdgeLength/2-Floater_Alpha]            
                    Floater.PosBottom[2]=-Floater.Height-MGThickness
                    Draft1=Floater.Draft
                    Floater.MassOb.CoG[2]=Floater.Height-4.4578
    
                    #### Set mooring lengths
                    if Mooring100t.VLink:
                        MoorLength = np.array([np.sum(Mooring100t.M2Lengths)])
                        # MoorLength=np.array([RefMoorLen])-1#,RefMoorLen-50, RefMoorLen-100])-1
                    else:
                        MoorLength = np.array([np.sum(Mooring100t.M1Lengths)])
                        # MoorLength=np.array([RefMoorLen])-1#,RefMoorLen-50, RefMoorLen-100])-1

                    MoorList=[Mooring100t]#, Mooring80t, Mooring55t]
                    for iM, Mooring in enumerate(MoorList):
                        Mooring.set_top_chain_length(Floater, BridalMooring, TotalLength=MoorLength[iM]) #TotalLength=192.313          

                        # if MG:
                        #     Mooring.set_top_chain_length(Floater, BridalMooring, TotalLength=MoorLength[iM]-LengthChange) #TotalLength=192.313
                        # else:
                        #     Mooring.set_top_chain_length(Floater, BridalMooring, TotalLength=MoorLength[iM]) #TotalLength=192.313          
                        # PTTool.BuoyLocation=np.array([Mooring.M1Lengths[0]])
                        
                    Sheet='UpdateForCylinderVolumes'
                    Coupling.set_non_linear_stiffness(wbCylinders, Sheet)
                    Coupling.NLForce=Coupling.NLForce*k 
                    Coupling.Stiffness['In']=(Coupling.NLForce[-1]-Coupling.NLForce[0])/(Coupling.NLDisplacement[-1]-Coupling.NLDisplacement[0]) 
                    Coupling.Stiffness['Out']= Coupling.Stiffness['In']
            
                    Coupling.set_critical_damping(Platform, Floater, Ratio=DampingRatio)
                    Coupling.set_damping_power(0.5)
                    Coupling.calc_names_and_locations(Layout, Platform, Compose)
                    Coupling.LoadCells=True #### This turns on the option to include constraints for structural loading 
                    
                   
                    #### Pick up a clean base model to create new OF model in                
                    model=ofx.Model(Mooring.LocationLineTypes)
                    Gen=OrcaFlexConstructionClasses.Gen(model)
                    Compose=OrcaFlexConstructionClasses.Compose(model)  
                    if Mooring100t.VLink:
                        Mooring100t.calc_vlink_mooring_details(Floater,VLinkMidLength)
                        Mooring100t.VLinkLength = 20
                    Compose.whole_model(Gen, General, Environment, Platform,PlatformEdge, Layout, Floater, DragPlate, Mooring100t, Coupling, MooringBuoy1)#         
    
                    # MooringPTList=[1,2,1,2,1] #### List of platform type and next line direction combinations for the mooring lines 
                    # MooringDirectionList=[0,120,240,0,120] 
                    # MooringTypeList=[Mooring55t,Mooring55t,Mooring80t,Mooring100t,Mooring100t] 
                    # MooringTypeList=[Mooring55t,Mooring55t,Mooring80t,Mooring100t,Mooring100t] 
                    # MooringBuoyList=[MooringBuoy3,MooringBuoy3,MooringBuoy2,MooringBuoy1,MooringBuoy1]
                    
                    # for iM, Mooring in enumerate(MooringTypeList):
                    #     Mooring.ConnectedDirection=[MooringDirectionList[iM]]
                    #     Mooring.ConnectedType=[MooringPTList[iM]]
                    #     Compose.connect_mooring_automatic(Layout, Platform, Mooring, Gen, MooringBuoyList[iM], Coupling)
    
                    # Compose.connect_extra_mooringbuoy(MooringBuoy2, Gen)
                    # Compose.connect_extra_mooringbuoy(PTTool, Gen)
                    
                    #### Set larger diameter on outer floaters
                    if LargerFloterOnEdgePlf:
                        Floater1 = copy.deepcopy(Floater)
                        Dia1 = EdgeFltD
                        DiaRatio = Dia1/Dia
                        Floater1.MassBottom=9.182/3
                        Floater1.Mass=Floater1.MassTop+Floater1.MassBottom+Floater1.MassTempex
                        PlatformEdge.MassNoFloatersMean = (PlatformEdge.MassNoFloaters_T1+PlatformEdge.MassNoFloaters_T2)/2
                        Floater1.VBottomCone=1.489*DiaRatio**2 ### if HeightBottomCone changed accordingly, the value should be times power3
                        Floater1.apply_double_taper(2.1*DiaRatio, Dia1, Dia1 ,PlatformEdge.Height, Floater.Height-Floater.HeightBottomCone, 3,8) #### Diameters, Heigths, Segements (top to bottom)
                        Floater1.calc_taper_draft_and_inertia(PlatformEdge.MassNoFloatersMean, PlatformEdge)
                        Floater1.apply_taper_marine_growth(PlatformEdge,DensityGrowth=1.325, Thickness=MGThickness, Cd=CdBelowWL)
                        Compose.connect_edgeplatform_floaters(PlatformEdge, Floater1, DragPlate, Mooring, Gen)
                        Floater2=copy.deepcopy(Floater1)
                    else:
                    #### Add tempex to outer floaters
                        Floater2=copy.deepcopy(Floater)
                    if EstiMass or Platform.AdjustMassOfEdgePlatforms:
                        Floater2.MassTempex=1.12/2
                    else:
                        Floater2.MassTempex=0.46
                    Floater2.calc_taper_draft_and_inertia(Platform.MassNoFloatersMean, Platform)   
                    Floater2.apply_taper_marine_growth(Platform,DensityGrowth=1.325, Thickness=MGThickness, Cd=CdBelowWL)
                    Floater2.MassOb.CoG[2]=Floater.Height-4.85                                                                 
                    Compose.connect_outer_floaters(Platform, Floater2, DragPlate, Mooring, Gen)
                    
                    #### set larger cylinder stiffness on edged platforms
                    # Coupling1 = copy.deepcopy(Coupling)
                    # k_CyldOnEdgePlats = 1.5
                    # Coupling1.NLForce=Coupling.NLForce*k_CyldOnEdgePlats/k
                    # Coupling1.Stiffness['In']=(Coupling1.NLForce[-1]-Coupling1.NLForce[0])/(Coupling1.NLDisplacement[-1]-Coupling1.NLDisplacement[0]) 
                    # Coupling1.Stiffness['Out']= Coupling1.Stiffness['In']
                    # Compose.connect_EdegePlatform_constraint_links(Platform, Coupling1, Gen)

                    #### Hardcode the extra components variant1
                    # model.DestroyObject('M1_10_0_1') #### destroy the mooring line at boatlander
                    
                    # MassBoatLander=4.955#### 6t minus the fairlead
                    # OFBoatLander=Gen.mass_buoy('BoatLander', MassBoatLander, [0.23,1.45,8])
                    # OFBoatLander.Connection='Platform_10_0'
                    # Compose.set_object_position(OFBoatLander, [-0.828,10.395,-4.659])
                    # OFBoatLander.Connection='Floater10_0'
                    
                    # MassGangWay=1.309#### 2t minus the fairlead
                    # OFGangWay=Gen.mass_buoy('GangWayAndCrane', MassGangWay, [1,1,1])
                    # OFGangWay.Connection='Platform_2_0'
                    # Compose.set_object_position(OFGangWay, [0.235,8.915,4.329])
                    # # OFBoatLander.Connection='Floater2_0'
                    # SkidMass = 35.305
                    # TransSkid=Input.MassObject(SkidMass, [24.7*np.sin(np.deg2rad(60))/3,0,5.7]) #[21.57,0,4.895]
                    # SkidSize=np.array([23.38, 8, 4.5])#4.5*np.cbrt(SkidMass/6) #### [2.35, 5.9, 2.39] was for Merganser
                    # OFTransformer=Gen.mass_buoy('TransSkid_35t', TransSkid.Mass, SkidSize)
                    # Gen.set_inertia(OFTransformer,np.diag(TransSkid.Mass/12*np.array([SkidSize[1]**2+SkidSize[2]**2,SkidSize[0]**2+SkidSize[2]**2,SkidSize[0]**2+SkidSize[1]**2])))
                    # Gen.set_drag_area(OFTransformer, [SkidSize[1]*SkidSize[2],SkidSize[0]*SkidSize[2],SkidSize[0]*SkidSize[1]])
                    # Gen.set_drag_Cd_6Dbuoy(OFTransformer,[1,1,1])
                    # OFTransformer.Connection='Platform_'+str(Platform.SkidPlatform)+'_240'
                    # Compose.set_object_position(OFTransformer, TransSkid.CoG) #### This is the position of the connected floater
                    # TransSkid=Input.MassObject(20, [7,0,4.895]) #[21.57,0,4.895]
                    # SkidSize=np.array([2.35, 5.9, 2.39])*np.cbrt(20/7) #### mass increases fron 7t to 20t
                    # OFTransformer=Gen.mass_buoy('TransSkid_20t', TransSkid.Mass, SkidSize)
                    # Gen.set_drag_area(OFTransformer, [SkidSize[1]*SkidSize[2],SkidSize[0]*SkidSize[2],SkidSize[0]*SkidSize[1]])
                    # Gen.set_drag_Cd_6Dbuoy(OFTransformer,[1,1,1])
                    # OFTransformer.Connection='Platform_26_240'
                    # Compose.set_object_position(OFTransformer, TransSkid.CoG) #### This is the position of the connected floater
    
                    #### change platform mass    
                    # for deg in [0,120,240]:
                    #     for pfn in [5,6,7,8,9,26]:
                    #         PF = model['Platform_'+str(pfn)+'_'+str(deg)]
                    #         if pfn == 5:
                    #             PF.Mass = PF.Mass + 2
                    #         elif pfn == 26:
                    #             PF.Mass = PF.Mass - 10/3
                    #         elif pfn == 9:
                    #             PF.Mass = PF.Mass + 10/3
                    #         else:
                    #             PF.Mass = PF.Mass + 3
                    #### add braces between floaters
                    # Compose.floater_braces(model,Layout,Floater,Gen)
                        
                    ### End of hardcode extra compoenents
                               
                ### Lines optmization
    
                # Environment.Hs=1 #### Neccesary when running irreg waves 
                Ur=[]
                Environment.NumberOfWaveComponentPerDirection=200
                for iW, WDirOF in enumerate(WDirsOF):      
                    # if WDirOF == 150:
                    #### Irregular Seastate description
                    # Environment.Tp=dfEnvBins['Tp'][it]
                    # Environment.Hs=dfEnvBins['Hs'][it]
                    # Environment.set_gamma()
                    # Environment.HsOrHmax='Hs'
                    # General.DurationSimulation=1800
                    #### End Irregular Seastate description 
                    
                    Environment.WaveHeading=WDirOF
                    Environment.WindDirection=WDirOF
                    # CurrentDirections=[45,225]
                    # if 225-90<WDirsOF[iW]<225+90:    
                    #     Environment.CurrentDirection=225
                    # else: 
                    #     Environment.CurrentDirection=45
                    Environment.CurrentDirection=WDirOF
                    # Hmax=dfEnv['Estimate for Agnes'][iW]
                    # Hmax=np.max(dfEnv['Estimate for Agnes'])
                    Environment.HsOrHmax='Hmax'
                    
                    for c in [0.8]:#[0,0.8]:
                        Environment.CurrentSpeed=c
                        # Environment.CurrentPowerLawExponent=1.8
                        if c > 0:
                            Environment.WindSpeed=27.8
                        else:
                            Environment.WindSpeed=27.8                                            
                        # Environment.WindSpeed=26.6
                        for iT,T in enumerate(Ts):
                            #### Regular seastate description
                            # Environment.Hmax=H[ii]
                            Environment.Tass=T
                            Environment.Hmax=1
                            Environment.calc_details_Hmax_Hb()
                            if Environment.Hb<Hmax:    
                                Environment.Hmax=Environment.Hb
                            else:
                                Environment.Hmax=Hmax
                            # if it==1:
                            #     while Environment.Hb<Hmax: 
                            #         Environment.Tass=Environment.Tass+0.01
                            #         Environment.calc_details_Hmax_Hb()
                            #     Environment.Hmax=Hmax
                            # print(Environment.Tass)
                            # print(f"Hmax={Environment.Hmax}")
                            Environment.calc_details_Hmax_Hb()
    
                            # Ur.append(Environment.UrHmax)
                            General.DurationSimulation=6*Environment.Tass
                            Environment.select_wave_type()
                            # print(f"T={T}, Hmax={Environment.Hmax}, WaveType={Environment.WaveType}")
                            # print(f"urhmax={Environment.UrHmax}")
                            #### End regular wave selection                  
                            
                            # if MG and FloaterMG:
                            #     MGEx='-MoorFltMG'
                            # elif MG:
                            #     MGEx='-MoorMG'
                            # elif FloaterMG:
                            #     MGEx='-FltMG'
                            # else:
                            #     MGEx=''
                                
                                                    
                            Gen.enviroment_non_linear_hmax(Environment)
                            
                            if IrregularEnv:
                                Environment.calc_TpHs_from_TassHmax()
                                Gen.environment_jonswap(Environment)
                                # General.DurationSimulation=3600*3
                                
                            Environment.name()
                            General.StaticsMinDamping=10
                            
                            Gen.general(General, Environment)
                            
                            filename = 'Tri_' + str(N[0]) + '-' + str(N[1]) + '_Tag_Alt' + str(Alt) + MGEx + Environment.Name +'.yml'
                            Input.CheckSummary(filename, LocSims, model, Chk_Summary=False, model_save=True)
                            # model.SaveData(os.path.join(LocSims, filename))
        # if not MG:
            ChkInstance = Input.CheckSummary(filename, LocSims, model=False, Chk_Summary=True, model_save=False)
            self = ChkInstance
        # globals()[f"ChkInstance{iMG}"]
    print(f'Model building complete.\n You can check model in {LocSims}\n ---')
    
    if False: #RunSims:
        YmlFiles = glob.glob(os.path.join(LocSims,'*yml'))
        ExtractData = DataExtractionClasses.DataExtraction(LayoutShape,LocSims)
        with multiprocessing.Pool(processes=Threads) as pool:
            #### Using a 'with' statement, don't need to explicitly call pool.close() and pool.join()
            #### The 'with' statement automatically takes care of closing the pool and waiting for all tasks to finish when it exits the block
            print(f' --- Simulation running with {Threads} threads.')
            
            #### Run model simulations
            pool.map(OrcaFlexConstructionClasses.run_simulation,YmlFiles)
            print(' ---\n Simulation finished.\n ---\n Data Extracting')

            #### Run data extraction
            SimFiles = glob.glob(os.path.join(LocSims,'*sim'))
            
            pool.map(ExtractData.extract_and_process_ofdata,SimFiles)
            
        ExtractData.postprocess_propkl()
    # winsound.Beep(1000, 1000)
    