# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 09:52:51 2021

@author: Sander
"""
# import collections
# try:
#     from collections import abc
#     collections.MutableMapping = abc.MutableMapping
# except:
#     pass
# import sys
# sys.path.insert(1, 'Classes')
## Script to extract details
import OrcFxAPI as ofx
import numpy as np
import numpy.matlib
import glob
import os
import pickle
from openpyxl import load_workbook
import InputClasses 
import OrcaFlexConstructionClasses 
import PostProcessingClasses as PPC
import matplotlib.pyplot as plt
import multiprocessing
import time

                                
class DataExtraction:
    def __init__(self, LayoutShape=None, SimsPath=None):
        self.SimsPath=SimsPath
        self.LayoutShape = LayoutShape

    def check_time(self, old_time, command_no) :
        cur_time = time.time()
        print(f"Time taken for command {command_no}: {(cur_time - old_time):6.2f}")
        return cur_time
    
    def postprocess_propkl(self):
        # LocSims=self.SimsPath
        # simfiles=glob.glob(os.path.join(LocSims,'*sim'))
        
        # for sim in simfiles:
        #     result=self.extract_and_process_ofdata(sim)
            
        # # # # ## Locate processed pkl files to summerize in excel for pivot charts and for plot loop
        Sims=glob.glob(os.path.join(self.SimsPath,'*_Pro.pkl'))
        # #### Activate Excel class 
        
        ToExcel=PPC.ToExcel(Sims)
        # #### Run specific extractions 
        # dfAirGap=ToExcel.air_gap(os.path.join(self.SimsPath,'ResultsAirGap.xlsx'), SaveExcel=False) 
        # dfRelAngels=ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsRelativeAngles.xlsx'), 'RelativeAngle', SaveExcel=False)
        # dfMooring=ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsMooring.xlsx'), 'Mooring' , SaveExcel=True)
        # dfCoupling=ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsCoupling.xlsx'), 'Coupling', SaveExcel=True)
        # dfFloater =ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsFloater.xlsx'), 'Floater', SaveExcel=True)
        # dfClash=ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsClash.xlsx'), 'Clash', SaveExcel=False)
        # dfPlatform=ToExcel.platform(os.path.join(self.SimsPath,'ResultsPlatform.xlsx'), 'Platform', SaveExcel=False)
        # # # dfGap=ToExcel.coupling_gap(os.path.join(self.SimsPath,'ResultsGap.xlsx'), SaveExcel=False)
        # dfBeams=ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsBeams.xlsx'), 'Beams', SaveExcel=False)
        # dfCons=ToExcel.floater_cons_coupling_clash(os.path.join(self.SimsPath,'ResultsConstraints.xlsx'), 'Cons', SaveExcel=False)
    # # # ## Check data and load or create pickle files 
    def extract_and_process_ofdata(self,sim):
        start_time = time.time()   
        General, Environment, Platform, Floater, DragPlate, Coupling, Mooring, Layout = OrcaFlexConstructionClasses.Activate.all_input(InputClasses)  
        model=1   
        
        # Gen=OrcaFlexConstructionClasses.Gen(model)
        Compose=OrcaFlexConstructionClasses.Compose(model)   
        # Coupling.calc_names_and_locations(Layout, Platform, Compose)
        
        OverwriteRaw=False
        OverwritePro=False
        
        ## Extract the platform configuration for extraction of data
    
        SimBase=os.path.basename(sim)
        Tri=SimBase.split('_')[1]
        # N1=Tri.split('-')[0]
        # N2=Tri.split('-')[1]
        Layout.calc_triangle(self.LayoutShape[0],self.LayoutShape[1])
        N=Tri
        try:
            Flip=Tri.split('-')[2]
            Layout.flip_matrix()
        except:
            pass         
        # Layout.Matrix=np.array([[0,2,1]])
        
        # ########## Towing Condition #########    
        # N='Tow'
        # Layout.Matrix=np.array([[1,2,1,2,1,2,1,2]])
        # ########## Towing Condition #########
        # if SimBase.split('_')[3]=='LrgFlt-NMG-Dmp0.25':
        Platform.EdgeLength=23.989
        Platform.EdgeLengthOuter=37
        Platform.calc_details_basic()
        Platform.add_airgap_probes_xy_outer(Platform.EdgeLengthOuter)
        # elif SimBase.split('_')[3]=='SmlPlt-NMG-Dmp0.25':
        #     Platform.EdgeLength=20.57
        #     Platform.EdgeLengthOuter=30
        #     Platform.calc_details_basic()
        #     Platform.add_airgap_probes_xy_outer(Platform.EdgeLengthOuter)
        
        Layout.calc_details(Platform, Coupling)
        Coupling.calc_names_and_locations(Layout, Platform, Compose)
        Mooring.calc_connected_platforms(Layout)
        # Platform.add_airgap_probes_xy_outer(34)
    
        ## Set baisis to sim complete check 
        Complete=True
        if Complete: 
            ##  Check that pkl file exists or run get data function
            if os.path.isfile(sim[0:-4]+'.pkl') and not OverwriteRaw:
                print(os.path.basename(sim) + ' Exists')
            else:
                model=ofx.Model(sim)
        ## ActivateGetData
                GetData=PPC.GetData(model, Compose ) #,SpecifiedTimePeriod=2
                if model.simulationComplete==False:
                    GetData.move_incomplete()
                    Complete=False
                else:
                    old_time = self.check_time(start_time, 1) 
                    GetData.environment(N, LastTwoWaves=False)
                    old_time = self.check_time(old_time, 2) 
                    GetData.platform_raw(Layout, Platform) #, MotionsFromStaticPosPOI=[6.6455,0,1.5],SkidCoG=[24.7*np.sin(np.deg2rad(60))/3,0,5.7])
                    ## Activate the next line when working with constraints in the simulations
                    old_time = self.check_time(old_time, 3)
                    GetData.constraints_raw(Coupling)
                    ##
                    # GetData.coupling_raw(Coupling)
                    old_time = self.check_time(old_time, 4)
                    GetData.floater_raw()
                    # if GetData.Split:
                        # GetData.beams_raw()
                    # GetData.mooring_raw(Mooring)
                    old_time = self.check_time(old_time, 5)
                    RawData=GetData.RawData
                    pickle.dump(GetData.RawData,open(sim[0:-4]+'.pkl','wb'))
                    print(os.path.basename(sim) + ' Extracted')
    
        if Complete:
            ##  Check that pkl file exists or run process data function 
            if os.path.isfile(sim[0:-4]+'_Pro.pkl') and not OverwritePro:
                print(os.path.basename(sim)+'_Pro' + ' Exists')
            else:
                RawData=pickle.load(open(sim[0:-4]+'.pkl','rb'))
                old_time = self.check_time(old_time, 6)
                ProcessData=PPC.ProcessData(RawData, Compose)    
                ProcessData.environment()
                #### Include_Accelrations option has been added as it take a lot of time to post process, for irregular wave it is suggested to put this on false
                # ProcessData.platform(Platform, Include_Accelerations=True, StructuralOutputPOI_Acc=[[21.627,0,5.047]], StructuralOutputPOI_Name=['SkidLocation'],PositionPOIs=[[18.1486,-0.8,	-10.752],[21.627,0,5.047]], PositionPOIs_Name=['UmbilicalHangOff', 'SkidLocation'])
                ProcessData.platform(Platform, Include_Accelerations=False) 
                ProcessData.floater_beams(RawData, Mooring, 'Floater')
                # if RawData['Environment']['PlatformSplit']:
                #     ProcessData.floater_beams(RawData, Mooring, 'Beams')
                # ProcessData.airgap(RawData)
                # ProcessData.coupling(RawData) 
                ## Activate the next line when working with constraints in the simulations
                old_time = self.check_time(old_time, 7)
                ProcessData.constraints()
                ##
                # ProcessData.gap_details(RawData)
                # ProcessData.relative_angles()
                # ProcessData.mooring()
                # ProcessData.clash_details(Layout)
                # ProcessData.calc_bottom_clash(0.409)
                old_time = self.check_time(old_time, 8)
                ProData=ProcessData.ProData
                # RawData=ProcessData.RawData
                pickle.dump(ProData,open(sim[0:-4]+'_Pro.pkl','wb'))
                # pickle.dump(RawData,open(sim[0:-4]+'.pkl','wb'))
                print(os.path.basename(sim) + ' Processed')
                old_time = self.check_time(start_time, 9)

        
        return sim
    

# def winapi_path(dos_path, encoding=None):
#     if (not isinstance(dos_path, str) and encoding is not None): 
#         dos_path = dos_path.decode(encoding)
#     path = os.path.abspath(dos_path)
#     if path.startswith(u"\\\\"):
#         return u"\\\\?\\UNC\\" + path[2:]
#     return u"\\\\?\\" + path

# SimsPath=r'C:\Users\Sander\Documents\03_LocalHydrodynamics\01_OrcaFlex\25_20221110_UpdFinal\ALS_EnvSelection\CylinderALS_0.5mm'
# SimsPath=r'C:\Users\Sander\Documents\03_LocalHydrodynamics\01_OrcaFlex\25_20221110_UpdFinal\ALS_EnvSelection\MooringALS'
# SimsPath=r'Z:\01_Hydro\46_230616_Hex_ConceptUpd\Hex3\FltD2.25_MBUpd0_BUpd1_FltHPlus2_Cplk1.5_SkidNearCoG1'
# SimsPath=r'C:\Users\Sander\Documents\03_LocalHydrodynamics\01_OrcaFlex\25_20221110_UpdFinal\FLS_EnvSelection'
# LocSims=winapi_path(SimsPath)
# LocSims=SimsPath
# 
# simfiles=glob.glob(os.path.join(LocSims,'*sim'))

# if __name__ == '__main__':    
#     pool = multiprocessing.Pool(processes=2)
#     results=pool.map(extract_and_process_ofdata, simfiles)
#     pool.close()
    
# for sim in simfiles:
#     result=extract_and_process_ofdata(sim)

# # ### This is for to check prodata content
# sim=simfiles[-1]
# ProData= pickle.load(open(sim[0:-4]+'_Pro.pkl','rb'))   
# RawData= pickle.load(open(sim[0:-4]+'.pkl','rb'))   


# ProData=pickle.load(open(os.path.join(LocSims,'Tri_3-1_Tag_39.2t-T1Inv-Mr14-MG_Wdp22.5_Hmx12.42_T12.20_Wdr120_Cv0.80_Cdr120_Wnv27.00_Wndr120_Pro.pkl'),'rb'))  

# # # # ## Locate processed pkl files to summerize in excel for pivot charts and for plot loop
# Sims=glob.glob(os.path.join(LocSims,'*_Pro.pkl'))


# # #### Activate Excel class 
# ToExcel=PPC.ToExcel(Sims)
# # #### Run specific extractions 
# dfAirGap=ToExcel.air_gap(os.path.join(LocSims,'ResultsAirGap.xlsx'), SaveExcel=False) 
# dfRelAngels=ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsRelativeAngles.xlsx'), 'RelativeAngle', SaveExcel=False)
# dfMooring=ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsMooring.xlsx'), 'Mooring' , SaveExcel=True)
# dfCoupling=ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsCoupling.xlsx'), 'Coupling', SaveExcel=False)
# dfFloater =ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsFloater.xlsx'), 'Floater', SaveExcel=False)
# dfClash=ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsClash.xlsx'), 'Clash', SaveExcel=False)
# dfPlatform=ToExcel.platform(os.path.join(LocSims,'ResultsPlatform.xlsx'), 'Platform', SaveExcel=False)
# # # # dfGap=ToExcel.coupling_gap(os.path.join(LocSims,'ResultsGap.xlsx'), SaveExcel=False)
# dfBeams=ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsBeams.xlsx'), 'Beams', SaveExcel=False)
# dfCons=ToExcel.floater_cons_coupling_clash(os.path.join(LocSims,'ResultsConstraints.xlsx'), 'Cons', SaveExcel=False)


# # Compose=OrcaFlexConstructionClasses.Compose()  
# Plot=PPC.Plot(Compose, LocSims) 

# for Sim1 in Sims:
   
#     Sim=os.path.join(LocSims,Sim1)    
#     ProData=pickle.load(open(Sim,'rb'))

#     Plot.air_gap(ProData, Sim,'Min', Type='Platform')
#     Plot.coupling_clash_platform(ProData, Sim,'y bend moment', 'AbsMax', Type='Beams') 
    
#     Plot.coupling_clash_platform(ProData, Sim,'Clearance', 'Min', Type='Clash') 
#     Plot.coupling_clash_platform(ProData, Sim,'Accelerations', 'AbsMax', Type='Platform', index=2) 
             
        
            
        
##### Check to see if relative yaw makes sense

# a=RawData['Platform']['9']['MotionsFromStaticPos']['Origin']['TimeHistory'][:,5]
# b=RawData['Platform']['10']['MotionsFromStaticPos']['Origin']['TimeHistory'][:,5]
# c=a+b

# plt.plot(c)

# a1=RawData['Platform']['10']['Yaw120']['TimeHistory']
# plt.plot(a1)

# a2=RawData['Platform']['10']['Yaw0']['TimeHistory']
# plt.plot(a2)

# a3=RawData['Platform']['10']['Yaw240']['TimeHistory']
# plt.plot(a3)

# plt.grid()

# b1=RawData['Platform']['10']['Yaw120']['TimeHistory']
# c1=180-np.abs(a1-b1)

# plt.plot(c1)

# b2=RawData['Platform']['10']['Yaw0']['TimeHistory']
# c2=180-np.abs(a2-b2)

# plt.plot(c2)
