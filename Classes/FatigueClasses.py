# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 13:09:55 2022

@author: Sander
"""

#### Set of function / classes used for the fatigue calculations
import numpy.matlib
import glob
import os
import pickle
from openpyxl import load_workbook
import InputClasses 
import OrcaFlexConstructionClasses 
import PostProcessingClasses as PPC
import matplotlib.pyplot as plt
import pandas as pd
import time
import numpy as np
import scipy.special as sc
import scipy
import re
import random

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

class read_struct_input:
    def __init__(self, LocSims):
        self.LocSims=LocSims
        self.StructInput={}
        # self.VersionNameInput1='v30c_FAT'
        # self.VersionNameInput2='_nodes_FAT_set'
                
    def read_stress_transfer_input_t1(self):
        #### For this function the version name of the files always had to be manually adjusted here 
        
        #### Load stress transfer matrix
        StructuralInputFile=os.path.join(self.LocSims,'FatigueInput','T1')
        
        ### shape of this matrix is hardcoded and should be adjusted 
        data=pd.DataFrame()
        # for j in [1,2]:
            # i=0
            # v27b_FAT1_nodes_b_set1
            
        fileObject = open(os.path.join(StructuralInputFile, 'v32g_FAT1_nodes_b_set1.txt'), "r")
        tmp = pd.read_csv(fileObject, skiprows=(0,1,2), delim_whitespace=True)
        data=pd.concat([data,tmp])   
                
        self.StructInput['Nodes_T1']=data['NODE']
        
        Out=np.zeros([6,54+12,len(data['NODE'])])
        
        for i in range(54+12):
            data=pd.DataFrame()
            # for j in [1,2]:
            if i==8:
                a=1
            if i+1 <= 18:
                fileObject = open(os.path.join(StructuralInputFile, 'v32g_FAT1_nodes_b_set'+ str(i+1) + '.txt'), "r")
            elif i+1<=36:
                fileObject = open(os.path.join(StructuralInputFile, 'v32g_FAT2_nodes_b_set'+ str(i+1-18) + '.txt'), "r")
            elif i+1<=54:
                fileObject = open(os.path.join(StructuralInputFile, 'v32g_FAT3_nodes_b_set'+ str(i+1-36) + '.txt'), "r")
            else:
                fileObject = open(os.path.join(StructuralInputFile, 'v32j_FAT4_nodes_b_set'+ str(i+1-54) + '.txt'), "r")
            
            tmp = pd.read_csv(fileObject, skiprows=(0,1,2), delim_whitespace=True)
            ### RawData is not save for space reasoning
            data=pd.concat([data,tmp])
            Out[:,i,:]=data.iloc[:,1:7].transpose()
            print('v32g_FAT_nodes_b_set' + str(i+1) + '.txt')
                
        self.StructInput['StressTransferMatrix_T1']=Out
        
    def read_pin_ltm_t1(self):
        #### This function is to read out and set-up the pin load transfer matrix
        import csv
        #### Load stress transfer matrix
        StructuralInputFile=os.path.join(self.LocSims,'')
                
        #### shape of this matrix is hardcoded and should be adjusted 
        def read_pin_csv(FileName,StructuralInputFile):
            #### Create a read data function here
            fileObject=os.path.join(StructuralInputFile, FileName)
            data = pd.read_csv(fileObject,sep='\\t',header=(1), engine='python')
            # data=data.rename(columns={'"Probe Name': 'Probe Name'})
            # data['Probe Name']=data['Probe Name'].str.replace('"','')
            # data=data.rename(columns={'Time"': 'Time'})
            # data['Time']=data['Time'].str.replace('"','')
            return data    
        
        data1=read_pin_csv('Static matrix FAT1 - probe results -pin4 cleaned.csv',StructuralInputFile)
        data2=read_pin_csv('Static matrix FAT1 - probe results -pin4 cleaned.csv',StructuralInputFile)
        data3=read_pin_csv('Static matrix FAT1 - probe results -pin4 cleaned.csv',StructuralInputFile)
        
        self.StructInput['Nodes_T1']=data1['Probe Name'].unique()
        
        #### Out matrix [xzy, unit loads, pins,]
        Out=np.zeros([3,54,int(len(data1['Probe Name'].unique()))])
        
        for i, PName in enumerate(self.StructInput['Nodes_T1']):
            
            # for j in [1,2]:
            ind=data1['Probe Name']==PName
            #### Transfer of N output to kN output
            Out[:,0:18,i]=np.array(data1.loc[ind,['X','Y','Z']]).transpose()/1000
            Out[:,18:36,i]=np.array(data2.loc[ind,['X','Y','Z']]).transpose()/1000
            Out[:,36:,i]=np.array(data3.loc[ind,['X','Y','Z']]).transpose()/1000
                
        self.StructInput['StressTransferMatrix_T1']=Out
    
        
    def read_stress_transfer_input_t2(self):
        #### This is currently outdated
        #### Load stress transfer matrix
        StructuralInputFile=os.path.join(self.LocSims,'FatigueInput','T2')
                
        ### shape of this matrix is hardcoded and should be adjusted 
        data=pd.DataFrame()
        for j in [0,120,240]:
            i=0
            fileObject = open(os.path.join(StructuralInputFile, 'v25c_FAT1_nodes_'+str(j)+'b_set' + str(i+1) + '.txt'), "r")
            tmp = pd.read_csv(fileObject, skiprows=(0,1), delim_whitespace=True)
            data=pd.concat([data,tmp])   
        self.StructInput['Nodes_T2']=data['NODE']
        
        Out=np.zeros([6,90,len(data['NODE'])])
               
        for i in range(42):
            data=pd.DataFrame()
            for j in [0,120,240]:
                fileObject = open(os.path.join(StructuralInputFile, 'v25c_FAT1_nodes_'+str(j)+'b_set' + str(i+1) + '.txt'), "r")
                tmp = pd.read_csv(fileObject, skiprows=(0,1), delim_whitespace=True)
                ### RawData is not save for space reasoning
                data=pd.concat([data,tmp])
            Out[:,i,:]=data.iloc[:,1:7].transpose()
            print('v25c_FAT1_nodes_0-120-240-b_set' + str(i+1) + '.txt')
            
        for i in range(48):
            data=pd.DataFrame()
            for j in [0,120,240]:
                fileObject = open(os.path.join(StructuralInputFile, 'v25c_FAT2_nodes_'+str(j)+'b_set' + str(i+1) + '.txt'), "r")
                tmp = pd.read_csv(fileObject, skiprows=(0,1), delim_whitespace=True)
                ### RawData is not save for space reasoning
                data=pd.concat([data,tmp])
            Out[:,i+42,:]=data.iloc[:,1:7].transpose()
            print('v25c_FAT2_nodes_0-120-240-b_set' + str(i+1) + '.txt')
                
        self.StructInput['StressTransferMatrix_T2']=Out
        
    def select_nodes_ULS_Lim(self, ULSStressLimit, Type='T1'):
       #### Function to select criticalnodes
       if Type=='T1':
           NodeStressInpt=os.path.join(self.LocSims,'FatigueInput','T1','v30c_blockD_seqv_maxovertime_allbodies.txt')
       elif Type=='T2':
           NodeStressInpt=os.path.join(self.LocSims,'FatigueInput','T2','v25c_blockD_seqv_maxovertime_allbodies.txt')
       fileObject=open(NodeStressInpt,"r")
       NodesStress=pd.read_csv(fileObject, sep='\t')
       NodesStress=NodesStress[NodesStress['Equivalent (von-Mises) Stress (MPa)']>ULSStressLimit]
       NodesStress=NodesStress.set_index('Node Number')
       # NodesSelection=list(NodesStress['Node Number'])
       # StructInput['SelectedNodes']=NodesStress
       NodeIndex=[]
       # StessSelection=[]
       print('StartNodeSelection '+ Type)
       #### Select only the surface nodes based reported nodes in the data output
       for i in self.StructInput['Nodes'+'_'+Type]:
           NodeIndex.append(i in NodesStress.index)
          
       self.StructInput['All Nodes Stresses Selection '+ Type]=NodesStress   #### All node slection includes non-surface nodes
       self.StructInput['Nodes Selection '+Type]=self.StructInput['Nodes'+'_'+Type][NodeIndex]
       self.StructInput['StressTransferMatrix Selection '+ Type]=self.StructInput['StressTransferMatrix'+'_'+Type][:,:,NodeIndex]
       self.StructInput['StressLimit Selection [MPa]' +Type]=ULSStressLimit
       print('End Nodes selection '+ Type)
       
    def random_node_selection(self, NNodes, Type='T1'):
        Nodes=self.StructInput['Nodes'+'_'+Type]
        NodeIndex=np.random.choice(range(len(Nodes)),size=NNodes, replace=False)                                       
        self.StructInput['Nodes Selection '+Type]=self.StructInput['Nodes'+'_'+Type][NodeIndex]
        self.StructInput['StressTransferMatrix Selection '+ Type]=self.StructInput['StressTransferMatrix'+'_'+Type][:,:,NodeIndex]
        print('Rndm Nodes selected '+ Type)
       
    def select_nodes_from_file(self, StructInput, SelectionFile, SRFFile=None, MaxNodes=None ):
        #### Function that can be used for nodes selection based on an input file
        fileObject = open(SelectionFile, "r")
        tmp = pd.read_csv(fileObject, skiprows=(0), delim_whitespace=True)
        FilterNodes=np.array(tmp['Node'])
        BaseNodes=np.array(StructInput['Nodes_T1'])
        STM=StructInput['StressTransferMatrix_T1']
        FilterInds=np.where(np.isin(BaseNodes,FilterNodes))[0]
        
        # NodesSelect=BaseNodes[np.isin(BaseNodes,FilterNodes)]
        if SRFFile:
            BaseNodes=np.take(BaseNodes, FilterInds)
            STM=np.take(STM, FilterInds, axis=-1)
            
            fileObject = open(SRFFile, "r")
            dfSRF = pd.read_csv(fileObject)
            dfSRF = dfSRF[np.isin(dfSRF['Node'],FilterNodes)]
            
            if MaxNodes: 
                #### Select max x nodes from SRF df -> Check ChatGTP
                FilterNodes=dfSRF.nsmallest(MaxNodes, 'SRF - min of three platforms')['Node']
            else: 
                FilterNodes=dfSRF['Node']
                
            FilterInds=np.where(np.isin(BaseNodes,FilterNodes))[0]
     
        StructInput['SelectionFile']=SelectionFile
        StructInput['Nodes Selection T1']=np.take(BaseNodes, FilterInds)
        StructInput['StressTransferMatrix Selection T1']=np.take(STM, FilterInds, axis=-1)
        StructInput['SRFSelectionFile']=SRFFile
        StructInput['NoMaxNodes']=MaxNodes
        del StructInput['StressTransferMatrix_T1']
        del StructInput['Nodes_T1']
        return StructInput
        
class FatCalc(GeneralFunctions):
    
    def __init__(self,LocSims):
        self.LocSims=LocSims
        self.FatDict={}
    
    def set_sn_properties(self, Dict=None):
        if Dict== None:
            self.K1=10**10.9699
            self.K2=10**13.6166
            self.m1=3
            self.m2=5
            self.S1=36
            self.Sq=21.05
            self.Nd=10**7
            self.Nc=2*10**6
            self.FATName='FAT 36'
        else:
            self.K1=Dict['K1']
            self.K2=Dict['K2']
            self.m1=Dict['m1']
            self.m2=Dict['m2']
            self.S1=Dict['S1']
            self.Sq=Dict['Sq']
            self.Nd=Dict['Nd']
            self.Nc=Dict['Nc']
            self.FATName=Dict['Name']
        ### Add a print line at some point here so that variable are shown in console
        # print('m1='+str(self.m1)+' m2='+str(self.m2))
            
    def plot_sn_curve(self):
        plt.figure('SN Curve')
        N1=10**(np.array([4,5,np.log10(self.Nd)]))
        N2=10**(np.array([np.log10(self.Nd),9]))
        S1=(self.K1/N1)**(1/self.m1)
        S2=(self.K2/N2)**(1/self.m2)
        plt.loglog(N1,S1)
        plt.loglog(N2,S2)
        plt.title(self.FATName)
        plt.xlabel('Cylcles Log(N)')
        plt.ylabel('Stress range (MPa)')    
        plt.grid()
       
        
    def calc_yearly_load_cycles(self,df, Type, Variable,  ScatterDict, DOF='', GeneralSelector=None, GeneralSelection=None):
        df=self.filterdf(df,['Variable'], [Variable])
        if GeneralSelector:
            df=self.filterdf(df,GeneralSelector, GeneralSelection)
        DOFs=[' X',' Y', ' Z', ' rX', ' rY', ' rZ', 'Horz']
        if DOF=='':
            Index=None
        else:
            Index=DOFs.index(DOF)
        MaxLoad=df['Max'+DOF].max()
        MinLoad=df['Min'+DOF].min()
        LoadEdges=np.linspace(0,(MaxLoad-MinLoad),100)
        LoadCens=[(LoadEdges[i]+LoadEdges[i+1])/2 for i in range(len(LoadEdges[:-1]))]
        
        for UN in df['UniqueName'].unique(): 
            if UN not in self.FatDict:
                self.FatDict[UN]={}
            self.FatDict[UN][DOF]={}
            self.FatDict[UN][DOF]['HourlyLoadHist']=np.zeros([len(self.Hss), len(self.Tps), len(self.WDirs),  len(LoadCens)])
            self.FatDict[UN][DOF]['LoadCens']=LoadCens
               
        for isim, sim in enumerate(self.Sims):
            ProData=pickle.load(open(sim,'rb'))
            #### Find the index of the index of bin
            iHs=np.where(self.Hss==np.round(ProData['Environment']['Hs'],2))[0][0]
            iTp=np.where(self.Tps==np.round(ProData['Environment']['Tp'],2))[0][0]
            iWDir=np.where(self.WDirsOF==np.round(ProData['Environment']['WaveDirection'],2))[0][0]
            
            for UN in df['UniqueName'].unique(): 
                # plt.figure(str(Hs)+str(Tp)+str(i), figsize=[15,10])
                #### Plot load spectrum based on rainflow count 
                #### This give the load range spectrum, no longer amplitude
                #### If statement added for platform variables as pro and raw data structure is different
                if UN.find('Platform')==-1: 
                    if Index:
                        Hist1= np.histogram(ProData[Type][UN][Variable]['HalfCycles'][Index], bins=LoadEdges)
                    else:
                        Hist1= np.histogram(ProData[Type][UN][Variable]['HalfCycles'], bins=LoadEdges)
                else:
                    Splt=UN.split('_')
                    Con='_'
                    UN1=Con.join(Splt[0:2])
                    UN1=UN1.replace('Platform','')
                    Hist1= np.histogram(ProData[Type][UN1][Variable][Splt[-1]]['HalfCycles'][Index], bins=LoadEdges)
                #### Go from half cylces to full cylces
                Hist2=Hist1[0]/2*(3600/self.SimulationDuration) ### Cycles per hour per sea-state
                # Hist3=plt.hist(LoadEdges[:-1], LoadEdges, weights=Hist2, alpha=0.3)
                # plt.yscale('log')
                self.FatDict[UN][DOF]['HourlyLoadHist'][iHs,iTp, iWDir,:]=Hist2
                
        for UN in df['UniqueName'].unique(): 
            #### Next step, multiply out the load matrix get load cylcles for the pins
            self.FatDict[UN][DOF]['YearlyLoadHistScatter']=self.FatDict[UN][DOF]['HourlyLoadHist']*ScatterDict['ScatterNorm'][:,:,:,None]*365*24/100 #### Set the amount of yearly cylces
            self.FatDict[UN][DOF]['YearlyLoadHist']=np.sum( self.FatDict[UN][DOF]['YearlyLoadHistScatter'], axis=(0,1,2))
    
    def calc_DELs(self):
        for i in self.FatDict.keys(): 
            for j in self.FatDict[i].keys():
                self.FatDict[i][j]['DEL']=(np.sum(self.FatDict[i][j]['YearlyLoadHist']*self.DesignLife*np.array(self.FatDict[i][j]['LoadCens'])**self.m1)/self.Nc)**(1/self.m1)
    
    def calc_yearly_stress_cylces(self, df, ScatterDict, SelectionName):
        FatDict={}
        df=self.filterdf(df, ['FATName'], SelectionName)
        MaxLoad=df['Max'].max()
        # MinLoad=df['Min'].min()
        LoadEdges=np.linspace(0,(MaxLoad-0),100)
        LoadCens=[(LoadEdges[i]+LoadEdges[i+1])/2 for i in range(len(LoadEdges[:-1]))]
        StressData=pickle.load(open(self.Sims[0],'rb'))
        Nodes=StressData[SelectionName]['General']['Nodes']
        HourlyStressHist=np.zeros([len(self.Hss), len(self.Tps), len(self.WDirs),  len(Nodes), len(StressData[SelectionName])-1,  len(LoadCens)])
        for isim, sim in enumerate(self.Sims):
            StressData=pickle.load(open(sim,'rb'))
            #### Find the index of the index of bin
            iHs=np.where(self.Hss==np.round(StressData['Environment']['Hs'],2))[0][0]
            iTp=np.where(self.Tps==np.round(StressData['Environment']['Tp'],2))[0][0]
            iWDir=np.where(self.WDirsOF==np.round(StressData['Environment']['WaveDirection'],2))[0][0]  
            for iPlat, Plat in enumerate(StressData[SelectionName].keys()):
                if Plat != 'General':
                    for iN, HalfCycles in enumerate(StressData[SelectionName][Plat]['HalfCyclesY']):
                        Hist1= np.histogram(HalfCycles, bins=LoadEdges) #### HalfCyclesY should be HalfCycles, this should renamed in a new version 
                        #### Go from half cylces to full cylces
                        Hist2=Hist1[0]/2*(3600/self.SimulationDuration) ### Cycles per hour per sea-state
                        # Hist3=plt.hist(LoadEdges[:-1], LoadEdges, weights=Hist2, alpha=0.3)
                        # plt.yscale('log')
                        HourlyStressHist[iHs,iTp, iWDir, iN, iPlat-1, :]=Hist2
                    # print(Plat)
            print(sim)
        # CycleHistYearly=HourlyStressHist*ScatterDict['ScatterNorm'][:,:,:,None,None, None]*365*24/100 #### Set the amount of yearly cylces
        # CylceHistTotalYearly=np.sum(CycleHistYearly, axis=(0,1,2))
        FatDict['General']=StressData[SelectionName]['General']
        FatDict['General']['LoadCens']=LoadCens
        FatDict['YearlyStressHistScatter']=HourlyStressHist*ScatterDict['ScatterNorm'][:,:,:,None,None, None]*365*24/100 #### Set the amount of yearly cylces
        FatDict['YearlyStressHist']=np.sum(FatDict['YearlyStressHistScatter'], axis=(0,1,2))
        return FatDict
        
    def calc_node_damage_and_SRF(self, FatDict):
        YearlyStressHist=FatDict['YearlyStressHist']
        StressCens=np.array(FatDict['General']['LoadCens'])
        
        #### Check SN curve Sq
        Sq = self.S1 * (self.Nc/self.Nd)**(1/self.m1)
        if round(Sq-self.Sq,2)==0:
            print('Sq matches')
        else:
            print('Sq does not match')
        
        Slope=np.where(StressCens<Sq, self.m2, self.m1)
        NMax= self.Nd * np.power(Sq/StressCens, Slope)
        
        DamageSpectrum=YearlyStressHist*self.DesignLife/NMax[None,None,:]
        Damage=np.sum(DamageSpectrum, axis=-1)
        
        SRFs = np.where(Damage <= self.D_ADL, (self.D_ADL/Damage)**(1/self.m2), (self.D_ADL/Damage)**(1/self.m1))
        
        return Damage, SRFs
            
        
        # fatigue_damage = calculate_damage(load_cycles, slope_low, slope_high, endurance_limit)
        # print(fatigue_damage)
        
    
    def get_DELs_T1(self, Plat):
        #### This fuction appends DELs in a order specified in load matrix given by Niels
        #### calc_DELs needs to be run first
        Out=[]
        for j in [0,120,240]:
            for DOF in [' X',' Y', ' Z', ' rX', ' rY', ' rZ']:
                Out.append(self.FatDict['Floater'+Plat+'_'+str(j)][DOF]['DEL'])
        for j in [0,120,240]:
            for DOF in [' X',' Y', ' Z']:
                try:
                    Out.append(self.FatDict['Con_Cen_'+Plat+'_'+str(j)+'_1'][DOF]['DEL'])
                except:     
                    Out.append(0)                               
            for DOF in [' X',' Y', ' Z']:
                try:
                    Out.append(self.FatDict['Con_Cen_'+Plat+'_'+str(j)+'_-1'][DOF]['DEL'])
                except:     
                    Out.append(0)                               
        for j in [0,120,240]:
            for DOF in [' X',' Y', ' Z']:
                try:
                    Out.append(self.FatDict['Con_Frld_'+Plat+'_'+str(j)+'_1'][DOF]['DEL'])
                except:     
                    Out.append(0)    
            for DOF in [' X',' Y', ' Z']:
                try:
                    Out.append(self.FatDict['Con_Frld_'+Plat+'_'+str(j)+'_-1'][DOF]['DEL'])
                except:     
                    Out.append(0)            
            
        for j in [0,120,240]:                           
            for DOF in [' rX',' rZ']:
                try:
                    Out.append(self.FatDict['Con_Cen_'+Plat+'_'+str(j)+'_1'][DOF]['DEL'])
                except:     
                    Out.append(0) 
            for DOF in [' rX',' rZ']:
                try:
                    Out.append(self.FatDict['Con_Cen_'+Plat+'_'+str(j)+'_-1'][DOF]['DEL'])
                except:     
                    Out.append(0) 
                    
        return np.array(Out)
        
    
    def extract_raos(self, dfStressRAOs,  TagName='',Overwrite=True,):
        #### This function extracts the RAOs from a dataframe and places them into a matrix it alos interpolates the matrix to the desired wavefrequencies
        SLS=PPC.StucturalLoadSelection()       
        
        #### Set fine range of wave periods to that work well for all spectra
        WavePeriods=np.concatenate([np.arange(2.5,14.5,0.25),np.arange(15,20,0.5), np.arange(20,30,3)])#, np.arange(20,30,3)
        
        df1m= dfStressRAOs[dfStressRAOs['HMax']==1] 
        dfBr= dfStressRAOs[dfStressRAOs['HMax']>1] 
        Nodes=dfStressRAOs.columns[16:-1]
        
        WaveDirections=np.sort(dfBr['WaveDirection'].unique())
        WavePeriodsOrcaFlex=np.sort(dfBr['WavePeriod'].unique())
        
        Overwrite=False
        
        if os.path.isfile(os.path.join(self.LocSims,'StressRAOs_'+TagName+'.pkl')) and not Overwrite:
            print(TagName+ ' Truss RAOs Loaded')
            DictFatigue=pickle.load(open(os.path.join(self.LocSims,'StressRAOs_'+TagName+'.pkl'),'rb'))
        else:
            RAO=np.zeros([2, len(dfBr['Platform'].unique()), len(Nodes), len(WaveDirections),len(WavePeriodsOrcaFlex)])            
            for idf, df in enumerate([dfBr,df1m]):#, 
                for idirec, direc in enumerate(WaveDirections):
                    print(str(direc))
                    for iplt, plat in enumerate(dfBr['Platform'].unique()):
                        print(str(plat))
                        for it, t in enumerate(WavePeriodsOrcaFlex):
                            # print(time.perf_counter())
                            tmp1=SLS.filterdf(df, ['WavePeriod', 'WaveDirection', 'Platform'], [t,direc,plat], Type='Exact')
                            # print(time.perf_counter())
                            RAO[idf,iplt,:, idirec,it]=tmp1.iloc[0,16:-1]
                            # print(time.perf_counter())
            
            print(time.perf_counter())   
            
            Platforms = [float(i) for i in dfBr['Platform'].unique()]
            
            #### Interpolate RAOS 
            #### Interpolation has been included again (in other version it was removed)
            #### Interpolation took the most time of this script
            Points=(np.array(Platforms) ,range(len(Nodes)),dfBr['WaveDirection'].unique(),WavePeriodsOrcaFlex)
                
            PointsOut=np.stack(np.meshgrid(Platforms,range(len(Nodes)),dfBr['WaveDirection'].unique(),WavePeriods))
            PointsOut1=np.moveaxis(PointsOut,[0,2],[4,0]) 
            
            DictFatigue={}    
            print(time.perf_counter())   
            DictFatigue['RAO_Br']=scipy.interpolate.interpn(Points, np.array(RAO[0,:,:,:,:]),PointsOut1)     
            # DictFatigue['RAO_Br']=np.array(RAO[0,:,:,:,:])
            print(time.perf_counter())   
            DictFatigue['RAO_1m']=scipy.interpolate.interpn(Points, np.array(RAO[1,:,:,:,:]),PointsOut1)     
            # DictFatigue['RAO_1m']=np.array(RAO[1,:,:,:,:])  
            print(time.perf_counter())   
             
            # HMaxBr=pd.pivot_table(dfBr, values=['HMax'], index=['WavePeriod'], aggfunc=np.mean)    
            # HMax1m=pd.pivot_table(dfBr, values=['HMax'], index=['WavePeriod'], aggfunc=np.mean)                        
            # WavePeriods=WavePeriodsOrcaFlex
            
            DictFatigue['General']={}
            DictFatigue['General']['WaveDirections']=WaveDirections
            DictFatigue['General']['WavePeriods_OrcaFlex']=WavePeriodsOrcaFlex
            DictFatigue['General']['WavePeriods']=WavePeriods
            DictFatigue['General']['WaveFrequencies']=2*np.pi/WavePeriods
            DictFatigue['General']['WaveFrequencies_OrcaFlex']=2*np.pi/WavePeriodsOrcaFlex
            DictFatigue['General']['WaveHeights_OrcaFlexBr']=np.array(pd.pivot_table(dfBr, values=['HMax'], index=['WavePeriod'], aggfunc=np.mean)['HMax'])
            DictFatigue['General']['WaveHeights_Br']=np.interp(WavePeriods,WavePeriodsOrcaFlex,DictFatigue['General']['WaveHeights_OrcaFlexBr'])
            DictFatigue['General']['WaveHeights_1m']=np.ones(np.shape(DictFatigue['General']['WaveHeights_Br']))
            DictFatigue['General']['Platforms']=Platforms
            DictFatigue['General']['Nodes']=Nodes
            pickle.dump(DictFatigue,open(os.path.join(self.LocSims,'StressRAOs_'+TagName+'.pkl'),'wb'))
            print('StressRAOs_'+TagName+'.pkl saved')
        
        print(time.perf_counter())    
        return DictFatigue
        
    ### for read scatter please use the PPC.FDcalculations class 
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
        Scatter[Scatter=='-']=0
        Scatter=np.asarray(Scatter,dtype=float)
        # self.Scatter=Scatter
        
        return Hss, Tps, Scatter
                  
    def gen_spect_dict(self,Hss,Tps,Scatter,WaveFrequencies):
        import spectral
        SpecDict={}
        for iTp, Tp in enumerate(Tps):
            for iHs, Hs in enumerate(Hss):
                if Scatter[iHs,iTp]>0:
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]=spectral.Spectrum('JONSWAP', Hs, 'Tp', Tp, WaveFrequencies, gamma_val=3.3)
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['Hs']=Hs
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['Tp']=Tp
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['TpInd']=iTp
                    SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['HsInd']=iHs
        #             plt.plot(SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['frequencies'], SpecDict['Hs_'+str(Hs)+'_Tp_'+str(Tp)]['density'])
        # plt.xlabel('Wave frequency [rad/s]')
        # plt.ylabel('Spectral density [m^2 s/rad]')
        # plt.grid()
        print('SpecDict Generated')
        SpecDict['General']={}
        SpecDict['General']['Hss']=Hss
        SpecDict['General']['Tps']=Tps
        SpecDict['General']['Scatter']=Scatter
        SpecDict['General']['WaveFrequecies']=WaveFrequencies
        return SpecDict    
            
    def calc_fd_responses_and_damage(self, SpecDict, DictFatigue, TagName=''):
        #### This function calculates the frequency domain responses and damages for all sea-states in the spectdict
        #### Responses and damages are summarized in a large excel
        WavePeriods=DictFatigue['General']['WavePeriods']
        WaveFrequencies=DictFatigue['General']['WaveFrequencies']
        WaveDirections=DictFatigue['General']['WaveDirections']
        Platforms=DictFatigue['General']['Platforms']
        Tps=SpecDict['General']['Tps']
        Hss=SpecDict['General']['Hss']
        Scatter=SpecDict['General']['Scatter']
        Nodes=DictFatigue['General']['Nodes']
                
        WH1m=DictFatigue['General']['WaveHeights_1m']
        WHBrLm=DictFatigue['General']['WaveHeights_Br']
        WHBrLmTpss=np.interp(0.9*Tps,  WavePeriods,WHBrLm)
                       
        DictFatigue['HsTp']={}
        DictFatigue['SigSummary']=np.zeros([len(Hss),len(Tps), len(Platforms), len(Nodes),len(WaveDirections)])
        DictFatigue['TzSummary']=np.zeros(np.shape(DictFatigue['SigSummary']))
        # DictFatigue['CyclesPerHourSummary']=np.zeros([len(Tps),len(Hss),len(Platforms), len(Nodes),len(WaveDirections)])
        DictFatigue['DamagePerHour']=np.zeros([len(Hss),len(Tps), len(Platforms), len(Nodes),len(WaveDirections)])
        
        #### Test
        DictFatigue['RAO_Br']=np.array(DictFatigue['RAO_Br'],dtype=float)
        DictFatigue['RAO_1m']=np.array(DictFatigue['RAO_1m'],dtype=float)
        
        # for HsTp in ['Hs_0.25_Tp_10.5']: # SpecDict.keys():
        for HsTp in SpecDict.keys():
            if HsTp!='General':
                print(HsTp)
                # tmp=SpecDict[HsTp]
                #### Interp Hs is based on the ratio 2(Hs-H1m)/(HBrLm-H1m) H1m is hardcoded as 1
                Ratio=(SpecDict[HsTp]['Hs']*2-1)/(WHBrLmTpss[SpecDict[HsTp]['TpInd']]-1)
                Ratio=np.max([Ratio,0])
                Ratio=np.min([Ratio,1])
                
                RAOintp=(DictFatigue['RAO_Br']-DictFatigue['RAO_1m'])*Ratio+DictFatigue['RAO_1m']
                # tmp['RAOintp']=(DictFatigue['RAO_Br']-DictFatigue['RAO_1m'])*Ratio+DictFatigue['RAO_1m']
                
                RespDens=RAOintp**2*SpecDict[HsTp]['density']                
                # tmp['RespDens']=tmp['RAOintp']**2*SpecDict[HsTp]['density']
                
                # Nnode=22300
                # iplatform=0
                # fig,ax=plt.subplots()
                # ax.set_title('RAOs Platform'+ str(round(Platforms[iplatform]))+ ' Node '+str(Nnode) + ' '+ HsTp)
                # ax.plot(WaveFrequencies, RAOintp[iplatform,Nnode,:,:].transpose())
                # ax.legend(WaveDirections)
                # ax.set_ylabel('Stress range amp / meter wave amplitudue [MPa/m]')
                # ax.set_xlabel('Wave frequency [rad/s]')
                # ax.grid()
                # # ax.legend()
                # ax2=ax.twinx()
                # ax2.plot(WaveFrequencies,SpecDict[HsTp]['density'],'--')
                # ax2.set_ylabel('Spectral density [m^2 s/rad]')
                # plt.figure(2)
                # plt.plot(WaveFrequencies,RespDens[iplatform,Nnode,:,:].transpose())
                # plt.xlabel('Wave frequency [rad/s]')
                # plt.ylabel('Response density [MPa^2 s/rad]')
                # plt.legend(WaveDirections)
                # plt.title('Response Platform'+ str(round(Platforms[iplatform]))+ ' Node '+str(Nnode)+ ' '+ HsTp)
                # plt.grid()
                
                
                # tmp['m2']=np.trapz(tmp['RespDens']*WaveFrequencies**2, WaveFrequencies)*-1
                m2=np.trapz(RespDens*WaveFrequencies**2, WaveFrequencies)*-1
                
                # tmp['m0']=np.trapz(tmp['RespDens'],WaveFrequencies)*-1
                m0=np.trapz(RespDens,WaveFrequencies)*-1
                
                # tmp['m0']=np.trapz(np.flip(tmp['RespDens'],axis=1),np.flip(WaveFrequencies))
                
                # tmp['Tz']=2*np.pi*(tmp['m0']/tmp['m2'])**0.5
                Tz=2*np.pi*(m0/m2)**0.5
                
                # plt.plot(tmp['m0'][0,:,0])
                # tmp['CyclesPerHour']=3600/tmp['Tz']
                
                # tmp['Sig']=2*(tmp['m0']**0.5)
                Sig=2*(m0**0.5)
                
                # tmp['DamagePer3H']=3600*3/tmp['Tz']*((1/self.K1)*((2*(2*tmp['m0'])**0.5)**self.m1)*sc.gamma(1+self.m1/2)*sc.gammaincc(1+self.m1/2, self.Sq**2/(8*tmp['m0'])) + (1/self.K2)*((2*(2*tmp['m0'])**0.5)**self.m2)*sc.gamma(1+self.m2/2)*sc.gammainc(1+self.m2/2, self.Sq**2/(8*tmp['m0'])))
                DamagePer3H=3600*3/Tz*((1/self.K1)*((2*(2*m0)**0.5)**self.m1)*sc.gamma(1+self.m1/2)*sc.gammaincc(1+self.m1/2, self.Sq**2/(8*m0)) + (1/self.K2)*((2*(2*m0)**0.5)**self.m2)*sc.gamma(1+self.m2/2)*sc.gammainc(1+self.m2/2, self.Sq**2/(8*m0)))
            
                # tmp['DamagePer3H']=3600*3/tmp['Tz']*((1/self.K1)*((2*(2*tmp['m0'])**0.5)**self.m1)*sc.gamma(1+self.m1/2)*sc.gammaincc(1+self.m1/2, self.Sq**2/(8*tmp['m0'])))
                # tmp['DamagePer3H']=3600*3/tmp['Tz']*((1/self.K2)*((2*(2*tmp['m0'])**0.5)**self.m2)*sc.gamma(1+self.m2/2)*sc.gammainc(1+self.m2/2, self.Sq**2/(8*tmp['m0'])))
                                
                # DictFatigue['HsTp'][HsTp]=tmp.copy()
                # DictFatigue['HsTp'][HsTp]={}
                # DictFatigue['SigSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['Sig']
                # DictFatigue['TzSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['Tz']
                # DictFatigue['CyclesPerHourSummary'][SpecDict[HsTp]['TpInd'],SpecDict[HsTp]['HsInd'],:]=tmp['CyclesPerHour']
                # DictFatigue['DamagePerHour'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=tmp['DamagePer3H']/3
               # 
                DictFatigue['SigSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=Sig
                DictFatigue['TzSummary'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=Tz
                # DictFatigue['CyclesPerHourSummary'][SpecDict[HsTp]['TpInd'],SpecDict[HsTp]['HsInd'],:]=tmp['CyclesPerHour']
                DictFatigue['DamagePerHour'][SpecDict[HsTp]['HsInd'],SpecDict[HsTp]['TpInd'],:]=DamagePer3H/3    
                
        print(time.perf_counter())                    
        
        #### Complete the fatigue dictionary to calculate the short term damages #####
        tmp1=np.moveaxis(DictFatigue['DamagePerHour'].mean(axis=4), [0,1], [2,3])
        
        ScatterNormalized=Scatter/np.sum(Scatter)
        DictFatigue['DamagePerYearScatter']=tmp1*ScatterNormalized*365*24
        DictFatigue['TotalDamagePerYear']=np.sum(np.sum(DictFatigue['DamagePerYearScatter'],axis=3),axis=2)
        DictFatigue['LifeTime']=1/DictFatigue['TotalDamagePerYear']
        
        pickle.dump(DictFatigue,open(os.path.join(self.LocSims,'FatigueDict'+TagName+ '.pkl'),'wb'))
        
        return DictFatigue
        
    def to_csv(self, DictFatigue, TagName):
        #### Old version
        #### Send the max damage of all platforms to csv
        DamageOut=DictFatigue['TotalDamagePerYear'].max(axis=0)
        Out=np.stack([DictFatigue['General']['Nodes'],DamageOut.transpose()])
        df=pd.DataFrame(Out.transpose(), columns=['Nodes','DamagePerYearMax'])
        
        df.to_csv(os.path.join(self.LocSims,'FatigueOutput_'+TagName+'.csv'), index=False)
        pickle.dump(df,open(os.path.join(self.LocSims,'FatigueOutput_'+TagName+ '.pkl'),'wb'))
    
    def to_csv_damage_and_SRF(self, DictFatigue, Damage, SRFs, SelectionName):
        #### Old version
        #### Send the max damage of all platforms to csv
        DamageOut=Damage.max(axis=-1)
        SRFOut=SRFs.min(axis=-1)
        Out=np.stack([DictFatigue['General']['Nodes'],DamageOut.transpose(), SRFOut.transpose()])
        df=pd.DataFrame(Out.transpose(), columns=['Node','Damage', 'SRF'])
        df.to_csv(os.path.join(self.LocSims,'FatigueOutput_'+SelectionName+'.csv'), index=False)
        # pickle.dump(df,open(os.path.join(self.LocSims,'FatigueOutput_'+TagName+ '.pkl'),'wb'))
        
    def append_excel_sheet(self, DataFrame, ExcelFile, SheetName='Sheet1'):
        if os.path.exists(ExcelFile):
            book = load_workbook(ExcelFile)
            writer = pd.ExcelWriter(ExcelFile,  engine = 'openpyxl')
            writer.book = book
               
            # writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
            
            DataFrame.to_excel(writer, SheetName)
            writer.save()
        else: 
            DataFrame.to_excel(ExcelFile, sheet_name=SheetName)
    

class Plot:
    
    def __init__(self,LocSims):
        self.LocSims=LocSims
        
    def get_xyz(self, dfFatigueOut):
        fileObject = open(os.path.join(self.LocSims,'FatigueInput','T1','v27b_FAT1_nodes_b_nodeXYZ.txt'), "r")
        tmp = pd.read_csv(fileObject, skiprows=(0,1,2,3,4), delim_whitespace=True)

        tmp2=tmp.set_index('NODE')

        Nodes=list(dfFatigueOut['Nodes'])
        X=list(tmp2['X'][Nodes])
        Y=list(tmp2['Y'][Nodes])
        Z=list(tmp2['Z'][Nodes])
        
        dfFatigueOut['X']=X
        dfFatigueOut['Y']=Y
        dfFatigueOut['Z']=Z
                
        return dfFatigueOut
    
    def scatter(self, df, Key):
        fig=plt.figure()
        ax = plt.axes(projection ="3d")    
        p=ax.scatter3D(df['X'], df['Y'], df['Z'], c=df[Key],vmin=0, vmax=1,s=0.05)
        # ax.set_zlim([-7,0])
        # ax.set_xlim([-3.5,3.5])
        # ax.set_ylim([-3.5,3.5])
        fig.colorbar(p)
 
class SN_Library:
    def __init__(self):
        self.Dict={}
        #### FAT 36
        Dict={}
        Dict['K1']=10**10.9699
        Dict['K2']=10**13.6166
        Dict['m1']=3
        Dict['m2']=5
        Dict['S1']=36
        Dict['Sq']=21.05
        Dict['Nd']=10**7
        Dict['Nc']=2*10**6 #### Nref used for DEL calculation
        Dict['Name']='FAT 36'
        self.Dict[Dict['Name']]=Dict
        #### FAT 71
        Dict1={}
        Dict1['K1']=10**11.855
        Dict1['K2']=10**15.091
        Dict1['m1']=5
        Dict1['m2']=5
        Dict1['S1']=71
        Dict1['Sq']=41.52
        Dict1['Nd']=10**7
        Dict1['Nc']=2*10**6 #### Nref used for DEL calculation
        Dict1['Name']='FAT 71'
        self.Dict[Dict1['Name']]=Dict1
        #### FAT 91
        Dict2={}
        Dict2['K1']=10**12.1818
        Dict2['K2']=10**15.6363
        Dict2['m1']=3
        Dict2['m2']=5
        Dict2['S1']=91.25
        Dict2['Sq']=53.36
        Dict2['Nd']=10**7
        Dict2['Nc']=2*10**6 #### Nref used for DEL calculation
        Dict2['Name']='FAT 91.25'
        self.Dict[Dict2['Name']]=Dict2
