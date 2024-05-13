# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 09:02:19 2023

@author: SanderMorshuisSolarD
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import os
import pickle
import datetime

from scipy.stats import norm
from scipy.signal import find_peaks
import OrcFxAPI as ofx



class ProcessWaveData:
    def __init__(self):
        self.CMap=plt.get_cmap("tab10")
        self.ColorIndex=0
        self.Distribution=ofx.evdWeibull
        self.Tail=ofx.exUpperTail
        self.TresholdReturnPeriod=0.5 ### Years
        self.DeclusterPeriodHour=3
    
    def calc_weibull_return_periods(self,ReturnPeriodsYear, TimeTrace, dTimeStepHour ):    
        #### Return
        #### Set threshold limit at 90%    
        def calc_probability(dtTotal, Extremes):
            return (np.flip(np.arange(len(Extremes))+1))/(dtTotal) ### x time in ... hours
        
        # Sorted=np.sort(TimeTrace)
        # Occurance=np.array(range(len(Sorted)))/len(Sorted)
        # ReturnPeriod=1/(1-Occurance)/(365*24)
        # Probability=calc_probability(dTimeStepHour*len(TimeTrace), Sorted)
        # ReturnPeriod1=1/Probability/(365*24)
        
        #### Find independent extreme events 
        PeaksInd, Peaks =find_peaks(np.array(TimeTrace), distance=self.DeclusterPeriodHour/dTimeStepHour)
        # plt.plot(df['DateTime'], TimeTrace)
        # plt.scatter(df['DateTime'][PeaksInd], TimeTrace[PeaksInd])
        PeaksSort=np.sort(TimeTrace[PeaksInd])
        ProbabilityPeaks=calc_probability(dTimeStepHour*len(TimeTrace), PeaksSort)
        ReturnPeriodPeaks=1/ProbabilityPeaks/(365*24)
        # Occurance=(np.array(range(len(PeaksSort)))+len(TimeTrace)-len(PeaksSort))/len(TimeTrace)
        # ReturnPeriod1=1/(1-Occurance)/(365*24)
        
        # plt.scatter(df['DateTime'],HsPeak.peaks)
        PeakSelect=PeaksSort[ReturnPeriodPeaks>self.TresholdReturnPeriod] #### only top 5 event per year same assumption as is made by DHI #### 
        Threshold=PeakSelect[0]
        print(Threshold)
        
        stats = ofx.ExtremeStatistics(TimeTrace, dTimeStepHour*3600)
        distribution = self.Distribution
        declusterPeriod = self.DeclusterPeriodHour*3600 ### Something is still going wrong with the declusting period of find peaks vs OrcaFlex
        #### Need to sort out which data OrcaFlex uses to fit extremes
        confidenceLevel = 95 # percentage
        spec = ofx.LikelihoodStatisticsSpecification(distribution, Threshold, declusterPeriod, self.Tail)
        # query results
        stats.Fit(spec)
        StormDurations=np.array(ReturnPeriodsYear)*365*24 
        ReturnLevel=[]
        ReturnLevelRawData=[] #Hs1yrRtPDirectional.
        for isD, stormDuration in enumerate(StormDurations):
            query = ofx.LikelihoodStatisticsQuery(stormDuration, confidenceLevel)
            stats.Fit(spec)    
            extremes = stats.Query(query)
            print('(1) {} hour return level = {}'.format(stormDuration,
            extremes.ReturnLevel))
            # print('(2) {}% confidence limits = [{}, {}]'.format(confidenceLevel,
            # extremes.ConfidenceInterval.Lower, extremes.ConfidenceInterval.Upper))
            # print('(3) sigma = {} (se = {})'.format(extremes.Sigma, extremes.SigmaStdError))
            # print('(4) xi = {} (se = {})'.format(extremes.Xi, extremes.XiStdError))
            ReturnLevel.append(extremes.ReturnLevel)
            ReturnLevelRawData.append(PeaksSort[np.argmin((ReturnPeriodPeaks-ReturnPeriodsYear[isD])**2)])
        
        return ReturnLevel, PeaksSort, ReturnPeriodPeaks, Threshold, ReturnLevelRawData
    
class Plot:
    def __init__(self):
        self.CMap=plt.get_cmap("tab10")
        self.PlotPath=r'X:\Development\SD02_Merganser\01_Hydromechanics\Reference Info\06_WeatherReferences\10_HKZ-Scheveningen'
    
    def plt_scatter(self, df, xlabel, ylabel, title = "", save_fig=False):
        fig,ax = plt.subplots()
        plt.scatter(df[xlabel],df[ylabel], s=0.1)
        plt.grid()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        if save_fig:
            plt.savefig(os.path.join(self.PlotPath, ylabel +' vs '+ xlabel+ '.png'))
        return fig, ax
    
    def plt_desinty_scatter(df, xlabel, ylabel ):
        cmap = LinearSegmentedColormap.from_list("mycmap", [(0, "blue"), (0.9/len(df[xlabel]), "green"), (1, "red")])
        plt.hexbin(df[xlabel],df[ylabel], cmap=cmap)
        plt.colorbar()
        plt.grid()
    
    def plt_polar_scatter(df, xlabel, ylabel):
        fig = plt.figure()
        ax = fig.add_subplot(projection='polar')
        ax.scatter(df[xlabel]*np.pi/180,df[ylabel], s=0.1)
        plt.grid()
        # plt.xlabel(xlabel)
        # plt.ylabel(ylabel)
        
        # PlotPath=r'X:\Development\SD02_Merganser\01_Hydromechanics\Reference Info\06_WeatherReferences\10_HKZ-Scheveningen'

class Wave_Filter:
    def __call__(self, np_hist, Tp_upper_filter_col, Tp_lower_filter_col, do_period_filter, Wave_min_Percentage = 0, Direction_Threshold = 0 ) :

        print(f"\n{'Total number of cases before filtering =':>50} {np.count_nonzero(np_hist)}")
        print(f"{'Before filtering sum check =':>50} {np_hist.sum().round(1)}")

        filtered_hist = self.Tp_Filter(np_hist, Tp_upper_filter_col, Tp_lower_filter_col, do_period_filter)
        filtered_hist = self.Dir_Filter(filtered_hist, Wave_min_Percentage, Direction_Threshold)
        
        return filtered_hist
    
    def Tp_Filter(self, np_hist, Tp_upper_filter_col, Tp_lower_filter_col, do_period_filter = False):
        
        filtered_hist = np_hist.copy()
        if do_period_filter:
            for iDir in range(filtered_hist.shape[2]):     
                filtered_hist[:,  Tp_upper_filter_col-1, iDir] = filtered_hist[:, Tp_upper_filter_col-1:, iDir].sum(axis=1)
                filtered_hist[:,  Tp_upper_filter_col: , iDir] = 0
                filtered_hist[:,  Tp_lower_filter_col+1, iDir] = filtered_hist[:, :Tp_lower_filter_col+2, iDir].sum(axis=1)
                filtered_hist[:, :Tp_lower_filter_col+1, iDir] = 0

            print(f"{'Total number of cases after filtering Tp =':>50} {np.count_nonzero(filtered_hist)}")
            print(f"{'Tp filter sum check =':>50} {filtered_hist.sum().round(1)}")

        return filtered_hist
    
    def Dir_Filter(self, filtered_hist, Wave_min_Percentage = 0, Direction_Threshold = 0):
        if Direction_Threshold > 0:
            sum_hist_filtered = np.sum(filtered_hist, axis=(0, 1))
            int_Dir = np.argmin(sum_hist_filtered)
            do_while = True
            while do_while:
                # print(int_Dir, sum_hist_filtered)

                if int_Dir == len(sum_hist_filtered):
                    change_Dir_right = 0
                else:
                    change_Dir_right = int_Dir+1
                change_Dir_left = int_Dir-1

                if   sum_hist_filtered[change_Dir_left]  != 0 and sum_hist_filtered[change_Dir_right] != 0:
                    filtered_hist[:, :, change_Dir_left]  += 0.5 * filtered_hist[:, :, int_Dir]
                    filtered_hist[:, :, change_Dir_right] += 0.5 * filtered_hist[:, :, int_Dir]                
                elif sum_hist_filtered[change_Dir_left]  == 0 and sum_hist_filtered[change_Dir_right] != 0:
                    filtered_hist[:, :, change_Dir_right] += 1.0 * filtered_hist[:, :, int_Dir]
                elif sum_hist_filtered[change_Dir_right] == 0 and sum_hist_filtered[change_Dir_left]  != 0:
                    filtered_hist[:, :, change_Dir_left]  += 1.0 * filtered_hist[:, :, int_Dir]
                elif sum_hist_filtered[change_Dir_left]  == 0 and sum_hist_filtered[change_Dir_right] == 0:
                    filtered_hist[:, :, change_Dir_left]  += 0.5 * filtered_hist[:, :, int_Dir]
                    filtered_hist[:, :, change_Dir_right] += 0.5 * filtered_hist[:, :, int_Dir]                
                    # print('Both sides 0 percent!')
                    
                filtered_hist[:,:,int_Dir]       = 0
                
                sum_hist_filtered[change_Dir_left]  = filtered_hist[:,:,change_Dir_left ].sum().round(2)
                sum_hist_filtered[change_Dir_right] = filtered_hist[:,:,change_Dir_right].sum().round(2)

                sum_hist_filtered[int_Dir]          = 0
                
                tmp_sum_hist = np.array([101 if x == 0 else x for x in sum_hist_filtered])
                int_Dir  = np.argmin(tmp_sum_hist)
                do_while = min(tmp_sum_hist) < Direction_Threshold
                # print(int_Dir, sum_hist_filtered)
        
        print(f"{'Total number of cases after filtering Dir =':>50} {np.count_nonzero(filtered_hist)}")
        print(f"{'Dir filter sum check =':>50} {filtered_hist.sum().round(1)}")

        Wave_min_Ratio = Wave_min_Percentage/100
        while_count = 0
        while np.any((0.0 < filtered_hist) & (filtered_hist < Wave_min_Ratio)):
            while_count += 1
            # print('while_count', while_count)
            for iDir in range(filtered_hist.shape[2]):  
                for iTp in range(filtered_hist.shape[1]):  
                    for iHs in range(filtered_hist.shape[0]):  
                        # print('Checking:', filtered_hist[iHs, iTp, iDir].round(2), iHs, iTp, iDir)
                        if 0.0 < filtered_hist[iHs, iTp, iDir] < Wave_min_Ratio:
                            # print('modifynig:', filtered_hist[iHs, iTp, iDir].round(2), iHs, iTp, iDir)
                            if iHs == 0:
                                if filtered_hist[iHs, iTp+1, iDir] == 0:
                                    filtered_hist[iHs, iTp-1, iDir] += filtered_hist[iHs, iTp, iDir]  
                                else:
                                    filtered_hist[iHs, iTp+1, iDir] += filtered_hist[iHs, iTp, iDir]  
                            else:
                                filtered_hist[iHs-1, iTp, iDir] += filtered_hist[iHs, iTp, iDir]  
                                
                            filtered_hist[iHs, iTp, iDir]    = 0.0
        
        print(f"{'Total number of cases after filtering low % wave =':>50} {np.count_nonzero(filtered_hist)}")
        print(f"{'Low % wave filter sum check =':>50} {filtered_hist.sum().round(1)}")
        return filtered_hist

class Create_Yml_Files():
    
    def Generate(self, Base_File_Name, Base_File_Path, np_hist, wave_direction_bins, wave_height_bins, wave_period_bins, do_Yml = False, do_Wave_Spreading = False):
    
        if do_Yml:
            print('\n')
            BaseFile  = os.path.join(Base_File_Path, Base_File_Name)
            BaseModel = ofx.Model(BaseFile)

            Fatigue_File_Start = '_'.join(Base_File_Name.split('_')[0:5])
            Case_No = 0
            
            Wdp = BaseModel["Environment"].WaterDepth
            BaseModel["General"].ImplicitConstantTimeStep           = 0.05
            BaseModel["General"].TargetLogSampleInterval            = 0.1
            BaseModel["General"].StageDuration[0]                   = 20
            BaseModel["General"].StageDuration[1]                   = 1800
            BaseModel["Environment"].WaveType                       = 'JONSWAP'
            BaseModel["Environment"].WaveJONSWAPParameters          = 'Partially specified'
            BaseModel["Environment"].WaveGamma                      = 3.3
            # BaseModel["Environment"].WaveKinematicsCutoffDepth    = Wdp-1
            
            BaseModel["Environment"].CurrentSpeedAtSurface          = 0 
            BaseModel["Environment"].CurrentSpeedAtSeabed           = 0
            BaseModel["Environment"].RefCurrentDirection            = 0
            BaseModel["Environment"].WindSpeed                      = 0
            BaseModel["Environment"].WindDirection                  = 0
            
            if do_Wave_Spreading:
                BaseModel["Environment"].WaveNumberOfComponents         = 100
                BaseModel["Environment"].WaveNumberOfSpectralDirections = 7
                BaseModel["Environment"].WaveDirectionSpreadingExponent = 8
                Spreading_Text = '_sp'
            else:
                BaseModel["Environment"].WaveNumberOfComponents         = 200
                BaseModel["Environment"].WaveNumberOfSpectralDirections = 1
                Spreading_Text = ''
                            
            for iDir in range(len(wave_direction_bins[:-1])):
                Dir = 0.5 * (wave_direction_bins[iDir] + wave_direction_bins[iDir+1])
                
                for iHs in range(len(wave_height_bins[:-1])):
                    Hs = 0.5 * (wave_height_bins[iHs] + wave_height_bins[iHs+1])
                    
                    for iTp in range(len(wave_period_bins[:-1])):
                        Tp = 0.5 * (wave_period_bins[iTp] + wave_period_bins[iTp+1])
                        
                        if np_hist[iHs, iTp, iDir] != 0:
                            Case_No += 1
                            Ofx_Dir = np.mod((360-Dir)+90,360)
                            Fatigue_File_Name = os.path.join(Base_File_Path, Fatigue_File_Start + "_Hs"+str(Hs)+"_Tp"+str(Tp)+"_Dir"+str(Ofx_Dir)+".yml")
                            
                            BaseModel["Environment"].WaveDirection                  = Ofx_Dir
                            BaseModel["Environment"].WaveHs                         = Hs
                            BaseModel["Environment"].WaveTp                         = Tp

                            BaseModel.SaveData(Fatigue_File_Name)
                            print('Saved: ',Case_No, Dir, Hs, Tp, Ofx_Dir)
                                
    