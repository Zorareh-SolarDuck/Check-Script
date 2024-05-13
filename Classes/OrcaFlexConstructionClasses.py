# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 15:54:11 2021

@author: Sander
"""
## Classes that are use for construction or setup of OrcaFlex models
import OrcFxAPI as ofx
import numpy as np
import sys
sys.path.insert(1, 'Classes')
import InputClasses as Input
from multiprocessing import Pool
import os
import glob

class Activate:
    
    def all_input(InputClasses):
        
        # Activate Environment
       Environment=InputClasses.Environment()
       Environment.calc_details_steep_seastate()
       
       General=InputClasses.General()
    
       # Activate Platform and floaters
       Platform=InputClasses.Platform() 
       Platform.calc_details_basic()
       Floater=InputClasses.Floater() 
       Floater.calc_details(Platform.MassT1)
       
       DragPlate=InputClasses.DragPlate()
       DragPlate.calc_details(Platform, Floater)
       
       # Activate Coupling and mooring + extra calcualtions 
       Coupling=InputClasses.Coupling()      
       Coupling.calc_details()    
       
       ## Class used to calculate the farm layout and associated geometrical functions
       Layout=InputClasses.Layout()
       Layout.calc_details(Platform, Coupling) 
     
       Mooring=InputClasses.Mooring(Platform, Coupling, Environment)
       Mooring.calc_details(Floater, Layout)
       Mooring.calc_connected_platforms(Layout)
           
       return General, Environment, Platform, Floater, DragPlate, Coupling, Mooring, Layout
  

## Class with functions that generate the components(lego blocks) that are needed to compose the OrcaFlex model
class Gen: 
    
    def __init__(self, model):
        self.model=model
        
    def set_CoG(self, Object, CoG):
        Object.CentreOfMassX=CoG[0]
        Object.CentreOfMassY=CoG[1]
        Object.CentreOfMassZ=CoG[2]
        
    def set_inertia(self,Object, I):
        ## I is the 6 dof inertia matrix 
        Object.MassMomentOfInertiaX=I[0,0]
        Object.MassMomentOfInertiaY=I[1,1]
        Object.MassMomentOfInertiaZ=I[2,2]
        
    def set_drag_area(self, Object, Array):
        Object.DragAreaX=Array[0]
        Object.DragAreaY=Array[1]
        Object.DragAreaZ=Array[2]
    def set_drag_Cd_6Dbuoy(self, Object, Array):
        Object.DragForceCoefficientX=Array[0]
        Object.DragForceCoefficientY=Array[1]
        Object.DragForceCoefficientZ=Array[2]
    
    def platform(self, Platform, name, Body):               
        OFPlatform=self.model.CreateObject(ofx.ot6DBuoy)
        OFPlatform.Name='Platform_' + name
        ## Set Drawing 
        OFPlatform.VertexX=Platform.DrawingVertexX
        OFPlatform.VertexY=Platform.DrawingVertexY
        OFPlatform.VertexZ=Platform.DrawingVertexZ
        OFPlatform.EdgeFrom=[1,1,1,2,2,5,3,4,3]
        OFPlatform.EdgeTo=[2,3,5,6,4,6,5,6,4]
        OFPlatform.PenColour=Platform.Colour

        OFPlatform.Mass=Body.Mass
        self.set_CoG(OFPlatform,Body.CoG)
        self.set_inertia(OFPlatform,Body.I)

        ## Adjust volume
        OFPlatform.Volume=0
        return OFPlatform

    def layer_4(self, Platform, name, Body):               
        OFPlatform=self.model.CreateObject(ofx.ot6DBuoy)
        OFPlatform.Name='Layer4_' + name
        ## Set Drawing 
        OFPlatform.VertexX=Platform.DrawingVertexX_L4
        OFPlatform.VertexY=Platform.DrawingVertexY_L4
        OFPlatform.VertexZ=Platform.DrawingVertexZ_L4
        OFPlatform.EdgeFrom=[1,1,1,2,2,5,3,4,3]
        OFPlatform.EdgeTo=[2,3,5,6,4,6,5,6,4]
        OFPlatform.PenColour=Platform.Colour_L4

        OFPlatform.Mass=Body.Mass
        self.set_CoG(OFPlatform,Body.CoG)
        self.set_inertia(OFPlatform,Body.I)

        ## Adjust volume
        OFPlatform.Volume=0
        return OFPlatform
    
    def node(self,Name,Platform):
        OFNode=self.model.CreateObject(ofx.ot6DBuoy)
        OFNode.Name='Node'+Name
        self.set_neglegible_buoy_properties(OFNode)
        OFNode.Mass=Platform.MassNodesTotal/3
        OFNode.VertexX=np.array(OFNode.VertexX)/10
        OFNode.VertexY=np.array(OFNode.VertexY)/10
        OFNode.VertexZ=np.array(OFNode.VertexZ)/10
        return OFNode
    
    def floater(self, Floater, name):
        try: 
            self.model.DestroyObject('Floater'+ name)    
        except: 
            pass
        OFFloater=self.model.CreateObject(ofx.ot6DBuoy)
        OFFloater.Name='Floater'+ name
        OFFloater.PenColour=Floater.Colour
        ## Set Mass Properties
        OFFloater.Mass=Floater.MassOb.Mass
        self.set_CoG(OFFloater,Floater.MassOb.CoG+[0,0,Floater.PosBottom[2]])
        self.set_inertia(OFFloater, Floater.MassOb.I)
        ## Set geometry and spar discretizaion
        OFFloater.BuoyType='Spar buoy'
        OFFloater.NumberOfCylinders=Floater.NumCylinders

        OFFloater.CylinderLength=Floater.CylinderLength
        OFFloater.CylinderOuterDiameter=Floater.CylinderDiameter
        OFFloater.CylinderNormalDragForceCoefficient= Floater.CdNormal

        ## Set drag properties
        OFFloater.NormalDragAreaCalculatedFromGeometry='Yes'
        OFFloater.CylinderAxialDragArea= [Floater.AreaAxial for i in range(Floater.NumCylinders)]
        
        ## Axial properties have only bee applied to the floater end as they are only relevant at that location
        temp=np.zeros(Floater.NumCylinders)
        temp[-1]=Floater.CdAxial
        OFFloater.CylinderAxialDragForceCoefficient= temp
        
        SegementLength=OFFloater.CylinderLength[-1]
        VolumeSegment=Floater.CylinderDiameter[-1]**2*np.pi/4*SegementLength 
        Displacement=Floater.DisplacedMassPlatform/3/1.025
        temp[-1]=Floater.CaAxial*Displacement/(VolumeSegment)
        OFFloater.CylinderAxialAddedMassForceCoefficient= temp
        
        temp[-1]=Floater.WaveDampAxial
        OFFloater.CylinderUnitAxialDampingForce= temp
        
        ## Set added mass and damping 
        OFFloater.CylinderNormalAddedMassForceCoefficient= [Floater.CaNormal for i in range(Floater.NumCylinders)]
        OFFloater.CylinderUnitNormalDampingForce= [Floater.WaveDampNormalSegment for i in range(Floater.NumCylinders)]
                
        return OFFloater
    
    def mass_buoy(self, Name, Mass, Size):
        Buoy=self.model.CreateObject(ofx.ot6DBuoy)
        Buoy.Name=Name
        self.set_neglegible_buoy_properties(Buoy)
        Buoy.Mass=Mass
        if np.size(Size)==3: ### is an array is put in this will be the edge length of the rectange
            Buoy.VertexX=np.array(Buoy.VertexX)/3*Size[0]/2
            Buoy.VertexY=np.array(Buoy.VertexY)/3*Size[1]/2
            Buoy.VertexZ=np.array(Buoy.VertexZ)/3*Size[2]/2
        else:
            Buoy.VertexX=np.array(Buoy.VertexX)/3*Size/2
            Buoy.VertexY=np.array(Buoy.VertexY)/3*Size/2
            Buoy.VertexZ=np.array(Buoy.VertexZ)/3*Size/2
            
        return Buoy
        
    
    def set_neglegible_buoy_properties(self, Buoy):
        Buoy.Volume=0
        Buoy.Mass=0
        Buoy.MassMomentOfInertiaX=0
        Buoy.MassMomentOfInertiaY=0
        Buoy.MassMomentOfInertiaZ=0
        Buoy.FluidAccelerationForceCoefficientX=0
        Buoy.FluidAccelerationForceCoefficientY=0
        Buoy.FluidAccelerationForceCoefficientZ=0
        Buoy.HydrodynamicMassX=0        
        Buoy.HydrodynamicMassY=0
        Buoy.HydrodynamicMassY=0
        
    def drag_plate(self, DragPlate, Name):
        try: 
            self.model.DestroyObject('DragPlate'+Name)    
        except: 
            pass
        OFDragPlate=self.model.CreateObject(ofx.ot6DBuoy)
        OFDragPlate.Name='DragPlate'+Name
        self.set_neglegible_buoy_properties(OFDragPlate)
        i=np.arange(0,2*np.pi+2*np.pi/12,(2*np.pi/12))
        OFDragPlate.VertexX=DragPlate.Radius*np.cos(i)
        OFDragPlate.VertexY=DragPlate.Radius*np.sin(i)
        OFDragPlate.VertexZ=np.zeros(i.shape)
        OFDragPlate.EdgeFrom=np.arange(1,len(i)+1)
        OFDragPlate.EdgeTo=np.concatenate((np.arange(2,len(i)),np.array([1])))
        OFDragPlate.Height=0.1
        OFDragPlate.FluidAccelerationForceCoefficientZ=DragPlate.Ca
        OFDragPlate.AddedMassCoefficientZ=DragPlate.Ca
        OFDragPlate.HydrodynamicMassZ=DragPlate.HydrodynamicMass
        OFDragPlate.DragAreaZ=DragPlate.DragArea
        OFDragPlate.DragForceCoefficientZ=DragPlate.Cd
        return OFDragPlate
           
    def coupling(self, Coupling, Name, Type):
        if Name.find('Mid')!=-1 and Coupling.DetailedModel :
            try:
                self.model.DestroyObject(Name) 
            except:     
                Link=self.line(['CouplingBar'],[Coupling.Len['Mid']],[2],Name) #### Watch out the segement lengths of the coupling is hardcoded here 
        else:
            Link=self.model.CreateObject(ofx.otLink)
            Link.name=Name
            Link.LinkType='Spring/damper'
            Link.LinearSpring='Yes'
            Link.LinearDamper='Yes'
            Link.UnstretchedLength=Coupling.Len[Type]
            Link.Stiffness=Coupling.Stiffness[Type]
            Link.DampingConstant=Coupling.Damping[Type]
            if Coupling.LinearSpring[Type]=='No':
                Link.LinearSpring='No'
                ## Added to make unstreched length correct hardcoded that excel unstretched length is 1.95
                ExtraLength=Coupling.Len[Type]-Coupling.NLLengthZeroIntercept
                Link.NumberOfSpringTableEntries=len(Coupling.NLForce)
                Link.SpringLength=Coupling.NLDisplacement + ExtraLength
                Link.SpringTension=Coupling.NLForce  
            if Coupling.LinearDamping[Type]=='No':
                Link.LinearDamper='No'
                Link.NumberOfDamperTableEntries=len(Coupling.NLDampingForce)
                Link.DamperTension=Coupling.NLDampingForce
                Link.DamperVelocity=Coupling.NLDampingVelocity  
        return Link
    
    def edge_coupling(self, Coupling, Name, Type):
        if Name.find('Mid')!=-1 and Coupling.DetailedModel :
            try:
                self.model.DestroyObject(Name) 
            except:     
                Link=self.line(['CouplingBar'],[Coupling.Len['Mid']],[2],Name) #### Watch out the segement lengths of the coupling is hardcoded here 
        else:
            # Link=self.model.CreateObject(ofx.otLink)
            Link = self.model[Name]
            # Link.name=Name
            Link.LinkType='Spring/damper'
            Link.LinearSpring='Yes'
            Link.LinearDamper='Yes'
            Link.UnstretchedLength=Coupling.Len[Type]
            Link.Stiffness=Coupling.Stiffness[Type]
            Link.DampingConstant=Coupling.Damping[Type]
            if Coupling.LinearSpring[Type]=='No':
                Link.LinearSpring='No'
                ## Added to make unstreched length correct hardcoded that excel unstretched length is 1.95
                ExtraLength=Coupling.Len[Type]-Coupling.NLLengthZeroIntercept
                Link.NumberOfSpringTableEntries=len(Coupling.NLForce)
                Link.SpringLength=Coupling.NLDisplacement + ExtraLength
                Link.SpringTension=Coupling.NLForce  
            if Coupling.LinearDamping[Type]=='No':
                Link.LinearDamper='No'
                Link.NumberOfDamperTableEntries=len(Coupling.NLDampingForce)
                Link.DamperTension=Coupling.NLDampingForce
                Link.DamperVelocity=Coupling.NLDampingVelocity  
        return Link
    
    def constraint(self, Name):
        try: 
            self.model.DestroyObject(Name)    
        except: 
            pass
        Constraint=self.model.CreateObject(ofx.otConstraint)
        Constraint.Name= Name
        return Constraint       
       
                      
    def line(self, LineTypes, LineLengths, SegementLengths, Name):
        Line=self.model.CreateObject(ofx.otLine)
        Line.Name=Name
        Line.NumberOfSections=len(LineTypes)
        Line.LineType=LineTypes
        Line.Length=LineLengths
        Line.TargetSegmentLength=SegementLengths      
        return Line
    
    def mooring_tether(self, LineTypes, LineLengths, SegementLengths, Name):
        Link=self.tether(Name)
        # Set stiffness to the equivalent
        Link.UnstretchedLength=np.sum(LineLengths)
        k_inv=[]
        for i,LineType in enumerate(LineTypes):
            LineType=self.model[LineType]
            k_inv.append(1/(LineType.EA/LineLengths[i]))
        k=1/(np.sum(k_inv))
        Link.Stiffness=k*np.sum(LineLengths)
        return Link
    
    def tether(self, Name):
        Link=self.model.CreateObject(ofx.otLink)
        Link.name= Name
        Link.LinkType='Tether'
        return Link
    
    def set_buoy_properties(self, OFBuoy, Buoy):
        #### 3D buoy
        OFBuoy.Mass=Buoy.Mass
        OFBuoy.Volume=Buoy.Volume
        OFBuoy.Height=Buoy.Height
        OFBuoy.DragAreaX=Buoy.DragAreaSurge
        OFBuoy.DragAreaY=Buoy.DragAreaSurge
        OFBuoy.DragAreaZ=Buoy.DragAreaHeave
        OFBuoy.CdX=Buoy.CdSurge
        OFBuoy.CdY=Buoy.CdSurge
        OFBuoy.CdZ=Buoy.CdHeave
        OFBuoy.CaX=Buoy.CaSurge
        OFBuoy.CaY=Buoy.CaSurge
        OFBuoy.CaZ=Buoy.CaHeave
        return OFBuoy
    
    def mooring_3dbuoy(self, Name, MooringBuoy): ## This can be for a mooring bouy or for a shakel, shakel properties are hardcoded here
        OFBuoy=self.model.CreateObject(ofx.ot3DBuoy)
        OFBuoy.Name=Name
        if MooringBuoy==0:
            #### This was an override for the VLink mooring connection component
            OFBuoy.Mass=0.3
            OFBuoy.Volume=0
            OFBuoy.Height=0.1
        else:
            self.set_buoy_properties(OFBuoy, MooringBuoy)
        return OFBuoy
    
    def vlink_buoy(self,Name, VLinkBuoy=0):
        # return self.mooring_3dbuoy(Name, MooringBouy)
        OFBuoy=self.model.CreateObject(ofx.ot3DBuoy)
        OFBuoy.Name=Name
        if VLinkBuoy==0:
            #### This was an override for the VLink mooring connection component
            OFBuoy.Mass=0.001
            OFBuoy.Volume=0
            OFBuoy.Height=0.1
        else:
            self.set_buoy_properties(OFBuoy, VLinkBuoy)
        return OFBuoy
    
    def clump_weight(self, Name, MooringBuoy): 
        try:
            OFClump=self.model[Name]
        except: 
            OFClump=self.model.CreateObject(ofx.otClumpType)
        OFClump.Name=Name
        self.set_buoy_properties(OFClump, MooringBuoy)
        OFClump.Offset=MooringBuoy.Offset
        return OFClump
     
    def pendulum_buoy(self, Name, Pendulum, Type): ## This can be for a mooring bouy or for a shakel, shakel properties are hardcoded here
        Buoy=self.model.CreateObject(ofx.ot3DBuoy)
        Buoy.Name=Name
        Buoy.Mass=Pendulum.Mass[Type]
        Buoy.Volume=Pendulum.Volume[Type]
        Buoy.Height=Pendulum.Height[Type]
        Buoy.DragAreaX=Pendulum.DragArea[Type]
        Buoy.DragAreaY=Pendulum.DragArea[Type]
        Buoy.DragAreaZ=Pendulum.DragArea[Type]
        Buoy.CdX=Pendulum.Cd
        Buoy.CdY=Pendulum.Cd
        Buoy.CdZ=Pendulum.Cd
        Buoy.CaX=Pendulum.Cd
        Buoy.CaY=Pendulum.Cd
        Buoy.CaZ=Pendulum.Cd
        return Buoy
    
    def wingtype(self, Wing, Name):
        OFWing=self.model.CreateObject(ofx.otWingType)
        OFWing.Name=Name
        OFWing.NumberOfAngles=len(Wing.Angles)
        OFWing.Angle=Wing.Angles 
        OFWing.Lift=Wing.Cl[Name]
        OFWing.Drag=Wing.Cd[Name]
        OFWing.Moment=Wing.Cm[Name]
        OFWing.PenColour=Wing.Colour[Name]
    
    def environment_basics(self, OFEnv, Environment):
        OFEnv.WaterDepth=Environment.Depth
        OFEnv.CurrentMethod='Power law'
        OFEnv.CurrentRamped='Yes'
        OFEnv.CurrentExponent=Environment.CurrentPowerLawExponent
        OFEnv.CurrentSpeedAtSurface=Environment.CurrentSpeed
        OFEnv.RefCurrentDirection=Environment.CurrentDirection 
        OFEnv.WindSpeed=Environment.WindSpeed
        OFEnv.WindDirection=Environment.WindDirection
    
    def environment_jonswap(self, Environment):
        OFEnv=self.model.environment
        OFEnv.WaveType='JONSWAP'
        OFEnv.WaveNumberOfSpectralDirections=1
        OFEnv.WaveGamma=Environment.Gamma
        OFEnv.WaveHs=Environment.Hs
        OFEnv.WaveTp=Environment.Tp
        OFEnv.WaveDirection=Environment.WaveHeading
        OFEnv.WaveNumberOfComponents=Environment.NumberOfWaveComponentPerDirection
        if Environment.NumberOfWaveComponentPerDirection<200 and Environment.SpreadingCoefficient == None:
            print('Env number of wave compenents likely too low')
        self.environment_basics(OFEnv, Environment)
        OFEnv.KinematicStretchingMethod=Environment.WaveStretching
        
        if Environment.SpreadingCoefficient != None :
            OFEnv.WaveNumberOfSpectralDirections=Environment.NumberOfWaveDirections
            OFEnv.WaveDirectionSpreadingExponent=Environment.SpreadingCoefficient
            OFEnv.WaveNumberOfComponentsPerDirection=Environment.NumberOfWaveComponentPerDirection
            if Environment.NumberOfWaveComponentPerDirection>200:
                print('Likely to many wave components per direction')
            
        if Environment.Tp > 5*np.sqrt(Environment.Hs) or Environment.Tp < 3.6*np.sqrt(Environment.Hs):
            print('**************************************************')
            print('Warning: Hm0,Tp is outside the JONSWAP range')
            print('The validity of the spectral density is questionable')
            print('Hs: ' + str(Environment.Hs) +' Tp:'+ str(Environment.Tp) + 'Tp/(sqrt(Hs))=' + str(Environment.Tp/np.sqrt(Environment.Hs)))
    
        if Environment.Gamma > 7 or Environment.Gamma < 1:
            print('**************************************************')
            print('Warning: gamma_val is outside the valid range')
            print('The validity of the spectral density is questionable')
            print('Hs: ' + str(Environment.Hs) +' Tp:'+ str(Environment.Tp)+ 'Tp/(sqrt(Hs))=' + str(Environment.Tp/np.sqrt(Environment.Hs)))
            
    def enviroment_non_linear_hmax(self, Environment, Order=11):
        ## Type must be "Stokes' 5th" or "Cnoidal"
        OFEnv=self.model.environment
        OFEnv.WaveType=Environment.WaveType
        OFEnv.WaveDirection=Environment.WaveHeading
        OFEnv.WaveHeight=Environment.Hmax
        OFEnv.WavePeriod=Environment.Tass  
        if Environment.WaveType=='Dean stream':
            OFEnv.WaveStreamFunctionOrder=Order
        self.environment_basics(OFEnv, Environment)
           
    def general(self, General, Environment):
        if Environment.HsOrHmax=='Hs':
            self.model.general.StageDuration=[2*Environment.Tp,General.DurationSimulation]
        else:
            self.model.general.StageDuration=[np.max([10,2*Environment.Tass]),General.DurationSimulation]
        self.model.general.StaticsMaxDamping=General.StaticsMaxDamping
        self.model.general.StaticsMinDamping=General.StaticsMinDamping
        self.model.general.StaticsMaxIterations=General.StaticsMaxIterations
        # self.model.general.ImplicitVariableMaxTimeStep=General.ImplicitConstantTimeStep
        self.model.general.ImplicitConstantTimeStep=General.ImplicitConstantTimeStep
        self.model.general.TargetLogSampleInterval=General.TargetLogSampleInterval
        
    
## Clase with functions that position and connect everything        
class Compose:
    
    def __init__(self, model=None):
        self.model=model
        
    def whole_model(self, Gen, General, Environment, Platform, PlatformEdge, Layout, Floater, DragPlate, Mooring, Coupling, MooringBuoy=0):
        Gen.general(General, Environment)
        # Gen.environment_jonswap(Environment)
        self.position_platforms(Platform, Layout, Floater, Gen, Mooring, PlatformEdge)      
        self.connect_floaters(Platform, Layout, Floater, DragPlate, Gen) 
        if Coupling.LoadCells:
    	    if Coupling.DetailedModel:
    	    	self.connect_constraint_links_with6Dbuoy(Platform, Coupling, Gen)
    	    else:
    	    	self.connect_constraint_links(Platform, Coupling, Gen)
        else:
            self.connect_links(Platform, Coupling, Gen)
        self.connect_mooring_automatic(Layout, Platform, Mooring, Gen, MooringBuoy, Coupling)
       
    def rotate_coordinates_around_z(self, Pos, Phi,):
        #Input phi is in degrees
        Phi=Phi*np.pi/180
        Pos=np.matrix(Pos)
        RotMat=np.matrix([[np.cos(Phi), -np.sin(Phi), 0], [np.sin(Phi), np.cos(Phi),0],[0,0,1]])
        NewPos=np.squeeze(np.asarray(RotMat*Pos.transpose()))
        NewPos=NewPos.transpose()
        return NewPos      
        
    def translate_coordinates(self, Platform, Pos, Phi):
        #Translates from input axis system to plaform axis system with certain rotation       
        PosFromCoG=np.array(Pos)-np.array([Platform.GeoCoG[0],0,0])
        PosNewFromCoG=self.rotate_coordinates_around_z(PosFromCoG, Phi)
        PosNew=PosNewFromCoG+np.array([Platform.GeoCoG[0],0,0])
        return PosNew
    
    def local_to_global(self, Layout, N, Pos):
        #N is plaformform number
        PosGlobalOrientation=self.rotate_coordinates_around_z(Pos, Layout.Orientation[N])
        PosGlobal=[Layout.PlatformsX[N], Layout.PlatformsY[N], 0] + PosGlobalOrientation
        return PosGlobal
    
    def set_object_position(self, Object, XYZ):
        Object.InitialX=XYZ[0]
        Object.InitialY=XYZ[1]
        Object.InitialZ=XYZ[2]
        
    def set_object_EndA_EndB(self, Object, XYZA, XYZB):
        self.set_object_EndA(Object,XYZA)
        self.set_object_EndB(Object,XYZB)
        
    def set_object_EndA(self,Object, XYZA):
        Object.EndAX=XYZA[0]
        Object.EndAY=XYZA[1]
        Object.EndAZ=XYZA[2]
        
    def set_object_EndB(self,Object, XYZB):
        Object.EndBX=XYZB[0]
        Object.EndBY=XYZB[1]
        Object.EndBZ=XYZB[2]
    
    def set_platform_position(self, OFPlatform, Layout, Floater, i):
        OFPlatform.InitialX=Layout.PlatformsX[i]
        OFPlatform.InitialY=Layout.PlatformsY[i]
        OFPlatform.InitialRotation3=Layout.Orientation[i]
        OFPlatform.InitialZ=-Floater.PosBottom[2]-Floater.Draft
    
    def position_nodes(self, Platform, i,j, Gen, Layout):
        if Layout.Types[i]==1:
            Platform.NodesXYZ=Platform.NodesXYZ_T1
        elif Layout.Types[i]==2:
            Platform.NodesXYZ=Platform.NodesXYZ_T2
        for inn, n in enumerate(Platform.NodesXYZ):
            OFNode=Gen.node('_'+str(i)+'_'+str(j)+'_'+str(inn),Platform)
            OFNode.Connection='Platform_'+str(i)+'_0'
            XYZ=self.translate_coordinates(Platform, n, j)
            OFNode.InitialRotation3=0
            self.set_object_position(OFNode, XYZ)
            OFNode.Connection='Free'
    
    def connect_beams(self, Platform, i, j, Gen, Layout):
        ### set the correct masses
        if Platform.PlaceSkid and i == Platform.SkidPlatform:
            for key in Platform.LineMassPerLengthSkid.keys():
                self.model[key].MassPerUnitLength=Platform.LineMassPerLengthSkid[key]
        else:
            for key in Platform.LineMassPerLength.keys():
                self.model[key].MassPerUnitLength=Platform.LineMassPerLength[key]
                if Platform.BoxTruss:
                    self.model[key].EIx = Platform.BoxTrussEIx
                    self.model[key].EIy = Platform.BoxTrussEIy
                    self.model[key].EA = Platform.BoxTrussEA
                    self.model[key].GJ = Platform.BoxTrussGJ
                    self.model[key].NormalDragLiftDiameter = Platform.TrussNormalDragLiftDiameter
                    self.model[key].AxialDragLiftDiameter = Platform.TrussAxialDragLiftDiameter
                self.model[key].CGx=-Platform.BeamCoGZOffset[key]
            
                
        if Layout.Types[i]==1:
            Platform.BeamConLocA=Platform.BeamConLocA_T1
            Platform.BeamConLocB=Platform.BeamConLocB_T1
            Platform.BeamLengths=Platform.BeamLengths_T1
        else:
            Platform.BeamConLocA=Platform.BeamConLocA_T2
            Platform.BeamConLocB=Platform.BeamConLocB_T2
            Platform.BeamLengths=Platform.BeamLengths_T2
            
        for j in [0,120,240]:
            for ib, BeamType in enumerate(Platform.BeamTypes):
                # print(i)
                if Platform.WindShielding:
                    if i in Platform.EdgePlatforms:
                        if Platform.AdjustMassOfEdgePlatforms:
                            BeamName = BeamType+str(Layout.Types[i])+'Edge'
                        else:
                            BeamName = BeamType+str(Layout.Types[i])
                    else:
                        BeamName = BeamType+str(Layout.Types[i])+'WindShielding'+str(Platform.CoeReduceRatio)
                else:
                    if i in Platform.EdgePlatforms:
                        if Platform.AdjustMassOfEdgePlatforms:
                            BeamName = BeamType+str(Layout.Types[i])+'Edge'
                        else:
                            BeamName = BeamType+str(Layout.Types[i])
                    else:
                        if i == Platform.SkidPlatform:
                            BeamName = BeamType+str(Layout.Types[i]) + 'Skid'
                        else:
                            BeamName = BeamType+str(Layout.Types[i])

                ###
                # print(BeamName)
                OFBeam=Gen.line([BeamName],[Platform.BeamLengths[ib]],[Platform.BeamLengths[ib]/Platform.BeamSegments],'Beam_' +str(i)+'_'+str(j)+'_'+str(ib))
                
                if Platform.BoxTruss:
                    OFBeam.IncludeTorsion = 'Yes'
                    OFBeam.EndATwistingStiffness = 1e+307
                    OFBeam.EndBTwistingStiffness = 1e+307
                OFBeam.EndADeclination=90
                OFBeam.EndBDeclination=90
                OFBeam.EndAGamma=0
                OFBeam.EndBGamma=0
                OFBeam.EndBxBendingStiffness=1e+307
                OFBeam.EndAyBendingStiffness=1e+307
                OFBeam.EndAxBendingStiffness=1e+307
                OFBeam.EndByBendingStiffness=1e+307
                EndAName=Platform.BeamConnectionsA[ib].split('_')
                EndAName[1]=str(i)
                EndAName[2]=str(round(np.mod(float(EndAName[2])+j,360)))
                OFBeam.EndAConnection='_'.join(EndAName)
                EndBName=Platform.BeamConnectionsB[ib].split('_')
                EndBName[1]=str(i)
                EndBName[2]=str(round(np.mod(float(EndBName[2])+j,360)))
                OFBeam.EndBConnection='_'.join(EndBName)
                if EndAName[0]=='Platform':
                    XYZA=self.translate_coordinates(Platform,Platform.BeamConLocA[ib], j)
                else:
                    XYZA=[0,0,0]
                if EndBName[0]=='Platform':
                    XYZB=self.translate_coordinates(Platform,Platform.BeamConLocB[ib], j)
                else:
                    XYZB=[0,0,0]
                OFBeam.EndAAzimuth=np.mod(Platform.BeamAzimuth[ib]+j,360)
                OFBeam.EndBAzimuth=OFBeam.EndAAzimuth
                self.set_object_EndA_EndB(OFBeam, XYZA, XYZB)
    
    def position_platforms(self, PlatformInner, Layout, Floater, Gen, Mooring, PlatformOuter):
        
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0: 
                if i in PlatformInner.EdgePlatforms: 
                    Platform=PlatformOuter
                    # print(i)
                else:
                    Platform=PlatformInner
                if Platform.Split:
                    
                    for j in [0,120,240]:
                       if Layout.Types[i]==1:
                           DrawingPos=np.array([Platform.DrawingVertexX0_T1, Platform.DrawingVertexY0_T1, Platform.DrawingVertexZ0])
                       elif Layout.Types[i]==2:
                           DrawingPos=np.array([Platform.DrawingVertexX0_T2, Platform.DrawingVertexY0_T2, Platform.DrawingVertexZ0])
                       DrawingPos=np.transpose(DrawingPos)
                       NewDrawingPos=self.translate_coordinates(Platform, DrawingPos, j)
                       Platform.DrawingVertexX=NewDrawingPos[:,0]
                       Platform.DrawingVertexY=NewDrawingPos[:,1]
                       Platform.DrawingVertexZ=NewDrawingPos[:,2]
                       
                       #### Section that adjusts masses of platforms at the edges
                       #### This section automatically regonizes an edge platform based on the mooring selection
                       Position=Platform.CouplingCoG
                       Ex_T1=Input.MassObject(0,[0,0,0])
                       Ex_T2=Input.MassObject(0,[0,0,0])
                       MFrld=Input.MassObject(0,[0,0,0])
                       
                       #### to adjust edge platform mass before adjust outer platforms mass with mooring
                       # if Platform.AdjustMassOfEdgePlatforms:
                       #     Platform.MassCouplingsTotal_T1 = Platform.MassCouplingsTotal_T1Edge
                       #     Platform.MassCouplingsTotal_T2 = Platform.MassCouplingsTotal_T2Edge
                       #     Platform.T1CornerOut = Platform.T1CornerOutEdge
                       #     Platform.T2CornerOut = Platform.T2CornerOutEdge
                           
                       if Platform.AdjustMassOfOuterPlatforms: 
                           if any((np.array(Mooring.ConnectedPlatformsTotal)== i) & (np.array(Mooring.OrientationLocalTotal)==j)):
                               Ex_T1=Input.MassObject(-Platform.MassCouplingsTotal_T1/6,self.translate_coordinates(Platform, Position, j))
                               Ex_T2=Input.MassObject(-Platform.MassCouplingsTotal_T2/6,self.translate_coordinates(Platform, Position, j))
                               MFrld=Input.MassObject(Platform.MassFairlead, self.translate_coordinates(Platform, Mooring.Position, j))
                               
                           if any((np.array(Mooring.ConnectedPlatformsTotal)== i) & (np.array(Mooring.OrientationLocalTotal)==np.mod(j-120,360))):
                               Ex_T1=Input.MassObject(-Platform.MassCouplingsTotal_T1/6, self.translate_coordinates(Platform, Position*np.array([1,-1,1]), j+240))
                               Ex_T2=Input.MassObject(-Platform.MassCouplingsTotal_T2/6, self.translate_coordinates(Platform, Position*np.array([1,-1,1]), j+240))                                               
                               MFrld=Input.MassObject(Platform.MassFairlead, self.translate_coordinates(Platform, Mooring.Position*np.array([1,-1,1]), j+240)) 
                           
                       if Layout.Types[i]==1:
                           T1Corner=Platform.add_mass_objects([Platform.T1CornerOut[str(j)], Ex_T1, MFrld])
                           OFPlatform=Gen.platform(Platform, str(i)+'_'+str(j),T1Corner)
                       elif Layout.Types[i]==2:
                           T2Corner=Platform.add_mass_objects([Platform.T2CornerOut[str(j)], Ex_T2, MFrld])
                           OFPlatform=Gen.platform(Platform, str(i)+'_'+str(j),T2Corner)
                           

                       self.set_platform_position(OFPlatform, Layout, Floater, i)
                       if not Platform.BoxTruss:
                           self.position_nodes(Platform,i,j,Gen, Layout)
                    self.connect_beams(Platform,i,j,Gen, Layout)
                else:
                    if Layout.Types[i]==1:
                        OFPlatform=Gen.platform(Platform, str(i), Platform.MassObj_T1_AllMinL4Flt)
                    elif Layout.Types[i]==2:
                       OFPlatform=Gen.platform(Platform, str(i), Platform.MassObj_T2_AllMinL4Flt)
                    self.set_platform_position(OFPlatform, Layout, Floater, i)
                    OFL_4=Gen.layer_4(Platform, str(i), Platform.MassObj_Layer4)
                    OFL_4.connection=OFPlatform.Name
                    self.set_object_position(OFL_4, [0,0,0])
                    OFL_4.InitialRotation3=0
                    
    def connect_reference_constraints(self, Layout, Platform, ReferenceConstraint, Gen):
        #### This does not yet work for split platforms
        if Platform.Split:
            print('This function has not been added for split platforms')
        else:    
            for i in range(len(Layout.Types)):
                if Layout.Types[i] != 0: 
                    Cons=Gen.constraint('Con_Ref_'+str(i))
                    Cons.DOFFree=['Yes' for i in range(6)]
                    Cons.Connection='Platform_'+str(i)
                    self.set_object_position(Cons, ReferenceConstraint.Position)
                    Cons.Connection='Fixed'
                    Cons.InitialGamma=ReferenceConstraint.Orientation
                    self.model['Platform_'+str(i)].connection='Con_Ref_'+str(i)                  
                                                                                    
    def connect_drag_plate(self, DragPlate, Gen, PlatNum, Dir, ExtraFloaterName=''):
        if DragPlate.Type !=0:
            OFDragPlate=Gen.drag_plate(DragPlate, str(PlatNum)+'_' +str(Dir) + ExtraFloaterName)
            OFDragPlate.Connection='Floater'+str(PlatNum)+'_'+str(Dir) + ExtraFloaterName
            self.set_object_position(OFDragPlate, [0,0,DragPlate.Z])
            
            
    def connect_floaters(self, Platform, Layout, Floater, DragPlate, Gen, ExtraFloaterName='', Type=[1,2]):
        for i in range(len(Layout.Types)):
            if Layout.Types[i] in Type:
                for j in [0, 120, 240]:
                    OFFloater=Gen.floater(Floater, str(i)+'_'+str(j) + ExtraFloaterName)
                    OFFloaterPos=self.translate_coordinates(Platform, Floater.PosBottom, j)
                    if Platform.Split:
                        OFFloater.Connection='Platform_'+str(i)+'_'+str(j)
                    else:
                        OFFloater.Connection='Platform_'+str(i)
                    OFFloater.InitialX=OFFloaterPos[0]
                    OFFloater.InitialY=OFFloaterPos[1]
                    OFFloater.InitialZ=Floater.ZConnectionPoint
                    OFFloater.StackBaseCentreZ=OFFloaterPos[2]
                    OFFloater.InitialRotation3=0
                    if Platform.WindShielding:
                        if i in Platform.InnerPlatforms:
                            OFFloater.CylinderNormalDragForceCoefficient= Floater.CdNormalIn
                        else:
                            OFFloater.CylinderNormalDragForceCoefficient= Floater.CdNormal
                    else:
                        OFFloater.CylinderNormalDragForceCoefficient= Floater.CdNormal
                    self.connect_drag_plate(DragPlate, Gen, i, j, ExtraFloaterName=ExtraFloaterName)   

    def connect_outer_floaters(self, Platform, Floater2, DragPlate, Mooring, Gen, Type=[1,2], DeclinationAngle=False):
        for i, plat in enumerate(Mooring.ConnectedPlatformsTotal):
            if Mooring.ConnectedPlatformType[i] in Type:
                ## These directions can be adjusted to add floaters on all sides
                for j in [0,120]:
                    Orient=np.mod(Mooring.OrientationLocalTotal[i]+j,360)
                    OFFloater=Gen.floater(Floater2, str(plat)+'_'+str(Orient))
                    OFFloaterPos=self.translate_coordinates(Platform, Floater2.PosBottom, Orient)
                    if Platform.Split:
                        OFFloater.Connection='Platform_'+str(plat)+'_'+str(Orient)
                    else:
                        OFFloater.Connection='Platform_'+str(plat)
                    OFFloater.InitialX=OFFloaterPos[0]
                    OFFloater.InitialY=OFFloaterPos[1]
                    OFFloater.InitialZ=0 ## This values determines where floater loads are read out
                    OFFloater.StackBaseCentreZ=OFFloaterPos[2]
                    OFFloater.InitialRotation3=0           
                    if DeclinationAngle:
                        #### Add some extra length to the floataer to compensate the rotation, such that only an extra part is added
                        OFFloater.StackBaseCenterZ=OFFloater.StackBaseCenterZ-np.max(Floater2.CylinderDiameter)/2*np.tan(DeclinationAngle*np.pi/180)
                        #### Connect to the fairlead for and easy floater rotation
                        OFFloater.Connection='Con_Frld_'+str(plat)+'_'+str(Mooring.OrientationLocalTotal[i])+'_'+str(Mooring.LocTotal[i])
                        OFFloater.InitialRotation2=DeclinationAngle
                        OFFloater.Connection='Platform_'+str(plat)+'_'+str(Orient)
                    self.connect_drag_plate(DragPlate, Gen, plat, Orient)
                      
    def connect_edgeplatform_floaters(self, Platform, Floater, DragPlate, Mooring, Gen, Type=[1,2], DeclinationAngle=False):
        
        ### this to change the floater properties at the edge platforms, 3 floater on each platform included
        for ip, plat in enumerate(Platform.EdgePlatforms):
            for j in [0,120,240]:
                OFFloater = Gen.floater(Floater, str(plat)+'_'+str(j))
                OFFloaterPos = self.translate_coordinates(Platform, Floater.PosBottom, j)
                if Platform.Split:
                    OFFloater.Connection='Platform_'+str(plat)+'_'+str(j)
                else:
                    OFFloater.Connection='Platform_'+str(plat)
                OFFloater.InitialX=OFFloaterPos[0]
                OFFloater.InitialY=OFFloaterPos[1]
                OFFloater.InitialZ=0 ## This values determines where floater loads are read out
                OFFloater.StackBaseCentreZ=OFFloaterPos[2]
                OFFloater.InitialRotation3=0     
                self.connect_drag_plate(DragPlate, Gen, plat, j)
                    
    def connect_floater_chain(self, Platform, Layout,  Floater, Gen, PreMoment):
        #### function to connect chains, the chain is selected in the line-types
        for i in range(len(Layout.Types)):
            if Layout.Types[i] in [1,2]:
                for j in [0, 120, 240]:
                    OFFloaterPos1=self.translate_coordinates(Platform, Floater.PosBottom, j)
                    OFFloaterPos2=self.translate_coordinates(Platform, Floater.PosBottom, np.mod(j+120,360))
                    Length=np.linalg.norm(OFFloaterPos1-OFFloaterPos2)-Floater.CylinderDiameter[-2]
                    Arm=Platform.BeamHeight - Floater.PosBottom[2] - Floater.CylinderLength[-1] #### Adjust for the platforms connection locations to the truss 
                    #### Unstreched length should be determined by the sitffness of the structure and the length of the arm
                    EI=self.model['TrussOutT1'].EIx
                    Phi=PreMoment*Platform.BeamLengths_T1[0]/EI
                    UnstretchtedLength=Length-Arm*np.sin(Phi)*2
                    OFChain=Gen.line(['FloaterChain40mm'],[UnstretchtedLength],[UnstretchtedLength/3],'Chain_'+str(i)+'_'+str(j))
                    OFChain.EndAConnection='Floater'+str(i)+'_'+str(j)
                    OFChain.EndBConnection='Floater'+str(i)+'_'+str(np.mod(j+120,360))
                    Normal=self.rotate_coordinates_around_z([0,Floater.CylinderDiameter[-2]/2,Floater.PosBottom[2]+Floater.CylinderLength[-1]], j)
                    self.set_object_EndB(OFChain,Normal)
                    Normal1=np.zeros(3)
                    Normal1[:2]=Normal[:2]*-1
                    Normal1[2]=Normal[2]
                    self.set_object_EndA(OFChain,Normal1)

    def connect_links(self, Platform , Coupling, Gen):
         for i in range(len(Coupling.Names)):
             Link=Gen.coupling(Coupling, Coupling.Names[i], Coupling.Types[i])
             Link.EndAConnection=Coupling.EndAConnections[i]
             Link.EndBConnection=Coupling.EndBConnections[i]
             self.set_object_EndA_EndB(Link, Coupling.PosEndAs[i], Coupling.PosEndBs[i])
           
    def connect_constraint_links(self, Platform , Coupling, Gen):
         for i in range(len(Coupling.Names)):
             Link=Gen.coupling(Coupling, Coupling.Names[i], Coupling.Types[i])
             Constraint=Gen.constraint('Con_' + Coupling.Names[i])
             
             Constraint.Connection=Coupling.EndAConnections[i]
             Constraint.InitialX=Coupling.PosEndAs[i][0]
             Constraint.InitialY=Coupling.PosEndAs[i][1]
             Constraint.InitialZ=Coupling.PosEndAs[i][2]
             Constraint.InitialAzimuth=Coupling.Directions[i]
             Constraint.InitialDeclination=0
             Constraint.InitialGamma=0

             Link.EndAConnection='Con_'+ Coupling.Names[i]
             ConCenterName='Con_Cen_'+ Coupling.EndBConnections[i][9:].split('_')[0]+ '_'+ str(Coupling.Directions[i])+'_'+ str(Coupling.Locs[i]*-1)
             try: 
                 Link.EndBConnection=ConCenterName
             except:
                 ConstraintCenter=Gen.constraint(ConCenterName)
                 ConstraintCenter.Connection=Coupling.EndBConnections[i]
                 ConstraintCenter.InitialX=Coupling.PosEndBs[i][0]
                 ConstraintCenter.InitialY=Coupling.PosEndBs[i][1]
                 ConstraintCenter.InitialZ=Coupling.PosEndBs[i][2]
                 ConstraintCenter.InitialAzimuth=Coupling.Directions[i]
                 ConstraintCenter.InitialDeclination=0
                 ConstraintCenter.InitialGamma=0
                 Link.EndBConnection=ConCenterName  
             self.set_object_EndA_EndB(Link, [0,0,0],[0,0,0])
             
    
    def connect_constraint_links_with6Dbuoy(self, Platform , Coupling, Gen):
        for i in range(len(Coupling.Names)):    
            Link=Gen.coupling(Coupling, Coupling.Names[i], Coupling.Types[i])
            Constraint=Gen.constraint('Con_' + Coupling.Names[i])
            
            ConCenterName='Con_Cen_'+ Coupling.EndBConnections[i][9:].split('_')[0]+ '_'+ str(Coupling.Directions[i])+'_'+ str(Coupling.Locs[i]*-1)
            ConLinkName = 'Con_Link_'+ Coupling.EndBConnections[i][9:].split('_')[0]+ '_'+ str(Coupling.Directions[i])+'_'+ str(Coupling.Locs[i]*-1)
            # BuoyCenterName='Buoy_Cen_'+ Coupling.EndBConnections[i][9:].split('_')[0]+ '_'+ str(Coupling.Directions[i])+'_'+ str(Coupling.Locs[i]*-1)
            ConBarName='Con_Bar_'+ Coupling.EndBConnections[i][9:].split('_')[0]+ '_'+ str(Coupling.Directions[i])+'_'+ str(Coupling.Locs[i]*-1)
            Constraint.Connection=Coupling.EndAConnections[i]
            Link.EndAConnection='Con_'+ Coupling.Names[i]
            # Link.EndAAzimuth=180
            # Link.EndADeclination=50
            Constraint.InitialX=Coupling.PosEndAs[i][0]
            Constraint.InitialY=Coupling.PosEndAs[i][1]
            Constraint.InitialZ=Coupling.PosEndAs[i][2]
            Constraint.InitialAzimuth=Coupling.Directions[i]
            Constraint.InitialDeclination=0
            Constraint.InitialGamma=0
            try: 
                # Link.EndBConnection=BuoyCenterName
                if 'Mid' in Coupling.Names[i]:
                    Link.EndBConnection=ConBarName
                else:
                    Link.EndBConnection=ConLinkName
                    # Link.EndBConnection=BuoyCenterName
            except:
                ConBar = Gen.constraint(ConBarName)
                ConLink = Gen.constraint(ConLinkName)
                ConstraintCenter=Gen.constraint(ConCenterName)
                # BuoyCenter = Gen.of_6dbuoy(BuoyCenterName,0.001,0.5)
                ConBar.Connection = ConLinkName
                ConLink.Connection = ConCenterName
                # BuoyCenter.Connection = ConCenterName
                ConBar.InitialX=-0.2218
                ConBar.InitialY=0
                ConBar.InitialZ=-0.16817
                ConBar.InitialDeclination=np.arctan(ConBar.InitialX/ConBar.InitialZ)*180/np.pi
                ConBar.DOFFree[3] = 'Yes'
                ConLink.DOFFree[4] = 'Yes'
                # ConBouy.DOFFree[5] = 'Yes'
                
                ConstraintCenter.Connection=Coupling.EndBConnections[i]
                ConstraintCenter.InitialX=Coupling.PosEndBs[i][0]
                ConstraintCenter.InitialY=Coupling.PosEndBs[i][1]
                ConstraintCenter.InitialZ=Coupling.PosEndBs[i][2]
                ConstraintCenter.InitialAzimuth=Coupling.Directions[i]
                ConstraintCenter.InitialDeclination=0
                ConstraintCenter.InitialGamma=0
                # for dof in range(6):
                #     ConstraintCenter.DOFFree[dof] = 'Yes'
                # ConstraintCenter.DOFFree[2] = 'Yes'
                # ConstraintCenter.DOFFree[4] = 'Yes'
                
                # BuoyCenter.Connection = ConCenterName
                # BuoyCenter.InitialX=0
                # BuoyCenter.InitialY=0
                # BuoyCenter.InitialZ=0
                # BuoyCenter.InitialRotation3=0
                # Link.EndBConnection=BuoyCenterName
                if 'Mid' in Coupling.Names[i]:
                    Link.EndBConnection=ConBarName
                else:
                    # Link.EndBConnection=BuoyCenterName
                    Link.EndBConnection=ConLinkName
                
            if 'Mid' in Coupling.Names[i]:
                # Link.IncludeTorsion='Yes'
                Link.TopEnd='End B'
                Link.EndAAzimuth=180
                Link.EndBAzimuth=0
                Link.EndADeclination=ConBar.InitialDeclination
                Link.EndBDeclination=0 #ConBar.InitialDeclination
                Link.EndAGamma=0
                Link.EndBGamma=0
                Link.EndBxBendingStiffness=1e+307
                Link.EndByBendingStiffness=1e+307
                # Link.EndBTwistingStiffness=100e6
                # self.set_object_EndA_EndB(Link, [0,0,0],[-0.2218,0,-0.16817])
                self.set_object_EndA_EndB(Link, [0,0,0],[0,0,0])
            elif 'In' in Coupling.Names[i]:
                if '-1' in str(Coupling.Locs[i]):
                    self.set_object_EndA_EndB(Link, [0,0,0],[0,-0.125,0])
                else:
                    self.set_object_EndA_EndB(Link, [0,0,0],[0,0.125,0])
            else:
                if '-1' in str(Coupling.Locs[i]):
                    self.set_object_EndA_EndB(Link, [0,0,0],[0,0.125,0])

                else:
                    self.set_object_EndA_EndB(Link, [0,0,0],[0,-0.125,0])
    
    def connect_EdegePlatform_constraint_links(self, Platform , Coupling, Gen):
        for i in range(len(Coupling.Names)):
            if int(Coupling.Names[i].split('_')[1]) in Platform.EdgePlatforms:
                Link=Gen.edge_coupling(Coupling, Coupling.Names[i], Coupling.Types[i])
                Constraint=Gen.constraint('Con_' + Coupling.Names[i])
                
                Constraint.Connection=Coupling.EndAConnections[i]
                Constraint.InitialX=Coupling.PosEndAs[i][0]
                Constraint.InitialY=Coupling.PosEndAs[i][1]
                Constraint.InitialZ=Coupling.PosEndAs[i][2]
                Constraint.InitialAzimuth=Coupling.Directions[i]
                Constraint.InitialDeclination=0
                Constraint.InitialGamma=0
       
                Link.EndAConnection='Con_'+ Coupling.Names[i]
                ConCenterName='Con_Cen_'+ Coupling.EndBConnections[i][9:].split('_')[0]+ '_'+ str(Coupling.Directions[i])+'_'+ str(Coupling.Locs[i]*-1)
                try: 
                    Link.EndBConnection=ConCenterName
                except:
                    ConstraintCenter=Gen.constraint(ConCenterName)
                    ConstraintCenter.Connection=Coupling.EndBConnections[i]
                    ConstraintCenter.InitialX=Coupling.PosEndBs[i][0]
                    ConstraintCenter.InitialY=Coupling.PosEndBs[i][1]
                    ConstraintCenter.InitialZ=Coupling.PosEndBs[i][2]
                    ConstraintCenter.InitialAzimuth=Coupling.Directions[i]
                    ConstraintCenter.InitialDeclination=0
                    ConstraintCenter.InitialGamma=0
                    Link.EndBConnection=ConCenterName  
                self.set_object_EndA_EndB(Link, [0,0,0],[0,0,0])
                             
    def connect_mooring(self, Platform, Layout, Mooring, Gen, Coupling, MooringBuoy):
        ## posability of mooring and tether connection
        for i in range(len(Mooring.ConnectedPlatform)):
            UniqueName=str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.OrientationLocal[i])+  '_' +str(round(Mooring.Loc[i]))
            if Mooring.Type=='Tether':
                Line=Gen.mooring_tether(Mooring.M1LineTypes, Mooring.M1Lengths, Mooring.M1SegementLengths, 'M1_' + UniqueName )
                # Line=Gen.mooring_tether(Mooring.M1LineTypes, Mooring.M1Lengths, Mooring.M1SegementLengths, 'M1_' + str(Mooring.ConnectedPlatform[i])+'_' +str(Mooring.Loc[i]))
            elif Mooring.Type=='Line':
                Line=Gen.line(Mooring.M1LineTypes, Mooring.M1Lengths, Mooring.M1SegementLengths, 'M1_' +UniqueName)
                Line.FullStaticsConvergenceControlMethod=Mooring.FullStaticsConvergenceControlMethod
                Line.StaticsStep2=Mooring.StaticsStep2
                # Line.LayAzimuth=np.mod(Mooring.OrientationGlobal[i]+180, 360) #OrcaFlex regonizes layazimuth from anchor point, therefore the 180 deg shift
                Line.EndAAzimuth = Mooring.OrientationGlobal[i]
                Line.EndADeclination = 90
                Line.EndBAzimuth = Mooring.OrientationGlobal[i]
                Line.EndBDeclination = 90
                Line.ContentsMethod=Mooring.ContentsMethod
                Line.IncludeAxialContentsInertia=Mooring.IncludeAxialContentsInertia
                Line.StaticsSeabedFrictionPolicy=Mooring.StaticsSeabedFrictionPolicy
                Line.Representation=Mooring.Representation

            MooringPos=self.translate_coordinates(Platform, Mooring.Position*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])               
            if Platform.Split:
                ConnectedPlatform=str(Mooring.ConnectedPlatform[i])+'_'+ str(round(np.mod(Mooring.OrientationLocal[i]+60-60*Mooring.Loc[i],360)))
            else:
                ConnectedPlatform=str(Mooring.ConnectedPlatform[i])

            if Coupling.LoadCells:
                OFCon=Gen.constraint('Con_Frld_' + UniqueName)
                OFCon.Connection='Platform_'+ConnectedPlatform
                OFCon.InitialAzimuth=Mooring.OrientationLocal[i]
                OFCon.InitialDeclination=0
                OFCon.InitialGamma=0
                self.set_object_position(OFCon, MooringPos)
                Line.EndAConnection='Con_Frld_' + UniqueName
                self.set_object_EndA(Line,[0,0,0])                
            else:
                Line.EndAConnection='Platform_'+ConnectedPlatform
                self.set_object_EndA(Line, MooringPos)
                
            Line.EndBConnection='Platform_'+ConnectedPlatform
            Line.EndBX=MooringPos[0]
            Line.EndBY=MooringPos[1]
            Line.EndBConnection='Anchored' #'Fixed'
            Line.EndBxBendingStiffness=0.2
            Line.EndBX=Line.EndBX+Mooring.RadiusAnchor*np.cos((Mooring.OrientationGlobal[i]-Mooring.Angle*Mooring.Loc[i])*np.pi/180)
            Line.EndBY=Line.EndBY+Mooring.RadiusAnchor*np.sin((Mooring.OrientationGlobal[i]-Mooring.Angle*Mooring.Loc[i])*np.pi/180)
            Line.EndBHeightAboveSeabed=-0.1
            # Line.EndBZ=-1*Mooring.Depth
            
            if Mooring.MooringBuoy=='Yes':
                Line.NumberOfAttachments=len(MooringBuoy.BuoyLocation)
                Line.AttachmentType=[MooringBuoy.Name for i in range(len(MooringBuoy.BuoyLocation))]
                Line.Attachmentz=MooringBuoy.BuoyLocation
                Line.AttachmentzRelativeTo=[MooringBuoy.ReferenceEnd for i in range(len(MooringBuoy.BuoyLocation))]
    
    def place_mooring_line(self,Name,Gen,Mooring,i,VLink=False):
        if VLink:
            LineTypes = Mooring.VLinkLineTypes
            LineLengths = Mooring.M2Lengths
            SegementLengths = Mooring.M2SegementLengths
            LineName = 'M2_' + Name
        else:
            LineTypes = Mooring.M1LineTypes
            LineLengths = Mooring.M1Lengths
            SegementLengths = Mooring.M1SegementLengths
            LineName = 'M1_' + Name
        
        if Mooring.Type=='Tether':
            Line=Gen.mooring_tether(LineTypes, LineLengths, SegementLengths, LineName )
        elif Mooring.Type=='Line':
            Line=Gen.line(LineTypes, LineLengths, SegementLengths, LineName)
            Line.FullStaticsConvergenceControlMethod=Mooring.FullStaticsConvergenceControlMethod
            Line.StaticsStep2=Mooring.StaticsStep2
            # Line.LayAzimuth=np.mod(Mooring.OrientationGlobal[i]+180, 360) #OrcaFlex regonizes layazimuth from anchor point, therefore the 180 deg shift
            Line.EndAAzimuth = Mooring.OrientationGlobal[i]
            Line.EndADeclination = 90
            Line.EndBAzimuth = Mooring.OrientationGlobal[i]
            Line.EndBDeclination = 90
            Line.ContentsMethod=Mooring.ContentsMethod
            Line.IncludeAxialContentsInertia=Mooring.IncludeAxialContentsInertia
            Line.StaticsSeabedFrictionPolicy=Mooring.StaticsSeabedFrictionPolicy
        if not VLink:
            Line.Representation=Mooring.Representation
        return Line
        
    def connect_vlink_mooring_sameplatform(self,Platform, Layout, Mooring, Gen, Coupling, MooringBuoy):
        #### Posibility of VLink mooring of the same platform
        for i in range(len(Mooring.ConnectedPlatform)):
            MooringPos=self.translate_coordinates(Platform, Mooring.Position*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
            MooringPos_anchor=self.translate_coordinates(Platform, np.array([Mooring.Position[0],0,Mooring.Position[2]])*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
            # GlobalPos=self.local_to_global(Layout, Mooring.ConnectedPlatform[i], MooringPos) 
            if Mooring.Loc[i]==1:
                #### place mooring vlink line
                VLinkName=str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.OrientationLocal[i])+  '_' +str(round(Mooring.Loc[i]))
                Line1 = self.place_mooring_line(VLinkName,Gen,Mooring,i,VLink=True)
    
                if Platform.Split:
                    ConnectedPlatform=str(Mooring.ConnectedPlatform[i])+'_'+ str(round(np.mod(Mooring.OrientationLocal[i]+60-60*Mooring.Loc[i],360)))
                else:
                    ConnectedPlatform=str(Mooring.ConnectedPlatform[i])
    
                if Coupling.LoadCells:
                    OFCon=Gen.constraint('Con_Frld_' + VLinkName)
                    OFCon.Connection='Platform_'+ConnectedPlatform
                    OFCon.InitialAzimuth=Mooring.OrientationLocal[i]
                    OFCon.InitialDeclination=0
                    OFCon.InitialGamma=0
                    self.set_object_position(OFCon, MooringPos)
                    Line1.EndAConnection='Con_Frld_' + VLinkName
                    self.set_object_EndA(Line1,[0,0,0])                
                else:
                    Line1.EndAConnection='Platform_'+ConnectedPlatform
                    self.set_object_EndA(Line1, MooringPos)
            
                BuoyName = 'VLinkBuoy'+str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.OrientationLocal[i])
                Buoy=Gen.vlink_buoy(BuoyName)
                VLinkBuoyPos = self.translate_coordinates(Platform, Mooring.VLinkBuoyPositon*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
                # GlobalPosVLinkBuoy = self.local_to_global(Layout, Mooring.ConnectedPlatform[i], Mooring.VLinkBuoyPositon) 
                Buoy.Connection = 'Platform_'+ConnectedPlatform
                self.set_object_position(Buoy, VLinkBuoyPos) 
                Buoy.Connection = 'Free'                           
                
                #### place mooring line M2
                MainLineName = str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.OrientationLocal[i]) + '_' +str(round(Mooring.Loc[i]))
                Line2 = self.place_mooring_line(MainLineName,Gen,Mooring,i,VLink=False)
    
                if Platform.Split:
                    ConnectedPlatform=str(Mooring.ConnectedPlatform[i])+'_'+ str(round(np.mod(Mooring.OrientationLocal[i]+60-60*Mooring.Loc[i],360)))
                else:
                    ConnectedPlatform=str(Mooring.ConnectedPlatform[i])
                
                Line2.EndAConnection = BuoyName
                self.set_object_EndA(Line2,[0,0,0]) 
                Line2.EndBConnection='Platform_'+ConnectedPlatform
                # Line2.EndBX=MooringPos[0]
                # Line2.EndBY=MooringPos_anchor[1] 
                # Line2.EndBX=Line2.EndBX+Mooring.RadiusAnchor*np.cos((Mooring.OrientationGlobal[i]-Mooring.Angle*Mooring.Loc[i])*np.pi/180)
                # Line2.EndBY=Line2.EndBY+Mooring.RadiusAnchor*np.sin((Mooring.OrientationGlobal[i]-Mooring.Angle*Mooring.Loc[i])*np.pi/180)
                if not Mooring.AdjustBridalCornerAnchorPos:
                    VLinkAnchorPos = self.translate_coordinates(Platform, Mooring.VlinkMidAnchorPos*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
                else:
                    #### hardcode platform numbers, later need update to autochange for Hex1,2,3,4
                    if Mooring.ConnectedPlatform[i] in Platform.CornerPlatforms:
                        if  Mooring.ConnectedPlatform[i] in [Platform.InCorners[0][0],Platform.InCorners[0][2],Platform.InCorners[1][1],Platform.InCorners[2][0],Platform.InCorners[3][1],Platform.InCorners[3][3]]:#[10,47,96,91,54,5]:
                            VLinkAnchorPos = self.translate_coordinates(Platform, Mooring.VlinkCornerAnchorPos1*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
                        elif Mooring.ConnectedPlatform[i] in [Platform.InCorners[0][1],Platform.InCorners[0][3],Platform.InCorners[1][0],Platform.InCorners[2][1],Platform.InCorners[3][0],Platform.InCorners[3][2]]:#[6,11,64,95,90,37]:
                            VLinkAnchorPos = self.translate_coordinates(Platform, Mooring.VlinkCornerAnchorPos2*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
                    else:
                        VLinkAnchorPos = self.translate_coordinates(Platform, Mooring.VlinkMidAnchorPos*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
                Line2.EndBX = VLinkAnchorPos[0]
                Line2.EndBY = VLinkAnchorPos[1]
                Line2.EndBConnection='Anchored' #'Fixed'
                Line2.EndBxBendingStiffness=0.2
                Line2.EndBHeightAboveSeabed=-0.1
                
                Line1.EndBConnection = BuoyName
                
                self.set_object_EndB(Line1,[0,0,0]) 
                if Mooring.MooringBuoy=='Yes':
                    Line2.NumberOfAttachments=len(MooringBuoy.BuoyLocation)
                    Line2.AttachmentType=[MooringBuoy.Name for i in range(len(MooringBuoy.BuoyLocation))]
                    Line2.Attachmentz=MooringBuoy.BuoyLocation
                    Line2.AttachmentzRelativeTo=[MooringBuoy.ReferenceEnd for i in range(len(MooringBuoy.BuoyLocation))]
            elif Mooring.Loc[i]==-1:
                # MooringPos=self.translate_coordinates(Platform, Mooring.Position*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i])
                # GlobalPos=self.local_to_global(Layout, Mooring.ConnectedPlatform[i], MooringPos) 
                
                #### place mooring vlink line
                VLinkName=str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.OrientationLocal[i])+  '_' +str(round(Mooring.Loc[i]))
                Line1 = self.place_mooring_line(VLinkName,Gen,Mooring,i,VLink=True)
    
                if Platform.Split:
                    ConnectedPlatform=str(Mooring.ConnectedPlatform[i])+'_'+ str(round(np.mod(Mooring.OrientationLocal[i]+60-60*Mooring.Loc[i],360)))
                else:
                    ConnectedPlatform=str(Mooring.ConnectedPlatform[i])
    
                if Coupling.LoadCells:
                    OFCon=Gen.constraint('Con_Frld_' + VLinkName)
                    OFCon.Connection='Platform_'+ConnectedPlatform
                    OFCon.InitialAzimuth=Mooring.OrientationLocal[i]
                    OFCon.InitialDeclination=0
                    OFCon.InitialGamma=0
                    self.set_object_position(OFCon, MooringPos)
                    Line1.EndAConnection='Con_Frld_' + VLinkName
                    self.set_object_EndA(Line1,[0,0,0])                
                else:
                    Line1.EndAConnection='Platform_'+ConnectedPlatform
                    self.set_object_EndA(Line1, MooringPos)
                Line1.EndBConnection = 'VLinkBuoy'+str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.OrientationLocal[i])
                self.set_object_EndB(Line1,[0,0,0])
                
    def connect_mooring_with_buoy(self, Platform, Layout, Mooring, Gen, MooringBuoy):
        #### This function is outdated, it however may be used when wanting to connect a mooring buoy to a tether 
        #### For this function seperate buoys are created, the more optimal version is to add attatchements to a line 
        #### The above is now included in the connect_mooring_function 
        #### This function is activated with mooring.mooringbuoy='OldFunction'
        for i in range(len(Mooring.ConnectedPlatform)):
            if Mooring.Type=='Tether':
                # Line=Gen.mooring_tether(Mooring.M2LineTypes, Mooring.M2Lengths, Mooring.M2SegementLengths, 'M2_' + str(Mooring.ConnectedPlatform[i])+'_' +str(Mooring.Loc[i]))
                Line=Gen.mooring_tether(Mooring.M2LineTypes, Mooring.M2Lengths, Mooring.M2SegementLengths, 'M2_' + str(Mooring.ConnectedPlatform[i]))
            elif Mooring.Type=='Line':
                Line=Gen.mooring(Mooring.M2LineTypes, Mooring.M2Lengths, Mooring.M2SegementLengths, 'M2_' + str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.Loc[i]))
                Line.LayAzimuth=np.mod(Mooring.OrientationLocal[i], 360) #OrcaFlex regonizes layazimuth from anchor point, therefore the 180 deg shift
                #### This has been added for the springs in the Cork model test
                Line.ContentsMethod=Mooring.ContentsMethod
                Line.IncludeAxialContensInertia='No'
                
            MooringPos=self.translate_coordinates(Platform, Mooring.Position*[1,Mooring.Loc[i],1], Mooring.OrientationLocal[i]) ## Local MooringConncetion
            GlobalPos=self.local_to_global(Platform, Layout, Mooring.ConnectedPlatform[i], MooringPos)
            Buoy=Gen.vlink_buoy('3DBuoy'+str(Mooring.ConnectedPlatform[i])+'a', MooringBuoy)
            Buoy.InitialX=GlobalPos[0]+ MooringBuoy.MooringRadius*np.cos(Mooring.OrientationGlobal[i]*np.pi/180)
            Buoy.InitialY=GlobalPos[1]+ MooringBuoy.MooringRadius*np.sin(Mooring.OrientationGlobal[i]*np.pi/180)
            Buoy.InitialZ=0
            ## Connect mooring line M2 
            Line.EndAConnection='3DBuoy' +str(Mooring.ConnectedPlatform[i])+'a'
            Line.EndAX=0
            Line.EndAY=0
            Line.EndAZ=0
            Line.EndBConnection='Platform'+str(Mooring.ConnectedPlatform[i])
            Line.EndBConnection='Anchored' #'Fixed'
            Line.EndBxBendingStiffness=0.2
            Line.EndBX=GlobalPos[0]+Mooring.RadiusAnchor*np.cos(Mooring.OrientationGlobal[i]*np.pi/180)
            Line.EndBY=GlobalPos[1]+Mooring.RadiusAnchor*np.sin(Mooring.OrientationGlobal[i]*np.pi/180)
            Line.EndBHeightAboveSeabed=-0.1
            #Line.EndBZ=-1*Mooring.Depth
            ## Conncet mooring line M1a
            Line=Gen.mooring(Mooring.M1aLineTypes, Mooring.M1aLengths, Mooring.M1aSegementLengths,'Link_' + str(Mooring.ConnectedPlatform[i])+'_' + str(Mooring.Loc[i]))
                # Line.LayAzimuth=np.mod(Mooring.Orientation[i]+180, 360) #OrcaFlex regonizes layazimuth from anchor point, therefore the 180 deg shift
            Line.EndAConnection='Platform'+str(Mooring.ConnectedPlatform[i])
            Line.EndAX=MooringPos[0]
            Line.EndAY=MooringPos[1]
            Line.EndAZ=MooringPos[2]
            Line.EndBConnection='3DBuoy' +str(Mooring.ConnectedPlatform[i])+'a'
            Line.EndBX=0
            Line.EndBY=0
            Line.EndBZ=0
    
    def connect_vlink_mooring(self, Platform, Mooring, Layout, Gen, MooringBuoy):
        for i in range(0,np.size(Mooring.ConnectedPlatform)-1):
            Position=[Mooring.Position[0], Layout.EdgeLength/2, Mooring.Position[2]]
            MooringPosCenter=self.translate_coordinates(Platform, Position*np.array([1,Mooring.Loc[i],1]), Mooring.OrientationLocal[i])
            GlobalPos=self.local_to_global(Layout, Mooring.ConnectedPlatform[i], MooringPosCenter)          
            # Place Buoy                
            Buoy=Gen.vlink_buoy('3DBuoy'+str(Mooring.ConnectedPlatform[i]), MooringBuoy)
            Buoy.InitialX=GlobalPos[0]+ Mooring.VLinkMidLength*np.cos(Mooring.OrientationGlobal[i]*np.pi/180)
            Buoy.InitialY=GlobalPos[1]+ Mooring.VLinkMidLength*np.sin(Mooring.OrientationGlobal[i]*np.pi/180)
            Buoy.InitialZ=GlobalPos[2]
            # Place Mooring line M2
            if Mooring.Type== 'Line':
                Line=Gen.mooring(Mooring.M2LineTypes, Mooring.M2Lengths, Mooring.M2SegementLengths, 'M2_'+ str(Mooring.ConnectedPlatform[i])+'_' +str(Mooring.Loc[i]))
                Line.LayAzimuth=np.mod(Mooring.OrientationLocal[i], 360)
            if Mooring.Type=='Tether':
                Line=Gen.mooring_tether(Mooring.M2LineTypes, Mooring.M2Lengths, Mooring.M2SegementLengths, 'M2_'+ str(Mooring.ConnectedPlatform[i])+'_' +str(Mooring.Loc[i]))
            Line.EndAConnection='3DBuoy' +str(Mooring.ConnectedPlatform[i])
            Line.EndAX=0
            Line.EndAY=0
            Line.EndAZ=0
            Line.EndBConnection='Anchored' #'Fixed'
            Line.EndBX=GlobalPos[0]+Mooring.RadiusAnchor*np.cos(Mooring.OrientationGlobal[i]*np.pi/180)
            Line.EndBY=GlobalPos[1]+Mooring.RadiusAnchor*np.sin(Mooring.OrientationGlobal[i]*np.pi/180)
            # Line.EndBZ=-1*Mooring.Depth
            Line.EndBHeightAboveSeabed=-0.1
            # Place VLinks
            for j in range(2):
                ## Commented lines below are to make a tether tether combination, senstivities have show this does not make simulations faster
                # if Type== 'Line':
                VLink=Gen.mooring(Mooring.VLinkLineType, Mooring.VLinkLength, Mooring.VLinkSegementLength, 'VLink_'+ str(Mooring.ConnectedPlatform[i])+'_'+str(j))
                # if Type== 'Tether':
                    # VLink=Gen.mooring_tether([Mooring.VLinkLineType], [Mooring.VLinkLength], [Mooring.VLinkSegementLength], 'VLink_'+ str(i)+'_'+str(j))
                VLink.EndBConnection='3DBuoy' +str(Mooring.ConnectedPlatform[i])
                VLink.EndBX=0
                VLink.EndBY=0
                VLink.EndBZ=0
                Loc=1
                if j==1:
                    Loc=-1                
                MooringPos=self.translate_coordinates(Platform, Mooring.Position*[1,Mooring.Loc[i]*Loc,1], Mooring.OrientationLocal[i])
                VLink.EndAConnection='Platform' +str(Mooring.ConnectedPlatform[i+j])
                VLink.EndAX=MooringPos[0]
                VLink.EndAY=MooringPos[1]
                VLink.EndAZ=MooringPos[2]
                
    def connect_mooring_automatic(self, Layout, Platform, Mooring, Gen, MooringBuoy, Coupling):
        ## this function can be used to autormatically generate the mooring is all sides are the same           
        ## sort the connections to be able to make Vlinks
            ConPlat=np.asarray(Mooring.ConnectedPlatformsTotal)
            OrientLocalTot=np.asarray(Mooring.OrientationLocalTotal)
            ConType=np.asarray(Mooring.ConnectedPlatformType)
            LocTot=np.asarray(Mooring.LocTotal)
            
            if Mooring.MooringBuoy== 'Yes':
                try:
                    OFMooringBuoy=self.model[MooringBuoy.Name]
                except: 
                    OFMooringBuoy=Gen.clump_weight(MooringBuoy.Name, MooringBuoy)

            for i,t in enumerate(Mooring.ConnectedType):
                for j,direc in enumerate(Mooring.ConnectedDirection):
                    ## Connect V Links -> order of connections is important as they are linked
                    ## That is sorted out in the next few lines
                    A=np.array([ConType==t,OrientLocalTot==direc])
                    Plats=ConPlat[A.all(0)]
                    # Locs=LocTot[A.all(0)]
                    Mooring.ConnectedPlatform=Plats
                    Mooring.OrientationLocal=OrientLocalTot[A.all(0)]
                    Mooring.OrientationGlobal=np.mod(OrientLocalTot[A.all(0)]+180*(t-1)+180,360)
                    Mooring.Loc=LocTot[A.all(0)]
                    
                    if Mooring.VLink:
                        ## Create Vlink type mooring system
                        # self.connect_vlink_mooring(Platform,Mooring, Layout, Gen, MooringBuoy)  
                        self.set_parrallel_mooring_locations(Mooring)
                        self.connect_vlink_mooring_sameplatform(Platform, Layout, Mooring, Gen, Coupling, MooringBuoy)
                        ## Connect single mooring lines at the corners
                        if Mooring.CornerMooring:
                            Mooring.ConnectedPlatform=[Mooring.ConnectedPlatform[0], Mooring.ConnectedPlatform[-1]]
                            Mooring.OrientationLocal=[Mooring.OrientationLocal[0], Mooring.OrientationLocal[-1]]
                            Mooring.OrientationGlobal=[Mooring.OrientationGlobal[0],Mooring.OrientationGlobal[-1]]
                            Mooring.Loc=np.array([Mooring.Loc[0],Mooring.Loc[-1]])*[-1,1]
                            if MooringBuoy==0:
                                self.connect_mooring(Platform, Layout, Mooring, Gen)
                            else: 
                                self.connect_mooring_with_buoy(Platform, Layout, Mooring, Gen, MooringBuoy)
                    else: 
                        ## Create distributed mooring system
                        if not Mooring.Center:
                            self.set_parrallel_mooring_locations(Mooring)
                            
                        elif Mooring.Center:
                            pass
                        
                        if Mooring.MooringBuoy=='OldFunction':
                            self.connect_mooring_with_buoy(Platform, Layout, Mooring, Gen, MooringBuoy)
                        else:
                            self.connect_mooring(Platform, Layout, Mooring, Gen, Coupling, MooringBuoy)    
                            
    def set_parrallel_mooring_locations(self,Mooring):
        Mooring.Loc=np.concatenate((np.ones(Mooring.ConnectedPlatform.shape),np.ones(Mooring.ConnectedPlatform.shape)*-1))
        Mooring.OrientationLocal=np.concatenate((Mooring.OrientationLocal,  Mooring.OrientationLocal))
        Mooring.OrientationGlobal= np.concatenate((Mooring.OrientationGlobal,Mooring.OrientationGlobal))
        Mooring.ConnectedPlatform=np.concatenate((Mooring.ConnectedPlatform, Mooring.ConnectedPlatform))
                            
    def connect_extra_mooringbuoy(self, MooringBuoy, Gen):
        OFMooringBuoy=Gen.clump_weight(MooringBuoy.Name, MooringBuoy)
        for obj in self.model.objects:
            if obj.Name.find('M1_')!=-1:
                Line=obj
                # Line.NumberOfAttachments=Line.NumberOfAttachments+len(MooringBuoy.BuoyLocation)
                AT=list(Line.AttachmentType)
                Az=list(Line.Attachmentz)
                ARt=list(Line.AttachmentzRelativeTo)
                Line.AttachmentType=AT+[MooringBuoy.Name for i in range(len(MooringBuoy.BuoyLocation))]
                Line.Attachmentz=Az+list(MooringBuoy.BuoyLocation)
                Line.AttachmentzRelativeTo=ARt+[MooringBuoy.ReferenceEnd for i in range(len(MooringBuoy.BuoyLocation))]
                       
    def connect_pendulum(self, Layout, Pendulum, Gen):
        ## Must be connected after the platform and floaters
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0:   
                c=0.2
                for j in Pendulum.Mass.keys():
                    Buoy=Gen.PendulumBuoy('Buoy_'+str(j)+'_'+str(i), Pendulum, j)
                    Buoy.Connection='Platform'+str(i)
                    Buoy.InitialX=Pendulum.Pos[j][0]
                    Buoy.InitialY=Pendulum.Pos[j][1]
                    Buoy.InitialZ=Pendulum.Pos[j][2]-c
                    Buoy.Connection='Free'
                    c=c+0.2
                for d in [0,120,240]:
                    Link=Gen.tether('Link_Center_'+str(i)+'_'+str(d))
                    Link.EndAConnection='Floater'+str(i)+'_'+str(d)
                    Link.EndBConnection='Buoy_Center_'+str(i)
                    self.set_EndA_EndB(Link, [0,0,0,0,0,0])
                    Link.EndAConnection='Platform'+str(i)
                    Link.UnstretchedLength=Pendulum.LenLink
                    Link.Stiffness=Pendulum.StiffnessRope
                    
                Link=Gen.tether('Link_Pend_'+str(i)+'_'+str(d))
                Link.EndAConnection='Buoy_Center_'+str(i)
                Link.EndBConnection='Buoy_Pend_'+str(i)
                self.set_EndA_EndB(Link, [0,0,0,0,0,0])
                Link.UnstretchedLength=Pendulum.LenPedulum
                Link.Stiffness=Pendulum.StiffnessRope
            
    def connect_wings_old(self,Wing, Layout, Pendulum, Gen):
        #### This is an onlder function that has been kept in this script as and example 
        for Type in Wing.WingTypes:
            Gen.wingtype(Wing, Type)
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0:   
                OFPlat=self.model['Platform'+str(i)]
                OFPlat.NumberOfWings=1
                Row=np.floor_divide(i,Layout.Matrix.shape[1])
                if Row>=len(Wing.RowReductionFactors):
                    RF=Wing.RowReductionFactors[-1]
                else:
                    RF=Wing.RowReductionFactors[Row]
                OFPlat.WingSpan=[Wing.Span*RF**0.5]
                OFPlat.WingChord=[Wing.Chord*RF**0.5]
                OFPlat.WingAzimuth=[Wing.Azimuth]
                OFPlat.WingDeclination=[Wing.Declination]
                OFPlat.WingGamma=[Wing.Gamma]
                OFPlat.WingType=[Wing.WingTypes[Layout.Types[i]-1]]
                OFPlat.WingCenterX=[Wing.ConnectionPosition[0]]
                OFPlat.WingCenterY=[Wing.ConnectionPosition[1]]
                OFPlat.WingCenterZ=[Wing.ConnectionPosition[2]]   
                
    def connection_wings(self, Wing, Layout, Platform, Env):
        #### This function assumes that the wing types have been pre-defined, similar to the line types
        #### This fuctions should be applied after the envriomental direction is set
        for i in range(len(Layout.Types)):
            if Layout.Types[i] != 0:   
                if Platform.split: 
                    OFPlat=self.model['Platform'+str(i)+'_0']
                else:
                    OFPlat=self.model['Platform'+str(i)]
                
                OFPlat.NumberOfWings=1
                # The code below only works fora single direction
                # Row=np.floor_divide(i,Layout.Matrix.shape[1])
                # if Row>=len(Wing.RowReductionFactors):
                #     RF=Wing.RowReductionFactors[-1]
                # else:
                #     RF=Wing.RowReductionFactors[Row]
                OFPlat.WingSpan=[Wing.Span]
                OFPlat.WingChord=[Wing.Chord]
                OFPlat.WingAzimuth=[Wing.Azimuth]
                OFPlat.WingDeclination=[Wing.Declination]
                OFPlat.WingGamma=[Wing.Gamma]
                
                dOrientation=(OFPlat.InitialRotation3-Env.Heading)
                dOrientation1=np.mod(dOrientation+30,120)
                if dOrientation1<60:
                    OFPlat.WingType=[Wing.WingTypes[Layout.Types[i]-1]]
                
                
                OFPlat.WingCenterX=[Wing.ConnectionPosition[0]]
                OFPlat.WingCenterY=[Wing.ConnectionPosition[1]]
                OFPlat.WingCenterZ=[Wing.ConnectionPosition[2]]   
                
    def floater_braces(self,model,Layout,Floater,Gen):
        self.BraceLength = Floater.PosBottom[1]*2 - Floater.Diameter
        self.BraceZpos = -9 # wrt floater axes
        for obj in model.objects:
            if obj.name.find('Floater') != -1:
                # print(obj.name)
                if Layout.Types[int(obj.name.split('_')[0][-1])]==1:
                    LineType = ['TrussOutT1']
                else:
                    LineType = ['TrussOutT2']
                
                BraceName = 'Brace_'+obj.name.split('Floater')[1]
                OFLine = Gen.line(LineType, [self.BraceLength], [self.BraceLength/4], BraceName)
                OFLine.EndBConnection = obj.name
                deg = int(obj.name.split('_')[1])
                if deg == 0:
                    OFLine.EndAConnection = obj.name.split('_')[0]+'_'+str(deg+120)
                    self.set_object_EndA_EndB(OFLine, [0,Floater.Diameter/2,self.BraceZpos], [0,-Floater.Diameter/2,self.BraceZpos])
                elif deg == 120:
                    OFLine.EndAConnection = obj.name.split('_')[0]+'_'+str(deg+120)
                    self.set_object_EndA_EndB(OFLine, [-Floater.Diameter/2*np.sin(np.radians(60)),-Floater.Diameter/2*np.cos(np.radians(60)),self.BraceZpos], [Floater.Diameter/2*np.sin(np.radians(60)),Floater.Diameter/2*np.cos(np.radians(60)),self.BraceZpos])
                else:
                    OFLine.EndAConnection = obj.name.split('_')[0]+'_'+str(0)
                    self.set_object_EndA_EndB(OFLine, [Floater.Diameter/2*np.sin(np.radians(60)),-Floater.Diameter/2*np.cos(np.radians(60)),self.BraceZpos], [-Floater.Diameter/2*np.sin(np.radians(60)),Floater.Diameter/2*np.cos(np.radians(60)),self.BraceZpos])
                OFLine.EndAAzimuth = deg + 90
                OFLine.EndBAzimuth = deg + 90 
                # self.set_object_EndA_EndB(OFLine, [0,0,-9], [0,0,-9])
                OFLine.EndADeclination = 90
                OFLine.EndBDeclination = 90
                OFLine.EndAGamma = 0
                OFLine.EndBGamma = 0
                OFLine.EndAxBendingStiffness = 1e+307
                OFLine.EndBxBendingStiffness = 1e+307
                OFLine.EndAyBendingStiffness = 1e+307
                OFLine.EndByBendingStiffness = 1e+307
            
            
# class ModelRunning:
    #### this function is to achieve OrcaFlex automation of running models in batch processing
    #### Yamls path and threads are as input
    #### simulation will be automatically save at yams path with the same name as Yamls but with .sim extension
    #### If plan to use single thread, then give Threads=1
    # def __init__(self):#,Threads):
    #     # self.Threads = Threads
    #     self.LocYamls = ModelPath
    #     self.YmlFiles = glob.glob(os.path.join(self.LocYamls,'*yml'))
        
def run_simulation(yaml_file):
    #### function for auto run OrcaFlex simulations using python parallel
    oFmodel = ofx.Model(yaml_file)
    oFmodel.RunSimulation()
    return oFmodel.SaveSimulation(yaml_file[:-4]+'.sim')

    # def parallel_running(self): 
        
        # if __name__ == '__main__':
        # with Pool(processes = self.Threads) as pool: 
        #     pool.map(self.run_simulation,self.YmlFiles)            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
