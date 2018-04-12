#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	IGD
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
YFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	3363.250000
XWID	249.500000
YPTS	800
YMIN	0.112500
YWID	179.775000
*
*	Documentational Text:
*
TITL	'imaging_CNRS3'
IRNAM	'Intensity'
XNAM	'Field'
YNAM	'Angle'
IRUNI	''
XUNI	'G'
YUNI	''
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    herve
DATE    06/11/11
TIME    13:54:50
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CWIMG
OXS1    IADC
AXS1    B0VL
AXS2    GRAD
AXS3    
A1CT    0.3488
A1SW    0.02495
MWFQ    9.785e+09
IMGP    0.000
IMGT    -90.000
GRAD    175.000
IMTY    2D
MWPW    0.04004
A2RS    800
A2LL    1
A2SW    1
AVGS    1
RESO    TMHS
SPTP    0.08192
RCAG    77
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    1024
RCTC    0.04096
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        3488.00 G
Delay              1.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3363.250 G
SweepDirection     Up
SweepWidth         249.5 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.785 GHz
QMonitBridge       On

.DVC     grdUnit, 1.0

AnglePhi           0.000 deg
AnglePsi           175.0000 G/cm
AngleTheta         -90.000 deg
BandWidth          0.342 G
CurrentCheck       Wait stable
FieldOfView        10.00 mm
FirstAlpha         0.113 deg
ImageOrient        ZY
ImageType          2D
MaxPsi             0.00 deg
MinSW              249.49 G
NPixel             512
NrOfAlpha          800
NrOfBeta           1
NrOfProfs          800
NrOfPsi            1
PixelSize          0.0195 mm
SWTime             83.890 s
SleepSecs          1.0 s
SpecWidth          1.000 G
SpectralRes        0.00195 G/p
TotalTime          1163.0 min

.DVC     mwBridge, 1.0

AcqFineTuning      Each Slice Scan
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        100.0
Power              40.04 mW
PowerAtten         7 dB
QValue             4400

.DVC     ramp2, 1.0

Flyback            False
NbPoints           800
StartValue         1
SweepData          
SweepWidth         1
sweepType          Ramp up

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        False
AllowShortCt       False
Calibrated         True
ConvTime           81.92 ms
DModAFCTrap        True
DModAmp            1.00 G
DModCalibrated     True
DModDetectSCT      First
DModEliDelay       1.0 us
DModExtLockIn      False
DModExtTrigger     False
DModFieldMod       First
DModGain           60 dB
DModHighPass       True
DModIntegrator     True
DModModOutput      Internal
DModSignalInput    Internal
DModTimeConst      1.28 ms
DoubleModFreq      5.00 kHz
DoubleModPhase     0.0
DoubleMode         False
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               77 dB
Harmonic           1
HighPass           True
Integrator         True
IsCalibExp         False
ModAmp             1.000 G
ModFreq            100.00 kHz
ModInput           Internal
ModOutput          Internal
ModPhase           0.0
Offset             0.0 %
QuadMode           False
QuadPhase          90.0
Resolution         1024
Resonator          1
SctNorm            False
SctRevision        SCT
SignalInput        Internal
SpuExtension       False
SweepTime          83.89 s
TimeConst          40.96 ms
TimeExp            False
TuneCaps           39

*
************************************************************
