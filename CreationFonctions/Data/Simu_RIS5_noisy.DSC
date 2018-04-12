#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1760
XMIN	11979.000000
XWID	199.886364
*
*	Documentational Text:
*
TITL	'Simu_RIS5_noisy'
IRNAM	'Intensity'
IINAM	''
XNAM	'Field'
IRUNI	''
IIUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    01/12/18
TIME    07:43:27
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    1.2079
A1SW    0.02
MWFQ    3.400362e+10
MWPW    0.0005557
AVGS    2
RESO    ER5106QT
SPTP    0.02273
RCAG    20
RCHM    2
B0MA    0.0005
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    1024
RCTC    1e-05
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     ezcw, 1.0

MetalIonsDef       False
OrgRadicalsDef     False

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        12079.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     11979.000 G
SweepDirection     Up
SweepWidth         200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       34.003624 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              0.5557 mW
PowerAtten         23.0 dB
QValue             3800

.DVC     recorder, 1.0

BaselineCorr       On
NbScansAcc         2
NbScansDone        2
NbScansToDo        2
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       0

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
Calibrated         True
DModDetectSCT      First
DualDetect         1./2. Harm.
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               20 dB
GateTime           22.73 ms
Harmonic           2
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             5.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           44
NbPoints           1760
Offset             0.0 %
QuadMode           True
Resonator          1
SamplingTime       22.73 ms
ScanTime           40.00 s
SctNorm            True
SctRevision        Allegro
TimeConst          0.01 ms
TimeExp            False
TuneCaps           46

*
************************************************************
