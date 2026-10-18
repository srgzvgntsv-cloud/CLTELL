# Microsoft Developer Studio Project File - Name="CLTELL" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=CLTELL - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "CLTELL.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "CLTELL.mak" CFG="CLTELL - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "CLTELL - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "CLTELL - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "CLTELL - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "CLTELL - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "CLTELL - Win32 Release"
# Name "CLTELL - Win32 Debug"
# Begin Source File

SOURCE=.\CLTELL.f90
DEP_F90_CLTEL=\
	".\Debug\LINES.MOD"\
	".\Debug\TELL_ABU_FIT_MOL.MOD"\
	
# End Source File
# Begin Source File

SOURCE=.\CONVOLVE.f90
DEP_F90_CONVO=\
	".\Debug\LINES.MOD"\
	".\Debug\SLD.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\DATA_MODULES.f90
# End Source File
# Begin Source File

SOURCE=.\DEALLOCATE_ARRAY.f90
DEP_F90_DEALL=\
	".\Debug\INSTRUMENT.MOD"\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FBASE_PLINE.f90
DEP_F90_FBASE=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\LINES.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FIND_INSTRUMENTAL.f90
DEP_F90_FIND_=\
	".\Debug\INSTRUMENT.MOD"\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FIND_RV.f90
DEP_F90_FIND_R=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FORM_TELLURIC.f90
DEP_F90_FORM_=\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GET_COEF_LAYER.f90
DEP_F90_GET_C=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\LINES.MOD"\
	".\Debug\MOLECULE_ISOTOP.mod"\
	".\Debug\PART_FUNCTION.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ILINE_CORRECT.f90
DEP_F90_ILINE=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\INSTRUMENT.MOD"\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\INIT_SLD.f90
DEP_F90_INIT_=\
	".\Debug\LINES.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\LINE_ABSORBTION.f90
DEP_F90_LINE_=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\INSTRUMENT.MOD"\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	".\Debug\TELLURICA.mod"\
	".\Debug\WATER_VAROP_CONT.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\LINE_CORRECT.f90
DEP_F90_LINE_C=\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\LOAD_OBSERVED.f90
DEP_F90_LOAD_=\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\LOAD_STAR_LINE.f90
DEP_F90_LOAD_S=\
	".\Debug\LINES.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MATH_SUBROUTIES.f90
DEP_F90_MATH_=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\INSTRUMENT.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\PART_FUNCTION.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MERGING_LIST_LIN.f90
DEP_F90_MERGI=\
	".\Debug\LINES.MOD"\
	".\Debug\MOLECULE_ISOTOP.mod"\
	".\Debug\SLD.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MODEL_ATMOSPHERE.f90
DEP_F90_MODEL=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\EARTHLY_MODEL_ATM.mod"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELL_ABU_FIT_MOL.MOD"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\MODULES.f90
# End Source File
# Begin Source File

SOURCE=.\RAYLEIGH_SCATTERING.f90
DEP_F90_RAYLE=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\READ_CONFIG.f90
DEP_F90_READ_=\
	".\Debug\LINES.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\READ_HITRAN.f90
DEP_F90_READ_H=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\LINES.MOD"\
	".\Debug\MOLECULE_ISOTOP.mod"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\READ_SLD_GROUP_LINE.f90
DEP_F90_READ_S=\
	".\Debug\OBSERVED.mod"\
	".\Debug\SLD.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\TELL_ABU_CORRECT.f90
DEP_F90_TELL_=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\LINES.MOD"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELLURICA.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\vLINE_CORRECT.f90
DEP_F90_VLINE=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\INSTRUMENT.MOD"\
	".\Debug\LINES.MOD"\
	".\Debug\MODELS.MOD"\
	".\Debug\OBSERVED.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\WATER_VAROP_CONT.f90
DEP_F90_WATER=\
	".\Debug\CONSTANTS.mod"\
	".\Debug\OBSERVED.mod"\
	".\Debug\TELLURICA.mod"\
	".\Debug\WATER_VAROP_CONT.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\CXML.LIB
# End Source File
# End Target
# End Project
