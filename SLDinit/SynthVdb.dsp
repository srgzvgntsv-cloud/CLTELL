# Microsoft Developer Studio Project File - Name="SynthVdb" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=SynthVdb - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "SynthVdb.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "SynthVdb.mak" CFG="SynthVdb - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "SynthVdb - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "SynthVdb - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "SynthVdb - Win32 Release"

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

!ELSEIF  "$(CFG)" == "SynthVdb - Win32 Debug"

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

# Name "SynthVdb - Win32 Release"
# Name "SynthVdb - Win32 Debug"
# Begin Source File

SOURCE=.\arrays_allocation.f90
NODEP_F90_ARRAY=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\Hetab.mod"\
	".\Debug\HlineTab.mod"\
	".\Debug\lines_data.mod"\
	".\Debug\mol_dat.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\calculate_spectrum.f90
NODEP_F90_CALCU=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\lines_data.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\continuum_opacity.f90
NODEP_F90_CONTI=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\continuum_wavelengths.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\data_modules.f90
# End Source File
# Begin Source File

SOURCE=.\effective_depths.f90
NODEP_F90_EFFEC=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\first_record.f90
NODEP_F90_FIRST=\
	".\Debug\synth_data.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\HeI_prof.f90
NODEP_F90_HEI_P=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\Hetab.mod"\
	".\Debug\lines_data.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Hydrogen.f90
NODEP_F90_HYDRO=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\HlineTab.mod"\
	".\Debug\Lemkesdat.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Josh.f90
NODEP_F90_JOSH_=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\JoshTab.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\kapp.f90
NODEP_F90_KAPP_=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\KarsTab.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Lemkesdat.f90
# End Source File
# Begin Source File

SOURCE=.\lines_opacity.f90
NODEP_F90_LINES=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\chemical_elements.mod"\
	".\Debug\Hetab.mod"\
	".\Debug\HlineTab.mod"\
	".\Debug\lines_data.mod"\
	".\Debug\mol_dat.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\load_model.f90
NODEP_F90_LOAD_=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\chemical_elements.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\Math_subroutines.f90
# End Source File
# Begin Source File

SOURCE=.\modules.f90
# End Source File
# Begin Source File

SOURCE=.\number_density.f90
NODEP_F90_NUMBE=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\chemical_elements.mod"\
	".\Debug\mol_dat.mod"\
	".\Debug\PfSahaTab.mod"\
	".\Debug\synth_data.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\profiles.f90
DEP_F90_PROFI=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\Hetab.mod"\
	".\Debug\HlineTab.mod"\
	".\Debug\lines_data.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	".\Debug\temporary_arrays.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\synthV_input.f90
DEP_F90_SYNTH=\
	".\Debug\atmosphere_model.mod"\
	".\Debug\chemical_elements.mod"\
	".\Debug\mol_dat.mod"\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SynthVdb.f90
NODEP_F90_SYNTHV=\
	".\Debug\radiation_field.mod"\
	".\Debug\synth_data.mod"\
	
# End Source File
# End Target
# End Project
