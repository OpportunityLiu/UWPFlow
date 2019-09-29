                                      UWPFLOW
           Continuation and Direct Methods to Locate Fold Bifurcations
                              in AC/DC/FACTS Power Systems


                               Claudio A. Canizares
                              University of Waterloo
                            Waterloo, Ontario   N2L-3G1
                                      CANADA
                              ccanizar@uwaterloo.ca
                           http://www.power.uwaterloo.ca


                               Fernando L. Alvarado
                           University of Wisconsin-Madison
                            Madison, Wisconsin   53706
                                       USA
                               alvarado@ece.wisc.edu


                                   Shu Zhang
                              University of Waterloo
                            Waterloo, Ontario   N2L-3G1
                                      CANADA
                          syzhang@engmail.uwaterloo.ca 


                                   Mario Watson
                              University of Waterloo
                            Waterloo, Ontario   N2L-3G1
                                      CANADA
                          mawatson@engmail.uwaterloo.ca 

                                   April 9, 2010
                          (First Version: December 20, 1996)


This program is provided without charge for research purposes only.  The program 
or any of its parts may not be used for any commercial applications. 

The authors would appreciate any comments and suggestions on how to improve the 
program.  Any reports of problems should be directed to the authors, who reserve 
their right to modify the program at any time without previous notification.

DISCLAIMER:  THE AUTHORS DO NOT GUARANTEE THE ACCURACY OF THE RESULTS OBTAINED 
WITH THIS PROGRAM, NOR ITS PERFORMANCE.

________________________________________________________________________________

BRIEF PROGRAM DESCRIPTION:

UWPFLOW is a research tool that has been designed to calculate local 
bifurcations related to system limits or singularities in the system Jacobian.
The program also generates a series of output files that allow further analyses, 
such as tangent vectors, left and right eigenvectors at a singular bifurcation 
point, Jacobians, power flow solutions at different loading levels, voltage 
stability indices, etc.

The program reads ac/dc power flow data in WSCC/BPA/EPRI formats [1,3] or IEEE 
common format [2]; FACTS devices data in a special format described in the on-
line help file (WINDOWS) and using the models described in [5]; steady-state 
load model data in OH format [4]; and steady-state generator data in a simple 
free format as explained in the on-line help.  The program also reads ac data in 
other formats of interest to only some particular users (see on-line help).  

Additional unformatted data is required for bifurcation analysis, such as the 
direction of generation change, direction of load change, and maximum generator 
powers.  The program assumes that one parameter, the "loading factor," is 
allowed to change.  All steady state system controls remain operational unless 
otherwise specified through the program options.

The program has been developed in C and C++ and runs under WINDOWS 7 (and other
previous versions of WINDOWS) and UNIX environments.  It has no limitations on 
system size, other than those imposed by memory limitations in the corresponding 
environment, i.e., RAM and swap space in UNIX and WINDOWS.  The program has been 
successfully used to study a real 3000+ bus system in PCs and a variety of UNIX 
servers.

For more details about the program capabilities, models and the techniques used 
refer to [5,6,7,8,9,10,11,12].  For a PDF copy of most of these documents access 
the WEB server URL:  http://www.power.uwaterloo.ca

________________________________________________________________________________

FILES: 

To get a hold of the program, look under "downloads" at the following URL:

               http://www.power.uwaterloo.ca

The following directories and files form part of the distribution package:

uwpflow/-> Executable files needed to run the program.  In UNIX these files
           must be built from the source code using a C compiler, by
           running the "makefile" with the command "make all".  Script 
           (batch) AWK files for post-processing of output files are also 
           included in this directory.

           The executable "maxim" file is used to determine the maximum 
           entries in the output vectors (tangent, eigenvectors, mismatches)
           A DOS executable public domain version of AWK is also included, 
           so that the AWK file "tomatlab.awk" can be run in DOS and UNIX
           to transform the output Jacobians to MATLAB format for matrix and
           eigenvalue analyses.  

           Add this directory to the path to have access to all the files
           needed to run the program, especially when running the script 
           (batch) files.

uwpflow/source/->  UWPFLOW C (*.c) and C++ code files (*.cpp), headers (*.h), 
                   and other necessary files, including make and project files.

uwpflow/examples/-> Input data for 2 examples, a 173 bus ac/dc system (WSCC-
                    EPRI format) with a SVC, a TCSC and a STATCOM (the format 
                    is described on the data files and the on-line help), and 
                    the IEEE 300 bus ac test system (IEEE common format).  
                    Additional input files needed for bifurcation studies are 
                    also included in this directory.

                    The script (batch) files "173sys" and "ieee300", which 
                    explain some test cases and the files needed for the 
                    studies, are included in this directory also.  These files 
                    are used to automatically run the given examples (in 
                    WINDOWS XP and UNIX).  The script files are shown below.

uwpflow/examples/ -> The following files are needed for Cygwin AWK: 
                     awk.exe, cyggcc_s-1.dll, cygiconv-2.dll, cygintl-8.dll, 
                     cygwin1.dll. These are needed to run the examples 
                     (not really needed for the program). 

________________________________________________________________________________

INSTALLATION:

To install the program:

* WINDOWS 7 (Also tested on some WINDOWS XP & VISTA
 machines):
  

  1.- Download "UWPflow_setup.msi".

  
  
  2.- Run UWPflow_setup.msi.  Follow the steps in the setup file.
      The setup will allow you to create a start menu and a desktop shortcut.
  
  3.- Open the desktop shortcut "UWPflow".
     
  4.- Now you are ready to run the examples in the EXAMPLES subdirectory.
      The batch files IEEE300.BAT and 173SYS.BAT are provided to demonstrate 
      the program features, and to test the program in your system.  
      To run these files use the Batch/Script File option (F6) under the Execute
      Menu option.

  5.- To uninstall it, please use the software management files found 
      under the Control Panel.

* UNIX: 

  1.- Copy the compressed "tar.Z" file "uwpflow.tar.Z".

  2.- Uncompress the file by typing 2 commands:    uncompress uwpflow.tar
                                                   tar xf uwpflow.tar
      This creates a "uwpflow" directory with the corresponding subdirectories 
      "source" and "examples", and other UNIX and AWK utilities. 
     
  3.- Compile UWPFLOW and other routines by running "make all" in the 
      subdirectory "source".  If you have problems, check the "makefile" file 
      to see whether the compiler name and options are appropriate for your 
      system.
 
  4.- Make sure that the uwpflow directory is added to the PATH so that the 
      files needed for running the program can be accessed.

  5.- Run the script files "ieee300" and "173sys" in the "examples" 
      subdirectory to see whether the program operates properly in your system, 
      and to familiarize yourself with the program. 
      To run these files type the file name preceded by the "source" command, 
      e.g., "source ieee300".  
  
Study and run the script (batch) files to familiarize yourself with the 
program. 
________________________________________________________________________________

REFERENCES:

[1]  "Extended Transient-Midterm Stability Package: User's Manual for the Power
     Flow Program," EPRI computer code manual EL-2002-CCM, January 1987.

[2]  "Common Format for Exchange of Solved Load Flow Data," IEEE Trans. Power
     Apparatus and Systems, Vol. 92, No. 6, Nov./Dec. 1973, pp. 1916-1925. 
     Working Group report.

[3]  "Methodology for the Integration of HVDC Links in Large AC Systems-Phase 2:
     Advanced Concepts," Vol. 1, EPRI technical report EL-4365, April 1987.

[4]  "Small Signal Stability Analysis Program Package," Version 2, EPRI user
     manual EL-6678, January 1990.

[5]  C. A. Canizares, M. Pozzi, S. Corsi, and E. Uzunovic, "STATCOM Modeling 
     for Voltage and Angle Stability Studies," International Journal of
     Electrical Power & Energy Systems, Vol. 25, No. 6, June 2003, pp. 431-441.

[6]  C. A. Canizares and F. L. Alvarado, "Point of Collapse and Continuation
     Methods for Large AC/DC Systems," IEEE Trans. Power Systems, Vol. 8, No. 1,
     February 1993, pp. 1-8.

[7]  C. A. Canizares, F. L. Alvarado, C. L. DeMarco, I. Dobson, W. F. Long,
     "Point of Collapse Methods Applied to AC/DC Power Systems," IEEE Trans.
     Power Systems, Vol. 7, No. 2, May 1992, pp. 673-683.

[8]  C. A. Canizares, "On Bifurcations, Voltage Collapse and Load Modeling,"
     IEEE Trans. Power Systems, Vol. 10, No. 1, February 1995, pp. 512-522.

[9]  C. A. Canizares, A. Z. de Souza and  V. H. Quintana, "Improving
     Continuation Methods for Tracing Bifurcation Diagrams in Power Systems,"
     Bulk Power System Voltage Phenomena-III Seminar, ECC Inc., Davos,
     Switzerland, August 1994.

[10] A. Z. de Souza, C. A. Canizares and V. H. Quintana, "New Techniques to
     Speed Up Voltage Collapse Computations Using Tangent Vectors," IEEE 
     Trans. Power Systems, Vol. 12, No. 3, August 1997, pp. 1380-1387.

[11] C. A. Canizares and Z. Faur, "Analysis of SVC and TCSC Controllers in 
     Voltage Collapse," IEEE Trans. Power Systems, Vol. 14, No. 1, February 
     1999, pp. 158-165.

[12] C. A. Canizares, editor, "Voltage Stability Assessment: Concepts, Practices
     and Tools," IEEE-PES Power Systems Stability Subcommittee Special Publication, 
     SP101PSS, August 2002. 


