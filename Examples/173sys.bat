@echo off
echo                  173 BUS AC/DC/FACTS TEST SYSTEM:
echo                             Examples

echo For this tutorial to run, the directory where 'uwpflow' is located must
echo be added to the path. 

echo This tutorial briefly illustrates several possible applications of the 
echo program 'uwpflow', including the options and files needed to run
echo different cases for an ac/dc/FACTS system in WSCC (EPRI) format.

echo --------------------------------------------------------------------------
echo 1)  Run a standard ac/dc/FACTS power flow and write in ASCII the output
echo     file '173sys.pf0'.  Also create the IEEE Common Format file 
echo     '173sys.cf' with the final solution (-W option), and write the 
echo     Jacobian (-j option) in '173sys.jac' with the corresponding 
echo     mismatch and solution vectors in '173sys.mis' and '173sys.var', 
echo     respectively:
pause 

@echo on
uwpflow 173sys.wsc 173sys.pf0 -j173sys -W173sys.cf
@echo off

echo     Observe that the IEEE Common Format file created, '173sys.cf',
echo     has the HVDC and FACTS data in the same format as the original input
echo     file '173sys.wsc'.  However, the generator maximum power information
echo     is not available in the IEEE format.
pause
echo     The Jacobian can be transformed into MATLAB matrix format using
echo     the AWK filter 'tomatlab', i.e., 

@echo on
call tomatlab 173sys.jac sysjac.m
@echo off

echo     The MATLAB data is stored in this case in the output file 'sysjac.m'.
pause
echo --------------------------------------------------------------------------
echo 2)  Run a power flow with a load increase of 0.03 (-L option).  For
echo     this study one needs the generation and load direction information 
echo     in the file '173sys.k' (-K option).  Notice that this example
echo     illustrates the use of a distributed slack bus to solve the power
echo     flow problem.   The IEEE Common Format file created in the previous
echo     run, '173sys.cf', is used as the input data in this case 
echo     (-I option).  Also, the right and left eigenvectors of the smallest
echo     eigenvalue of the system Jacobian are stored in the files '173sys.v'
echo     (-Y option) and '173sys.w' (-y option), respectively.  
echo     The ASCII output is stored in the file '173sys.pf1'.
pause 

@echo on
uwpflow -I 173sys.cf 173sys.pf1 -L0.03 -K173sys.k -Y173sys.v -y173sys.w
@echo off

echo     The maximum normalized entries in '173sys.v' and '173sys.w' can 
echo     then be obtained by running 'maxim', e.g.,
pause

@echo on
maxim 173sys.w 173sys_w.max
@echo off

echo     For on-line help on 'maxim' type:   maxim -h
echo --------------------------------------------------------------------------
echo 3)  Run a continuation power flow (-c option).  To do this the 
echo     generation and load direction information stored in file 
echo     '173sys.k' is needed (-K option).  The last solution point is written
echo     in ASCII in the output file '173sys.pf2', and the log information is 
echo     stored in the file '173sys.lg1' (-l option).
echo     The tracing of the bifurcation diagram is stopped after the 
echo     bifurcation point at 0.8 of the maximum value of the loading factor 
echo     (-S option):
pause 

@echo on
uwpflow 173sys.wsc -K173sys.k -c 173sys.pf2 -l173sys.lg1 -S0.8 
@echo off

echo --------------------------------------------------------------------------
echo 4)  Run another continuation power flow (-c option), similar to
echo     to the previous case 3, but ignoring maximum active power
echo     generation limits (-X option).  The log information is stored in the 
echo     file '173sys.lg2' (-l option), and the tracing of the bifurcation 
echo     diagram is again stopped after the bifurcation point at 0.8 of the
echo     maximum value of the loading factor (-S option).  Also, in this
echo     case the ASCII output is suppressed (-s option).
pause 

@echo on
uwpflow 173sys.wsc -K173sys.k -c -X -s -l173sys.lg2 -S0.8
@echo off

echo --------------------------------------------------------------------------
echo 5)  Run a parameterized continuation power flow (-H option);
echo     once again the direction file '173sys.k' is needed (-K option).
echo     In this example the final ASCII output is suppressed (-s option), 
echo     and the convergence information is stored in the log file 
echo     '173sys.lg3'  (-l option).  The nose curve information is stored in 
echo     the file 'sysvp.m' in MATLAB format (-m option).  All ac limits are 
echo     turned off (-n option), and the -e option is used to print out some 
echo     HVDC variable information.  Finally, the option -k is used to 
echo     reduce the step size in the continuation method, i.e., obtain more
echo     points in the bifurcation diagram, and the -S option is used to stop
echo     the continuation process after the maximum loading point is reached.

echo     The computation takes a couple of minutes, please wait....

@echo on
uwpflow 173sys.wsc -K173sys.k -Hsysvp.m -m -n -l173sys.lg3  -s -k0.5 -e -S0.9
@echo off

echo --------------------------------------------------------------------------
echo 6)  Similar to examples 3, 4 and 5, but in this case the system load
echo     is modeled as a nonlinear function of the voltage, i.e., 
echo                            Pl = Pn*V^a + Pz*V^2
echo                            Ql = Qn*V^b + Qz*V^2
echo     The -D option is used in this case to define the values of Pn, Qn, 
echo     Pz, Qz, a, and b stored in OH (SMMS) format in '173sys.oh'. In this
echo     case, several voltage stability indices (SF, VSF, TV) are printed out 
echo     (-f option), together with the voltages and angles of the buses 
echo     defined in '173sys.vp' (-i option).
pause 

@echo on
uwpflow 173sys.wsc -K173sys.k -c -l173sys.lg4 -s -D173sys.oh -i173sys.vp -f -S0.8
@echo off

echo --------------------------------------------------------------------------
echo 7)  This example shows the use of the program to obtain a bifurcation
echo     manifold using a detailed generator steady state model, which
echo     considers Ra, Xd, Xq, and stator and rotor limits (in Ia and Eq).  
echo     The generator data is read from the file '173sys.gen' using the 
echo     -3 option.  Once more, the generation and load direction defined in 
echo     '173sys.k' are needed (-K option), and the convergence information 
echo     is stored in the file '173sys.lg5' using the -l option.
echo     The options -i and -e are used to obtain the desired output.
echo     This example also shows some convergence difficulties due to the 
echo     lack of adequate initial guesses for the internal generator variables
echo     (this problem can be significantly reduced by using an IEEE Common 
echo     Format input file; however, be aware that PgMax information must be
echo     defined then in a 'K' file).
pause 

@echo on
uwpflow 173sys.wsc -K173sys.k -3173sys.gen -c -i173sys.vp -e -s -l173sys.lg5
@echo off

echo --------------------------------------------------------------------------
echo 8) Finally, this example depicts how to use the program to try to 
echo    find the exact bifurcation point using a direct method (-C option), 
echo    but it fails due to the fact that the bifurcation point is not 
echo    associated with a singular Jacobian.  Again, the generation and 
echo    load direction defined in '173sys.k' are needed (-K option).  
echo    For the program to converge to the desired solution, an initial
echo    solution close to the bifurcation point is needed, thus the -L 
echo    option is used here to move the system closer to the maximum point
echo    detected in case 7.  Hence, the -3 option is used again to read 
echo    the steady state generator data stored in file '173sys.gen'.  
echo    Notice that the program fails to obtain a solution in this case.
pause 

@echo on
uwpflow 173sys.wsc -K173sys.k -C -s -L0.03 -3173sys.gen
@echo off


