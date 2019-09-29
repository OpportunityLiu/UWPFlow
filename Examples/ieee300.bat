@echo off
echo                    IEEE 300 BUS AC TEST SYSTEM:
echo                             Examples

echo For this tutorial to run, the directory where 'uwpflow' is located must
echo be added to the path. 

echo The tutorial briefly illustrates several possible applications of the 
echo program 'uwpflow', including the options and files needed to run
echo the different test cases.  An on-line description of the program
echo and all its options is available by running:   uwpflow -h | more
echo Thus:
pause 

@echo on
uwpflow -h | more
@echo off

echo In all cases the -I option (input data in IEEE Common Format) are used, 
echo since the input data is in IEEE Common Format.

echo --------------------------------------------------------------------------
echo 1)  Run a standard power flow and write ASCII output in the output
echo     file 'ieee300.pf0' (output redirection is optional):
pause 

@echo on
uwpflow -I ieee300.cf >ieee300.pf0
@echo off

echo --------------------------------------------------------------------------
echo 2)  Run a standard power flow suppressing the ASCII output (-s option),
echo     and write convergence information in the log file 'ieee300.lg0'
echo     (-l option):
pause 

@echo on
uwpflow -I ieee300.cf -s -lieee300.lg0
@echo off

echo --------------------------------------------------------------------------
echo 3)  Run a standard power flow suppressing ASCII output (-s option),
echo     and write the solution in TAPE IEEE Common Format in the file 
echo     'ieee300.tap' (-W option):
pause 

@echo on
uwpflow -I ieee300.cf -s -Wieee300.tap
@echo off

echo --------------------------------------------------------------------------
echo 4)  Run a standard power flow suppressing ASCII output (-s option),
echo     and write Jacobian information (-j option) in the files 'ieee300.jac' 
echo     (Jacobian), 'ieee300.var' (variables, i.e., V's, P's, Q's, etc.),
echo     and 'ieee300.mis' (mismatches, i.e., deltaP's, deltaQ's, etc.):
pause 

@echo on
uwpflow -I ieee300.cf -s -jieee300
@echo off

echo     The Jacobian can be transformed to MATLAB matrix format using
echo     the AWK filter 'tomatlab', i.e., 

@echo on
call tomatlab ieee300.jac ieeejac.m
@echo off

echo      The file 'ieeejac.m' contains the Jacobian matrix in MATLAB format
echo      to be used directly by this program.
pause

echo --------------------------------------------------------------------------
echo 5)  Run a power flow with a distributed slack bus by using the
echo     generation direction information defined in the file 'ieee300.k' 
echo     (-K option), and write the results in the output file 'ieee300.pf1':
pause 

@echo on
uwpflow -I ieee300.cf -Kieee300.k ieee300.pf1
@echo off

echo --------------------------------------------------------------------------
echo 6)  Run a power flow with a load increase of 0.01 p.u. (-L option).  For
echo     this study one needs the generation and load direction information 
echo     in the file 'ieee300.k' (-K option).  The results are written in CARD 
echo     IEEE Common Format in the file 'ieee300.crd', and the option -Y and -y
echo     are used to write the right and left e-vectors associated to the 
echo     smallest Jacobian eigenvalue in the files 'ieee300.v' and 'ieee300.w',
echo     respectively:
pause 

@echo on
uwpflow -I ieee300.cf -L0.01 -Kieee300.k -wieee300.crd -s -Yieee300.v -yieee300.w
@echo off

echo --------------------------------------------------------------------------
echo 7)  Run a continuation power flow (-c option).  To do this the 
echo     generation and load direction information in the file 'ieee300.k' 
echo     is needed (-K option).  An approximation of the right eigenvector
echo     at the bifurcation point is written in the file 'ieee300.rgt'
echo     (-E option).  The solution information is logged into the file 
echo     'ieee300.lg1' (-l option), and the final output is suppressed using
echo     the -s option.  The -N option is used to turn off all system controls.
pause 

@echo on
uwpflow -I ieee300.cf -Kieee300.k -c -N -Eieee300.rgt -s -lieee300.lg1
@echo off

echo     The maximum entries in 'ieee300.rgt' can then be obtained by running 
echo     'maxim', e.g., 

@echo on
maxim ieee300.rgt ieee300_rgt.max
@echo off

echo     For on-line help on 'maxim' type:   maxim -h
pause

echo --------------------------------------------------------------------------
echo 8)  Run a continuation power flow (-c option) similar to case 7.  Hence,
echo     the generation and load direction information in the file 'ieee300.k' 
echo     is needed (-K option).  However, in this example the continuation 
echo     process is stopped after 5 steps (-z option), and the corresponding
echo     tangent vector at that point is written in the file 'ieee300.tg'
echo     (-Z option).  The final output is suppressed using the -s option.  
pause 

@echo on
uwpflow -I ieee300.cf -Kieee300.k -c -s -z5 -Zieee300.tg 
@echo off

echo     The maximum entries in 'ieee300.tg' can then be obtained by running 
echo     'maxim'.

echo --------------------------------------------------------------------------
echo 9)  Run a continuation power flow similar to case 7, but now the -0 option
echo     is used to generate a series of MATLAB files ('sys*.m') that are used 
echo     to compute several voltage stability indices, namely, eigenvalues, 
echo     singular values, test functions, and reduced determinants.  The 
echo     solution information is logged in the same file 'ieee300.lg2':
pause

@echo on
uwpflow -I ieee300.cf -Kieee300.k -c -N -s -lieee300.lg2 -0sys
@echo off

echo     One of the MATLAB files created is 'sys.m', which can be used to 
echo     compute and plot the desired indices together with the generated file 
echo     'inviter.m'.

echo --------------------------------------------------------------------------
echo 10) Run a parameterized continuation power flow (-H option);
echo     once again the direction file 'ieee300.k' is needed (-K option).
echo     In this example the final ASCII output is stored in 'ieee300.pf3', and
echo     the convergence information is stored in 'ieee300.lg3' (-l option).
echo     The -S option stops the program after the maximum loading point
echo     is computed, at about 80% of the maximum value. The file 
echo     'ieee300.vp' is used to define the desired output buses for the 
echo     voltage profiles (-i option).  
pause 

@echo on
uwpflow -I ieee300.cf -Kieee300.k -H -lieee300.lg3 -S0.8 ieee300.pf3 -iieee300.vp
@echo off


echo --------------------------------------------------------------------------
echo 11) Similar to examples 7 and 8, but in this case the nose curve
echo     and Transient Energy Function (TEF) information are stored in 
echo     the file 'ieeevp.m' in MATLAB format (-m and -O options).  All limits 
echo     are turned off (-n option), so that the TEF can be correctly 
echo     processed in MATLAB (read program BUGS in README.TXT).  
echo     The option -k is used to reduce the step size in the continuation 
echo     method, so that more points on the bifurcation diagram can be
echo     computed.

echo     The computation takes a couple of minutes, please wait....

@echo on
uwpflow -I ieee300.cf -Kieee300.k -cieeevp.m -m -O -n -lieee300.lg4 -s -k0.5
@echo off

echo     The program generates the required MATLAB '*.m' files needed in this 
echo     case to run 'ieeevp.m' in MATLAB.

echo --------------------------------------------------------------------------
echo 12) Finally, this example depicts how to use the program to find the 
echo     exact bifurcation point detected above (case 11) using a direct 
echo     method (-C option).  Again, the generation and load direction defined 
echo     in 'ieee300.k' are needed (-K option); system limits are also turned
echo     off using the -n option.  The right eigenvector at the bifurcation 
echo     point is written in the file 'ieee300.lft', and the bifurcation point
echo     is stored in TAPE IEEE Common Format in 'ieee300.poc' (-W option).
echo     For the program to converge, an initial solution close to the
echo     bifurcation point is needed, thus the -L option is used here
echo     to move the system closer to the singularity point; however, observe
echo     the slow convergence as the system is not close enough to the desired
echo     solution point. 
pause 

@echo on
uwpflow -I ieee300.cf -Kieee300.k -Cieee300.lft -s -Wieee300.poc -L0.3 -n
@echo off

echo     The maximum entries in 'ieee300.lft' can then be obtained by 
echo     running 'maxim'.


