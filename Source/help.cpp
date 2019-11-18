/* AC/DC Power Flow Help. */

#include "param.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ErrorStop(char *Msg);

/* ---------------------- ErrorStop ------------------------------------ */
void ErrorStop(char *Msg) {
  printf("UW Continuation Power Flow (c)1992,1996,1999, 2006 C. Canizares, F. "
         "Alvarado and S. Zhang.\n");
  if (!NullName(Msg))
    printf("Error: %s\n", Msg);
  printf("\n");
  printf("Usage:\n");
  printf("       Like any other UNIX program, i.e., command-line options "
         "(-option)\n");
  printf("       with redirection of output (>) from screen into files:\n");
  printf("\n");
  printf("              uwpflow  [-options]  input_file   [[>]output_file]\n");
  printf("\n");
  printf("Input file:\n");
  printf("       The input_file could be in WSCC/BPA/EPRI format (or Anarede's "
         "variation)\n");
  printf("       or IEEE common format (or Electrocon's variation) for the ac "
         "system.\n");
  printf("       Other AC formats may be used (see options below).\n");
  printf("       The dc data could be on either WSCC/BPA or multiterminal EPRI "
         "format,\n");
  printf("       although it just deals with the standard two-terminal HVDC "
         "problem.\n");
  printf("       The FACTS devices format was specifically designed for this "
         "program and\n");
  printf("       is explained in the help file and test systems provided with "
         "the program.\n");
  printf("\n");
  printf("Output files:\n");
  printf("       The program writes the solution into the output_file in "
         "ASCII.  It\n");
  printf("       can also write the solved case in a file in IEEE common "
         "format\n");
  printf("       using the -W or -w option (HVDC links are written in EPRI's "
         "ETMSP format\n");
  printf("       and FACTS devices are written in their own format).\n");
  printf("       Additional files can be created for post-processing analyses "
         "(see\n");
  printf("       options below), such as the bifurcation diagram (nose curve) "
         "in column\n");
  printf("       form for plotting with MATLAB or Octave, Jacobians for SMMS "
         "or\n");
  printf("       MATLAB/Octave studies, etc.\n");
  printf("\n");
  printf("Solution Technique:\n");
  printf("       The power flows are solved with simultaneous N-R, allowing "
         "for\n");
  printf("       asynchronous systems, area interchange, remote voltage "
         "control, and\n");
  printf("       local and remote regulating transformers (LTCs and phase "
         "shifters\n");
  printf("       controlling voltages, angles, and/or active and reactive "
         "power flows).\n");
  printf("\n");
  printf("Options:\n");
  printf("         -a     Turns off tap and angle limits in regulating "
         "transformers.\n");
  printf("\n");
  printf("         -A     Turns off interchange area control.\n");
  printf("\n");
  printf("         -b     Solve base case before changing the loading factor "
         "lambda.\n");
  printf("\n");
  printf("         -Bnum  PQ bus number 'num' where the voltage is fixed in "
         "order to\n");
  printf("                find the loading factor (lambda) for voltage "
         "collapse studies.\n");
  printf("                Must be used with -K and -v options.\n");
  printf("\n");
  printf("         -cfile Increases the loading factor lambda using a "
         "continuation\n");
  printf("                method for finding voltage profiles.\n");
  printf("                The output (optional 'file') is a list of max. 8 ac "
         "voltages\n");
  printf("                that change the most, plus 3 additional variables "
         "for each dc\n");
  printf(
      "                bus and for each generator (see -e and -3 options).\n");
  printf("                Must be used with the -K option.\n");
  printf("\n");
  printf("         -Cfile Direct method studies, i.e., find the max. loading "
         "factor\n");
  printf("                lambda for a given generation and load direction. "
         "The base case\n");
  printf("                loading can be initialized using the -L option; "
         "however, the\n");
  printf("                program calculates an initial loading of the system "
         "before the\n");
  printf("                direct method is applied.  The left e-vector is "
         "written in\n");
  printf(
      "                'file' (optional).  Must be used with the -K option.\n");
  printf("\n");
  printf("         -d     Generates some debug output.\n");
  printf("\n");
  printf("         -Dfile Read load model data from 'file', using Ontario "
         "Hydro (OH)\n");
  printf("                format, which is based on the load model: "
         "Pl=Pn*V^a+Pz*V^2\n");
  printf("                                                          "
         "Ql=Qn*V^b+Qz*V^2\n");
  printf("                If a bus is not defined in the list, Pn=Qn=0 is "
         "assumed.\n");
  printf("\n");
  printf("         -e     Output 3 dc variables per dc link, and 3 internal "
         "variables for\n");
  printf("                generators in the output list during the "
         "continuation process\n");
  printf("                (see the -i option).\n");
  printf("\n");
  printf("         -Efile Print in 'file' the continuation method direction "
         "vector at the\n");
  printf("                maximum loading factor (PoC right e-vector).\n");
  printf("\n");
  printf("         -fnum  Output Sensitivity Factors (SF) and Voltage "
         "Sensitivity Factors\n");
  printf("                (VSF and tangent vector) during continuation "
         "computations.\n");
  printf("                The number 'num' defines the bus for which the "
         "voltage entry\n");
  printf("                and rank in the tangent vector are printed out "
         "(Tangent Vector\n");
  printf("                Index); if this number is not defined, the program "
         "chooses the\n");
  printf("                bus with the maximum initial voltage entry in the "
         "tangent\n");
  printf("                vector.\n");
  printf("\n");
  printf("         -Fval  Stability/sparsity value 'val' for factorization "
         "(def. 0.01).\n");
  printf("                A value of 0 means choose a pivot based on sparsity "
         "only;\n");
  printf("                a value of 1 means choose a pivot based on stability "
         "only.\n");
  printf("\n");
  printf("         -g     Force Q in generators to zero when reading data in "
         "IEEE\n");
  printf(
      "                common format, since sometimes a value of Qg creates\n");
  printf("                convergence problems.\n");
  printf("\n");
  printf("         -G     Turns off recovery from some ac device limits in the "
         "program.\n");
  printf("                For example, the program allows to recover voltage "
         "control\n");
  printf("                after a Q-limit is reached by monitoring the "
         "voltage; this\n");
  printf("                option eliminates that possibility.\n");
  printf("\n");
  printf("         -h     Prints this message in standard output.\n");
  printf("\n");
  printf("         -Hfile Increases the loading factor lambda using a "
         "parameterized\n");
  printf("                continuation method for finding voltage profiles.\n");
  printf("                The output (optional 'file') is a list of max. 8 ac "
         "voltages\n");
  printf("                that change the most, plus 3 additional variables "
         "for each dc\n");
  printf("                bus (see -e option).  Must be used with the -K "
         "option.\n");
  printf("\n");
  printf("         -ifile List of numbers and names in 'file' for printing "
         "variable\n");
  printf("                profiles with the -c and -H options. The input "
         "format is:\n");
  printf("                Number Bus/AreaName [VarType].\n");
  printf("                Use zero when either the number or the name are "
         "unknown.\n");
  printf("                If Name has spaces, wrap it in double or single "
         "quotes.\n");
  printf("                VarType is optional and can be: V for voltage "
         "(default),\n");
  printf("                D for angle, PG for MW generated, QG for Mvar "
         "generated,\n");
  printf("                PL for MW load, QL for Mvar load, or PA for MW area "
         "flow.\n");
  printf("                If Name and Number are both equal to 0 and VarType "
         "is either \n");
  printf("                PL, QL, PG or QG, the program will print the "
         "corresponding \n");
  printf("                total load or generation in MW or Mvar.\n");
  printf("\n");
  printf("         -I     AC input data in IEEE common format.\n");
  printf("\n");
  printf(
      "         -Ip    AC input data in IEEE common format with Power World\n");
  printf("                modifications.\n");
  printf("\n");
  printf("         -jname Write the Jacobian of the solved case in I J VALUE "
         "format in\n");
  printf("                'name.jac'.  The equation mismatches and the system "
         "variables\n");
  printf("                are also written in 'name.mis' and 'name.var', "
         "respectively.\n");
  printf("\n");
  printf("         -Jname With the -B option, it generates the Jacobian, "
         "mismatches and\n");
  printf("                variables corresponding to the system without the "
         "loading\n");
  printf("                factor as a variable.  For PoC studies (-C option) "
         "it generates\n");
  printf("                the nxn system Jacobian, mismatches and variables; "
         "for the\n");
  printf("                complete (2n+1)x(2n+1) PoC Jacobian, use the -j "
         "option.  The\n");
  printf("                corresponding Jacobian, mismatches and variables are "
         "written\n");
  printf("                in I J VALUE format in 'name.jac', 'name.mis' and "
         "'name.var',\n");
  printf("                respectively.\n");
  printf("\n");
  printf("         -kval  Factor 'val' used in the homotopy continuation "
         "method for\n");
  printf("                finding the increments in the loading factor lambda "
         "(def. 1).\n");
  printf("                Must be used with the -c and -H options.\n");
  printf("\n");
  printf("         -Kfile Read generation and load distribution factors from "
         "'file'.\n");
  printf("                All data is assumed p.u. and must be separated by "
         "spaces: \n");
  printf("                BusNumber BusName DPg Pnl Qnl PgMax [Smax Vmax Vmin "
         "Pzl Qzl].\n");
  printf("                If the input variables DPg, Pnl, Qnl or PgMax are "
         "unknown, give\n");
  printf("                them a value of zero; Smax, Vmax, Vmin, Pzl and Qzl "
         "are\n");
  printf("                optional.\n");
  printf("                The generation factors DPg are normalized for each "
         "area, i.e.,\n");
  printf("                ||DPg||=1 per area.\n");
  printf("                The load is represented by:\n");
  printf("                     Pl=(Pn+Pnl*lambda)*V^a+(Pz+Pzl*lambda)*V^2\n");
  printf("                     Ql=(Qn+Qnl*lambda)*V^b+(Qz+Qzl*lambda)*V^2\n");
  printf("                where Pn, Qn, Pz, Qz, a, and b are defined with the "
         "-D option,\n");
  printf("                and lambda corresponds to the loading factor.  If "
         "the -D option\n");
  printf("                is not used, the load model default values are: "
         "a=b=0, Pz=Qz=0.\n");
  printf("                Buses not in the list are assumed to have zero "
         "distribution\n");
  printf("                factors.   If BusName has spaces, wrap it in double "
         "or single\n");
  printf("                quotes.\n");
  printf("\n");
  printf("         -lfile Write standard error output to 'file' (log file).\n");
  printf("\n");
  printf("         -Lval  Loading factor 'val' (def. 0).  Simulates load "
         "changes in\n");
  printf("                conjunction with the load distribution factors (-K "
         "option).\n");
  printf("\n");
  printf("         -m     Output continuation profiles in MATLAB/Octave "
         "format.  If TEF profiles\n");
  printf("                are needed, use the -O option.\n");
  printf("\n");
  printf("         -Mnum  Number 'num' of max. N-R iterations, overriding "
         "input data\n");
  printf("                (default 50).\n");
  printf("\n");
  printf("         -n     Turns off all ac system limits.\n");
  printf("\n");
  printf("         -N     Turns off all ac system controls.\n");
  printf("\n");
  printf("         -otol  The tolerance 'tol' controls the application of "
         "limits during\n");
  printf("                the continuation process.  The smaller this value, "
         "the more\n");
  printf("                steps required of the continuation method (default "
         "0.000001).\n");
  printf("\n");
  printf("         -Onum  This option is used together with -m to output ac/dc "
         "TEF\n");
  printf("                information during the continuation method to "
         "determine\n");
  printf("                the energy profiles with the help of MATLAB or "
         "Octave. The integer\n");
  printf("                'num' corresponds to the number of significant "
         "digits in TEF\n");
  printf("                (default and min. 6; max. 10).\n");
  printf("                If the system has HVDC links, the program defaults "
         "the PI\n");
  printf("                controller gains (Kp and Ki) of all HVDC converters "
         "for dc \n");
  printf("                computations to 'typical' values of Ki=75 and Kp=1.  "
         "These\n");
  printf("                values can later be changed in the MATLAB/Octave "
         "output file.\n");
  printf("                The program also generates two MATLAB/Octave '.m' "
         "files needed for\n");
  printf("                plotting and computation of the ac/dc energy "
         "function, namely,\n");
  printf("                'addtotef.m' (only if dc lines present) and "
         "'tefprof.m'.\n");
  printf("\n");
  printf(
      "         -p     Turns off P and Q limits in regulating transformers.\n");
  printf("\n");
  printf("         -P     Turns off P and Q control by regulating "
         "transformers.\n");
  printf("\n");
  printf("         -q     Turns off Q limits in PV buses.\n");
  printf("\n");
  printf("         -qx    Turns off V limits in reactance-controlled (BX) "
         "buses.\n");
  printf("                These buses are defined in the ETMSP/EPRI/BPA input "
         "data file.\n");
  printf("\n");
  printf("         -qz    Turns off Q limits in reactance-controlled (BZ) "
         "buses.\n");
  printf("                These buses are defined in the Italian ADD input "
         "data file.\n");
  printf("\n");
  printf("         -Q     Turns off remote voltage generator control.  "
         "Generators\n");
  printf("                will control their terminal voltages at their given "
         "values.\n");
  printf("\n");
  printf("         -QX    Turns off remote voltage control in "
         "reactance-controlled (BX)\n");
  printf("                buses.  The local bus voltage will be used for the "
         "control.\n");
  printf("\n");
  printf("         -r     Turns off V limits in regulating transformers and PV "
         "buses.\n");
  printf("\n");
  printf("         -R     Turns off V control by regulating transformers.\n");
  printf("\n");
  printf("         -s     Suppress ASCII output_file.\n");
  printf("\n");
  printf("         -Sval  Stop value 'val' for the loading factor lambda in "
         "the continua-\n");
  printf("                tion method (-c and -H options), in p.u. of the "
         "maximum lambda.\n");
  printf("                The default is 0 for a complete trace of the "
         "bifurcation\n");
  printf("                diagram (min. 0; max. 1, i.e., lambda maximum).\n");
  printf("\n");
  printf("         -ttol  If the relative error of two consecutive iteration "
         "mismatches\n");
  printf("                is larger than 'tol', voltage limits and regulating "
         "transformer\n");
  printf("                limits are applied (default 0.1).\n");
  printf("\n");
  printf(
      "         -Ttol  P.U. tolerance 'tol' for N-R method (default 1e-4).\n");
  printf("\n");
  printf("         -uval  Value 'val' of the tolerance used to reduce the "
         "number of equa-\n");
  printf("                tions in the continuation method (-c and -H "
         "options), based on\n");
  printf("                the tangent vector technique.\n");
  printf("                The default is 1e-3 (min. 0 or no reduction, max. "
         "0.2).\n");
  printf("                WARNING: This option might produce cycling, "
         "back-tracking, or\n");
  printf("                         singularity problems. If this happens, "
         "increase the\n");
  printf("                         number of steps in the -U option and/or "
         "reduce the\n");
  printf("                         value of 'val' in -u. \n");
  printf("\n");
  printf("         -Unum  Number 'num' of steps of the continuation method (-c "
         "and -H\n");
  printf("                options) between system reductions.  Used with the "
         "-u option.\n");
  printf("                The default number is 10 (min. 2, max. 100).\n");
  printf("\n");
  printf("         -vmag  Voltage magnitude 'mag' at the first PQ bus (unless "
         "otherwise\n");
  printf("                specified by -B option) to find the corresponding "
         "lambda\n");
  printf("                for voltage collapse studies.  Must be used with -K "
         "option.\n");
  printf("\n");
  printf("         -Vfile Read initial guesses for ac and dc variables from "
         "'file'.\n");
  printf("                The data must be separated by spaces, i.e.,\n");
  printf("                BusNumber BusName V_mag V_ang(deg), for ac buses, "
         "and\n");
  printf("                BusNumber BusName Variables Values, for dc buses.\n");
  printf("                For defining the dc variable, use the same EPRI "
         "format used\n");
  printf("                to define the control modes in the dc lines (e.g., "
         "ALGAPA,\n");
  printf("                would represent the dc variables Alpha, Gamma and "
         "Power);\n");
  printf("                the values must be given in standard units, i.e., "
         "MW, MVAR, KV,\n");
  printf("                Amp, deg.  If the input variables are unknown give "
         "them a value\n");
  printf("                of zero.  If BusName has spaces, wrap it in double "
         "or single\n");
  printf("                quotes.\n");
  printf("                If no 'file' name is given, the system is given a "
         "'flat' start.\n");
  printf("\n");
  printf("         -wfile Write solved case into 'file' using IEEE CARD common "
         "format.\n");
  printf("\n");
  printf("         -Wfile Similar to -w option, but the solved case is written "
         "in 'file'\n");
  printf("                using IEEE TAPE common format.\n");
  printf("\n");
  printf("         -x     Do not use the distributed slack bus concept during "
         "the solution\n");
  printf("                process, i.e., use only one slack bus.  The "
         "generator powers\n");
  printf("                change based only on the direction defined in the K "
         "file and\n");
  printf("                the load level, as defined by the loading factor "
         "lambda.\n");
  printf("\n");
  printf("         -X     Do not enforce maximum active generation limits "
         "(PgMax).\n");
  printf("                See -K option.\n");
  printf("\n");
  printf("         -yfile Print in 'file' an approximation of the left "
         "e-vector of the\n");
  printf("                smallest real |e-value| at the current operating "
         "point.\n");
  printf("                Works with the -c/-H and -z options as well.\n");
  printf("\n");
  printf("         -Yfile Print in 'file' an approximation of the right "
         "e-vector of the\n");
  printf("                smallest real |e-value| at the current operating "
         "point.\n");
  printf("                Works with the -c/-H and -z options as well.\n");
  printf("\n");
  printf("         -znum  Stop continuation method after 'num' steps.  "
         "Together with the\n");
  printf("                -Z option, it can be used to print a tangent vector "
         "for a\n");
  printf("                particular value of the loading parameter lambda.\n");
  printf("\n");
  printf(
      "         -Zfile Print in 'file' the normalized tangent vector to the\n");
  printf("                bifurcation manifold at step 'num' of the "
         "continuation method\n");
  printf("                (as defined by -z).  It is used with the -c or -H "
         "options.\n");
  printf("\n");
  printf("         -0name Print in MATLAB/Octave format the Jacobian matrices "
         "needed to compute\n");
  printf("                several voltage stability indices. These matrices "
         "are printed\n");
  printf("                in the files 'name#.m', where # stands for the step "
         "number in\n");
  printf("                the continuation method (-c or -H options).\n");
  printf("                A file 'name.m' is also created with all the "
         "MATLAB/Octave\n");
  printf(
      "                instructions to compute and plot 6 distinct indices:\n");
  printf("                   * Minimum |e-value| and singular value for full "
         "matrix.\n");
  printf("                   * Minimum |e-value| and singular value for "
         "reduced Q matrix.\n");
  printf("                   * Test function and reduced determinant for a "
         "given bus\n");
  printf("                     (-1 option).\n");
  printf("                The program generates the MATLAB/Octave 'inviter.m' "
         "file needed for\n");
  printf("                the computation of these indices.\n");
  printf("                WARNING: This option generates a lot of output files "
         "that might\n");
  printf("                         clutter your system; it must be used with "
         "caution.\n");
  printf("\n");
  printf("         -1num  Used with -0 option, and defines the load bus at "
         "which the test\n");
  printf("                functions are computed.  If this option is not used, "
         "or if the\n");
  printf("                bus 'num' does not correspond to a load bus, the "
         "program\n");
  printf("                chooses the bus with the maximum voltage entry in "
         "the initial\n");
  printf("                tangent vector.\n");
  printf("\n");
  printf("         -2num  Define the number 'num' of steps used to determine "
         "the\n");
  printf("                change of direction of the loading parameter lambda "
         "in\n");
  printf("                the continuation method (-c or -H options) due to "
         "voltages\n");
  printf("                increasing after Q-limits are encountered (default "
         "5).\n");
  printf("\n");
  printf("         -3file Read generator steady-state data from 'file' using "
         "free format:\n");
  printf("                BusNumber BusName Ra Xd Xq Ia_max Eq_max Eq_min\n");
  printf("                For BusNames with spaces, wrap the word in double or "
         "single\n");
  printf("                quotes.  The program uses either the BusNumber or "
         "the BusName\n");
  printf("                to identify the bus; if one of this is not known, "
         "give it a 0\n");
  printf("                value.  For round-rotor machines, make Xd=Xq, or "
         "define Xq=0.\n");
  printf("                The program assumes the following default values: "
         "Ra=0, Xq=<Xd,\n");
  printf("                Ia_max=large, Eq_max=large, Eq_min=0.\n");
  printf("\n");
  printf("         -4     Turns off Eq limits in all generators.  See -3 "
         "option.\n");
  printf("\n");
  printf("         -5     Turns off Ia limits in all generators.  See -3 "
         "option.\n");
  printf("\n");
  printf("         -6file AC input data in ITALIAN format. If 'file' name is "
         "given,\n");
  printf("                a COLAS ADD file is read, which defines: new bus kV "
         "levels; \n");
  printf("                min. and max. bus voltages; load voltage "
         "coefficients 'a' and \n");
  printf("                and 'b', i.e.,\n");
  printf("                     Pl=(Pn+Pnl*lambda)*V^a \n");
  printf("                     Ql=(Qn+Qnl*lambda)*V^b \n");
  printf("                pilot nodes and generators for secondary voltage "
         "control \n");
  printf("                (-# option); and generator and load directions. \n");
  printf("                This file may be used instead of the -K and -D "
         "files; however,\n");
  printf("                these files take precedence over the ADD file in "
         "defining \n");
  printf("                similar variables for collapse studies.\n");
  printf("\n");
  printf("         -7     Enforce Vmax and Vmin on system buses during the "
         "continuation\n");
  printf("                process (-c or -H options).\n");
  printf("\n");
  printf("         -8     Enforce Imax limits on transmission elements during "
         "the \n");
  printf("                continuation process (-c or -H options).\n");
  printf("\n");
  printf("         -9     Do not enforce maximum power generation limits "
         "(Smax).\n");
  printf("                See -K option.\n");
  printf("\n");
  printf("         -$val  Define the base power value 'val', overriding the "
         "value given\n");
  printf("                in the input data file.\n");
  printf("\n");
  printf("         -#     Use secondary voltage control as defined by ENEL, "
         "i.e., remote\n");
  printf("                voltage control of pilot buses by generators with "
         "participation\n");
  printf("                factors defined as: \n");
  printf("                  * over-excited  -> q_i=Qmax_i/Sum Qmax of pilot "
         "bus gens.\n");
  printf("                  * under-excited -> q_i=Qmin_i/Sum Qmin of pilot "
         "bus gens.\n");
  printf("\n");
  if (!NullName(Msg))
    exit(ERROREXIT);
}
