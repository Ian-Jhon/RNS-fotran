# RNS
Rotating and Non-rotating Neutron Star models

How to use:

1. First you need to go to the folder EoS_prep, to prepare the eos
file. Once there you need to compile it: cc -O HnG.c -lm -o HnG.x
 
For that you are going to need to have a c compiler installed on the
machine. There will be a couple of warnings but an executable file
named HnG.x should be created.
 
 
2. Now you need to run the code with your EoS. make sure your EoS is
in the same directory as HnG.x. This file needs to have two columns:
energy density and pressure (in that order) IMPORTANT: THE EOS MUST BE
IN CGS!! (E.D. IN g/cm^3 AND P IN dyn/cm^2)
 
 
3. To run the code execute: ./HnG.x YOUR_EOS_FILE 105 > EOS.DAT 105
represents the number of lines in your file, in the example above that
is 105 lines. Note that the final eos is piped into EOS.DAT, which is
the file you will use to run the rotating program.
 
 
4. Now that you have EOS.DAT, you can run the rotating star
program. Go back to the directory of of rns.F, and compile it. To
compile you should type make
 
 
5. Once compiled a file named rns.x should be created. Make sure you
copy EOS.DAT to the directory of rns.x
 
 
6. You need to adjust the parameters.dat file. The file has a simple
explanation. For now you should only need option 1. Below is an
example.  The example below will run a rotating star whose ratio
between polar radius and eq. radius is 0.8, and whose central density
is 350.0 MeV/fm^3.
 
! INPUT DESIRED OPTION: 1 -> SINGLE STAR, NEEDS CENTRAL DENSITY AND RATIO r_p/r_e,
!                            ALSO NEEDS OPTION 1 OR 0:
!                            1->DETAILED OUTPUT IN "OUTPUT_STRUCT.DAT",
!			     0 -> SIMPLE OUTPUT ON TERMINAL
!                       2 -> KEPLER FREQUENCY SEQ., NEEDS CENTRAL DENSITY RANGE
!                              AND NUMBER OF STARS DESIRED
!                       3 -> CONSTANT REST MASS SEQ., NEEDS DESIRED REST MASS, 
!                            CENTRAL DENSITY RANGE AND NUMBERS OF STARS, ALSO NEEDS
!                            OPTION 1 OR 0:1->DETAILED OUTPUT IN "OUTPUT_STRUCT.DAT",
!			     0 -> SIMPLE OUTPUT ON TERMINAL
!                       4 -> SPHERICAL SEQUENCE, NEEDS DENSITY RANGE
!                       5 -> SINGLE KEPLER FREQUENCY STAR. ONLY NEEDS CENTRAL DENSITY.
!                             PRINTING OPTIONS SAME AS IN 1)
!                       6 -> CONSTANT CENTRAL DENSITY AND VARIABLE FREQUENCY. NEEDS
!                            VALUE FOR CENTRAL DENSITY,NUMBER OF STARS AND PRINTING
!			     OPTION (PRINTING OPTIONS SAME AS IN 1)
!                       7 -> Spin-up models
!                       8 -> !CONSTANT CENTRAL DENSITY AND VARIABLE FREQUENCY OPTION,
!                            NEEDS CENTRAL DENSITYcRANGE (IT WILL RUN A SEQUENCE OF
!                            VAR FREQUENCY FOR EACH DENSITY), AND NUMBER OF MODELS.
!                       9 -> !CONSTANT RATIO r_p/r_e. NEEDS CENTRAL DENSITY RANGE,
!                            DESIRED RATIO, NUMBER OF STARS AND PRINTING OPTIONS
!                            (SEE 1).
!INPUT OPTION IN NEXT LINE
1
!INPUT REQUIRED PARAMETERS IN NEXT LINE, SEE INSTRUCTIONS
350.0 0.8  1 
 
7. Once run, the code will output to the screen the general properties
of the star. The detailed structure will be output to
OUTPUT_STRUCT.DAT
 
8. The columns in the output file are as follow: Angle theta (from
z-axis), r, Pressure, alpha, gama, rho, omega, velocity_sq
