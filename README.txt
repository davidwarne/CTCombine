# Summary
#---------

  Combining files in the EGSPhant format.

  Currently, version 0.1

  Features:
   Will read two (2) egsphant files, combine them and write out.  No smarts implemented.
  

# Author
#--------

  Mark Dwyer  ( m2.dwyer@qut.edu.au )

# Institution
#-------------

 High Performance Computing and Research Support
 Queensland University of Technology

# Licenses
#----------

 Feel free to copy/modify all or any portion of this as so desired.

 No guarantee is made as to the accuracy of this code.  Probably best
 not to use in a nuclear reactor.

# Package structure :
#--------------------

  The directory Version_0.1/ is organized as follows :

  - EGSPhant.h                 : The single (header) file of the file format storage.
  - Main.cpp		       : My dodgy kick off file containing the main()
  - Makefile  		       : How I like to compile things
  - README.txt                 : You be readin' it

  - documentation/ : Ha Ha Ha haaa (wipes tear from eye) Good one!
            

# Getting started
#-----------------

  At the command line, just type "make".  This will work on *nix systems.

  To run, "./CTCombine ./Data/Second/cthead.egsphant ./Data/Second/EPID_for_cthead.egsphant output.egsphant"
   for example.
------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------
