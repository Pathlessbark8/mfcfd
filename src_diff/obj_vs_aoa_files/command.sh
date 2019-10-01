#!/bin/bash
alias tapenade="~/tapenade_3.14/bin/tapenade"

echo "differentiating meshfree solver"

#angle of attack
tapenade -b -fixinterface -O tapfiles/ -head "q_lskum_d" -outvars "cld cdd cmd" -vars "point%x point%y" *.F90 

rename 's/.f90$/.F90/' tapfiles/*.f90
