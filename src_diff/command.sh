#!/bin/bash
alias tapenade="~/tapenade_3.14/bin/tapenade"

echo "differentiating meshfree solver"

#explicit
tapenade -b -fixinterface -O tapfiles/ -head "q_lskum" -outvars "cl cd cm total_entropy total_enstrophy" -vars "point%x point%y" *.F90 explicit/*.F90

rename 's/.f90$/.F90/' tapfiles/*.f90
