alias tapenade="~/FORTRAN/tapenade_3.14/bin/tapenade"
rm -rf tangent_files
mkdir tangent_files
clear
echo '*********Running the TAPENADE script*************'

tapenade -forward -fixinterface -O tangent_files/ -head "q_lskum" -vars "point%x point%y" -outvars "vector_cost_func" \
q_lskum.F90 fpi_solver.F90 data_structure.F90 flux_residual.F90 state_update.F90 q_variables.F90 objective_function.F90 \
point_normals.F90 generate_connectivity.F90 split_fluxes.F90 interior_fluxes.F90 wall_fluxes.F90 outer_fluxes.F90 \
quadrant_fluxes.F90 limiters.F90 parameter.F90 compute_force_coeffs.F90 compute_entropy.F90 compute_enstrophy.F90 stagnation_values.F90

rename .f90 .F90 tangent_files/*.f90
