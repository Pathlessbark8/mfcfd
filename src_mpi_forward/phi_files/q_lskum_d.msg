1 Command: Command-line option -vars is deprecated: please use -head instead
2 Command: Command-line option -outvars is deprecated: please use -head instead
3 Command: Procedure q_lskum understood as q_lskum_mod.q_lskum
4 q_lskum_mod: (DD15) Module initial_conditions_mod is not defined
5 q_lskum_mod: (DD15) Module ieee_arithmetic is not defined
6 q_lskum: (TC20) Variable rank is not declared
7 q_lskum: (TC33) External procedure initial_conditions is not declared
8 q_lskum: (TC33) External procedure ieee_is_nan is not declared
9 q_lskum: (DF03) variable rank is used before initialized
10 q_lskum: (DF03) variable rank is used before initialized
11 q_lskum: (DF03) variable rank is used before initialized
12 q_lskum: (DF03) variable rank is used before initialized
13 q_lskum: (DF03) variable rank is used before initialized
14 q_lskum: (AD08) Head routine q_lskum has remanent variables lmin that one may want to declare as (in)dependents in case of repeated calls
15 q_lskum: (AD03) Varied variable residue written by I-O to file 301
16 q_lskum: (AD03) Varied variable dtg written by I-O to file 301
17 q_lskum: (AD03) Varied variable t written by I-O to file 301
18 q_lskum: (AD14) Differentiated procedure of q_lskum needs to initialize derivative of undeclared side-effect variable lmin
19 (DF05) External routine ieee_is_nan, called by (q_lskum), is assumed default In:.1. and Out:.11
20 (DF05) External routine initial_conditions, called by (q_lskum), is assumed default In:. and Out:.
21 fpi_solver_mod: (DD15) Module post_processing_mod is not defined
22 fpi_solver: (TC33) External procedure update_begin_dq_ghost is not declared
23 fpi_solver: (TC33) External procedure update_begin_qm_ghost is not declared
24 fpi_solver: (TC33) External procedure update_end_dq_ghost is not declared
25 fpi_solver: (TC33) External procedure update_end_qm_ghost is not declared
26 fpi_solver: (TC33) External procedure update_begin_ddq_ghost is not declared
27 fpi_solver: (TC33) External procedure update_end_ddq_ghost is not declared
28 fpi_solver: (TC33) External procedure update_begin_prim_ghost is not declared
29 fpi_solver: (TC33) External procedure update_end_prim_ghost is not declared
30 fpi_solver: (TC20) Variable mpi_double is not declared
31 fpi_solver: (TC20) Variable mpi_sum is not declared
32 fpi_solver: (TC20) Variable petsc_comm_world is not declared
33 fpi_solver: (TC20) Variable ierr is not declared
34 fpi_solver: (TC20) Variable rank is not declared
35 fpi_solver: (TC33) External procedure print_primal_output is not declared
36 fpi_solver: (DF03) variable rank is used before initialized
37 (DF05) External routine print_primal_output, called by (fpi_solver), is assumed default In:. and Out:.
38 (DF05) External routine update_end_prim_ghost, called by (fpi_solver), is assumed default In:. and Out:.
39 (DF05) External routine update_begin_prim_ghost, called by (fpi_solver), is assumed default In:. and Out:.
40 (DF05) External routine update_end_ddq_ghost, called by (fpi_solver), is assumed default In:. and Out:.
41 (DF05) External routine update_begin_ddq_ghost, called by (fpi_solver), is assumed default In:. and Out:.
42 (DF05) External routine update_end_qm_ghost, called by (fpi_solver), is assumed default In:. and Out:.
43 (DF05) External routine update_end_dq_ghost, called by (fpi_solver), is assumed default In:. and Out:.
44 (DF05) External routine update_begin_qm_ghost, called by (fpi_solver), is assumed default In:. and Out:.
45 (DF05) External routine update_begin_dq_ghost, called by (fpi_solver), is assumed default In:. and Out:.
46 state_update_mod: (DD15) Module petsc_data_structure_mod is not defined
47 func_delta: (TC20) Variable mpi_double is not declared
48 func_delta: (TC20) Variable mpi_min is not declared
49 func_delta: (TC20) Variable petsc_comm_world is not declared
50 func_delta: (TC20) Variable ierr is not declared
51 conserved_vector_ubar: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
52 conserved_vector_ubar: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
53 generate_connectivity_mod: (DD15) Module petsc_data_structure_mod is not defined
54 get_interior_neighbours: (TC20) Variable rank is not declared
55 get_interior_neighbours: (DF03) variable rank is used before initialized
56 get_interior_neighbours: (DF03) variable rank is used before initialized
57 get_interior_neighbours: (DF03) variable rank is used before initialized
58 get_interior_neighbours: (DF03) variable rank is used before initialized
59 get_wall_boundary_neighbours: (TC20) Variable rank is not declared
60 get_wall_boundary_neighbours: (DF03) variable rank is used before initialized
61 get_wall_boundary_neighbours: (DF03) variable rank is used before initialized
62 get_wall_boundary_neighbours: (DF03) variable rank is used before initialized
63 get_outer_boundary_neighbours: (TC20) Variable rank is not declared
64 get_outer_boundary_neighbours: (DF03) variable rank is used before initialized
65 get_outer_boundary_neighbours: (DF03) variable rank is used before initialized
66 get_outer_boundary_neighbours: (DF03) variable rank is used before initialized
67 flux_gxp: (TC33) External procedure derf is not declared
68 flux_gxp: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
69 flux_gxn: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
70 flux_gyp: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
71 flux_gyn: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
72 flux_gxp_d: (TC19) Equality test on reals is not reliable
73 flux_gxp_d: (TC19) Equality test on reals is not reliable
74 flux_gxn_d: (TC19) Equality test on reals is not reliable
75 flux_gxn_d: (TC19) Equality test on reals is not reliable
76 flux_gyp_d: (TC19) Equality test on reals is not reliable
77 flux_gyp_d: (TC19) Equality test on reals is not reliable
78 flux_gyn_d: (TC19) Equality test on reals is not reliable
79 flux_gyn_d: (TC19) Equality test on reals is not reliable
80 flux_quad_gxi: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
81 flux_quad_gxi: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
82 flux_quad_gxii: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
83 flux_quad_gxii: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
84 flux_quad_gxiii: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
85 flux_quad_gxiii: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
86 flux_quad_gxiv: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
87 flux_quad_gxiv: (AD09) Please provide a differential of function derf for arguments Arg1=(in) Result=(out)
88 flux_quad_gxi_d: (TC19) Equality test on reals is not reliable
89 flux_quad_gxi_d: (TC19) Equality test on reals is not reliable
90 flux_quad_gxi_d: (TC19) Equality test on reals is not reliable
91 flux_quad_gxi_d: (TC19) Equality test on reals is not reliable
92 flux_quad_gxii_d: (TC19) Equality test on reals is not reliable
93 flux_quad_gxii_d: (TC19) Equality test on reals is not reliable
94 flux_quad_gxii_d: (TC19) Equality test on reals is not reliable
95 flux_quad_gxii_d: (TC19) Equality test on reals is not reliable
96 flux_quad_gxiii_d: (TC19) Equality test on reals is not reliable
97 flux_quad_gxiii_d: (TC19) Equality test on reals is not reliable
98 flux_quad_gxiii_d: (TC19) Equality test on reals is not reliable
99 flux_quad_gxiii_d: (TC19) Equality test on reals is not reliable
100 flux_quad_gxiv_d: (TC19) Equality test on reals is not reliable
101 flux_quad_gxiv_d: (TC19) Equality test on reals is not reliable
102 flux_quad_gxiv_d: (TC19) Equality test on reals is not reliable
103 flux_quad_gxiv_d: (TC19) Equality test on reals is not reliable
104 (DF05) External routine derf, called by (flux_gxp flux_gxn flux_gyp flux_gyn flux_gxp_d flux_gxn_d flux_gyp_d flux_gyn_d conserved_vector_ubar flux_quad_gxi flux_quad_gxii flux_quad_gxiii flux_quad_gxiv flux_quad_gxi_d flux_quad_gxii_d flux_quad_gxiii_d flux_quad_gxiv_d), is assumed default In:.1. and Out:.11
105 compute_force_coeffs_mod: (DD15) Module petsc_data_structure_mod is not defined
106 compute_cl_cd_cm: (TC20) Variable proc is not declared
107 compute_cl_cd_cm: (TC24) Wrong number of dimensions for array itos: 1 expected, and 2 given
108 compute_cl_cd_cm: (TC20) Variable rank is not declared
109 compute_cl_cd_cm: (TC20) Variable mpi_double is not declared
110 compute_cl_cd_cm: (TC20) Variable mpi_sum is not declared
111 compute_cl_cd_cm: (TC20) Variable petsc_comm_world is not declared
112 compute_cl_cd_cm: (TC20) Variable ierr is not declared
113 compute_entropy_mod: (DD15) Module petsc_data_structure_mod is not defined
114 compute_entropy: (TC20) Variable mpi_double is not declared
115 compute_entropy: (TC20) Variable mpi_sum is not declared
116 compute_entropy: (TC20) Variable petsc_comm_world is not declared
117 compute_entropy: (TC20) Variable ierr is not declared
118 stagnation_values_mod: (DD15) Module petsc_data_structure_mod is not defined
119 data_structure_mod: (DD28) Save statement may save extra variables, please use "save::varname" 
120 objective_function_j: (TC20) Variable mpi_double is not declared
121 objective_function_j: (TC20) Variable mpi_sum is not declared
122 objective_function_j: (TC20) Variable petsc_comm_world is not declared
123 objective_function_j: (TC20) Variable ierr is not declared
124 objective_function_j: (TC20) Variable rank is not declared
125 objective_function_j: (DF03) variable rank is used before initialized
126 get_interior_neighbours: Tool: (TypeSpec tree regeneration) null wrapped type
127 get_wall_boundary_neighbours: Tool: (TypeSpec tree regeneration) null wrapped type
128 get_outer_boundary_neighbours: Tool: (TypeSpec tree regeneration) null wrapped type
129 compute_cl_cd_cm: Tool: (TypeSpec tree regeneration) null wrapped type
