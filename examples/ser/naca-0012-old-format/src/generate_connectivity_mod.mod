V33 :0x4 generate_connectivity_mod
25 generate_connectivity.F90 S624 0
01/23/2019  14:26:42
use parameter_mod public 0 indirect
use data_structure_mod public 0 direct
enduse
D 62 24 688 2704 686 7
D 236 20 7
D 238 20 7
D 240 20 7
D 242 20 7
D 244 20 7
D 246 20 7
D 248 20 7
D 250 20 7
D 252 20 7
D 254 20 7
D 256 20 7
D 258 20 7
D 260 20 7
D 262 20 7
D 264 20 7
D 266 20 7
D 268 20 7
D 270 20 7
D 272 20 7
D 274 20 7
D 276 20 7
D 278 20 7
D 280 20 7
D 282 20 7
D 284 20 7
D 286 20 7
D 288 20 7
D 290 20 7
S 624 24 0 0 0 8 1 0 5015 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 generate_connectivity_mod
S 668 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 674 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 678 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 686 25 5 data_structure_mod points
R 688 5 7 data_structure_mod x points
R 689 5 8 data_structure_mod x$sd points
R 690 5 9 data_structure_mod x$p points
R 691 5 10 data_structure_mod x$o points
R 693 5 12 data_structure_mod y points
R 695 5 14 data_structure_mod y$sd points
R 696 5 15 data_structure_mod y$p points
R 697 5 16 data_structure_mod y$o points
R 700 5 19 data_structure_mod local_id points
R 701 5 20 data_structure_mod local_id$sd points
R 702 5 21 data_structure_mod local_id$p points
R 703 5 22 data_structure_mod local_id$o points
R 706 5 25 data_structure_mod global_id points
R 707 5 26 data_structure_mod global_id$sd points
R 708 5 27 data_structure_mod global_id$p points
R 709 5 28 data_structure_mod global_id$o points
R 712 5 31 data_structure_mod left points
R 713 5 32 data_structure_mod left$sd points
R 714 5 33 data_structure_mod left$p points
R 715 5 34 data_structure_mod left$o points
R 717 5 36 data_structure_mod right points
R 719 5 38 data_structure_mod right$sd points
R 720 5 39 data_structure_mod right$p points
R 721 5 40 data_structure_mod right$o points
R 724 5 43 data_structure_mod flag_1 points
R 725 5 44 data_structure_mod flag_1$sd points
R 726 5 45 data_structure_mod flag_1$p points
R 727 5 46 data_structure_mod flag_1$o points
R 730 5 49 data_structure_mod flag_2 points
R 731 5 50 data_structure_mod flag_2$sd points
R 732 5 51 data_structure_mod flag_2$p points
R 733 5 52 data_structure_mod flag_2$o points
R 736 5 55 data_structure_mod nbhs points
R 737 5 56 data_structure_mod nbhs$sd points
R 738 5 57 data_structure_mod nbhs$p points
R 739 5 58 data_structure_mod nbhs$o points
R 743 5 62 data_structure_mod conn points
R 744 5 63 data_structure_mod conn$sd points
R 745 5 64 data_structure_mod conn$p points
R 746 5 65 data_structure_mod conn$o points
R 749 5 68 data_structure_mod nx points
R 750 5 69 data_structure_mod nx$sd points
R 751 5 70 data_structure_mod nx$p points
R 752 5 71 data_structure_mod nx$o points
R 754 5 73 data_structure_mod ny points
R 756 5 75 data_structure_mod ny$sd points
R 757 5 76 data_structure_mod ny$p points
R 758 5 77 data_structure_mod ny$o points
R 762 5 81 data_structure_mod prim points
R 763 5 82 data_structure_mod prim$sd points
R 764 5 83 data_structure_mod prim$p points
R 765 5 84 data_structure_mod prim$o points
R 769 5 88 data_structure_mod flux_res points
R 770 5 89 data_structure_mod flux_res$sd points
R 771 5 90 data_structure_mod flux_res$p points
R 772 5 91 data_structure_mod flux_res$o points
R 776 5 95 data_structure_mod q points
R 777 5 96 data_structure_mod q$sd points
R 778 5 97 data_structure_mod q$p points
R 779 5 98 data_structure_mod q$o points
R 784 5 103 data_structure_mod dq points
R 785 5 104 data_structure_mod dq$sd points
R 786 5 105 data_structure_mod dq$p points
R 787 5 106 data_structure_mod dq$o points
R 790 5 109 data_structure_mod entropy points
R 791 5 110 data_structure_mod entropy$sd points
R 792 5 111 data_structure_mod entropy$p points
R 793 5 112 data_structure_mod entropy$o points
R 795 5 114 data_structure_mod vorticity points
R 797 5 116 data_structure_mod vorticity$sd points
R 798 5 117 data_structure_mod vorticity$p points
R 799 5 118 data_structure_mod vorticity$o points
R 801 5 120 data_structure_mod vorticity_sqr points
R 803 5 122 data_structure_mod vorticity_sqr$sd points
R 804 5 123 data_structure_mod vorticity_sqr$p points
R 805 5 124 data_structure_mod vorticity_sqr$o points
R 808 5 127 data_structure_mod xpos_nbhs points
R 809 5 128 data_structure_mod xpos_nbhs$sd points
R 810 5 129 data_structure_mod xpos_nbhs$p points
R 811 5 130 data_structure_mod xpos_nbhs$o points
R 813 5 132 data_structure_mod xneg_nbhs points
R 815 5 134 data_structure_mod xneg_nbhs$sd points
R 816 5 135 data_structure_mod xneg_nbhs$p points
R 817 5 136 data_structure_mod xneg_nbhs$o points
R 819 5 138 data_structure_mod ypos_nbhs points
R 821 5 140 data_structure_mod ypos_nbhs$sd points
R 822 5 141 data_structure_mod ypos_nbhs$p points
R 823 5 142 data_structure_mod ypos_nbhs$o points
R 825 5 144 data_structure_mod yneg_nbhs points
R 827 5 146 data_structure_mod yneg_nbhs$sd points
R 828 5 147 data_structure_mod yneg_nbhs$p points
R 829 5 148 data_structure_mod yneg_nbhs$o points
R 833 5 152 data_structure_mod xpos_conn points
R 834 5 153 data_structure_mod xpos_conn$sd points
R 835 5 154 data_structure_mod xpos_conn$p points
R 836 5 155 data_structure_mod xpos_conn$o points
R 838 5 157 data_structure_mod xneg_conn points
R 841 5 160 data_structure_mod xneg_conn$sd points
R 842 5 161 data_structure_mod xneg_conn$p points
R 843 5 162 data_structure_mod xneg_conn$o points
R 847 5 166 data_structure_mod ypos_conn points
R 848 5 167 data_structure_mod ypos_conn$sd points
R 849 5 168 data_structure_mod ypos_conn$p points
R 850 5 169 data_structure_mod ypos_conn$o points
R 852 5 171 data_structure_mod yneg_conn points
R 855 5 174 data_structure_mod yneg_conn$sd points
R 856 5 175 data_structure_mod yneg_conn$p points
R 857 5 176 data_structure_mod yneg_conn$o points
R 860 5 179 data_structure_mod delta points
R 861 5 180 data_structure_mod delta$sd points
R 862 5 181 data_structure_mod delta$p points
R 863 5 182 data_structure_mod delta$o points
S 935 23 5 0 0 0 936 624 7994 0 0 A 0 0 0 0 B 0 35 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 generate_connectivity
S 936 14 5 0 0 0 1 935 7994 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 7 0 624 0 0 0 0 generate_connectivity
F 936 0
S 937 23 5 0 0 0 941 624 8016 0 0 A 0 0 0 0 B 0 107 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_interior_neighbours
S 938 1 3 0 0 6 1 937 8040 4 3000 A 0 0 0 0 B 0 107 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 939 1 3 0 0 9 1 937 5960 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nx
S 940 1 3 0 0 9 1 937 5989 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ny
S 941 14 5 0 0 0 1 937 8016 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 7 3 0 0 0 0 0 0 0 0 0 0 0 0 37 0 624 0 0 0 0 get_interior_neighbours
F 941 3 938 939 940
S 942 23 5 0 0 0 946 624 8042 0 0 A 0 0 0 0 B 0 167 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_wall_boundary_neighbours
S 943 1 3 0 0 6 1 942 8040 4 3000 A 0 0 0 0 B 0 167 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 944 1 3 0 0 9 1 942 5960 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nx
S 945 1 3 0 0 9 1 942 5989 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ny
S 946 14 5 0 0 0 1 942 8042 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 11 3 0 0 0 0 0 0 0 0 0 0 0 0 109 0 624 0 0 0 0 get_wall_boundary_neighbours
F 946 3 943 944 945
S 947 23 5 0 0 0 951 624 8071 0 0 A 0 0 0 0 B 0 230 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_outer_boundary_neighbours
S 948 1 3 0 0 6 1 947 8040 4 3000 A 0 0 0 0 B 0 230 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 949 1 3 0 0 9 1 947 5960 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nx
S 950 1 3 0 0 9 1 947 5989 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ny
S 951 14 5 0 0 0 1 947 8071 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 15 3 0 0 0 0 0 0 0 0 0 0 0 0 170 0 624 0 0 0 0 get_outer_boundary_neighbours
F 951 3 948 949 950
A 36 2 0 0 0 6 668 0 0 0 36 0 0 0 0 0 0 0 0 0 0 0
A 114 2 0 0 0 6 674 0 0 0 114 0 0 0 0 0 0 0 0 0 0 0
A 186 2 0 0 0 6 678 0 0 0 186 0 0 0 0 0 0 0 0 0 0 0
Z
T 686 62 0 0 0 0
A 690 7 236 0 1 2 1
A 689 6 0 36 1 2 1
A 696 7 238 0 1 2 1
A 695 6 0 36 1 2 1
A 702 7 240 0 1 2 1
A 701 6 0 36 1 2 1
A 708 7 242 0 1 2 1
A 707 6 0 36 1 2 1
A 714 7 244 0 1 2 1
A 713 6 0 36 1 2 1
A 720 7 246 0 1 2 1
A 719 6 0 36 1 2 1
A 726 7 248 0 1 2 1
A 725 6 0 36 1 2 1
A 732 7 250 0 1 2 1
A 731 6 0 36 1 2 1
A 738 7 252 0 1 2 1
A 737 6 0 36 1 2 1
A 745 7 254 0 1 2 1
A 744 6 0 114 1 2 1
A 751 7 256 0 1 2 1
A 750 6 0 36 1 2 1
A 757 7 258 0 1 2 1
A 756 6 0 36 1 2 1
A 764 7 260 0 1 2 1
A 763 6 0 114 1 2 1
A 771 7 262 0 1 2 1
A 770 6 0 114 1 2 1
A 778 7 264 0 1 2 1
A 777 6 0 114 1 2 1
A 786 7 266 0 1 2 1
A 785 6 0 186 1 2 1
A 792 7 268 0 1 2 1
A 791 6 0 36 1 2 1
A 798 7 270 0 1 2 1
A 797 6 0 36 1 2 1
A 804 7 272 0 1 2 1
A 803 6 0 36 1 2 1
A 810 7 274 0 1 2 1
A 809 6 0 36 1 2 1
A 816 7 276 0 1 2 1
A 815 6 0 36 1 2 1
A 822 7 278 0 1 2 1
A 821 6 0 36 1 2 1
A 828 7 280 0 1 2 1
A 827 6 0 36 1 2 1
A 835 7 282 0 1 2 1
A 834 6 0 114 1 2 1
A 842 7 284 0 1 2 1
A 841 6 0 114 1 2 1
A 849 7 286 0 1 2 1
A 848 6 0 114 1 2 1
A 856 7 288 0 1 2 1
A 855 6 0 114 1 2 1
A 862 7 290 0 1 2 1
A 861 6 0 36 1 2 0
Z
