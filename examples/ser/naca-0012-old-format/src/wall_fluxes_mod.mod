V33 :0x4 wall_fluxes_mod
15 wall_fluxes.F90 S624 0
01/23/2019  14:26:42
use quadrant_fluxes_mod public 0 direct
use split_fluxes_mod public 0 direct
use q_variables_mod public 0 direct
use parameter_mod public 0 indirect
use data_structure_mod public 0 direct
use limiters_mod public 0 direct
enduse
D 62 24 692 2704 690 7
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
D 506 21 9 1 3 35 0 0 0 0 0
 0 35 3 3 35 35
D 509 21 9 1 3 35 0 0 0 0 0
 0 35 3 3 35 35
D 512 21 9 1 3 35 0 0 0 0 0
 0 35 3 3 35 35
S 624 24 0 0 0 8 1 0 5015 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 wall_fluxes_mod
S 635 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 672 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 678 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 682 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 690 25 5 data_structure_mod points
R 692 5 7 data_structure_mod x points
R 693 5 8 data_structure_mod x$sd points
R 694 5 9 data_structure_mod x$p points
R 695 5 10 data_structure_mod x$o points
R 697 5 12 data_structure_mod y points
R 699 5 14 data_structure_mod y$sd points
R 700 5 15 data_structure_mod y$p points
R 701 5 16 data_structure_mod y$o points
R 704 5 19 data_structure_mod local_id points
R 705 5 20 data_structure_mod local_id$sd points
R 706 5 21 data_structure_mod local_id$p points
R 707 5 22 data_structure_mod local_id$o points
R 710 5 25 data_structure_mod global_id points
R 711 5 26 data_structure_mod global_id$sd points
R 712 5 27 data_structure_mod global_id$p points
R 713 5 28 data_structure_mod global_id$o points
R 716 5 31 data_structure_mod left points
R 717 5 32 data_structure_mod left$sd points
R 718 5 33 data_structure_mod left$p points
R 719 5 34 data_structure_mod left$o points
R 721 5 36 data_structure_mod right points
R 723 5 38 data_structure_mod right$sd points
R 724 5 39 data_structure_mod right$p points
R 725 5 40 data_structure_mod right$o points
R 728 5 43 data_structure_mod flag_1 points
R 729 5 44 data_structure_mod flag_1$sd points
R 730 5 45 data_structure_mod flag_1$p points
R 731 5 46 data_structure_mod flag_1$o points
R 734 5 49 data_structure_mod flag_2 points
R 735 5 50 data_structure_mod flag_2$sd points
R 736 5 51 data_structure_mod flag_2$p points
R 737 5 52 data_structure_mod flag_2$o points
R 740 5 55 data_structure_mod nbhs points
R 741 5 56 data_structure_mod nbhs$sd points
R 742 5 57 data_structure_mod nbhs$p points
R 743 5 58 data_structure_mod nbhs$o points
R 747 5 62 data_structure_mod conn points
R 748 5 63 data_structure_mod conn$sd points
R 749 5 64 data_structure_mod conn$p points
R 750 5 65 data_structure_mod conn$o points
R 753 5 68 data_structure_mod nx points
R 754 5 69 data_structure_mod nx$sd points
R 755 5 70 data_structure_mod nx$p points
R 756 5 71 data_structure_mod nx$o points
R 758 5 73 data_structure_mod ny points
R 760 5 75 data_structure_mod ny$sd points
R 761 5 76 data_structure_mod ny$p points
R 762 5 77 data_structure_mod ny$o points
R 766 5 81 data_structure_mod prim points
R 767 5 82 data_structure_mod prim$sd points
R 768 5 83 data_structure_mod prim$p points
R 769 5 84 data_structure_mod prim$o points
R 773 5 88 data_structure_mod flux_res points
R 774 5 89 data_structure_mod flux_res$sd points
R 775 5 90 data_structure_mod flux_res$p points
R 776 5 91 data_structure_mod flux_res$o points
R 780 5 95 data_structure_mod q points
R 781 5 96 data_structure_mod q$sd points
R 782 5 97 data_structure_mod q$p points
R 783 5 98 data_structure_mod q$o points
R 788 5 103 data_structure_mod dq points
R 789 5 104 data_structure_mod dq$sd points
R 790 5 105 data_structure_mod dq$p points
R 791 5 106 data_structure_mod dq$o points
R 794 5 109 data_structure_mod entropy points
R 795 5 110 data_structure_mod entropy$sd points
R 796 5 111 data_structure_mod entropy$p points
R 797 5 112 data_structure_mod entropy$o points
R 799 5 114 data_structure_mod vorticity points
R 801 5 116 data_structure_mod vorticity$sd points
R 802 5 117 data_structure_mod vorticity$p points
R 803 5 118 data_structure_mod vorticity$o points
R 805 5 120 data_structure_mod vorticity_sqr points
R 807 5 122 data_structure_mod vorticity_sqr$sd points
R 808 5 123 data_structure_mod vorticity_sqr$p points
R 809 5 124 data_structure_mod vorticity_sqr$o points
R 812 5 127 data_structure_mod xpos_nbhs points
R 813 5 128 data_structure_mod xpos_nbhs$sd points
R 814 5 129 data_structure_mod xpos_nbhs$p points
R 815 5 130 data_structure_mod xpos_nbhs$o points
R 817 5 132 data_structure_mod xneg_nbhs points
R 819 5 134 data_structure_mod xneg_nbhs$sd points
R 820 5 135 data_structure_mod xneg_nbhs$p points
R 821 5 136 data_structure_mod xneg_nbhs$o points
R 823 5 138 data_structure_mod ypos_nbhs points
R 825 5 140 data_structure_mod ypos_nbhs$sd points
R 826 5 141 data_structure_mod ypos_nbhs$p points
R 827 5 142 data_structure_mod ypos_nbhs$o points
R 829 5 144 data_structure_mod yneg_nbhs points
R 831 5 146 data_structure_mod yneg_nbhs$sd points
R 832 5 147 data_structure_mod yneg_nbhs$p points
R 833 5 148 data_structure_mod yneg_nbhs$o points
R 837 5 152 data_structure_mod xpos_conn points
R 838 5 153 data_structure_mod xpos_conn$sd points
R 839 5 154 data_structure_mod xpos_conn$p points
R 840 5 155 data_structure_mod xpos_conn$o points
R 842 5 157 data_structure_mod xneg_conn points
R 845 5 160 data_structure_mod xneg_conn$sd points
R 846 5 161 data_structure_mod xneg_conn$p points
R 847 5 162 data_structure_mod xneg_conn$o points
R 851 5 166 data_structure_mod ypos_conn points
R 852 5 167 data_structure_mod ypos_conn$sd points
R 853 5 168 data_structure_mod ypos_conn$p points
R 854 5 169 data_structure_mod ypos_conn$o points
R 856 5 171 data_structure_mod yneg_conn points
R 859 5 174 data_structure_mod yneg_conn$sd points
R 860 5 175 data_structure_mod yneg_conn$p points
R 861 5 176 data_structure_mod yneg_conn$o points
R 864 5 179 data_structure_mod delta points
R 865 5 180 data_structure_mod delta$sd points
R 866 5 181 data_structure_mod delta$p points
R 867 5 182 data_structure_mod delta$o points
S 1049 23 5 0 0 0 1052 624 8338 0 0 A 0 0 0 0 B 0 134 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 wall_dgx_pos
S 1050 7 3 0 0 506 1 1049 8064 800004 3000 A 0 0 0 0 B 0 134 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 g
S 1051 1 3 0 0 6 1 1049 8252 4 3000 A 0 0 0 0 B 0 134 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 1052 14 5 0 0 0 1 1049 8338 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 99 2 0 0 0 0 0 0 0 0 0 0 0 0 20 0 624 0 0 0 0 wall_dgx_pos
F 1052 2 1050 1051
S 1053 23 5 0 0 0 1056 624 8351 0 0 A 0 0 0 0 B 0 257 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 wall_dgx_neg
S 1054 7 3 0 0 509 1 1053 8064 800004 3000 A 0 0 0 0 B 0 257 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 g
S 1055 1 3 0 0 6 1 1053 8252 4 3000 A 0 0 0 0 B 0 257 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 1056 14 5 0 0 0 1 1053 8351 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 102 2 0 0 0 0 0 0 0 0 0 0 0 0 141 0 624 0 0 0 0 wall_dgx_neg
F 1056 2 1054 1055
S 1057 23 5 0 0 0 1060 624 8364 0 0 A 0 0 0 0 B 0 374 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 wall_dgy_neg
S 1058 7 3 0 0 512 1 1057 8064 800004 3000 A 0 0 0 0 B 0 374 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 g
S 1059 1 3 0 0 6 1 1057 8252 4 3000 A 0 0 0 0 B 0 374 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 1060 14 5 0 0 0 1 1057 8364 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 105 2 0 0 0 0 0 0 0 0 0 0 0 0 259 0 624 0 0 0 0 wall_dgy_neg
F 1060 2 1058 1059
A 35 2 0 0 0 6 635 0 0 0 35 0 0 0 0 0 0 0 0 0 0 0
A 36 2 0 0 0 6 672 0 0 0 36 0 0 0 0 0 0 0 0 0 0 0
A 114 2 0 0 0 6 678 0 0 0 114 0 0 0 0 0 0 0 0 0 0 0
A 186 2 0 0 0 6 682 0 0 0 186 0 0 0 0 0 0 0 0 0 0 0
Z
T 690 62 0 0 0 0
A 694 7 236 0 1 2 1
A 693 6 0 36 1 2 1
A 700 7 238 0 1 2 1
A 699 6 0 36 1 2 1
A 706 7 240 0 1 2 1
A 705 6 0 36 1 2 1
A 712 7 242 0 1 2 1
A 711 6 0 36 1 2 1
A 718 7 244 0 1 2 1
A 717 6 0 36 1 2 1
A 724 7 246 0 1 2 1
A 723 6 0 36 1 2 1
A 730 7 248 0 1 2 1
A 729 6 0 36 1 2 1
A 736 7 250 0 1 2 1
A 735 6 0 36 1 2 1
A 742 7 252 0 1 2 1
A 741 6 0 36 1 2 1
A 749 7 254 0 1 2 1
A 748 6 0 114 1 2 1
A 755 7 256 0 1 2 1
A 754 6 0 36 1 2 1
A 761 7 258 0 1 2 1
A 760 6 0 36 1 2 1
A 768 7 260 0 1 2 1
A 767 6 0 114 1 2 1
A 775 7 262 0 1 2 1
A 774 6 0 114 1 2 1
A 782 7 264 0 1 2 1
A 781 6 0 114 1 2 1
A 790 7 266 0 1 2 1
A 789 6 0 186 1 2 1
A 796 7 268 0 1 2 1
A 795 6 0 36 1 2 1
A 802 7 270 0 1 2 1
A 801 6 0 36 1 2 1
A 808 7 272 0 1 2 1
A 807 6 0 36 1 2 1
A 814 7 274 0 1 2 1
A 813 6 0 36 1 2 1
A 820 7 276 0 1 2 1
A 819 6 0 36 1 2 1
A 826 7 278 0 1 2 1
A 825 6 0 36 1 2 1
A 832 7 280 0 1 2 1
A 831 6 0 36 1 2 1
A 839 7 282 0 1 2 1
A 838 6 0 114 1 2 1
A 846 7 284 0 1 2 1
A 845 6 0 114 1 2 1
A 853 7 286 0 1 2 1
A 852 6 0 114 1 2 1
A 860 7 288 0 1 2 1
A 859 6 0 114 1 2 1
A 866 7 290 0 1 2 1
A 865 6 0 36 1 2 0
Z
