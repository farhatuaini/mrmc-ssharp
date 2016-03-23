help common
set method_path gauss_seidel
set method_path gauss_jacobi
P { >= 0.15 } [ tt U formula0 ]
print tree
$STATE[1]
$RESULT[1]
write_res_file 1 2 4
write_res_file 1 2 3 4 5 7 8 9 10 11 12
quit
