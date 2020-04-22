set terminal epslatex standalone size 4.5in,3.2in color

set output "plot-linscaling.tex"

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set border 31 linewidth 2

# linear fit
f(x) = m*x + b
fit f(x) 'SM.dat' using 2:3 via m,b

# logarithmic fit for linear curve
lf(x) = x + lb
fit lf(x) 'SM.dat' using (log($2)):(log($3)) via lb

set logscale xy
set xrange [500:100000]

set xlabel 'Number of atoms'
set ylabel 'Time (s)'

plot 'SM.dat' using 2:3 title 'Submatrix Method' with lp lw 5 lc rgb col1 pt 6, \
     exp(lf(log(x))) title 'Linear' with l lc rgb col2
#     f(x) title 'Linear fit' with l

set output "bla"
!latex plot-linscaling.tex
!dvips plot-linscaling.dvi
!ps2eps -f plot-linscaling.ps
!rm *.dvi *.aux *.ps *-inc.eps *.tex bla *.log
