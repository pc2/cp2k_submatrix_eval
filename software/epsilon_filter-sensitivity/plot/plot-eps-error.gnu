set terminal epslatex standalone size 4.5in,3.2in color

set output "plot-eps-error.tex"

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set border 31 linewidth 2

set key samplen 1 width -18 vertical maxrows 2 tmargin right

set logscale y 10
set xrange [-9:-2]
set yrange [1e-6:10]

set xlabel '$\epsilon_\textrm{filter}$'
set ylabel 'Abs. error in energy [meV/atom]'

set format x "$10^{%.0f}$"
set format y "$10^{%T}$"

set style line 1 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 6 pointsize 1
set style line 2 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 4 pointsize 1
set style line 3 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 2 pointsize 1
set style line 4 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 8 pointsize 1

prec = -12167.089585589709
atoms = 32 * 3 * 6 * 6 * 6 # NREP 6
hartree = 27.2116 # eV

plot 'SM.dat' using 1:(($2-prec)/atoms*hartree*1000) title 'Submatrix Method (pos. error)' with p ls 1, \
     'NS.dat' using 1:(($2-prec)/atoms*hartree*1000) title 'Newton-Schulz (pos. error)' with p ls 2, \
     'SM.dat' using 1:(-($2-prec)/atoms*hartree*1000) title '(neg. error)' with p ls 3, \
     'NS.dat' using 1:(-($2-prec)/atoms*hartree*1000) title '(neg. error)' with p ls 4, \
     'SM-shape.dat' using 1:(abs($2-prec)/atoms*hartree*1000) title '' with l ls 1, \
     'NS-shape.dat' using 1:(abs($2-prec)/atoms*hartree*1000) title '' with l ls 2

set output "bla"
!latex plot-eps-error.tex
!dvips plot-eps-error.dvi
!ps2eps -f plot-eps-error.ps
!ps2pdf -f plot-eps-error.ps
!rm *.dvi *.aux *.ps *-inc.eps *.tex bla
