set terminal epslatex standalone size 4.5in,3.2in color

set output "plot-strong-scaling.tex"

col1= 57*65536 + 106*256 + 177 #(57,106,177)
col2=218*65536 + 124*256 +  48 #(218,124,48)
col3= 62*65536 + 150*256 +  81 #(62,150,81)
col4=204*65536 +  37*256 +  41 #(204,37,41)
col5= 83*65536 +  81*256 +  84 #(83,81,84)
col6=107*65536 +  76*256 + 154 #(107,76,154)
col7=146*65536 +  36*256 +  40 #(146,36,40)
col8=148*65536 + 139*256 +  61 #(148,139,61)

set border 31 linewidth 2

set xlabel 'Number of CPU cores'
set ylabel 'Time (s)'
set yrange [2:13]
set ytics nomirror
set y2range [0.75:1.05]
set y2label 'Efficiency'
set y2tics (0.8, 0.9, "1.0" 1.0)


firsttime = 12.378
firstcores = 80

set style line 1 linecolor rgb col1 linewidth 5 dashtype 1 pointtype 6 pointsize default
set style line 2 linecolor rgb col2 linewidth 5 dashtype 1 pointtype 4 pointsize default
set style line 3 linecolor rgb col1 linewidth 5 dashtype 2 pointtype 6 pointsize default
set style line 4 linecolor rgb col2 linewidth 5 dashtype 2 pointtype 4 pointsize default
set style line 5 linecolor rgb "black" linewidth 5 dashtype 1
set style line 6 linecolor rgb "black" linewidth 5 dashtype 2

plot 'SM.dat' using 2:3 title 'Submatrix Method' with lp ls 1, \
     'SM.dat' using 2:(firsttime * firstcores / $2) title 'Perfect scaling' with lp ls 2, \
     'SM.dat' using 2:(firsttime * firstcores / $2 / $3) title '' with lp ls 3 axis x1y2, \
     NaN title 'Time' ls 5, \
     NaN title 'Efficiency' ls 6

set output "bla"
!latex plot-strong-scaling.tex
!dvips plot-strong-scaling.dvi
!ps2eps -f plot-strong-scaling.ps
!rm *.dvi *.aux *.ps *-inc.eps *.tex bla *.log
