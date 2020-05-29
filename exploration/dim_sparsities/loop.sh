set -x
wd=`pwd`
for b in SZV DZVP;
do
  #for N in `seq 1 5`; 
  for nodes in 14;
  do
    ntasks=`echo "$nodes*40" | bc`
    for N in 1 2 3 4 5 6 7; #NREP
    do
      for MODE in SUBMATRIX;
      do
        for tol in "0.00001";
        do
          cd $wd
          dir="H2OTEST_${b}_${N}_${MODE}_${tol}_${nodes}"
          rm -rf $dir
          mkdir $dir
          cp template.sh $dir
          cp H2O-dft-ls.inp $dir/H2O-dft-ls.inp

          cd $dir
          sed -i "s/_B_/$b/g" H2O-dft-ls.inp
          sed -i "s/_N_/$N/g" H2O-dft-ls.inp
          sed -i "s/_MODE_/$MODE/g" H2O-dft-ls.inp
          sed -i "s/_TOL_/$tol/g" H2O-dft-ls.inp
          sed -i "s/_NODES_/$nodes/g" template.sh
          sed -i "s/_n_/$ntasks/g" template.sh
          sbatch -J $dir template.sh
        done
      done
    done
  done
done
