
#./metis-5.1.0/build/Linux-x86_64/programs/gpmetis graph $1 -ptype=rb -niter=100 -seed=123984 > graph.${1}.log
./metis-5.1.0/build/Linux-x86_64/programs/gpmetis graph $1 -ptype=kway -objtype=vol -niter=100  -seed=123984 -minconn -ufactor=200 > graph.${1}.log
#./metis-5.1.0/build/Linux-x86_64/programs/gpmetis graph $1 -ptype=kway -seed=123984  > graph.${1}.log
mv graph.part.$1 graph.out
