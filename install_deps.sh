mkdir deps_temp
cd deps_temp

# Install hmmer
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz
cd hmmer-3.3.2
./configure
make
make install
cd ..

# Install prodigal
wget https://github.com/hyattpd/Prodigal/archive/refs/tags/v2.6.3.tar.gz
tar zxf v2.6.3.tar.gz
cd Prodigal-2.6.3
make
make install
cd ..

cd ..
rm -rf deps_temp

# Concatenate all .hmm profiles into a single database file
cd profiles && ./_build_hmm_db.sh