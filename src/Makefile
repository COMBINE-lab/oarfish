CC = g++

CFLAGS = -Wall -shared -fPIC

TARGET = mlpack_kde_wrapper

all: $(TARGET)

$(TARGET).so: $(TARGET).cpp
    $(CC) $(CFLAGS) -o $(TARGET).so $(TARGET).cpp -I /home/zzare/mlpack-4.2.0/build/install/include

g++ -c -Wall -fPIC mlpack_kde_wrapper.cpp -o mlpack_kde_wrapper.o -I /home/zzare/mlpack-4.2.0/build/install/include
g++ -shared -o mlpack_kde_wrapper.so mlpack_kde_wrapper.o

correct command:

    g++ -shared -o libmlpack_kde_wrapper.so mlpack_kde_wrapper.cpp -I /home/zzare/mlpack-4.2.0/build/install/include -fPIC
    RUSTFLAGS="-C link-arg=-Wl,-rpath,/mnt/scratch4/zahra/oarfish/src" cargo build --release
    ./target/release/oarfish --alignments ./scripts/output2/sirv/oarfish/Alignments/SRR6058583_c.bam --out-cov-prob 1D_kde_prob.tsv --stat-output 1D_kde_stat.tsv --coverage yes --prob kde_c --rate 1D --bins 1 > 1D_kde_quant.tsv


    ./target/release/oarfish --alignments ./scripts/output2/sirv/oarfish/Alignments/SRR6058583_c.bam --out-cov-prob 1D_kde_prob.tsv --stat-output 1D_kde_stat.tsv --coverage yes --prob multinomial --rate dr --bins 1 > 1D_kde_quant.tsv

    ./target/release/oarfish --alignments ./scripts/output2/sequin/oarfish/Alignments/ERR4352441.bam --out-cov-prob 1D_kde_prob.tsv --stat-output 1D_kde_stat.tsv --coverage yes --prob kde_c --rate 1D --bins 1 > 1D_kde_quant.tsv


    RUSTFLAGS="-C link-arg=-Wl,-rpath,/nfshomes/zzare/oarfish/src -L/nfshomes/zzare/oarfish/src" cargo build --release

    /opt/local/stow/gcc-11.2.0/bin/g++ -shared -o libmlpack_kde_wrapper.so mlpack_kde_wrapper.cpp -I /nfshomes/zzare /mlpack-4.2.0/build/install/include -fPIC -L /opt/local/stow/gcc-11.2.0/lib64

RUSTFLAGS="-C link-arg=-Wl,-rpath,/nfshomes/zzare/oarfish/src -L/opt/local/stow/gcc-11.2.0/lib64" cargo build --release