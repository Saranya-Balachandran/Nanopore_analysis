FROM ubuntu

RUN apt update
RUN apt install -y curl
RUN curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.0.2-linux-x64.tar.gz" -o dorado-1.0.2-linux-x64.tar.gz
RUN tar -xzf dorado-1.0.2-linux-x64.tar.gz
RUN dorado-1.0.2-linux-x64/bin/dorado --version
RUN mkdir dorado_models
RUN dorado-1.0.2-linux-x64/bin/dorado download --model dna_r10.4.1_e8.2_400bps_hac@v5.2.0 --directory dorado_models
