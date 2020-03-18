FROM ubuntu:18.04

MAINTAINER Ling-Hong Hung lhhunghimself@gmail.com

# Prepare R environment
ENV RHOME_DIR /usr/local/rhome
ENV PATH $RHOME_DIR/bin:$PATH
RUN mkdir -p $RHOME_DIR

# R pre-requisites
#To get R's blas and lapack must compile from source NOT from deb

RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils fonts-dejavu \
    build-essential xorg-dev gcc gcc-multilib gobjc++ gfortran libblas-dev libcairo2-dev liblzma-dev libreadline-dev aptitude \
    libbz2-dev libpcre3-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    software-properties-common wget texinfo texlive texlive-fonts-extra default-jdk && \
    		cd /tmp && wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz && \
    tar -xzvf R-3.6.1.tar.gz && \
    cd /tmp/R-* && ./configure && \
    cd /tmp/R-* && make -j 8 && \
    cd /tmp/R-* && make install rhome=$RHOME_DIR \
    && rm -rf /tmp/R-* \
    && apt-get remove -y build-essential wget\
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN echo R_UNZIPCMD="'"internal"'"  > /root/.Renviron


RUN apt-get update && apt-get install -y \
    t1-xfree86-nonfree ttf-xfree86-nonfree ttf-xfree86-nonfree-syriac \ 
    xfonts-100dpi xfonts-75dpi \
    libcairo2-dev libgtk2.0-dev \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get -y install firefox \ 
    build-essential python3-all python3-pip libncurses5-dev libncursesw5-dev libzmq3-dev \
    && pip3 install --upgrade pip \
    && pip install jupyter \
    && apt-get -y remove build-essential \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y git \
    && apt-get autoclean -y \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libx11-6 \
    libxss1 \
    libxt6 \
    libxext6 \
    libsm6 \
    libice6 \
    xdg-utils \
    libxext-dev \
    libxrender-dev \
    libxtst-dev \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("BiocInstaller",repos="http://bioconductor.org/packages/3.9/bioc",lib="/usr/local/rhome/lib/")'

RUN R -e "install.packages(c('IRdisplay', 'repr', 'devtools', 'evaluate', 'crayon','pbdZMQ', 'uuid', 'digest', 'BiocManager', 'stringr', 'RColorBrewer', 'Cairo', 'GlobalOptions', 'boot', 'foreign', 'mgcv', 'nlme'),repos = 'http://cran.us.r-project.org'); devtools::install_github('IRkernel/IRkernel'); IRkernel::installspec(); BiocManager::install(\"ComplexHeatmap\")"

RUN R -e "install.packages('config',repos = 'http://cran.us.r-project.org')"

#set mozilla preferences to not launch their homepage
COPY .mozilla /root/.mozilla

# Make a folder for AML config file and notebook
#ADD AML_Analysis.ipynb /AML
COPY . /AML

COPY /fonts/. /etc/fonts
COPY /X11/. /etc/X11

CMD ["/bin/sh"]