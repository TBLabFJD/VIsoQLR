# get shiny server and R from the rocker project
FROM rocker/shiny:4.0.5

## Documentation of the image
LABEL authors="Gonzalo Núñez Moreno" \
      description="VIsoQLR"

# Install libcurl4-gnutls-dev and libssl-dev
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev


# install R packages required 
RUN R -e 'install.packages(c( \
              "shiny", \
              "vctrs", \
              "plotly", \
              "DT", \
              "htmlwidgets", \
              "reticulate"\
            ), \
            repos="http://cran.rstudio.com/"\
          )'

RUN R -e "reticulate::install_miniconda()"
RUN R -e "reticulate::conda_install('r-reticulate', 'python-kaleido')"
RUN R -e "reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')"
RUN R -e "reticulate::use_miniconda('r-reticulate')"


# ENV CONDA_DIR /opt/conda
# RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# RUN /bin/bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH
# RUN conda install -c conda-forge python-kaleido




# copy the app directory into the image
COPY ./shiny-app/* /srv/shiny-server/


EXPOSE 8888

# run app
CMD R -e "shiny::runApp('/srv/shiny-server/VIsoQLR.R', host='0.0.0.0', port=8888)"


#docker build -t visoqlr -f Dockerfile .
#docker tag local-image:tagname new-repo:tagname
#docker push new-repo:tagname
