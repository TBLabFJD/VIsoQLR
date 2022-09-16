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
RUN R -e 'install.packages(c(\
              "shiny", \
              "vctrs", \
              "plotly", \
              "DT", \
              "htmlwidgets", \
              "shinyFiles", \
              "reticulate" \
            ), \
            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23"\
          )'


# copy the app directory into the image
COPY ./shiny-app/* /srv/shiny-server/


EXPOSE 8888

# run app
CMD R -e "shiny::runApp('/srv/shiny-server/VIsoQLR.R', host='0.0.0.0', port=8888)"
