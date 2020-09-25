FROM rocker/shiny-verse:latest

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libfftw3-dev

RUN install2.r --error \
  shiny \
  shinyFiles \
  shinythemes \
  shinyjs \
  stringr \
  tidyr
  
RUN R -e "BiocManager::install('EBImage')"

# copy the app to the image
COPY *.Rproj /srv/shiny-server/
COPY *.R /srv/shiny-server/
COPY src /srv/shiny-server/src

# select port
EXPOSE 3838

# allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]

# run app
CMD ["/usr/bin/shiny-server.sh"]
