# PMDetectionDemo

![Home screen.](../assets/pm-mixed.png)

## Overview
PMDetectionDemo is a visual tool to help users configure the CellQuant automated pipeline.  Demonstrates all possible membrane detection algorithms at all possible settings. 

## How to Run

1.  Install [Docker](https://www.docker.com/products/docker-desktop)

2.  Make sure Docker is running and pull the cellquant image from Docker Hub.  **You only need to do this once!**

`docker pull odonnelllab/pm-detection-demo`

3.  Run the following command to create a container and run on port 3838.  The -v flag maps the container directories with the local file system.

`docker run --rm -p 3838:3838 -v $HOME:/srv/shiny-server/home odonnelllab/pm-detection-demo`

The program will be running on localhost:3838.

A detailed installation tutorial for Docker and CellQuant is available [here](Tutorial/CellQuant-Installation-Instructions.pdf).


PMDetectionDemo is also included as part of the CellQuant [bundled release](https://github.com/sah129/CellQuant/releases/tag/v0.8-alpha).


## Dependencies
* Docker

CellQuant was built with R 3.6.3 with the following package dependencies: 
```
Bioconductor::EBImage
stringr
shiny
shinyFiles
shinyjs
shinythemes
tidyr
```
 
Created by Sarah Hawbaker for the O'Donnell Lab at the University of Pittsburgh, 2020.
