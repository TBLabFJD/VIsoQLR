# VIsoQLR: an interactive tool for the detection, quantification and fine-tuning of isoforms using long-read sequencing

VIsoQLR is an interactive analyzer, viewer and editor for the semi-automated identification and quantification of known and novel isoforms using long-read sequencing data. VIsoQLR is tailored to thoroughly analyze mRNA expression and maturation in low-throughput splicing assays. This tool takes sequences aligned to a reference, defines consensus splice sites, and quantifies isoforms. Users can edit splice sites through dynamic and interactive graphics and tables as part of their manual curation. Known transcripts, or isoforms detected by other methods, can also be imported as references for comparison.


## Developers
### Main developers
 - Gonzalo Núñez Moreno

### Contact
 - Gonzalo Núñez Moreno (gonzalo.nunezm@quironsalud.es)



## License
Mini-IsoQLR source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). Mini-IsoQLR includes several third party packages provided under other open source licenses, please check them for additional details.



## Installation
VIsoQLR has been wrapped into a Docker image. To set it up:

 1. Install Docker following the instructions in https://www.docker.com/. Docker is available for Linux, Windows and Mac users. For windows users, as the image has been developped for Linux kernel, it is neccesary to install Linux on Windows with WSL2 (Windows Subsystem for Linux) following the instructions in https://learn.microsoft.com/en-us/windows/wsl/install.
 2. Once installed type on your shell (or Windows PowerShell in Windows once Docker Desktop has been inicializated):
    ```
    docker pull tblabfjd/visoqlr:latest
    ```
    This will download the latest version of VIsoQLR image. This image has already all programs, packages and requirements needed to run VIsoQLR



## Run VIsoQLR
To run VIsoQLR type on your shell (or Windows PowerShell in Windows once Docker Desktop has been inicializated):
```
docker run -it -p 8888:8888 tblabfjd/visoqlr:latest
```
This will promp some text on the shell and once the line `Listening on http://0.0.0.0:8888` appears go to the browser and type `http://0.0.0.0:8888` on Linux systems or `http://localhost:8888` on Windows systems. If succesful you will see:

[![initialized_VIsoQLR](https://github.com/TBLabFJD/VIsoQLR/blob/main/images/initialized_VIsoQLR.png?raw=true)](https://github.com/TBLabFJD/VIsoQLR/tree/main/images/initialized_VIsoQLR.png)

Click on 'Browse...' button to navigate to your file system to upload your aligner reads in `GFF3` or `BED6` file format. Once selected the program will automatically run with its default parameters to identify and quantify isoforms. and you will see something like this:

[![default_analysis_VIsoQLR](https://github.com/TBLabFJD/VIsoQLR/blob/main/images/default_analysis_VIsoQLR.png?raw=true)](https://github.com/TBLabFJD/VIsoQLR/tree/main/images/default_analysis_VIsoQLR.png)


To know more about VIsoQLR have a look at the manuscript in bioRxiv: https://www.biorxiv.org/content/10.1101/2022.03.01.482488v2



