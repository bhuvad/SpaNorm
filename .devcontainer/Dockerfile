# syntax=docker/dockerfile:1

FROM bioconductor/bioconductor_docker:devel
LABEL Name=spanorm Version=0.0.1
RUN mkdir /root/rpkg
COPY . /root/rpkg
RUN R -e "devtools::install('/root/rpkg', dependencies=T)"

# VS Code R session attachment (see github.com/BhuvaLab/vscode_for_r)
COPY .devcontainer/.Rprofile /root/.Rprofile