FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    jq \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    pandoc libsz2 libhdf5-dev libcurl4-gnutls-dev libuv1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('shiny', 'renv'))"

COPY . /srv/shiny-server/
WORKDIR /srv/shiny-server

RUN R -e "renv::restore()"

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host = '0.0.0.0', port = 3838)"]
