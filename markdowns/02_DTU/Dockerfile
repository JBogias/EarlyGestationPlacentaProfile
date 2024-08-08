# Stage 1: Base image with R and renv
FROM bioconductor/r-ver:3.19-R-4.4.1 AS base

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    jq \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libsz2 \
    libhdf5-dev \
    pandoc \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables for Bioconductor
ENV BIOCONDUCTOR_VERSION=3.19
ENV R_BIOC_VERSION=$BIOCONDUCTOR_VERSION

# Set the Bioconductor repository URL
RUN echo 'options(BioC_mirror="https://bioconductor.org")' >> /usr/local/lib/R/etc/Rprofile.site
RUN echo 'options(repos = BiocManager::repositories())' >> /usr/local/lib/R/etc/Rprofile.site

# Install BiocManager package with retries
RUN R -e "options(repos = c(CRAN = 'https://cran.rstudio.com')); \
          tryCatch(install.packages('BiocManager', dependencies=TRUE, verbose=TRUE), error=function(e) { \
          Sys.sleep(5); \
          install.packages('BiocManager', dependencies=TRUE, verbose=TRUE) \
          })"

# Install renv package with retries
RUN R -e "options(repos = c(CRAN = 'https://cran.rstudio.com')); \
          tryCatch(install.packages('renv', dependencies=TRUE, verbose=TRUE), error=function(e) { \
          Sys.sleep(5); \
          install.packages('renv', dependencies=TRUE, verbose=TRUE) \
          })"

# Set working directory
WORKDIR /home/app

# Copy only the renv-related files first
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R

# Make sure cache directory exists
RUN mkdir -p /home/app/renv/cache

# Set the RENV_PATHS_CACHE environment variable to the shared cache
ENV RENV_PATHS_CACHE=/home/app/renv/cache

# Restore the renv environment
RUN R -e "renv::restore()"

# Now copy the rest of the application files
COPY . .

# Stage 2: Final image to run the Shiny app
FROM rocker/r-ver:4.4.1

RUN apt-get update && \
    apt-get install -y curl pandoc libsz2 libhdf5-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Ensure directory exists
RUN mkdir -p /usr/lib/R/etc

# Configure Shiny to use port 3838 and listen on all network interfaces
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site

# Create a new group and user for running the application
RUN groupadd --system app && useradd --system --gid app app

# Set the working directory
WORKDIR /home/app

# Copy the pre-loaded R packages from the base image
COPY --from=base /home/app .

# Change ownership of the app directory to the app user
RUN chown -R app:app /home/app

# Make the R markdown
RUN R -e "rmarkdown::render('/home/app/EarlyGestation_DTU.Rmd', output_file = '/home/app/www/EarlyGestation_DTU.html')"

# Switch to the app user
USER app

# Expose port 3838 for the Shiny app
EXPOSE 3838

# Command to run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app', port = 3838, host = '0.0.0.0')"]
