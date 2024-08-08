FROM rocker/r-ver:4.4.1

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
    libxml2-dev \
    pandoc \
    libsz2 \
    libhdf5-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install shiny package
RUN R -e "install.packages('shiny')"

# Ensure directory exists
RUN mkdir -p /usr/lib/R/etc

# Configure Shiny to use port 3838 and listen on all network interfaces
RUN echo "options(shiny.port = 3838, shiny.host = '0.0.0.0')" > /usr/lib/R/etc/Rprofile.site

# Create a new group and user for running the application
RUN groupadd --system app && useradd --system --gid app app

# Set the working directory
WORKDIR /home/app

# Copy application files
COPY . .

# Change ownership of the app directory to the app user
RUN chown -R app:app /home/app

# Switch to app user
USER app

# Expose port 3838
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/home/app/app', port = 3838, host = '0.0.0.0')"]

