library(httpuv)
library(here)

# Define the path to the HTML file
html_file <- here("EarlyGestation_RNAseq.html")

# Define a simple handler to serve the HTML file
app <- list(
  call = function(req) {
    list(
      status = 200L,
      headers = list(
        "Content-Type" = "text/html"
      ),
      body = readBin(html_file, "raw", n = file.info(html_file)$size)
    )
  }
)

# Start the server
port <- 8080
server <- startServer("0.0.0.0", port, app)

cat(sprintf("Serving HTML on port %d\n", port))

# Keep the server running
service()

# Keep the script running to keep the server alive
while (TRUE) {
  service()
  Sys.sleep(1)
}