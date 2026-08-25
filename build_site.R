# once
install.packages("pkgdown", repos = "https://cloud.r-project.org")
devtools::install()   # or: pak::local_install()

# build + open in browser
pkgdown::build_site()

browseURL("docs/index.html")
