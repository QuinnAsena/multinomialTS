# multinomialTS

This is the initial working version of `multinomialTS`. See [the vignette](https://quinnasena.github.io/multinomialTS/articles/multinomialTS-vignette.html) on how to use the model.

## Installation



### Windows

The best way to install at the moment is to compile from source:

- If you have Rtools44 or (Rtools45 for R 4.5) and devtools installed, you can build the latest version of the package directly from github using: `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`

Alternatively, if you are having issues with the build tools, install the binaries:

- `install.packages("https://github.com/QuinnAsena/multinomialTS/releases/download/v1.0.1-alpha/multinomialTS_1.0.1-alpha.zip", repos = NULL, type = "win.binary")`


### macOS

The best way to install at the moment is to compile from source:


- If you have xcode-select and devtools installed, you can build the latest version of the package directly from github using: `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`

- On apple it is nice and easy to download xcode-select by opening a terminal and copying this code: xcode-select --install. Then try and run: `devtools::install_github("https://github.com/QuinnAsena/multinomialTS")`


Alternatively, if you are having issues with the build tools, install the binaries:

For new M-series macs, use:

- `install.packages("https://github.com/QuinnAsena/multinomialTS/releases/download/v1.0.1-alpha/multinomialTS-arm64.tgz", repos = NULL)`

For Intel macs, use:

- `install.packages("https://github.com/QuinnAsena/multinomialTS/releases/download/v1.0.1-alpha/multinomialTS_1.0.0.tgz", repos = NULL)`

