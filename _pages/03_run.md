---
layout: archive
permalink: /run/
title: "Run!"
---

## How to use it?

The pipeline is separated in two sub-pipelines. _redccd_ and _redspec_.
The `--help` argument will print the argument plus some some description

```shell
$ redccd --help
  usage: redccd [-h] [-c] [--ignore-bias] [--auto-clean] [--saturation <Value>]
                [--raw-path raw_path] [--red-path red_path] [--debug]
                [--log-to-file] [--flat-normalize <Normalization Method>]
                [--flat-norm-order <Order>] [--dcr-par-dir <dcr.par directory>]
                [--keep-cosmic-files]

Goodman CCD Reduction - CCDreductions for Goodman spectroscopic data

optional arguments:
  -h, --help            Show this help message and exit
  -c, --cosmic          Clean cosmic rays from science data.
  --ignore-bias         Ignore bias correction
  --auto-clean          Automatically clean reduced data directory
  --saturation <Value>  Saturation limit. Default to 55.000 ADU (counts)
  --raw-path raw_path   Path to raw data.
  --red-path red_path   Path to reduced data.
  --debug               Show detailed information of the process.
  --log-to-file         Write log to a file.
  --flat-normalize <Normalization Method>
                        Choose a method to normalize the master flat
                        forspectroscoy. Choices are: mean, simple (model) and
                        full (fits model to each line).
  --flat-norm-order <Order>
                        Defines the order of the model to be fitted.
  --dcr-par-dir <dcr.par directory>
                        Directory of default dcr.par file.
  --keep-cosmic-files   After cleaning cosmic rays with dcr, do not remove the
                        input file and the cosmic rays file.

  ```

And for `redspec`:

```shell
$ redspec --help
  usage: redspec [-h] [--data-path <Source Path>]
                 [--proc-path <Destination Path>]
                 [--search-pattern <Search Pattern>]
                 [--output-prefix <Out Prefix>] [--extraction <Extraction Type>]
                 [--reference-files <Reference Dir>] [--interactive] [--debug]
                 [--log-to-file] [--save-plots] [--plot-results]

Extracts goodman spectra and does wavelength calibration.

optional arguments:
  -h, --help            show this help message and exit
  --data-path <Source Path>
                        Path for location of raw data. Default <./>
  --proc-path <Destination Path>
                        Path for destination of processed data. Default <./>
  --search-pattern <Search Pattern>
                        Pattern for matching the goodman's reduced data.
  --output-prefix <Out Prefix>
                        Prefix to add to calibrated spectrum.
  --extraction <Extraction Type>
                        Choose a which extraction to perform. Simple is a sum
                        across the spatial direction after the background has
                        been removed. Optimal is a more advanced method that
                        considers weights and profilefitting.
  --reference-files <Reference Dir>
                        Directory of Reference files location
  --interactive         Interactive wavelength solution.Disbled by default.
  --debug               Debugging Mode
  --log-to-file         Write log to a file
  --save-plots          Save all plots in a directory
  --plot-results        Show wavelength calibrated spectrum at the end.

   ```

You should always run `redccd` first and then `redspec`. There are certain
defaults values

### redccd Defaults

```shell
    --cosmic              False
    --ignore-bias         False
    --auto-clean          False
    --debug               False
    --log-to-file         False
    --keep-cosmic-files   False
    --saturation          55000
    --raw-path            ./
    --red-path            ./RED/
    --flat-normalize      simple
    --dcr-par-dir         files/
    --flat-norm-order     15
```
### redspec Defaults

```shell
    --data-path         ./
    --proc-path         ./
    --search-pattern    cfzsto
    --extraction        simple
    --reference-files   refdata
    --reference-lamp    (empty string)
    --output-prefix     g
    --interactive       False
    --debug             False
    --log-to-file       False
    --save-plots        False
    --plot-results      False
```
