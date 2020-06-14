# genome-map
Option to draw SVG of covered regions based on depth plots

## Authors
Jason Kwong (@kwongjc)  ::  [kwongj](https://github.com/kwongj)  

## Dependencies
* [Python 3.x](https://www.python.org/downloads/)
* [pandas](https://pypi.python.org/pypi/pandas/)
* [svgwrite](https://pypi.python.org/pypi/svgwrite/)

## Usage
```
$ genome-map.py -h
usage: 
  genome-map.py --svg map.svg <DEPTH-FILE>

Draw covered regions of genome

positional arguments:
  DEPTH-FILE           samtools depth output file (required)

optional arguments:
  -h, --help           show this help message and exit
  --mindepth INT       set threshold for min depth cutoff (default=0)
  --out FILE           save SVG output as specified file (default=map.svg)
  --size WIDExHIGH  specify width and height of SVG in pixels (default="800x600")
  --colour COLOUR   specify colour of recombination regions in HEX format (default=black)
  --version            show program's version number and exit
```

**Requires:**
* Depth file from samtools depth with all sites reported (-aa)

**Options:**
* Specify output file using `--out OUTFILE`
* Specify size of SVG in pixels eg. 800x600 `--size WIDTHxHEIGHT`
* Specify colour to show covered regions `--colour COLOUR`
* Specify threshold for displaying coverage `--mindepth`

## Bugs
Please submit via the GitHub issues page: [https://github.com/kwongj/genome-map/issues](https://github.com/kwongj/genome-map/issues)  

## Software Licence
GPLv3: [https://github.com/kwongj/genome-map/blob/master/LICENSE](https://github.com/kwongj/genome-map/blob/master/LICENSE)

## Other links
* [samtools depth](http://www.htslib.org/doc/samtools-depth.html)
