# rel2bed
usage: rel2bed [-h] -5 UTRREF5 -3 UTRREF3 -i INPUT

Translates relative 5u 3u positions into their absolute reference position in
BED format.

optional arguments:
  -h, --help            show this help message and exit
  -5 UTRREF5, --utrref5 UTRREF5
                        The path to the 5'utr reference file.
  -3 UTRREF3, --utrref3 UTRREF3
                        The path to the 3'utr reference file.
  -i INPUT, --input INPUT
                        The path to the file with the relative positions.
