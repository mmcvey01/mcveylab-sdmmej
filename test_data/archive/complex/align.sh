#!/bin/bash

muscle -in complex_i.fa -out complex_i.out.fa
cat complex_i.out.fa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > complex_i.out.single.fa

muscle -in complex_d.fa -out complex_d.out.fa
cat complex_d.out.fa | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > complex_d.out.single.fa
