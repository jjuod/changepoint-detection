#!/bin/bash

cd ~/Documents/kiwis/changepointdet/chipseq/BROAD/

~/soft/bigWigToBedGraph wgEncodeBroadHistoneGm12878ControlStdSig.bigWig gm12878Control.bedGraph
~/soft/bigWigToBedGraph wgEncodeBroadHistoneGm12878H3k36me3StdSig.bigWig gm12878H3k36me3.bedGraph

awk '$1=="chr1" && $2>=120000000 && $2<122000000' gm12878Control.bedGraph > filtered_control.bedGraph
awk '$1=="chr1" && $2>=120000000 && $2<122000000' gm12878H3k27ac.bedGraph > filtered_h3k27ac.bedGraph

