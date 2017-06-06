#!/bin/bash

awk '(NR-1)%(1225*200)<1225' merger.uv.vis > merger.uv.cropped.vis
awk '(NR-1)%(1225*200)<1225' merger.uv.true.vis > merger.uv.true.cropped.vis

wc -l merger.uv.cropped.vis merger.uv.true.cropped.vis

