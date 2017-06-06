#!/bin/bash

source ~/.bashrc > /dev/null

casa -c observe.py
casa -c get-data.py
./crop-data.sh
