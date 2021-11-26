#!/bin/bash

# create soft links to original data in $WORK
for f in $WORK/permafrost-pathogens-data/*; do ln -s $f; done
