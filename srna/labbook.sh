#!/bin/bash

# Small RNA-seq analyses for Cledi labbook bash script
#
# Copyright (C) 2018  Marcel Schilling
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#######################
# general information #
#######################

# file:        labbook.sh
# created:     2018-11-21
# last update: 2018-11-22
# author:      Marcel Schilling <marcel.schilling@mdc-berlin.de>
# purpose:     analyse sRNA-seq data for Cledi


######################################
# change log (reverse chronological) #
######################################

# 2018-11-22: added usable thread estimation to snakemake call
# 2018-11-21: initial version (instantiate Guix Manifest, run Snakefile,
#             change permissions)


##############
# parameters #
##############

# Define group to grant read access.
GROUP="AG_Rajewsky"


#########
# paths #
#########

# Define path of Guix manifest.
PROJECT_GUIX_MANIFEST=".guix-manifest.scm"

# Define path of Guix profile.
PROJECT_GUIX_PROFILE=".guix-profile"

# Define path of bash script used to estimate usable theads.
ESTIMATE_USABLE_THREADS_SH="./estimate_usable_threads.sh"


#################
# bash settings #
#################

# Exit with error as soon as any command exists with error.
set -e


############
# commands #
############

# Define command used to run the Guix package manager.
GUIX="guix"

# Define command used to run pipeline described in Snakefile.
SNAKEMAKE="snakemake"

# Define command used to estimate usable theads.
ESTIMATE_USABLE_THREADS="$ESTIMATE_USABLE_THREADS_SH"

# Define command used to change file permissions.
CHANGE_PERMISSIONS="chmod --recursive"

# Define command used to change file group ownership.
CHANGE_GROUP="chgrp --recursive"


##################
# setup software #
##################

# Install Guix packages.
$GUIX package \
  --profile="$PROJECT_GUIX_PROFILE" \
  --manifest="$PROJECT_GUIX_MANIFEST"

# Load Guix profile.
source "$PROJECT_GUIX_PROFILE"/etc/profile


################
# run pipeline #
################

# Run pipeline described in Snakefile estimating the number of usable threads.
$SNAKEMAKE \
  --cores $($ESTIMATE_USABLE_THREADS)


###################
# set permissions #
###################

# Revoke any permissions from other users.
$CHANGE_PERMISSIONS o-rwx .

# Revoke write permission from group members.
# Note: Not done here as both, Marcel and Seung, need to be able to modify stuff.
#$CHANGE_PERMISSIONS g-w .

# Change group ownership.
$CHANGE_GROUP "$GROUP" .

# Grant group members read permission.
# Note: Not done here as both, Marcel and Seung, need to be able to modify stuff.
#$CHANGE_PERMISSIONS g+r .

# Grant group members read and write permission.
# Note: Only done here as both, Marcel and Seung, need to be able to modify stuff.
$CHANGE_PERMISSIONS g+rw .
