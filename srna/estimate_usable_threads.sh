#!/bin/bash

# Small RNA-seq analyses for Cledi usable thread estimation bash script
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

# file:    estimate_usable_threads.sh
# created: 2018-11-22
# author:  Marcel Schilling <marcel.schilling@mdc-berlin.de>
# purpose: estimate usable threads for analyses of sRNA-seq data for
#          Cledi


######################################
# change log (reverse chronological) #
######################################

# 2018-11-22: initial version


##############
# parameters #
##############

# Define user group to bootstrap CPU usage for.
GROUP=AG_Rajewsky

# Define number of samples to bootstrap.
NBOOTSTRAPS=10000

# Define percentile of bootstrapped samples to assume for estimation.
PERCENTILE=95

# Define fraction of total available threads to leave free at any cost.
MARGIN=.1


#################
# bash settings #
#################

# Exit with error as soon as any command exists with error.
set -e


#################
# collecta data #
#################

# Get total number of available threads.
TOTAL=$(nproc)

# Get CPU usage per user.
USAGE=$(top -b -n1 \
        | awk -vORS='\0' 'NR > 7{print $2, $9}' \
       )

# Get users to bootstrap CPU usage for.
USERS=$(who \
        | awk '$0=$1' \
        | sort -u \
        | grep -v $(whoami) \
        | xargs -n1 id \
        | grep $GROUP \
        | sed -e 's/[^(]*(//' -e 's/).*//' \
       )


##################
# summarize data #
##################

# Sum up CPU usage for users not to bootstrap CPU usage for.
OTHER=$(tr ' ' '\n' <<< $USERS \
        | awk '
            NR==FNR{users[$0];next}
            !($1 in users){cpu+=$2}
            END{print cpu/100}
          ' - <(tr '\0' '\n' <<< $USAGE) \
       )

# Count users to bootstrap CPU usage for.
NUSERS=$(echo $USERS \
          | tr ' ' '\n' \
          | wc -l \
         )


#######################
# bootstrap CPU usage #
#######################

# Sum CPU usage for users to bootstrap CPU usage for (per user),
# bootstrap samples of users, sum up CPU usage per bootstrapped sample,
# and select specified percentile.
BOOTSTRAPPED=$(tr ' ' '\n' <<< $USAGE \
               | awk '
                   NR==FNR{cpu[$0];next}
                   $1 in cpu{cpu[$1]+=$2}
                   END{for(u in cpu)print u,cpu[u]/100}
                 ' - <(tr '\0' '\n' <<< $USAGE) \
               | shuf -r -n $(bc <<< $NBOOTSTRAPS*$NUSERS) \
               | awk '{$1=int((NR - 1)/'$NUSERS')}1' \
               | awk '
                   NR==1{i=$1;s=$2;next}
                   $1!=i{print s;i=$1;s=$2}
                   {s+=$2}
                   END{print s}
                 ' \
               | sort -n \
               | head -n $(bc <<< $NBOOTSTRAPS/100*$PERCENTILE) \
               | tail -n 1 \
              )


###########################
# estimate usable threads #
###########################

# Substract (bootstrapped) CPU usage from total number of available
# threads, minus the specified fraction to leave free at all cost, and
# floor the result.
printf '%.0f\n' $(bc <<< "(1-$MARGIN)*$(nproc) - $OTHER - $BOOTSTRAPPED")
