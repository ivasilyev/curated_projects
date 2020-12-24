#!/usr/bin/env bash


# Raspberry Pi 3B+ quickly reaches the 'soft throttle' point of 60 Celsius,
# designed to prevent the SoC hitting the hard-throttle maximum limit of 80 Celsius

# 1st window
# stress --cpu 4 --timeout 600
# 2nd window
rm -rf /tmp/*
i=1
while true
    do
        USAGE=$(top -bn1 | sed -n '/Cpu/p' | awk '{print $2}' | sed 's/..,//')
        CELSIUS=$(($(</sys/class/thermal/thermal_zone0/temp)/1000))
        echo "${i},${USAGE},${CELSIUS}" >> /tmp/temp_bench.csv
        ((i=i+1))
        sleep 1
    done
