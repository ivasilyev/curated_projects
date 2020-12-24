#!/usr/bin/env bash

# Raspberry Pi 3B+ quickly reaches the 'soft throttle' point of 60 Celsius,
# designed to prevent the SoC hitting the hard-throttle maximum limit of 80 Celsius

# 1st window
# stress --cpu 4

# 2nd window
OUT="/tmp/temp_bench.csv"
echo "time,usage,temperature,frequency" > "${OUT}"

TIME=1

while true
    do
        USAGE=$(top -bn1 | sed -n '/Cpu/p' | awk '{print $2}' | sed 's/..,//')
        CELSIUS=$(($(</sys/class/thermal/thermal_zone0/temp)/1000))
        MHZ=$(($(</sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_cur_freq)/1000))
        echo "${TIME},${USAGE},${CELSIUS},${MHZ}" >> "${OUT}"
        printf "%s °C; \t" "${CELSIUS}"
        ((TIME=TIME+1))
        sleep 1
    done
