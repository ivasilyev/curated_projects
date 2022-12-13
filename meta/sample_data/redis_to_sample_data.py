#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import cpu_count
from meta.utils.redis import redis_to_dicts
from meta.sample_data.sample_data_array import SampleDataArray


def parse_args():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Script to fetch sample data from Redis server",)
    p.add_argument("-s", "--sampledata", required=True, help="Output sample data file")
    p.add_argument("-z", "--size", default=cpu_count(), type=int,
                   help="Number of samples in the sample data file")
    p.add_argument("-t", "--host", default="localhost", help="Redis host")
    p.add_argument("-p", "--port", default=6379, type=int, help="Redis port")
    p.add_argument("-q", "--queue", default="sample_data", help="Redis queue name")
    p.add_argument("-u", "--pause", default=60, type=int, help="Polling pause")
    namespace = parse_args()
    return (
        namespace.sampledata,
        namespace.size,
        namespace.host,
        namespace.port,
        namespace.queue,
        namespace.pause,
    )


if __name__ == '__main__':
    (
        sampledata_file,
        sampledata_size,
        host_name,
        host_port,
        queue_name,
        pause,
    ) = parse_args()
    sampledata_dicts = redis_to_dicts(
        queue_name=queue_name,
        content_length=sampledata_size,
        pause=pause,
        host=host_name,
        port=host_port,
    )
    sampledata_array = SampleDataArray.import_from_dicts(sampledata_dicts)
    sampledata_array.dump_dict(sampledata_file)
