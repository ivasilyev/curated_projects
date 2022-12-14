#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from meta.utils.redis import dicts_to_redis
from meta.sample_data.sample_data_array import SampleDataArray


def parse_args():
    p = argparse.ArgumentParser(description="Script to upload sample data to Redis server",)
    p.add_argument("-s", "--sampledata", required=True, help="Input sample data file")
    p.add_argument("-t", "--host", default="localhost", help="Redis host")
    p.add_argument("-p", "--port", default=6379, type=int, help="Redis port")
    p.add_argument("-q", "--queue", default="sample_data", help="Redis queue name")
    p.add_argument("-f", "--flush", default=False, action="store_true", help="Flush all Redis queues")
    namespace = p.parse_args()
    return (
        namespace.sampledata,
        namespace.host,
        namespace.port,
        namespace.queue,
        namespace.flush
    )


if __name__ == '__main__':
    (
        sampledata_file,
        host_name,
        host_port,
        queue_name,
        is_flush,
    ) = parse_args()

    sampledata_array = SampleDataArray.load_dict(sampledata_file)
    dicts_to_redis(
        dicts=sampledata_array.export_to_dict().values(),
        is_flush=is_flush,
        queue_name=queue_name,
        host=host_name,
        port=host_port,
    )
