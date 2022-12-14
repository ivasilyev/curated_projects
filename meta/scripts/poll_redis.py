#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from meta.utils.redis import redis_to_strings
from meta.utils.io import dump_list


def parse_args():
    from argparse import ArgumentParser
    p = ArgumentParser(description="Script to fetch a string data from a Redis server",)
    p.add_argument("-o", "--output", required=True, help="Output file")
    p.add_argument("-z", "--size", default=1, type=int,
                   help="Number of samples in the sample data file")
    p.add_argument("-t", "--host", default="localhost", help="Redis host")
    p.add_argument("-p", "--port", default=6379, type=int, help="Redis port")
    p.add_argument("-q", "--queue", default="commands", help="Redis queue name")
    p.add_argument("-u", "--pause", default=60, type=int, help="Polling pause")
    namespace = p.parse_args()
    return (
        namespace.output,
        namespace.size,
        namespace.host,
        namespace.port,
        namespace.queue,
        namespace.pause,
    )


if __name__ == '__main__':
    (
        output_file,
        output_size,
        host_name,
        host_port,
        queue_name,
        pause,
    ) = parse_args()
    output_strings = redis_to_strings(
        queue_name=queue_name,
        content_length=output_size,
        pause=pause,
        host=host_name,
        port=host_port,
    )
    dump_list(output_strings, output_file)
