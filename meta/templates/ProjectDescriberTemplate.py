#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC, ABCMeta, abstractmethod


class ProjectDescriberTemplate(ABC):
    __metaclass__ = ABCMeta

    def __init__(self, owner: str, name: str, directory: str, groupdata: str, sampledata: str, mask: str):
        super().__init__()
        self.owner = owner
        self.name = name
        self.directory = directory
        self.groupdata = groupdata
        self.sampledata = sampledata
        self.mask = mask
