#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import ABC, ABCMeta, abstractmethod


class ProjectDescriberTemplate(ABC):
    __metaclass__ = ABCMeta
    owner = ""
    name = ""
    directory = ""
    groupdata = ""
    sampledata = ""
    mask = ""

    def __init__(self):
        super().__init__()
