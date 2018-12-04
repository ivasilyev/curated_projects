#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import matplotlib.pyplot as plt
import seaborn as sns

caps = """
2681
2702
2675
2675
2137
2672
2676
2650
2664
2688
2717
2716
2656
2617
2679
2656
2710
2654
2720
2643
2674
2656
2686
2696
2680
2731
2680
2654
2682
2701
2639
2688
2688
2661
2673
2602
2662
2677
2710
2639
2708
2693
"""

data = [int(i.strip()) for i in re.sub("[\r\n]+", "\n", caps).split("\n") if len(i.strip()) > 0]
print(len(data))
sns.set(rc={"figure.figsize": (10, 10)})
sns.distplot(data, bins=20, kde=True, rug=True)
