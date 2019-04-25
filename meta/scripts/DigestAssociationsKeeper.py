#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from meta.scripts.CounterWrapper import CounterWrapper


class DigestAssociationsKeeper:
    """
    This class is designed for digestion of annotated table.
    If some table row contain a word from the values list, it will be counted into the group corresponding to the key.
    Then the group values will be combined.
    Intersections (e.g. multiple selection of the same row into different groups) are allowed.
    """
    VIRULENCE_FACTORS = {"adhesion": ("adhesin", "adhesion", "laminin"),
                         "invasion": ("invasion", "invasin"),
                         "iron metabolism": (
                             "iron", "siderophore", "ferric", "ferrienterobactin", "ferrochelatase", "ferrichrome",
                             "aerobactin", "enterochelin", "enterobactin", "yersiniabactin", "yersinabactin",
                             "ferrienterobactin", "chrysobactin", "ornibactin", "precolibactin", "colibactin",
                             "ferripyoverdine"),
                         "pili-related": (
                             "pilin", "pilus", "prepilin", "fimbrial", "fimbriae", "fimbrillin", "fimbrin"),
                         "flagella-related": ("motor", "flagellin", "flagellar", "flagella", "flagellum"),
                         "regulation": ("regulator", "regulatory", "regulation", "receptor", "effector"),
                         "transport": ("transport", "transporter", "permease", "porin", "export"),
                         "toxin": (
                             "toxin", "enterotoxin", "cytotoxic", "necrotizing", "endotoxin", "leukotoxin", "exotoxin",
                             "cytolethal", "o-antigen", "lipooligosaccharide", "lipopolysaccharide"),
                         "chaperones": ("chaperone", "chaperonin"),
                         "cytolysins": ("hemolysin", "cytolysin"),
                         "cell wall related": ("cell wall",),
                         "cell membrane related": ("membrane",),
                         "capsule-related": ("capsular", "capsule"),
                         "lipoglycans": ("o-antigen", "lipooligosaccharide", "lipopolysaccharide"),
                         "secretion": ("secretion", "secreted", "secretory"),
                         "efflux": ()}
    DRUG_CLASSES = {"beta-lactam": ("cephalosporin", "penam", "penem"),
                    "aminoglycoside": (),
                    "fluoroquinolone": (),
                    "glycopeptide antibiotic": (),
                    "lincosamide": (),
                    "macrolide": (),
                    "nucleoside antibiotic": (),
                    "peptide antibiotic": (),
                    "phenicol": (),
                    "sulfonamide": (),
                    "tetracycline": (),
                    "triclosan": ()}
    RESISTANCE_MECHANISMS = {"efflux": ("efflux", "efflux pump"),
                             "inactivation": ("inactivation", "inhibition", "destruction"),
                             "reduced permeability": (
                                 "reduced permeability", "membrane modulation", "viscosity regulator"),
                             "target alteration": ("target alteration", "target methylation", "target glycosylation"),
                             "target protection": ("target protection", "target enhancement", "target binding"),
                             "target replacement": ("target replacement", )}

    @staticmethod
    def generate_keywords_dict(keywords: list, split_words: bool = False):
        from meta.scripts.Utilities import Utilities
        keywords = [i.strip() for i in keywords if isinstance(i, str)]
        if split_words:
            keywords = Utilities.flatten_2d_array([i.split(" ") for i in keywords])
        return {j: () for j in sorted([i for i in set(keywords)])}

    @staticmethod
    def get_n_majors_from_2d_array(arr: list, n: int = 10):
        """
        :param arr: [[key_1, value_1], ..., [last_key, last_value]]
        :param n: number of major keys to return without the "other" key
        :return: [<n keys with largest values>, ["other", the sum of all non-included values]]
        """
        arr = sorted([[str(i[0]), int(i[1])] for i in arr], key=lambda x: x[1], reverse=True)
        if n > len(arr):
            raise ValueError("The size of provided list is too small: {}".format(len(arr)))
        return arr[:n] + [["others", sum([i[1] for i in arr[n:]])]]

    @staticmethod
    def digest_df(df: pd.DataFrame, associations: dict, *columns_with_keywords):
        df_columns = list(df)
        columns_with_keywords = [i for i in columns_with_keywords if len(i) > 0]
        if len(columns_with_keywords) == 0:
            raise ValueError("No column for keyword search specified")
        try:
            df["lookup_column"] = df.loc[:, columns_with_keywords].astype(str).apply(
                lambda x: " ".join([CounterWrapper.prepare_string(i) for i in x]), axis=1)
        except KeyError:
            print(list(df), associations, columns_with_keywords)
        keywords_series = []
        for main_word in associations:
            key_words = associations.get(main_word)
            if not key_words or len(key_words) == 0:
                key_words = (main_word, )
            sub_df = df.loc[df["lookup_column"].apply(lambda x: any(i.lower() in str(x) for i in key_words)) == True,
                            [i for i in df_columns if i not in columns_with_keywords]]
            keywords_series.append(sub_df.sum().rename(main_word))
        out_df = pd.concat(keywords_series, axis=1, sort=False)
        out_df.index.name = "sample_path"
        out_df = out_df.transpose().fillna(0)
        out_df.index.name = "keyword"
        return out_df
