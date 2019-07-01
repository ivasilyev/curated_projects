#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from meta.scripts.Utilities import Utilities


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
    DRUG_CLASSES = {"beta-lactam": (
        "amoxiclav", "ampicillin", "aztreonam", "cephalosporin", "imipinem", "meropenem", "penam", "penem"),
                    "aminoglycoside": ("amikacin", "netilmicin", "gentamicin"),
                    "epoxide": ("fosfomycin",),
                    "fluoroquinolone": ("сiprofloxacin",),
                    "glycopeptide antibiotic": (),
                    "lincosamide": (),
                    "macrolide": (),
                    "nitrofuran": ("nitrofurantoin",),
                    "nucleoside antibiotic": (),
                    "peptide antibiotic": (),
                    "phenicol": ("сhloramphenicol",),
                    "sulfonamide": ("trimoxazol",),
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
        keywords = [j for j in [i.strip() for i in keywords if isinstance(i, str)] if len(j) > 0]
        if split_words:
            keywords = Utilities.flatten_2d_array([i.split(" ") for i in keywords])
        return {j: () for j in sorted([i.strip() for i in set(keywords)]) if len(j) > 0}

    @staticmethod
    def generate_genera_dict(keywords: list):
        return DigestAssociationsKeeper.generate_keywords_dict(
            [Utilities.safe_findall("([A-Z][a-z]{4,})", i) for i in keywords])

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
    def prepare_string(s: str, lowercase: bool = False):
        import re
        out = re.sub(" +", " ", re.sub("[^A-Za-z0-9\-]+", " ", s)).strip()
        if lowercase:
            return out.lower()
        return out

    @staticmethod
    def prepare_list(lst: list, lowercase: bool = False):
        return Utilities.remove_empty_values([DigestAssociationsKeeper.prepare_string(i, lowercase) for i in lst])

    def digest_df(self, df: pd.DataFrame, associations: dict, columns_with_keywords: list, include_key: bool = True,
                  all_in_lowercase: bool = False, strict: bool = False):
        """
        :param df: Pandas DataFrame object containing an index, keyword columns and value columns
        :param associations: Dictionary '{key: (keywords...)}'
        :param columns_with_keywords: List of columns to search keywords
        :param include_key: Should the key of keyword group be included?
        :param all_in_lowercase: Convert both strings to lowercase?
        :param strict: Only count full match
        :return: Pandas DataFrame object with keys as index and columns sums as values and dictionary with corresponding
                 intermediate grouped Pandas DataFrame objects
        """
        def __regular_search(s: str):
            return any(i in str(s) for i in key_words)

        def __strict_search(s: str):
            return any(i == str(s) for i in key_words)

        df = df.copy()
        df_columns = list(df)
        columns_with_keywords = Utilities.remove_empty_values([i for i in columns_with_keywords])
        columns_without_keywords = Utilities.remove_empty_values(
            [i for i in df_columns if i not in columns_with_keywords])
        if len(columns_with_keywords) == 0:
            print("No column for keyword search specified!")
            return
        try:
            # 'columns_with_keywords' might be more than 1
            df["lookup_column"] = df.loc[:, columns_with_keywords].astype(str).apply(
                lambda x: " ".join(self.prepare_list(x, lowercase=all_in_lowercase)), axis=1)
        except KeyError as e:
            print(e, list(df), associations, columns_with_keywords)
        keywords_series = []
        raw_values_ds = pd.DataFrame()
        for main_word in associations:
            key_words = associations.get(main_word)
            if not key_words or len(key_words) == 0:
                key_words = ()
            if include_key:
                key_words = list(key_words) + [main_word, ]
            key_words = sorted(set(self.prepare_list(key_words, lowercase=all_in_lowercase)))
            if len(key_words) == 0:
                raise ValueError("No values to search: '{}: {}'".format(main_word, key_words))
            if strict:
                df_with_keywords = df.loc[df["lookup_column"].apply(__strict_search) == True,
                                          columns_without_keywords]
            else:
                df_with_keywords = df.loc[df["lookup_column"].apply(__regular_search) == True,
                                          columns_without_keywords]
            keywords_series.append(df_with_keywords.sum().rename(main_word))
            # Reset index to avoid exceptions thrown by duplicates
            raw_values_df = df_with_keywords.reset_index()
            raw_values_df["keyword"] = main_word
            if raw_values_ds.shape[0] == 0:
                raw_values_ds = raw_values_df
            else:
                raw_values_ds = pd.concat([raw_values_ds, raw_values_df], axis=0, ignore_index=True)
        out_df = Utilities.merge_pd_series_list(keywords_series).fillna(0)
        out_df.columns.name = "value"
        out_df.index.name = "keyword"
        return out_df, raw_values_ds
