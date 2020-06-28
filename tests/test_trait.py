import pytest
import pandas as pd

from entity import Trait, UKBiobankTrait, Study, GTEXGWASTrait


class TestTrait:
    @pytest.mark.parametrize(
        "trait_code,trait_category",
        [
            ("50_raw", "Body size measures"),
            ("22617_1222", "Employment history"),
            ("G54", "Diseases (ICD10 main)"),
            ("M13_RHEUMATISM", "Diseases (FinnGen)"),
            ("20002_1473", "Diseases (cardiovascular)"),
            ("20002_1265", "Diseases (neurology/eye/psychiatry)"),
            ("6150_1", "Diseases (cardiovascular)"),
            ("6152_5", "Diseases (cardiovascular)"),
            ("6152_7", "Diseases (cardiovascular)"),
            ("3627_raw", "Diseases (cardiovascular)"),
            ("6152_6", "Diseases (respiratory/ent)"),
            ("6152_8", "Diseases (respiratory/ent)"),
            ("6151_4", "Diseases (musculoskeletal/trauma)"),
            ("6151_6", "Diseases (musculoskeletal/trauma)"),
            ("6152_9", "Diseases (allergies)"),
            ("20001_1068", "Cancer (other)"),
            ("20001_1020", "Cancer (gastrointestinal)"),
            # ("UKB_1160_Sleep_duration", "Psychiatric-neurologic"),
            # ("", ""),
        ],
    )
    def test_ukb_trait_category(self, trait_code, trait_category):
        trait = UKBiobankTrait(code=trait_code)
        assert trait.code == trait_code
        assert trait.category == trait_category

    @pytest.mark.parametrize(
        "trait_code,trait_category",
        [
            ("UKB_1160_Sleep_duration", "Psychiatric-neurologic"),
            (
                "UKB_20002_1154_self_reported_irritable_bowel_syndrome",
                "Digestive system disease",
            ),
        ],
    )
    def test_gtexgwas_trait_category(self, trait_code, trait_category):
        trait = GTEXGWASTrait(code=trait_code)
        assert trait.code == trait_code
        assert trait.category == trait_category

    def test_ukb_trait_no_cases_or_controls(self):
        pheno_from_code = UKBiobankTrait(code="50_raw")
        assert pheno_from_code is not None
        assert pheno_from_code.code == "50_raw"
        assert pheno_from_code.full_code == "50_raw-Standing_height"
        assert pheno_from_code.description == "Standing height"
        assert pheno_from_code.type == "continuous_raw"
        assert pheno_from_code.n == 360388
        assert pd.isnull(pheno_from_code.n_cases)
        assert pd.isnull(pheno_from_code.n_controls)
        assert pheno_from_code.source == "UK Biobank"
        assert pheno_from_code.study == Study.UK_BIOBANK
        assert pheno_from_code.get_plain_name() == "50_raw-Standing_height"

        pheno_from_full_code = UKBiobankTrait(
            full_code=pheno_from_code.get_plain_name()
        )
        assert pheno_from_code.code == pheno_from_full_code.code
        assert pheno_from_code.description == pheno_from_full_code.description
        assert pheno_from_code.type == pheno_from_full_code.type
        assert pheno_from_code.n == pheno_from_full_code.n
        assert pd.isnull(pheno_from_full_code.n_cases)
        assert pd.isnull(pheno_from_full_code.n_controls)
        assert pheno_from_code.source == pheno_from_full_code.source
        assert pheno_from_code.study == pheno_from_full_code.study
        assert pheno_from_code.get_plain_name() == pheno_from_full_code.get_plain_name()

    def test_ukb_trait_with_cases_and_controls(self):
        pheno_from_code = UKBiobankTrait(code="G54")
        assert pheno_from_code is not None
        assert pheno_from_code.code == "G54"
        assert (
            pheno_from_code.full_code
            == "G54-Diagnoses_main_ICD10_G54_Nerve_root_and_plexus_disorders"
        )
        assert (
            pheno_from_code.full_code
            == "G54-Diagnoses_main_ICD10_G54_Nerve_root_and_plexus_disorders"
        )
        assert (
            pheno_from_code.description
            == "Diagnoses - main ICD10: G54 Nerve root and plexus disorders"
        )
        assert pheno_from_code.type == "categorical"
        assert pheno_from_code.n == 361194
        assert pheno_from_code.n_cases == 143
        assert pheno_from_code.n_controls == 361051
        assert pheno_from_code.source == "UK Biobank"
        assert pheno_from_code.study == Study.UK_BIOBANK
        assert (
            pheno_from_code.get_plain_name()
            == "G54-Diagnoses_main_ICD10_G54_Nerve_root_and_plexus_disorders"
        )

        pheno_from_full_code = UKBiobankTrait(
            full_code=pheno_from_code.get_plain_name()
        )
        assert pheno_from_code.code == pheno_from_full_code.code
        assert pheno_from_code.description == pheno_from_full_code.description
        assert pheno_from_code.type == pheno_from_full_code.type
        assert pheno_from_code.n == pheno_from_full_code.n
        assert pheno_from_code.n_cases == pheno_from_full_code.n_cases
        assert pheno_from_code.n_controls == pheno_from_full_code.n_controls
        assert pheno_from_code.source == pheno_from_full_code.source
        assert pheno_from_code.study == pheno_from_full_code.study
        assert pheno_from_code.get_plain_name() == pheno_from_full_code.get_plain_name()

    def test_gtex_gwas_trait_no_cases_or_controls(self):
        pheno_from_code = GTEXGWASTrait(code="UKB_1160_Sleep_duration")
        assert pheno_from_code is not None
        assert pheno_from_code.code == "UKB_1160_Sleep_duration"
        assert pheno_from_code.description == "Sleep Duration UKB"
        assert pheno_from_code.type == "continuous_raw"
        assert pheno_from_code.n == 337119
        assert pd.isnull(pheno_from_code.n_cases)
        assert pd.isnull(pheno_from_code.n_controls)
        assert pheno_from_code.source == "UK Biobank"
        assert pheno_from_code.study == Study.GTEX_GWAS
        assert pheno_from_code.get_plain_name() == "UKB_1160_Sleep_duration"

    def test_gtex_gwas_trait_with_cases_and_controls(self):
        pheno_from_code = GTEXGWASTrait(
            code="UKB_20002_1094_self_reported_deep_venous_thrombosis_dvt"
        )
        assert pheno_from_code is not None
        assert (
            pheno_from_code.code
            == "UKB_20002_1094_self_reported_deep_venous_thrombosis_dvt"
        )
        assert pheno_from_code.description == "Deep Venous Thrombosis UKBS"
        assert pheno_from_code.type == "binary"
        assert pheno_from_code.n == 337119
        assert pheno_from_code.n_cases == 6767
        assert pheno_from_code.n_controls == 337119 - 6767
        assert pheno_from_code.source == "UK Biobank"
        assert pheno_from_code.study == Study.GTEX_GWAS
        assert (
            pheno_from_code.get_plain_name()
            == "UKB_20002_1094_self_reported_deep_venous_thrombosis_dvt"
        )

    def test_trait_get_trait_using_code_is_ukb(self):
        trait_code = "50_raw"

        trait = Trait.get_trait(trait_code)
        assert trait is not None
        assert isinstance(trait, UKBiobankTrait)
        assert trait.code == "50_raw"
        assert trait.get_plain_name() == "50_raw-Standing_height"

    def test_trait_get_trait_using_full_code_is_ukb(self):
        trait_code = "50_raw-Standing_height"

        trait = Trait.get_trait(full_code=trait_code)
        assert trait is not None
        assert isinstance(trait, UKBiobankTrait)
        assert trait.code == "50_raw"
        assert trait.get_plain_name() == "50_raw-Standing_height"

    def test_trait_get_trait_is_gtex_gwas(self):
        trait_code = "UKB_1160_Sleep_duration"

        trait = Trait.get_trait(trait_code)
        assert trait is not None
        assert isinstance(trait, GTEXGWASTrait)
        assert trait.code == "UKB_1160_Sleep_duration"
        assert trait.get_plain_name() == "UKB_1160_Sleep_duration"

    def test_trait_get_trait_using_full_code_is_gtex_gwas(self):
        trait_code = "UKB_1160_Sleep_duration"

        trait = Trait.get_trait(full_code=trait_code)
        assert trait is not None
        assert isinstance(trait, GTEXGWASTrait)
        assert trait.code == "UKB_1160_Sleep_duration"
        assert trait.get_plain_name() == "UKB_1160_Sleep_duration"
