from pydantic import BaseModel


class AssayTargetColumns(BaseModel):
    assay_design_area: str = "Assay Design Area"
    perc_tgt_match: str = "% target match"
    ratio_tgt_mismatch: str = "ratio of target mismatches"
    num_off_tgt_mismatch: str = "# off-target mismatches"
