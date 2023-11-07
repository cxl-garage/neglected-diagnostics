from pydantic import BaseModel


class AssayTargetColumns(BaseModel):
    target_id: str = "Target ID"
    assay_design_area: str = "Assay Design Area"
    perc_tgt_match: str = "% target match"
    ratio_tgt_mismatch: str = "Average of target mismatches"
    bpwise_error_percentage_tgt_mismatch: str = (
        "Bpwise error percentage of target mismatches"
    )
    num_off_tgt_mismatch: str = "Minimum # off-target mismatches"
