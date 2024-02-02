from pydantic import BaseModel


class AssayTargetColumns(BaseModel):
    target_id: str = "Assay Design Area ID"
    assay_design_area: str = "Assay Design Area"
    assay_design_area_start: str = "Assay Design Area Start Index"
    assay_design_area_end: str = "Assay Design Area End Index"
    assay_design_area_reference: str = "Assay Design Area Reference"
    perc_tgt_match: str = "% target match"
    ratio_tgt_mismatch: str = "Average of target mismatches"
    bpwise_error_percentage_tgt_mismatch: str = (
        "BP error percentage of target mismatches"
    )
    num_off_tgt_mismatch: str = "Minimum # off-target mismatches"
