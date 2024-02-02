from pydantic import BaseModel


class GroupSequenceColumns(BaseModel):
    id: str = "SequenceID"
    seq: str = "Sequence"
    count: str = "Count"
    group_number: str = "Group_Number"


class MultiSequenceAlginmentColumns(BaseModel):
    id: str = "SequenceID"
    seq: str = "Sequence"
