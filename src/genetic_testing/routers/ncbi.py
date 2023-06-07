from Bio import Entrez

# Entrez.api_key = ""  # Add your API Key
Entrez.email = "aniketf@uw.edu"  # Add your email address

# TODO: Add docstrings and comments
def get(operation, db, term):
    handle = None
    record = None
    if operation == "Search":
        handle = Entrez.esearch(db=db, term=term)

    if handle:
        record = Entrez.read(handle)
    
    return record
