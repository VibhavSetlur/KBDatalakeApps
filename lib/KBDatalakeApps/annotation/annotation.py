import os

def parse_kofam(data: str):
    annotation = {}
    for line in data.split('\n'):
        if line:
            _parts = line.strip().split()
            if len(_parts) == 2:
                feature_id, ko = _parts
                annotation[feature_id] = ko
    return annotation

def parse_bakta(data: str):
    pass

def parse_psortb(data: str):
    
    pass
