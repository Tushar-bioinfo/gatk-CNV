from pathlib import Path

dir = Path().cwd()

def cn(name:str):
    name = name.replace(".tumor",'')
    name = name.replace(".normal",'')
    return name

lst = [ (i.name,cn(i.name)) for i in dir.glob(".allelic.tsv")]

